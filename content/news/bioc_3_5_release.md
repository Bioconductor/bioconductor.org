April 25, 2017

Bioconductors:

We are pleased to announce Bioconductor 3.5, consisting of 1383
software packages, 316 experiment data packages, and 911 annotation
packages.

There are 88 new software packages, and many updates and improvements
to existing packages; Bioconductor 3.5 is compatible with R 3.4,
and is supported on Linux, 32- and 64-bit Windows, and Mac OS X.  This
release will include an updated Bioconductor [Amazon Machine Image][1]
and [Docker containers][2].

Visit [https://bioconductor.org][3]
for details and downloads.

[1]: https://bioconductor.org/help/bioconductor-cloud-ami/
[2]: https://bioconductor.org/help/docker/
[3]: https://bioconductor.org

Contents
--------

* [Getting Started with Bioconductor 3.5](#getting-started-with-bioconductor-35)
* [New Software Packages](#new-software-packages)
* [NEWS from new and existing packages](#news-from-new-and-existing-packages)
* [Deprecated and Defunct Packages](#deprecated-and-defunct-packages)

Getting Started with Bioconductor 3.5
======================================

To update to or install Bioconductor 3.5:

1. Install R 3.4.  Bioconductor 3.5 has been designed expressly for
   this version of R.

2. Follow the instructions at
   [http://bioconductor.org/install/](http://bioconductor.org/install/).

New Software Packages
=====================

There are 88 new software packages in this release of Bioconductor.

- [AnnotationFilter](https://bioconductor.org/packages/AnnotationFilter)
  This package provides class and other infrastructure to implement
  filters for manipulating Bioconductor annotation resources. The
  filters will be used by ensembldb, Organism.dplyr, and other
  packages.

- [ATACseqQC](https://bioconductor.org/packages/ATACseqQC) ATAC-seq,
  an assay for Transposase-Accessible Chromatin using sequencing, is
  a rapid and sensitive method for chromatin accessibility analysis.
  It was developed as an alternative method to MNase-seq, FAIRE-seq
  and DNAse-seq. Comparing to the other methods, ATAC-seq requires
  less amount of the biological samples and time to process. In the
  process of analyzing several ATAC-seq dataset produced in our labs,
  we learned some of the unique aspects of the quality assessment for
  ATAC-seq data.To help users to quickly assess whether their
  ATAC-seq experiment is successful, we developed ATACseqQC package
  partially following the guideline published in Nature Method 2013
  (Greenleaf et al.), including diagnostic plot of fragment size
  distribution, proportion of mitochondria reads, nucleosome
  positioning pattern, and CTCF or other Transcript Factor
  footprints.

- [banocc](https://bioconductor.org/packages/banocc) BAnOCC is a
  package designed for compositional data, where each sample sums to
  one. It infers the approximate covariance of the unconstrained data
  using a Bayesian model coded with `rstan`. It provides as output
  the `stanfit` object as well as posterior median and credible
  interval estimates for each correlation element.

- [basecallQC](https://bioconductor.org/packages/basecallQC) The
  basecallQC package provides tools to work with Illumina bcl2Fastq
  (versions >= 2.1.7) software.Prior to basecalling and
  demultiplexing using the bcl2Fastq software, basecallQC functions
  allow the user to update Illumina sample sheets from versions <=
  1.8.9 to >= 2.1.7 standards, clean sample sheets of common problems
  such as invalid sample names and IDs, create read and index
  basemasks and the bcl2Fastq command. Following the generation of
  basecalled and demultiplexed data, the basecallQC packages allows
  the user to generate HTML tables, plots and a self contained report
  of summary metrics from Illumina XML output files.

- [BiocFileCache](https://bioconductor.org/packages/BiocFileCache)
  This package creates a persistent on-disk cache of files that the
  user can add, update, and retrieve. It is useful for managing
  resources (such as custom Txdb objects) that are costly or
  difficult to create, web resources, and data files used across
  sessions.

- [BioCor](https://bioconductor.org/packages/BioCor) Calculates
  functional similarities based on the pathways described on KEGG and
  REACTOME or in gene sets. These similarities can be calculated for
  pathways or gene sets, genes, or clusters and combined with other
  similarities. They can be used to improve networks, gene selection,
  testing relationships...

- [BioMedR](https://bioconductor.org/packages/BioMedR) The BioMedR
  package offers an R/Bioconductor package generating various
  molecular representations for chemicals, proteins, DNAs/RNAs and
  their interactions.

- [biotmle](https://bioconductor.org/packages/biotmle) This package
  facilitates the discovery of biomarkers from biological sequencing
  data (e.g., microarrays, RNA-seq) based on the associations of
  potential biomarkers with exposure and outcome variables by
  implementing an estimation procedure that combines a generalization
  of the moderated t-statistic with asymptotically linear statistical
  parameters estimated via targeted minimum loss-based estimation
  (TMLE).

- [BLMA](https://bioconductor.org/packages/BLMA) Suit of tools for
  bi-level meta-analysis. The package can be used in a wide range of
  applications, including general hypothesis testings, differential
  expression analysis, functional analysis, and pathway analysis.

- [BPRMeth](https://bioconductor.org/packages/BPRMeth) BPRMeth
  package uses the Binomial Probit Regression likelihood to model
  methylation profiles and extract higher order features. These
  features quantitate precisely notions of shape of a methylation
  profile. Using these higher order features across promoter-proximal
  regions, we construct a powerful predictor of gene expression.
  Also, these features are used to cluster proximal-promoter regions
  using the EM algorithm.

- [branchpointer](https://bioconductor.org/packages/branchpointer)
  Predicts branchpoint probability for sites in intronic branchpoint
  windows. Queries can be supplied as intronic regions; or to
  evaluate the effects of mutations, SNPs.

- [BUMHMM](https://bioconductor.org/packages/BUMHMM) This is a
  probabilistic modelling pipeline for computing per- nucleotide
  posterior probabilities of modification from the data collected in
  structure probing experiments. The model supports multiple
  experimental replicates and empirically corrects coverage- and
  sequence-dependent biases. The model utilises the measure of a
  "drop-off rate" for each nucleotide, which is compared between
  replicates through a log-ratio (LDR). The LDRs between control
  replicates define a null distribution of variability in drop-off
  rate observed by chance and LDRs between treatment and control
  replicates gets compared to this distribution. Resulting empirical
  p-values (probability of being "drawn" from the null distribution)
  are used as observations in a Hidden Markov Model with a
  Beta-Uniform Mixture model used as an emission model. The resulting
  posterior probabilities indicate the probability of a nucleotide of
  having being modified in a structure probing experiment.

- [CATALYST](https://bioconductor.org/packages/CATALYST) Mass
  cytometry (CyTOF) uses heavy metal isotopes rather than fluorescent
  tags as reporters to label antibodies, thereby substantially
  decreasing spectral overlap and allowing for examination of over 50
  parameters at the single cell level. While spectral overlap is
  significantly less pronounced in CyTOF than flow cytometry,
  spillover due to detection sensitivity, isotopic impurities, and
  oxide formation can impede data interpretability. We designed
  CATALYST (Cytometry dATa anALYSis Tools) to provide a pipeline for
  preprocessing of cytometry data, including i) normalization using
  bead standards, ii) single-cell deconvolution, and iii) bead-based
  compensation.

- [cellbaseR](https://bioconductor.org/packages/cellbaseR) This R
  package makes use of the exhaustive RESTful Web service API that
  has been implemented for the Cellabase database. It enable
  researchers to query and obtain a wealth of biological information
  from a single database saving a lot of time. Another benefit is
  that researchers can easily make queries about different biological
  topics and link all this information together as all information is
  integrated.

- [cellscape](https://bioconductor.org/packages/cellscape) CellScape
  facilitates interactive browsing of single cell clonal evolution
  datasets. The tool requires two main inputs: (i) the genomic
  content of each single cell in the form of either copy number
  segments or targeted mutation values, and (ii) a single cell
  phylogeny. Phylogenetic formats can vary from dendrogram-like
  phylogenies with leaf nodes to evolutionary model-derived
  phylogenies with observed or latent internal nodes. The CellScape
  phylogeny is flexibly input as a table of source-target edges to
  support arbitrary representations, where each node may or may not
  have associated genomic data. The output of CellScape is an
  interactive interface displaying a single cell phylogeny and a
  cell-by-locus genomic heatmap representing the mutation status in
  each cell for each locus.

- [chimeraviz](https://bioconductor.org/packages/chimeraviz)
  chimeraviz manages data from fusion gene finders and provides
  useful visualization tools.

- [ChIPexoQual](https://bioconductor.org/packages/ChIPexoQual)
  Package with a quality control pipeline for ChIP-exo/nexus data.

- [clusterSeq](https://bioconductor.org/packages/clusterSeq)
  Identification of clusters of co-expressed genes based on their
  expression across multiple (replicated) biological samples.

- [coseq](https://bioconductor.org/packages/coseq) Co-expression
  analysis for expression profiles arising from high-throughput
  sequencing data. Feature (e.g., gene) profiles are clustered using
  adapted transformations and mixture models or a K-means algorithm,
  and model selection criteria (to choose an appropriate number of
  clusters) are provided.

- [cydar](https://bioconductor.org/packages/cydar) Identifies
  differentially abundant populations between samples and groups in
  mass cytometry data. Provides methods for counting cells into
  hyperspheres, controlling the spatial false discovery rate, and
  visualizing changes in abundance in the high-dimensional marker
  space.

- [DaMiRseq](https://bioconductor.org/packages/DaMiRseq) The DaMiRseq
  package offers a tidy pipeline of data mining procedures to
  identify transcriptional biomarkers and exploit them for
  classification purposes.. The package accepts any kind of data
  presented as a table of raw counts and allows including covariates
  that occur with the experimental setting. A series of functions
  enable the user to clean up the data by filtering genomic features
  and samples, to adjust data by identifying and removing the
  unwanted source of variation (i.e. batches and confounding factors)
  and to select the best predictors for modeling. Finally, a
  ``Stacking'' ensemble learning technique is applied to build a
  robust classification model. Every step includes a checkpoint that
  the user may exploit to assess the effects of data management by
  looking at diagnostic plots, such as clustering and heatmaps, RLE
  boxplots, MDS or correlation plot.

- [DelayedArray](https://bioconductor.org/packages/DelayedArray)
  Wrapping an array-like object (typically an on-disk object) in a
  DelayedArray object allows one to perform common array operations
  on it without loading the object in memory. In order to reduce
  memory usage and optimize performance, operations on the object are
  either delayed or executed using a block processing mechanism. Note
  that this also works on in-memory array-like objects like DataFrame
  objects (typically with Rle columns), Matrix objects, and ordinary
  arrays and data frames.

- [discordant](https://bioconductor.org/packages/discordant)
  Discordant is a method to determine differential correlation of
  molecular feature pairs from -omics data using mixture models.
  Algorithm is explained further in Siska et al.

- [DMRScan](https://bioconductor.org/packages/DMRScan) This package
  detects significant differentially methylated regions (for both
  qualitative and quantitative traits), using a scan statistic with
  underlying Poisson heuristics. The scan statistic will depend on a
  sequence of window sizes (# of CpGs within each window) and on a
  threshold for each window size. This threshold can be calculated by
  three different means: i) analytically using Siegmund et.al (2012)
  solution (preferred), ii) an important sampling as suggested by
  Zhang (2008), and a iii) full MCMC modeling of the data, choosing
  between a number of different options for modeling the dependency
  between each CpG.

- [epiNEM](https://bioconductor.org/packages/epiNEM) epiNEM is an
  extension of the original Nested Effects Models (NEM). EpiNEM is
  able to take into account double knockouts and infer more complex
  network signalling pathways.

- [EventPointer](https://bioconductor.org/packages/EventPointer)
  EventPointer is an R package to identify alternative splicing
  events that involve either simple (case-control experiment) or
  complex experimental designs such as time course experiments and
  studies including paired-samples. The algorithm can be used to
  analyze data from either junction arrays (Affymetrix Arrays) or
  sequencing data (RNA-Seq). The software returns a data.frame with
  the detected alternative splicing events: gene name, type of event
  (cassette, alternative 3',...,etc), genomic position, statistical
  significance and increment of the percent spliced in (Delta PSI)
  for all the events. The algorithm can generate a series of files to
  visualize the detected alternative splicing events in IGV. This
  eases the interpretation of results and the design of primers for
  standard PCR validation.

- [flowTime](https://bioconductor.org/packages/flowTime) This package
  was developed for analysis of both dynamic and steady state
  experiments examining the function of gene regulatory networks in
  yeast (strain W303) expressing fluorescent reporter proteins using
  a BD Accuri C6 and SORP cytometers. However, the functions are for
  the most part general and may be adapted for analysis of other
  organisms using other flow cytometers. Functions in this package
  facilitate the annotation of flow cytometry data with experimental
  metadata, as is requisite for dissemination and general
  ease-of-use. Functions for creating, saving and loading gate sets
  are also included. In the past, we have typically generated summary
  statistics for each flowset for each timepoint and then annotated
  and analyzed these summary statistics. This method loses a great
  deal of the power that comes from the large amounts of individual
  cell data generated in flow cytometry, by essentially collapsing
  this data into a bulk measurement after subsetting. In addition to
  these summary functions, this package also contains functions to
  facilitate annotation and analysis of steady-state or time-lapse
  data utilizing all of the data collected from the thousands of
  individual cells in each sample.

- [funtooNorm](https://bioconductor.org/packages/funtooNorm) Provides
  a function to normalize Illumina Infinium Human Methylation 450
  BeadChip (Illumina 450K), correcting for tissue and/or cell type.

- [GA4GHclient](https://bioconductor.org/packages/GA4GHclient)
  GA4GHclient provides an easy way to access public data servers
  through Global Alliance for Genomics and Health (GA4GH) genomics
  API. It provides low-level access to GA4GH API and translates
  response data into Bioconductor-based class objects.

- [gcapc](https://bioconductor.org/packages/gcapc) Peak calling for
  ChIP-seq data with consideration of potential GC bias in sequencing
  reads. GC bias is first estimated with generalized linear mixture
  models using weighted GC strategy, then applied into peak
  significance estimation.

-
  [geneClassifiers](https://bioconductor.org/packages/geneClassifiers)
  This packages aims for easy accessible application of classifiers
  which have been published in literature using an ExpressionSet as
  input.

- [GenomicDataCommons](https://bioconductor.org/packages/GenomicDataCommons)
  Programmatically access the NIH / NCI Genomic Data Commons RESTful
  service.

- [GenomicScores](https://bioconductor.org/packages/GenomicScores)
  Provide infrastructure to store and access genomewide
  position-specific scores within R and Bioconductor.

- [GISPA](https://bioconductor.org/packages/GISPA) GISPA is a method
  intended for the researchers who are interested in defining gene
  sets with similar, a priori specified molecular profile. GISPA
  method has been previously published in Nucleic Acid Research
  (Kowalski et al., 2016; PMID: 26826710).

- [goSTAG](https://bioconductor.org/packages/goSTAG) Gene lists
  derived from the results of genomic analyses are rich in biological
  information. For instance, differentially expressed genes (DEGs)
  from a microarray or RNA-Seq analysis are related functionally in
  terms of their response to a treatment or condition. Gene lists can
  vary in size, up to several thousand genes, depending on the
  robustness of the perturbations or how widely different the
  conditions are biologically. Having a way to associate biological
  relatedness between hundreds and thousands of genes systematically
  is impractical by manually curating the annotation and function of
  each gene. Over-representation analysis (ORA) of genes was
  developed to identify biological themes. Given a Gene Ontology (GO)
  and an annotation of genes that indicate the categories each one
  fits into, significance of the over-representation of the genes
  within the ontological categories is determined by a Fisher's exact
  test or modeling according to a hypergeometric distribution.
  Comparing a small number of enriched biological categories for a
  few samples is manageable using Venn diagrams or other means for
  assessing overlaps. However, with hundreds of enriched categories
  and many samples, the comparisons are laborious. Furthermore, if
  there are enriched categories that are shared between samples,
  trying to represent a common theme across them is highly
  subjective. goSTAG uses GO subtrees to tag and annotate genes
  within a set. goSTAG visualizes the similarities between the
  over-representation of DEGs by clustering the p-values from the
  enrichment statistical tests and labels clusters with the GO term
  that has the most paths to the root within the subtree generated
  from all the GO terms in the cluster.

- [GRridge](https://bioconductor.org/packages/GRridge) This package
  allows the use of multiple sources of co-data (e.g. external
  p-values, gene lists, annotation) to improve prediction of binary,
  continuous and survival response using (logistic, linear or Cox)
  group-regularized ridge regression. It also facilitates post-hoc
  variable selection and prediction diagnostics by cross-validation
  using ROC curves and AUC.

- [heatmaps](https://bioconductor.org/packages/heatmaps) This package
  provides functions for plotting heatmaps of genome-wide data across
  genomic intervals, such as ChIP-seq signals at peaks or across
  promoters. Many functions are also provided for investigating
  sequence features.

- [hicrep](https://bioconductor.org/packages/hicrep) Hi-C is a
  powerful technology for studying genome-wide chromatin
  interactions. However, current methods for assessing Hi-C data
  reproducibility can produce misleading results because they ignore
  spatial features in Hi-C data, such as domain structure and
  distance-dependence. We present a novel reproducibility measure
  that systematically takes these features into consideration. This
  measure can assess pairwise differences between Hi-C matrices under
  a wide range of settings, and can be used to determine optimal
  sequencing depth. Compared to existing approaches, it consistently
  shows higher accuracy in distinguishing subtle differences in
  reproducibility and depicting interrelationships of cell lineages
  than existing approaches. This R package `hicrep` implements our
  approach.

- [ideal](https://bioconductor.org/packages/ideal) This package
  provides functions for an Interactive Differential Expression
  AnaLysis of RNA-sequencing datasets, to extract quickly and
  effectively information downstream the step of differential
  expression. A Shiny application encapsulates the whole package.

- [IMAS](https://bioconductor.org/packages/IMAS) Integrative analysis
  of Multi-omics data for Alternative splicing.

- [ImpulseDE2](https://bioconductor.org/packages/ImpulseDE2)
  ImpulseDE2 is a differential expression algorithm for longitudinal
  count data sets which arise in sequencing experiments such as
  RNA-seq, ChIP-seq, ATAC-seq and DNaseI-seq. ImpulseDE2 is based on
  a negative binomial noise model with dispersion trend smoothing by
  DESeq2 and uses the impulse model to constrain the mean expression
  trajectory of each gene. The impulse model was empirically found to
  fit global expression changes in cells after environmental and
  developmental stimuli and is therefore appropriate in most cell
  biological scenarios. The constraint on the mean expression
  trajectory prevents overfitting to small expression fluctuations.
  Secondly, ImpulseDE2 has higher statistical testing power than
  generalized linear model-based differential expression algorithms
  which fit time as a categorial variable if more than six time
  points are sampled because of the fixed number of parameters.

- [IntEREst](https://bioconductor.org/packages/IntEREst) This package
  performs Intron-Exon Retention analysis on RNA-seq data (.bam
  files).

- [IWTomics](https://bioconductor.org/packages/IWTomics)
  Implementation of the Interval-Wise Testing (IWT) for omics data.
  This inferential procedure tests for differences in "Omics" data
  between two groups of genomic regions (or between a group of
  genomic regions and a reference center of symmetry), and does not
  require fixing location and scale at the outset.

- [karyoploteR](https://bioconductor.org/packages/karyoploteR)
  karyoploteR creates karyotype plots of arbitrary genomes and offers
  a complete set of functions to plot arbitrary data on them. It
  mimicks many R base graphics functions coupling them with a
  coordinate change function automatically mapping the chromosome and
  data coordinates into the plot coordinates. In addition to the
  provided data plotting functions, it is easy to add new ones.

- [Logolas](https://bioconductor.org/packages/Logolas) Produces logo
  plots of a variety of symbols and names comprising English
  alphabets, numerics and punctuations. Can be used for sequence
  motif generation, mutation pattern generation, protein amino acid
  geenration and symbol strength representation in any generic
  context.

- [mapscape](https://bioconductor.org/packages/mapscape) MapScape
  integrates clonal prevalence, clonal hierarchy, anatomic and
  mutational information to provide interactive visualization of
  spatial clonal evolution. There are four inputs to MapScape: (i)
  the clonal phylogeny, (ii) clonal prevalences, (iii) an image
  reference, which may be a medical image or drawing and (iv) pixel
  locations for each sample on the referenced image. Optionally,
  MapScape can accept a data table of mutations for each clone and
  their variant allele frequencies in each sample. The output of
  MapScape consists of a cropped anatomical image surrounded by two
  representations of each tumour sample. The first, a cellular
  aggregate, visually displays the prevalence of each clone. The
  second shows a skeleton of the clonal phylogeny while highlighting
  only those clones present in the sample. Together, these
  representations enable the analyst to visualize the distribution of
  clones throughout anatomic space.

- [MaxContrastProjection](https://bioconductor.org/packages/MaxContrastProjection)
  A problem when recording 3D fluorescent microscopy images is how to
  properly present these results in 2D. Maximum intensity projections
  are a popular method to determine the focal plane of each pixel in
  the image. The problem with this approach, however, is that
  out-of-focus elements will still be visible, making edges and fine
  structures difficult to detect. This package aims to resolve this
  problem by using the contrast around a given pixel to determine the
  focal plane, allowing for a much cleaner structure detection than
  would be otherwise possible. For convenience, this package also
  contains functions to perform various other types of projections,
  including a maximum intensity projection.

- [MCbiclust](https://bioconductor.org/packages/MCbiclust) Custom
  made algorithm and associated methods for finding, visualising and
  analysing biclusters in large gene expression data sets. Algorithm
  is based on with a supplied gene set of size n, finding the maximum
  strength correlation matrix containing m samples from the data set.

- [metavizr](https://bioconductor.org/packages/metavizr) This package
  provides Websocket communication to the metaviz web app
  (http://metaviz.cbcb.umd.edu) for interactive visualization of
  metagenomics data. Objects in R/bioc interactive sessions can be
  displayed in plots and data can be explored using a facetzoom
  visualization. Fundamental Bioconductor data structures are
  supported (e.g., MRexperiment objects), while providing an easy
  mechanism to support other data structures. Visualizations (using
  d3.js) can be easily added to the web app as well.

- [methylInheritance](https://bioconductor.org/packages/methylInheritance)
  Permutation analysis, based on Monte Carlo sampling, for testing
  the hypothesis that the number of conserved differentially
  methylated elements, between several generations, is associated to
  an effect inherited from a treatment and that stochastic effect can
  be dismissed.

- [MIGSA](https://bioconductor.org/packages/MIGSA) Massive and
  Integrative Gene Set Analysis. The MIGSA package allows to perform
  a massive and integrative gene set analysis over several expression
  and gene sets simultaneously. It provides a common gene expression
  analytic framework that grants a comprehensive and coherent
  analysis. Only a minimal user parameter setting is required to
  perform both singular and gene set enrichment analyses in an
  integrative manner by means of the best available methods, i.e.
  dEnricher and mGSZrespectively. The greatest strengths of this big
  omics data tool are the availability of several functions to
  explore, analyze and visualize its results in order to facilitate
  the data mining task over huge information sources. MIGSA package
  also provides several functions that allow to easily load the most
  updated gene sets from several repositories.

- [mimager](https://bioconductor.org/packages/mimager) Easily
  visualize and inspect microarrays for spatial artifacts.

- [motifcounter](https://bioconductor.org/packages/motifcounter)
  'motifcounter' provides functionality to compute the statistics
  related with motif matching and counting of motif matches in DNA
  sequences. As an input, 'motifcounter' requires a motif in terms of
  a position frequency matrix (PFM). Furthermore, a set of DNA
  sequences is required to estimated a higher-order background model
  (BGM). The package provides functions to investigate the the
  per-position and per strand log-likelihood scores between the PFM
  and the BGM across a given sequence of set of sequences.
  Furthermore, the package facilitates motif matching based on an
  automatically derived score threshold. To this end the distribution
  of scores is efficiently determined and the score threshold is
  chosen for a user-prescribed significance level. This allows to
  control for the false positive rate. Moreover, 'motifcounter'
  implements a motif match enrichment test based on two the number of
  motif matches that are expected in random DNA sequences. Motif
  enrichment is facilitated by either a compound Poisson
  approximation or a combinatorial approximation of the motif match
  counts. Both models take higher-order background models, the
  motif's self-similarity, and hits on both DNA strands into account.
  The package is in particular useful for long motifs and/or relaxed
  choices of score thresholds, because the implemented algorithms
  efficiently bypass the need for enumerating a (potentially huge)
  set of DNA words that can give rise to a motif match.

- [msgbsR](https://bioconductor.org/packages/msgbsR) Pipeline for the
  anaysis of a MS-GBS experiment.

- [multiOmicsViz](https://bioconductor.org/packages/multiOmicsViz)
  Calculate the spearman correlation between the source omics data
  and other target omics data, identify the significant correlations
  and plot the significant correlations on the heat map in which the
  x-axis and y-axis are ordered by the chromosomal location.

- [MWASTools](https://bioconductor.org/packages/MWASTools) MWAS
  provides a complete pipeline to perform metabolome-wide association
  studies. Key functionalities of the package include: quality
  control analysis of metabonomic data; MWAS using different
  association models (partial correlations; generalized linear
  models); model validation using non-parametric bootstrapping;
  visualization of MWAS results; NMR metabolite identification using
  STOCSY.

- [NADfinder](https://bioconductor.org/packages/NADfinder) Call peaks
  for two samples: target and control. It will count the reads for
  tiles of the genome and then convert it to ratios. The ratios will
  be corrected and smoothed. The z-scores is calculated for each
  counting windows over the background. The peaks will be detected
  based on z-scores.

- [netReg](https://bioconductor.org/packages/netReg) netReg fits
  linear regression models using network-penalization. Graph prior
  knowledge, in the form of biological networks, is being
  incorporated into the likelihood of the linear model. The networks
  describe biological relationships such as co-regulation or
  dependency of the same transcription factors/metabolites/etc.
  yielding a part sparse and part smooth solution for coefficient
  profiles.

- [Organism.dplyr](https://bioconductor.org/packages/Organism.dplyr)
  This package provides an alternative interface to Bioconductor
  'annotation' resources, in particular the gene identifier mapping
  functionality of the 'org' packages (e.g., org.Hs.eg.db) and the
  genome coordinate functionality of the 'TxDb' packages (e.g.,
  TxDb.Hsapiens.UCSC.hg38.knownGene).

- [pathprint](https://bioconductor.org/packages/pathprint) Algorithms
  to convert a gene expression array provided as an expression table
  or a GEO reference to a 'pathway fingerprint', a vector of discrete
  ternary scores representing high (1), low(-1) or insignificant (0)
  expression in a suite of pathways.

- [pgca](https://bioconductor.org/packages/pgca) Protein Group Code
  Algorithm (PGCA) is a computationally inexpensive algorithm to
  merge protein summaries from multiple experimental quantitative
  proteomics data. The algorithm connects two or more groups with
  overlapping accession numbers. In some cases, pairwise groups are
  mutually exclusive but they may still be connected by another group
  (or set of groups) with overlapping accession numbers. Thus, groups
  created by PGCA from multiple experimental runs (i.e., global
  groups) are called "connected" groups. These identified global
  protein groups enable the analysis of quantitative data available
  for protein groups instead of unique protein identifiers.

- [phosphonormalizer](https://bioconductor.org/packages/phosphonormalizer)
  It uses the overlap between enriched and non-enriched datasets to
  compensate for the bias introduced in global phosphorylation after
  applying median normalization.

- [POST](https://bioconductor.org/packages/POST) Perform orthogonal
  projection of high dimensional data of a set, and statistical
  modeling of phenotye with projected vectors as predictor.

- [PPInfer](https://bioconductor.org/packages/PPInfer) Interactions
  between proteins occur in many, if not most, biological processes.
  Most proteins perform their functions in networks associated with
  other proteins and other biomolecules. This fact has motivated the
  development of a variety of experimental methods for the
  identification of protein interactions. This variety has in turn
  urshered in the development of numerous different computational
  approaches for modeling and predicting protein interactions.
  Sometimes an experiment is aimed at identifying proteins closely
  related to some interesting proteins. A network based statistical
  learning method is used to infer the putative functions of proteins
  from the known functions of its neighboring proteins on a PPI
  network. This package identifies such proteins often involved in
  the same or similar biological functions.

- [RaggedExperiment](https://bioconductor.org/packages/RaggedExperiment)
  This package provides a flexible representation of copy number,
  mutation, and other data that fit into the ragged array schema for
  genomic location data. The basic representation of such data
  provides a rectangular flat table interface to the user with range
  information in the rows and samples/specimen in the columns.

- [ramwas](https://bioconductor.org/packages/ramwas) RaMWAS provides
  a complete toolset for methylome-wide association studies (MWAS).
  It is specifically designed for data from enrichment based
  methylation assays, but can be applied to other data as well. The
  analysis pipeline includes seven steps: (1) scanning aligned reads
  from BAM files, (2) calculation of quality control measures, (3)
  creation of methylation score (coverage) matrix, (4) principal
  component analysis for capturing batch effects and detection of
  outliers, (5) association analysis with respect to phenotypes of
  interest while correcting for top PCs and known covariates, (6)
  annotation of significant findings, and (7) multi-marker analysis
  (methylation risk score) using elastic net. Additionally, RaMWAS
  include tools for joint analysis of methlyation and genotype data.

- [REMP](https://bioconductor.org/packages/REMP) Machine
  learing-based tools to predict DNA methylation of locus-specific
  repetitive elements (RE) by learning surrounding genetic and
  epigenetic information. These tools provide genomewide and
  single-base resolution of DNA methylation prediction on RE that are
  difficult to measure using array-based or sequencing-based
  platforms, which enables epigenome-wide association study (EWAS)
  and differentially methylated region (DMR) analysis on RE.

- [RITAN](https://bioconductor.org/packages/RITAN) Tools for
  comprehensive gene set enrichment and extraction of multi-resource
  high confidence subnetworks.

- [RIVER](https://bioconductor.org/packages/RIVER) An implementation
  of a probabilistic modeling framework that jointly analyzes
  personal genome and transcriptome data to estimate the probability
  that a variant has regulatory impact in that individual. It is
  based on a generative model that assumes that genomic annotations,
  such as the location of a variant with respect to regulatory
  elements, determine the prior probability that variant is a
  functional regulatory variant, which is an unobserved variable. The
  functional regulatory variant status then influences whether nearby
  genes are likely to display outlier levels of gene expression in
  that person. See the RIVER website for more information,
  documentation and examples.

- [RJMCMCNucleosomes](https://bioconductor.org/packages/RJMCMCNucleosomes)
  This package does nucleosome positioning using informative
  Multinomial-Dirichlet prior in a t-mixture with reversible jump
  estimation of nucleosome positions for genome-wide profiling.

- [RnaSeqGeneEdgeRQL](https://bioconductor.org/packages/RnaSeqGeneEdgeRQL)
  A workflow package for RNA-Seq experiments

- [rqt](https://bioconductor.org/packages/rqt) Despite the recent
  advances of modern GWAS methods, it still remains an important
  problem of addressing calculation an effect size and corresponding
  p-value for the whole gene rather than for single variant. The R-
  package rqt offers gene-level GWAS meta-analysis. For more
  information, see: "Gene-set association tests for next-generation
  sequencing data" by Lee et al (2016), Bioinformatics, 32(17),
  i611-i619, <doi:10.1093/bioinformatics/btw429>.

- [RTNduals](https://bioconductor.org/packages/RTNduals) RTNduals is
  a tool that searches for possible co-regulatory loops between
  regulon pairs generated by the RTN package. It compares the shared
  targets in order to infer 'dual regulons', a new concept that tests
  whether regulon pairs agree on the predicted downstream effects.

- [samExploreR](https://bioconductor.org/packages/samExploreR) This R
  package is designed for subsampling procedure to simulate
  sequencing experiments with reduced sequencing depth. This package
  can be used to anlayze data generated from all major sequencing
  platforms such as Illumina GA, HiSeq, MiSeq, Roche GS-FLX, ABI
  SOLiD and LifeTech Ion PGM Proton sequencers. It supports multiple
  operating systems incluidng Linux, Mac OS X, FreeBSD and Solaris.
  Was developed with usage of Rsubread.

- [sampleClassifier](https://bioconductor.org/packages/sampleClassifier)
  The package is designed to classify gene expression profiles.

- [scDD](https://bioconductor.org/packages/scDD) This package
  implements a method to analyze single-cell RNA- seq Data utilizing
  flexible Dirichlet Process mixture models. Genes with differential
  distributions of expression are classified into several interesting
  patterns of differences between two conditions. The package also
  includes functions for simulating data with these patterns from
  negative binomial distributions.

- [scone](https://bioconductor.org/packages/scone) SCONE is an R
  package for comparing and ranking the performance of different
  normalization schemes for single-cell RNA-seq and other
  high-throughput analyses.

- [semisup](https://bioconductor.org/packages/semisup) This R
  packages moves away from testing interaction terms, and move
  towards testing whether an individual SNP is involved in any
  interaction. This reduces the multiple testing burden to one test
  per SNP, and allows for interactions with unobserved factors.
  Analysing one SNP at a time, it splits the individuals into two
  groups, based on the number of minor alleles. If the quantitative
  trait differs in mean between the two groups, the SNP has a main
  effect. If the quantitative trait differs in distribution between
  some individuals in one group and all other individuals, it
  possibly has an interactive effect. Implicitly, the membership
  probabilities may suggest potential interacting variables.

- [sparseDOSSA](https://bioconductor.org/packages/sparseDOSSA) The
  package is to provide a model based Bayesian method to characterize
  and simulate microbiome data. sparseDOSSA's model captures the
  marginal distribution of each microbial feature as a truncated,
  zero-inflated log-normal distribution, with parameters distributed
  as a parent log-normal distribution. The model can be effectively
  fit to reference microbial datasets in order to parameterize their
  microbes and communities, or to simulate synthetic datasets of
  similar population structure. Most importantly, it allows users to
  include both known feature-feature and feature-metadata correlation
  structures and thus provides a gold standard to enable benchmarking
  of statistical methods for metagenomic data analysis.

- [splatter](https://bioconductor.org/packages/splatter) Splatter is
  a package for the simulation of single-cell RNA sequencing count
  data. It provides a simple interface for creating complex
  simulations that are reproducible and well-documented. Parameters
  can be estimated from real data and functions are provided for
  comparing real and simulated datasets.

- [STROMA4](https://bioconductor.org/packages/STROMA4) This package
  estimates four stromal properties identified in TNBC patients in
  each patient of a gene expression datasets. These stromal property
  assignments can be combined to subtype patients. These four stromal
  properties were identified in Triple negative breast cancer (TNBC)
  patients and represent the presence of different cells in the
  stroma: T-cells (T), B-cells (B), stromal infiltrating epithelial
  cells (E), and desmoplasia (D). Additionally this package can also
  be used to estimate generative properties for the Lehmann subtypes,
  an alternative TNBC subtyping scheme (PMID: 21633166).

- [swfdr](https://bioconductor.org/packages/swfdr) This package
  allows users to estimate the science-wise false discovery rate from
  Jager and Leek, "Empirical estimates suggest most published medical
  research is true," 2013, Biostatistics, using an EM approach due to
  the presence of rounding and censoring. It also allows users to
  estimate the proportion of true null hypotheses in the presence of
  covariates, using a regression framework, as per Boca and Leek, "A
  regression framework for the proportion of true null hypotheses,"
  2015, bioRxiv preprint.

- [TCGAbiolinksGUI](https://bioconductor.org/packages/TCGAbiolinksGUI)
  "TCGAbiolinksGUI: A Graphical User Interface to analyze cancer
  molecular and clinical data. A demo version of GUI is found in
  https://tcgabiolinksgui.shinyapps.io/tcgabiolinks/"

- [TCseq](https://bioconductor.org/packages/TCseq) Quantitative and
  differential analysis of epigenomic and transcriptomic time course
  sequencing data, clustering analysis and visualization of temporal
  patterns of time course data.

- [timescape](https://bioconductor.org/packages/timescape) TimeScape
  is an automated tool for navigating temporal clonal evolution data.
  The key attributes of this implementation involve the enumeration
  of clones, their evolutionary relationships and their shifting
  dynamics over time. TimeScape requires two inputs: (i) the clonal
  phylogeny and (ii) the clonal prevalences. Optionally, TimeScape
  accepts a data table of targeted mutations observed in each clone
  and their allele prevalences over time. The output is the TimeScape
  plot showing clonal prevalence vertically, time horizontally, and
  the plot height optionally encoding tumour volume during
  tumour-shrinking events. At each sampling time point (denoted by a
  faint white line), the height of each clone accurately reflects its
  proportionate prevalence. These prevalences form the anchors for
  bezier curves that visually represent the dynamic transitions
  between time points.

- [treeio](https://bioconductor.org/packages/treeio) Base classes and
  functions for parsing and exporting phylogenetic trees.

- [TSRchitect](https://bioconductor.org/packages/TSRchitect) In
  recent years, large-scale transcriptional sequence data has yielded
  considerable insights into the nature of gene expression and
  regulation in eukaryotes. Techniques that identify the 5' end of
  mRNAs, most notably CAGE, have mapped the promoter landscape across
  a number of model organisms. Due to the variability of TSS
  distributions and the transcriptional noise present in datasets,
  precisely identifying the active promoter(s) for genes from these
  datasets is not straightforward. TSRchitect allows the user to
  efficiently identify the putative promoter (the transcription start
  region, or TSR) from a variety of TSS profiling data types,
  including both single-end (e.g. CAGE) as well as paired-end
  (RAMPAGE, PEAT). Along with the coordiantes of identified TSRs,
  TSRchitect also calculates the width, abundance and Shape Index,
  and handles biological replicates for expression profiling.
  Finally, TSRchitect imports annotation files, allowing the user to
  associate identified promoters with genes and other genomic
  features. Three detailed examples of TSRchitect's utility are
  provided in the User's Guide, included with this package.

- [twoddpcr](https://bioconductor.org/packages/twoddpcr) The twoddpcr
  package takes Droplet Digital PCR (ddPCR) droplet amplitude data
  from Bio-Rad's QuantaSoft and can classify the droplets. A summary
  of the positive/negative droplet counts can be generated, which can
  then be used to estimate the number of molecules using the Poisson
  distribution. This is the first open source package that
  facilitates the automatic classification of general two channel
  ddPCR data. Previous work includes 'definetherain' (Jones et al.,
  2014) and 'ddpcRquant' (Trypsteen et al., 2015) which both handle
  one channel ddPCR experiments only. The 'ddpcr' package available
  on CRAN (Attali et al., 2016) supports automatic gating of a
  specific class of two channel ddPCR experiments only.

- [wiggleplotr](https://bioconductor.org/packages/wiggleplotr) Tools
  to visualise read coverage from sequencing experiments together
  with genomic annotations (genes, transcripts, peaks). Introns of
  long transcripts can be rescaled to a fixed length for better
  visualisation of exonic read coverage.

NEWS from new and existing packages
===================================

[ABAEnrichment](https://bioconductor.org/packages/ABAEnrichment)
-------------

Changes in version 1.5.10:

NEW FEATURES

- added hg20 as reference genome for gene coordinates (used for
  genomic-regions-input and gene_len=TRUE) Previously only hg19 was
  available, which is still the default

USER-LEVEL CHANGES

- hgnc-symbols for gene coordinates are replaced by gene symbols which
  are available for a higher number of genes (in most cases they are
  identical)

Changes in version 1.5.9:

NEW FEATURES

- New function get_annotated_genes returns genes annotated to enriched
  or user-defined brain regions

- New function get_id returns the ID of a brain region given its name

USER-LEVEL CHANGES

- results from aba_enrich get sorted on times_FWER_under_0.05 followed
  by min_FWER and mean_FWER (order of min and mean switched)

- genes in aba_enrich(...)&#91;&#91;2&#93;&#93; are sorted alphabetically

- when genomic regions are provided as input, candidate regions are
  implicitly also part of the background for the randomsets (like it is
  for single-genes-input)

Changes in version 1.5.8:

IMPROVEMENTS

- Use dynamic tolerance for FWER calculation

- Added more tests and checks

- Improved checking of arguments when genomic regions are provided as
  input

[AMOUNTAIN](https://bioconductor.org/packages/AMOUNTAIN)
---------

Changes in version 1.1.1:

NEW FEATURES

- C implementations

Changes in version 1.0.1:

USER-LEVEL CHANGES

- R markdown vignette

[anamiR](https://bioconductor.org/packages/anamiR)
------

Changes in version 1.2.0:

- GSEA workflow : GSEA analysis for finding specific pathaways related
  miRNAs is now available.

[AneuFinder](https://bioconductor.org/packages/AneuFinder)
----------

Changes in version 1.3.4:

NEW FEATURES

- Proper print() methods for AneuFinder objects.

- plotProfile(..., normalize.counts='2-somy') option added for plotting
  of normalized counts.

- heatmapGenomewideCluster() added for convenient assessment of the
  clusterByQuality() result.

BUG FIXES

- Fixed option Aneufinder(..., strandseq = TRUE) for DNAcopy method.

- Fixed a bug for SCE plotting in heatmapGenomewide().

- Proper creation of variable-width bins for huge reference files.

- Better heatmap dimensions for few cells.

- Bugfix for Inf values in clusterByQuality().

SIGNIFICANT USER-LEVEL CHANGES

- Renamed plotPCA() to plot_pca() to avoid namespace conflict with
  BiocGenerics.

Changes in version 1.3.3:

NEW FEATURES

- New function plotPCA() to do principal component analysis.

- Introduced parameter 'exclude.regions' to karyotypeMeasures(),
  plotHeterogeneity() and heatmapGenomewide(). This should facilitate
  excluding artifact regions from the clustering and karyotype
  measures.

- Parameter 'regions' now also available in plotHeterogeneity() (only
  in karyotypeMeasures() before).

SIGNIFICANT USER-LEVEL CHANGES

- The method to compute the dendrogram in heatmapGenomewide() was
  changed to simple hierarchical clustering on the copy number at
  bin-level (was segment-level before).

Changes in version 1.3.1:

NEW FEATURES

- Parameter 'normalChromosomeNumbers' in karyotypeMeasures() can handle
  mixture samples now.

[AnnotationFilter](https://bioconductor.org/packages/AnnotationFilter)
----------------

Changes in version 0.99.5:

NEW FEATURES

- Add convertFilterExpressionQuoted function.

- Add field method.

[AnnotationForge](https://bioconductor.org/packages/AnnotationForge)
---------------

Changes in version 1.18.0:

MODIFICATIONS

- RSQLite deprecated dbGetPreparedQuery/dbSendPreparedQuery; Updated
  with dbGetQuery/dbSendQuery/dbBind/dbFetch

- update for building BioC 3.5

BUG FIXES

- Resolve tmp issue, change outputDir to NCBIFilesDir

- Fixed a bug in makeOrgPackageFromNCBI when there are no GO terms

[AnnotationHub](https://bioconductor.org/packages/AnnotationHub)
-------------

Changes in version 2.8.0:

NEW FEATURES

- add .get1,RDSResource-method

- add RdsResource class

- add EnsDb dispatch class

- expose rdatapath in metadata

MODIFICATIONS

- modify records exposed as metadata - expose records added <= snapshot
  date - expose a single OrgDb per organism per BioC version

- edits to .get1,GenomicScores-method and
  .get1,GenomicScoresResource-method

- work on biocVersion and snapshotDate relationship: - snapshotDate()
  must be <= biocVersion() release date - possibleDates() are now
  filtered by snapshotDate()

- remove GenomicScoresResource; Robert Castelo will handle loading
  these resources in his GenomicScores software package

- Changed show method for hub object - removed sourcelastmodifieddate -
  added rdatadateadded

BUG FIXES

- fix bug in ordering of output from .uid0()

- fix bugs in 'snapshotDate<-' method

[AnnotationHubData](https://bioconductor.org/packages/AnnotationHubData)
-----------------

Changes in version 1.6.0:

NEW FEATURES

- add makeStandardTxDbsToSqlite() recipe

- add 'ensembl' and 'MySQL' as possible SourceType values

- tidy and export makeStandard*ToAHMs and makeNCBIToOrgDbsToAHMs

MODIFICATIONS

- move currentMetadata

- tidy pushResources interface

- modified parsing of species name and genome in
  .ensemblMetadataFromUrl()

- modified standard OrgDb recipe

- enhance and clean vignette

- move 'Tags' check from readCsvFromMetadata() to
  makeAnnotationHubMetadata()

- remove dependency on xml2, curl, httr and probably other wheel
  reinventions, alter imports and suggests

- specify multiple 'Tags' as colon separated string instead of comma
  separated; avoids problems with read.csv()

- select data moved to GenomeInfoDbData package

- Added additional documentation instructions for core members to add
  contributed data to AnnotationHub

- rename files; remove old JSON test file no longer applicable

- pass 'install' argument down through recipe

- General code tidy; remove unused functions and comments; clarify
  checks

BUG FIXES

- readMetadataFromCsv() fills in DataProvider and Coordinate_1_based if
  missing

- fix bug introduced in checking 'release' in makeEnsemblTwoBit recipe

- makeAnnotationHubMetadata() now processes all inst/extdata/*.csv
  files

- fix subset and import bug in makeAnnotationHubMetadata()

- Fix bug in Rdatapath and sourceurl for makeEnsemblFasta.R

[annotationTools](https://bioconductor.org/packages/annotationTools)
---------------

Changes in version 1.49.1 (2017-03-15):

- Easy mining of NCBI's new ortholog database ('Orthologs from
  Annotation pipeline', or 'gene_group' database) using the function
  getHOMOLOG.

[annotatr](https://bioconductor.org/packages/annotatr)
--------

Changes in version 1.2.0:

NEW FEATURES

- Add support for CpG annotations for hg38, mm10, and rn6 via the UCSC
  goldenpath URLs.

- Add a function to build annotations from AnnotationHub resources,
  build_ah_annots().

- Add support for chromHMM tracks (chromatin state) from the UCSC
  Genome Browser.

- Users may annotate to chromatin states in multiple cell lines, if
  desired.

- Use rtracklayer::liftOver to lift hg19 and mm9 enhancers into hg38
  and mm10.

USER-FACING CHANGES

- Add minoverlaps parameter to annotate_regions() that is passed to
  GenomicRanges::findOverlaps().

- Change supported_annotations() and supported_genomes() into
  builtin_annotations() and builtin_genomes(). This enables more
  flexibility required for AnnotationHub annotations.

- Added documentation for coercing result of annotate_regions() to
  data.frame and subsetting based on gene symbol to the vignette.

BUGFIXES

- Fixed a bug in coercion of GRanges to data.frame where row.names
  could be duplicated. Thanks to @kdkorthauer.

- Require GenomeInfoDb >= 1.10.3 because of changes to NCBI servers.

- Change scale_fill_brewer() to scale_fill_hue() in plot_categorical()
  to enable more categories and avoid plotting abnormalities.

- Fixed bug that mistakenly displayed some supported annotations.

- Fixed a bug in lncRNA annotation building caused by incomplete
  reference.

[aroma.light](https://bioconductor.org/packages/aroma.light)
-----------

Changes in version 3.5.1 (2017-04-14):

SIGNIFICANT CHANGES

- robustSmoothSpline() uses a re-weighted re-iterative algorithm that
  fits a smooth spline using stats::smooth.spline(), calculates the
  residuals and which are used to fit a re-weighted smooth spline and
  so on until converence. Due to updates to stats::smooth.spline() in R
  (>= 3.4.0) it is no longer feasible to maintain a highly optimized
  version of the algorithm, because it was based on internal
  stats::smooth.spline() code that has no completely changed.  Instead
  the re-iterative algorithm calls stats::smooth.spline() as is, which
  slows it down.  More importantly, it will now give slightly different
  estimates.

SOFTWARE QUALITY

- In addition to continous integration (CI) tests and nightly
  Bioconductor tests, the package is now also tested regularly against
  all reverse package depencies available on CRAN and Bioconductor.

Changes in version 3.5.0 (2016-10-18):

- The version number was bumped for the Bioconductor devel version,
  which is now BioC v3.5 for R (>= 3.4.0).

[BaalChIP](https://bioconductor.org/packages/BaalChIP)
--------

Changes in version 1.1.1:

NEW FEATURES

- Added Citation file

- Added citation reference to documentation

- Corrected typo in documentation

[basecallQC](https://bioconductor.org/packages/basecallQC)
----------

Changes in version 0.99:

- Initial Submission.

[BgeeDB](https://bioconductor.org/packages/BgeeDB)
------

Changes in version 2.2.0:

- Instead of returning a topGO-compatible object, topAnat.R now returns
  an object from the topAnatData class, an extension of topGOdata class
  from the topGO package.

- Fixed small issue with management of data types given as input by the
  user (dataType argument when creating new Bgee class)

- Fixed bug in experiment Id check step. Now accomodates SRA Ids.

- Fixed data frames header names that included double dots.

- Removed dependency to biomaRt in vignette. Code is still detailed but
  not run, instead a pre-created gene list object is loaded from the
  data/ directory.

[bigmelon](https://bioconductor.org/packages/bigmelon)
--------

Version: 1.1.14
Text:

[bioCancer](https://bioconductor.org/packages/bioCancer)
---------

Version: 1.2.21
Category: Last commit was not completed
Text:

[BiocFileCache](https://bioconductor.org/packages/BiocFileCache)
-------------

Changes in version 0.99.0:

SIGNIFICANT NEW FEATURES

- First Bioconductor release.

[BiocStyle](https://bioconductor.org/packages/BiocStyle)
---------

Changes in version 2.4.0:

NEW FEATURES

- Vignette "Authoring R Markdown vignettes"

- R Markdown templates for 'pdf_document2' and 'html_document2'

- Standard way of specifying author affiliations

- Support for short title in R Markdown PDF output

- Argument 'relative.path' to 'latex2()'
  (https://support.bioconductor.org/p/90352/)

SIGNIFICANT USER-VISIBLE CHANGES

- Increase column width in order to accommodate 80 characters wide code
  chunks

- Separate caption title from description with newline

- Use canonical URL to link to CRAN packages
  (https://github.com/Bioconductor/BiocStyle/issues/24)

- Consistently number equations on right hand side across different
  output formats

- Numerous CSS tweaks

BUG FIXES AND IMPROVEMENTS

- Support for PDFs typeset with 9pt and 8pt font size

- Proper formatting of 'longtable' captions

- Fix to retain spaces in '\texttt'

- Replace carets "\^{}" by "\textasciicircum" to fix incompatibility
  with LaTeX 'soul' package used for inline code highlighting

- Patch to avoid overfull pages containing a float followed by a
  longtable

[biomaRt](https://bioconductor.org/packages/biomaRt)
-------

Changes in version 2.32.0:

BUG FIXES

- Fixed bug when columns were not returned in the order requested,
  which resulted in the wrong column names being added to the result.

[biosigner](https://bioconductor.org/packages/biosigner)
---------

Changes in version 1.3.6:

INTERNAL MODIFICATIONS

- minor internal modification for compatibility with the Galaxy module

Changes in version 1.3.4:

MINOR MODIFICATIONS

- minor modification in vignette

Changes in version 1.3.2:

MINOR MODIFICATIONS

- vignette now in pdf format

[BLMA](https://bioconductor.org/packages/BLMA)
----

Changes in version 1.0.0:

Initial release of package BLMA includes

- bilevelAnalysisGeneset: a function to perform a bi-level
  meta-analysis in conjunction with geneset enrichment methods
  (ORA/GSA/PADOG) to integrate multiple gene expression datasets.

- bilevelAnalysisPathway: a function to perform a bi-level
  meta-analysis conjunction with Impact Analysis to integrate multiple
  gene expression datasets.

- intraAnalysisClassic: a function to perform an intra-experiment
  analysis in conjunction with any of the classical hypothesis testing
  methods, such as t-test, Wilcoxon test, etc.

- bilevelAnalysisClassic: a function to perform a bi-level
  meta-analysis in conjunction with any of the classical hypothesis
  testing methods, such as t-test, Wilcoxon test, etc.

- intraAnalysisGene: a function to perform an intra-experiment analysis
  in conjunction with the moderated t-test (limma package) for the
  purpose of differential expression analysis of a gene expression
  dataset

- bilevelAnalysisGene: a function to perform a bi-level meta-analysis
  in conjunction with the moderate t-test (limma package) for the
  purpose of differential expression analysis of multiple gene
  expression datasets

- loadKEGGPathways: this function loads KEGG pathways and names

- addCLT: a function to combine independent studies using the average
  of p-values

- fisherMethod: a function to combine independent p-values using the
  minus log product

- stoufferMethod: a function to combine independent studies using the
  sum of p-values transformed into standard normal variables

[BrowserViz](https://bioconductor.org/packages/BrowserViz)
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

[bsseq](https://bioconductor.org/packages/bsseq)
-----

Changes in version 1.11:

- bsseq now uses DelayedMatrix objects from the DelayedArray package
  for all matrix-like data. This enables large data to be stored on
  disk rather than in memory.

- Serialized (saved) BSseq, BSseqTstat, and BSseqStat objects will need
  to be updated by invoking x <- updateObject(x).

[BUMHMM](https://bioconductor.org/packages/BUMHMM)
------

Changes in version 0.99.6:

- Updated version information.

Changes in version 0.99.5:

- Updated citations.

Changes in version 0.99.3:

- Addressed Notes and Warnings after build.

Changes in version 0.99.2:

- Implemented technical changes for better integration with existing
  Bioconductor classes.

Changes in version 0.99.1:

- Fixed minor build errors, e.g. subscribing to mailing list.

Changes in version 0.99.0:

- Submitted to Bioconductor.

[CAGEr](https://bioconductor.org/packages/CAGEr)
-----

Changes in version 1.18.0:

- Added Charles Plessy as co-maintainer.

- Remove warning by replacing deprecated ignoreSelf() with drop.self().

[Cardinal](https://bioconductor.org/packages/Cardinal)
--------

Changes in version 1.7.2:

SIGNIFICANT USER-VISIBLE CHANGES

- In 'image' method, 'superpose = TRUE' now supports multiple LHS
  arguments in the formula (e.g., formula = a + b ~ x * y) >>>>>>>
  e51eb75... v1.7.1 updates for matter compatibility

Changes in version 1.7.1 (2016-11-29):

NEW FEATURES

- PCA is now supported for larger-than-memory on-disk datasets

- External 'matter' matrices replace 'Binmat' matrices for on-disk
  support

- Added 'image3D' aliases for all 'ResultSet' subclasses

SIGNIFICANT USER-VISIBLE CHANGES

- Added 'matter' package to Depends list <<<<<<< HEAD >>>>>>>
  41b2923... v.1.7.2 bug fixes and cleanup >>>>>>> e51eb75... v1.7.1
  updates for matter compatibility

[CFAssay](https://bioconductor.org/packages/CFAssay)
-------

Version: 1.9.1
Category: Bugfixes in functions pes, sfpmean and plot.cellsurvLQfit:
        calculation of mean values at different doses for curve and
        point plotting when PEmethod = "fix
Text:

[ChAMP](https://bioconductor.org/packages/ChAMP)
-----

Changes in version 2.6.2:

- DMRcate pacakge get updated, Error like "Error in if (nsig == 0) { :
  missing value where TRUE/FALSE needed" has been solved.

- In champ.load(), instead of replacing all 0 and negative value into
  0.0001, we relplace them as smallest positive now.

- Fixed warnings() in GUI() functions.

- In champ.runCombat() function, removed restriction on factors like
  Sample_Group. Also, added "variable" parameter so that user may
  assign other variables other then "Sample_Group".

- Modified champ.DMR() function, for ProbeLasso, there is no need to
  input myDMP anymore, ProbeLasso function would calculate inside the
  function.

[chipenrich](https://bioconductor.org/packages/chipenrich)
----------

Changes in version 2.0.0:

NEW FEATURES

- A new method for enrichment, polyenrich() is designed for gene set
  enrichment of experiments where the presence of multiple peaks in a
  gene is accounted for in the model. Use the polyenrich() function for
  this method.

- New features resulting from chipenrich.data 2.0.0:

- New genomes in chipenrich.data: danRer10, dm6, hg38, rn5, and rn6.

- Reactome for fly in chipenrich.data.

- Added locus definitions, GO gene sets, and Reactome gene sets for
  zebrafish.

- All genomes have the following locus definitions: nearest_tss,
  nearest_gene, exon, intron, 1kb, 5kb, 10kb, 1kb_outside_upstream,
  5kb_outside_upstream, 10kb_outside_upstream, 1kb_outside,
  5kb_outside, and 10kb_outside.

IMPROVEMENTS

- The chipenrich method is now significantly faster. Chris Lee figured
  out that spline calculations in chipenrich are not required for each
  gene set. Now a spline is calculated as peak ~ s(log10_length) and
  used for all gene sets. The correlation between the resulting
  p-values is nearly always 1. Unfortunately, this approach cannot be
  used for broadenrich().

- The chipenrich(..., method='chipenrich', ...) function automatically
  uses this faster method.

- Clarified documentation for the supported_locusdefs() to give
  explanations for what each locus definition is.

- Use sys.call() to report options used in chipenrich() in opts.tab
  output. We previously used as.list(environment()) which would also
  output entire data.frames if peaks were loaded in as a data.frame.

- Various updates to the vignette to reflect new features.

SIGNIFICANT USER-LEVEL CHANGES

- As a result of updates to chipenrich.data, ENRICHMENT RESULTS MAY
  DIFFER between chipenrich 1.Y.Z and chipenrich 2.Y.Z. This is because
  revised versions of all genomes have been used to update
  LocusDefinitions, and GO and Reactome gene sets have been updated to
  more recent versions.

- The broadenrich method is now its own function, broadenrich(),
  instead of chipenrich(..., method = 'broadenrich', ...).

- User interface for mappability has been streamlined. 'mappability'
  parameter in broadenrich(), chipenrich(), and polyenrich() functions
  replaces the three parameters previously used: 'use_mappability',
  'mappa_file', and 'read_length'. The unified 'mappability' parameter
  can be 'NULL', a file path, or a string indicating the read length
  for mappability, e.g. '24'.

- A formerly hidden API for randomizations to assess Type I Error rates
  for data sets is now exposed to the user. Each of the enrich
  functions has a 'randomization' parameter. See documentation and
  vignette for details.

- Many functions with the 'genome' parameter had a default of 'hg19',
  which was not ideal. Now users must specify a genome and it is
  checked against supported_genomes().

- Input files are read according to their file extension. Supported
  extensions are bed, gff3, wig, bedGraph, narrowPeak, and broadPeak.
  Arbitrary extensions are also supported, but there can be no header,
  and the first three columns must be chr, start, and end.

SIGNIFICANT BACKEND CHANGES

- Harmonize all code touching LocusDefinition and tss objects to
  reflect changes in chipenrich.data 2.0.0.

- Alter setup_ldef() function to add symbol column. If a valid genome
  is used use orgDb to get eg2symbol mappings and fill in for the user.
  Users can give their own symbol column which will override using
  orgDb. Finally, if neither symbol column or valid genome is used,
  symbols are set to NA.

- Any instance of 'geneid' or 'names' to refer to Entrez Gene IDs are
  now 'gene_id' for consistency.

- Refactor read_bed() function as a wrapper for rtracklayer::import().

- Automatic extension handling of BED3-6, gff3, wig, or bedGraph.

- With some additional code, automatic extension handling of narrowPeak
  and broadPeak.

- Backwards compatible with arbitrary extensions: this still assumes
  that the first three columns are chr, start, end.

- The purpose of this refactor is to enable additional covariates for
  the peaks for possible use in future methods.

- Refactor load_peaks() to use
  GenomicRanges::makeGRangesFromDataFrame().

- Filtering gene sets is now based on the locus definition, and can be
  done from below (min) or above (max). Defaults are 15 and 2000,
  respectively.

- Randomizations are all done on the LocusDefinition object.

- Added lots of unit tests to increase test coverage.

- Make Travis builds use sartorlab/chipenrich.data version of data
  package for faster testing.

DEPRECATED AND DEFUNCT

- Calling the broadenrich method with chipenrich(..., method =
  'broadenrich', ...) is no longer valid. Instead, use broadenrich().

- Various utility functions that were used in the original development
  have been removed. Users never saw or used them.

BUG FIXES

- Fixed bug in randomization with length bins where artifactually,
  randomizations would sort genes on Entrez ID introducing problems in
  Type I error rate.

- Fixed a bug where the dependent variable used in the enrichment model
  was used to name the rows of the enrichment results. This could be
  confusing for users. Now, rownames are simply integers.

- Fixed a bug that expected the result of read_bed() to be a list of
  IRanges from initial development. Big speed bump.

Changes in version 1.12.1:

BUG FIXES

- Fixed a bug in the check for proper organism + geneset combinations.
  Prevented combinations that are actually valid from running.

[ChIPpeakAnno](https://bioconductor.org/packages/ChIPpeakAnno)
------------

Changes in version 3.9.19:

- add FAQs

Changes in version 3.9.18:

- import rowRanges() generic from DelayedArray instead of
  SummarizedExperiment

Changes in version 3.9.17:

- Improve the function featureAlignedExtendSignal

Changes in version 3.9.16:

- Fix the bug that unable to find an inherited method for function
  'reverseComplement' for signature '"BString"' in oligoSummary
  function.

Changes in version 3.9.14:

- re-order the signals for cumulativePercentage.

Changes in version 3.9.13:

- add new function cumulativePercentage.

Changes in version 3.9.12:

- remove unused code in oligoSummary.

Changes in version 3.9.11:

- update the documents of featureAlignedSignal.

Changes in version 3.9.10:

- fix the bug of featureAlignedSignal when there are seqnames in
  featues but not in cvglists.

Changes in version 3.9.9:

- add 2 more parameters to getEnrichedGO

Changes in version 3.9.8:

- update the colnames of toGRanges

Changes in version 3.9.7:

- fix the bug of findOverlapsOfPeaks when the peaks are exactly same

Changes in version 3.9.6:

- fix the bug of findOverlapsOfPeaks when there is no peak names for
  peaklist

Changes in version 3.9.5:

- update the documents of findOverlapsOfPeaks

Changes in version 3.9.4:

- update the documents of getEnrichedGO

- add more output for findOverlapsOfPeaks.

Changes in version 3.9.3:

- trim seqnames when using toGRanges

Changes in version 3.9.2:

- add adjustFragmentLength for featureAlignedExtendSignal.

Changes in version 3.9.1:

- fix the bug for pair end reads for featureAlignedExtendSignal.

[ChIPseeker](https://bioconductor.org/packages/ChIPseeker)
----------

Changes in version 1.11.4:

- bug fixed of intron rank <2017-04-19, Wed> +
  https://github.com/GuangchuangYu/ChIPseeker/issues/54

Changes in version 1.11.3:

- bug fixed of dropAnno <2017-04-10, Mon>

- bug fixed of peak width generated by shuffle <2017-03-31, Fri> +
  <https://github.com/GuangchuangYu/ChIPseeker/issues/51>

Changes in version 1.11.2:

- optimize getGeneAnno <2016-12-21, Wed>

- change plotAnnoBar and plotDistToTSS according to stacking bar order
  change in ggplot2 (v2.2.0) <2016-12-16, Fri> +
  https://github.com/GuangchuangYu/ChIPseeker/issues/47 +
  https://blog.rstudio.org/2016/11/14/ggplot2-2-2-0/

Changes in version 1.11.1:

- update startup message <2016-11-09, Wed>

[chromstaR](https://bioconductor.org/packages/chromstaR)
---------

Changes in version 1.1.4:

BUG FIXES

- Fixed a mistake in the calculation of differential scores from
  version 1.1.2

Changes in version 1.1.2:

NEW FEATURES

- Proper print() methods for all objects.

BUG FIXES

- Fixed a bug where chromosomes with a single bin were making problems.

Changes in version 1.1.1:

NEW FEATURES

- Selection of peaks can be done with 'changeFDR'.

- Peak calls are available in each chromstaR-object as list entry
  '$peaks'.

SIGNIFICANT USER-LEVEL CHANGES

- 'plotFoldEnrichment' renamed to 'plotEnrichment'.

- 'exportBinnedData', 'exportUnivariates', 'exportMultivariates',
  'exportCombinedMultivariates' replaced by 'exportCounts',
  'exportPeaks', 'exportCombinations'.

BUG FIXES

- Proper computation of fold enrichments, with < 1 indicating depletion
  and > 1 indicating enrichment.

DEPRECATED AND DEFUNCT

- 'changePostCutoff'.

- 'plotFoldEnrichment'.

- 'exportBinnedData', 'exportUnivariates', 'exportMultivariates',
  'exportCombinedMultivariates'.

[CHRONOS](https://bioconductor.org/packages/CHRONOS)
-------

Version: 1.3.1
Text:

Version: 1.3.0
Text:

[ClassifyR](https://bioconductor.org/packages/ClassifyR)
---------

Changes in version 1.10.0:

- errorMap replaced by samplesMetricMap. The plot can now show either
  error rate or accuracy.

[cleaver](https://bioconductor.org/packages/cleaver)
-------

Changes in version 1.13.2 (2017-02-14):

- Use `expect_equal` instead of `expect_identical` to avoid failing of
  the `.revsequence` unit test on "toluca2" (Mac OS X Mavericks
  (10.9.5) / x86_64).

- Add "importFrom(stats, setNames)" to NAMESPACE.

Changes in version 1.13.1 (2017-01-04):

- Remove GPL headers from *.R files.

- Remove useless return calls.

[clusterExperiment](https://bioconductor.org/packages/clusterExperiment)
-----------------

Changes in version 1.1.2 (2017-04-04):

Changes

- RSEC now has option `rerunClusterMany`, which if FALSE will not rerun
  the clusterMany step if RSEC is called on an existing
  clusterExperiment object (assuming of course, clusterMany has been
  run already on the object)

- setBreaks now has option `makeSymmetric` to force symmetric breaks
  around zero when using the quantile option.

- setBreaks now has a default for breaks (i.e. for minimal use, the
  user doesn't have to give the argument, just the data) in which case
  setBreaks will automatically find equal-spaced breaks of length 52
  filling the range of data compatible with aheatmap. The order of the
  arguments `data` and `breaks` has been switched, however, to better
  accomodate this usage.

- plotClusters can now handle NA values in the colData

- plotClusters for `clusterExperiment` object now allows for setting
  `sampleData=TRUE` to indicate the plotting all of the sampleData in
  the colData slot.

- nPCADims now allows values between 0,1 to allow for keeping
  *proportion* of variance explained.

- addClusters now allows for argument `clusterLabel` to assign a
  clusterLabel when the added cluster is a vector (if matrix, then
  clusterLabel is just the column names of the matrix of cluster
  assignments)

Bug fixes

- fixed bug in clusterExperiment subsetting to deal with orderSamples
  correctly.

- fixed bug in mergeClusters unable to plot when too big of edge
  lengths (same as plotDendrogram)

- fixed bug in subsetting, where unable to subset samples by character

- fixed bug in removeClusters so that correctly updates dendro_index
  and primary_index slots after cluster removed.

Changes in version 1.1.1 (2016-10-14):

Changes

- Inverted definition of contrast for one-versus-all so now is
  X-ave(all but X); this means logFC positive -> cluster X greater than
  average of rest

Bug fixes

- add check in clusterMany that non-zero dimensions

- changed 'warning' to 'note' in combineMany when no clusters
  specified.

- fixed bug in plotDendrogram unable to plot when makeDendrogram used
  dimReduce="mad"

- fixed bug in clusterMany where beta set to NA if clusteringFunction
  is type K. Added check that beta cannot be NA.

- Added check that alpha and beta in (0,1)

[clusterProfiler](https://bioconductor.org/packages/clusterProfiler)
---------------

Changes in version 3.3.6:

- update kegg_species information <2017-03-26, Sun>

- bug fixed of bitr_kegg for converting ID to kegg <2017-03-01, Wed> +
  https://support.bioconductor.org/p/93170/#93174

- better support of plotGOgraph for gseGO, use core enriched gene info
  <2017-02-27, Mon>

Changes in version 3.3.5:

- solve #81 <2017-02-12, Sun>

- bitr_kegg support converting Path/Module to geneID and vice versa
  <2017-01-03, Tue>

Changes in version 3.3.4:

- bug fixed of download_KEGG for supporting different keyType
  <2017-01-02, Mon>

- split 3 GO sub-ontologies for enrichGO <2016-12-12, Mon>

- dotplot for compareClusterResult to support 3 GO sub-ontologies
  <2016-12-8, Thu>

- bug fixed of determine `ont` from `expand.dots` <2016-12-06, Mon> +
  see https://github.com/GuangchuangYu/clusterProfiler/issues/72

Changes in version 3.3.3:

- switch from BiocStyle to prettydoc <2016-11-30, Wed>

- change summary to as.data.frame in internal calls to prevent warning
  message <2016-11-14, Mon>

Changes in version 3.3.2:

- define simplify generics as it was removed from IRanges <2016-11-11,
  Fri>

Changes in version 3.3.1:

- update startup message <2016-11-09, Wed>

- bug fixed in enrichDAVID <2016-11-06, Sun> +
  https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/#comment-2987574851

- bug fixed in head, tail and dim for compareCluster object
  <2016-10-18, Tue>

[ClusterSignificance](https://bioconductor.org/packages/ClusterSignificance)
-------------------

Changes in version 1.3.2:

- added p-value estimate confidence intervals via conf.int function

- Some vignette updates including adding the "common questions"
  section.

Changes in version 1.3.1:

- Added plotting of individual comparisons in classification and
  permutation plots.

- Added df argument (degrees of freedom, passed to smooth.spline) to
  projection and permute function. This allows some degree of control
  over how linear or crooked the principal curve is drawn. NOTE: you
  must (!) give the same df value to the projection and permute
  functions for your results to be valid and, at the moment, there is
  no automated checking that this is the case.

- The classification step now uses a mew method for separation scoring
  and fixes a previous bug which could occur when group sizes were not
  equal. This change is reflected in the vignette.

- Fixed bug in concatenation c(), when some iterations fail for a
  specific comparison thus causing a error when concatenating permute
  results.

[CNEr](https://bioconductor.org/packages/CNEr)
----

Changes in version 3.5:

NEW FEATURES

- Add function orgKEGGIds2EntrezIDs to fetch the mapping between KEGG
  IDs and Entrez IDs

- Add function makeAxtTracks

- Add function addAncestorGO

Changes in version 3.4:

NEW FEATURES

- Updated CNE class for storing all the information about running the
  pipeline.

- Add read.rmMask.GRanges to read RepeatMasker .out file.

- Add read.rmskFasta to read soft repeat masked fasta file.

- Add the distribution plot of axt alignment matches.

- Add the distribution plot of CNE length.

- Add the syntenicDotplot for axt alignment and GRangePairs.

- When readAxt and readBed, seqinfo is kept when available.

- New fixCoordinates function makes the coordinates of Axt alignments
  always relative to positive strands.

- Add parallel subAxt for Axt alignment.

- Add the plot of genomic distribution of CNE.

- Add the function to make bed and bigwig files of CNEs.

BUG FIXES

- Instead of an error, an empty GRangePairs is returned when no CNEs
  identified.

Changes in version 3.3:

BUG FIXES

- Fix a bug caused by "format" in the blat step.

NEW FEATURES

- Add the pairwise whole genome alignment pipeline

- Add a new class "GRangePairs"

- "Axt" class is now based on "GRangePairs" class.

- readAncora for reading Ancora format CNE files.

[ComplexHeatmap](https://bioconductor.org/packages/ComplexHeatmap)
--------------

Changes in version 1.13.2:

- add `packLegend()`

- Legend() supports to add txt labels on grids

Changes in version 1.13.1:

- `Heatmap()`: add `km_title` to set the format of row title when `km`
  is set

- `anno_link()`: add `extend` to extend the regions for the labels

- `anno_boxplot()`: for row annotation, outliers are in the correct in
  y-axis. Thanks @gtg602c for the fix

- `HeatmapAnnotation()`: gaps are included in the size of the
  annotations

- `anno_link()`: graphic parameters are correctly reordered

- `densityHeatmap()`: viewport is created with `clip = TRUE`

- `decorate_*()`: add `envir` option to control where to look for
  variables inside `code`

- `Legend()`: title supports expression

- `anno_*()`: if the input is a data frame, warn that users may convert
  it to matrix

[CrispRVariants](https://bioconductor.org/packages/CrispRVariants)
--------------

Changes in version 1.3.7:

- Updates to vignette

- Fix bug removing variants by name in variantCounts

- Fixed argument legend.symbol.size being ignored in
  plotAlignmenta,DNAString-method.

Changes in version 1.3.6:

- Fix in new function mergeChimeras when no chimeras present

Changes in version 1.3.5:

- New option "minoverlap" in readsToTarget allows reads that do not
  span the target region to be considered

- plotAlignments now works with character as well as DNAString objects

- Merging of long gaps mapped as chimeras now possible

Changes in version 1.3.4:

- New function refFromAlns infers the reference sequence from aligned
  reads

- Fixed bug causing an empty plot when plotting a single alignment with
  a large deletion

- Changed annotateGenePlot from panel.margin to panel.spacing in
  accordance with recent ggplot2 versions

- Added "create.plot" argument to plotAlignments for signature
  CrisprSet to make plot customisation easier.

- Fixed bug in argument names when all alignments are chimeric

- CrisprRun name now defaults to the coordinates when no name is
  provided

Changes in version 1.3.3:

- Fixed bug causing incorrect x-axis position in plotAlignments when
  strand unspecified

[csaw](https://bioconductor.org/packages/csaw)
----

Changes in version 1.10.0:

- Added calculation of dominant directionality in combineTests().
  Fixed out-of-array indexing bug in the C++ code.

- Supported factor input for ids argument in combineTests(),
  getBestTest().

- Added the empiricalFDR(), empiricalOverlaps() functions for
  controlling the empirical FDR.

- Added the mixedClusters(), mixedOverlaps() functions for testing for
  mixed clusters.

- Ensured that window-level FDR threshold chosen by controlClusterFDR()
  is not above the cluster-level FDR.

- Minor fix to scaledAverage() to avoid slightly inaccurate results.
  Also, zero or negative scale factors now return -Inf and NA,
  respectively.

- Switched to new scaleOffset() function for adding offsets in
  asDGEList(). Added option to specify the assay to be used.

- Added multi-TSS support in detailRanges().

- Modified paired-end machinery in windowCounts(), getPESizes() to be
  more accommodating of overruns.

- Ignored secondary and supplementary alignments in all functions.

- Added options to specify assay in SE objects in filterWindows().

- Replaced weighting with normalization options in profileSites().

- Updated user's guide.

[customProDB](https://bioconductor.org/packages/customProDB)
-----------

Changes in version 1.15.2:

UPDATED FUNCTIONS

- Update function SharedJunc.R

Changes in version 1.14.1:

UPDATED FUNCTIONS

- Update mouse dbsnp version part in function
  PrepareAnnotationEnsembl.R

- Update function OutputVarproseq.R, JunctionType.R and calculateRPKM.R

NEW FEATURES

- Update function OutputNovelJun.R to save '_jun_anno.RData',
  '_coding.RData' and '_junpep.RData', which could be used as input for
  proBAMr package

[cydar](https://bioconductor.org/packages/cydar)
-----

Changes in version 1.0.0:

- New package cydar, for detecting differential abundance in mass
  cytometry data.

[cytofkit](https://bioconductor.org/packages/cytofkit)
--------

Changes in version 1.6.6 (2017-04-11):

NEW FEATURES

- add fixedLogicle transformation option on the GUI, with window popup
  to allow specifing the w, t, m, a parameters for logicle
  transformation.

- add openShinyAPP (boolean parameter) option in cytofkit main
  function, which can open shinyAPP once the analysis was done and
  automatically load the result object into the shinyAPP for
  exploration.

- add cytofkitShinyAPP2 function which can take cytofkit
  analysis_results (either file name of R object) as input and
  automatically load to shinyAPP once launched.

Changes in version 1.6.5 (2017-03-27):

MODIFICATION

- debug that FlowSOM_k doesn't work in cytofkit main function

Changes in version 1.6.4 (2017-03-17):

MODIFICATION

- in the shiny server code change function call of c in do.call to
  base::c

Changes in version 1.6.3 (2017-03-16):

MODIFICATION

- debug the FSC|SSC channel processing error

Changes in version 1.6.2 (2017-03-08):

NEW FEATURES

- add default linear transformation to FSC and SSC channels

- add support for PDF figure download on shinyAPP, update the side
  panel to be tab dependent

- add new color palatte in heatmap (greenred and spectral) and level
  plot (spectral)

- Add cluster filter in cluster plot on shinyAPP

- Allow multiple annotation for same cluster (specify
  cluster_annotation name) on shinyAPP

- Allow color selection for each cluster on shinyAPP

- Allow modification of the marker name on shinyAPP

Changes in version 1.6.1 (2016-10-27):

MODIFICATION

- updated colorPalette options in cytof_colorPlot function, also the
  Shiny APP, added spectral.

- debugged rowname conflication when regroup the samples in shinyAPP,
  now only use global ID, discarded the local cell ID, which avoid the
  dumplicate rownames conflication but results in failure in saving new
  FCS files.

[CytoML](https://bioconductor.org/packages/CytoML)
------

Version: 1.1.2
Category: support diva workspace parsing
Text:

[DaMiRseq](https://bioconductor.org/packages/DaMiRseq)
--------

Changes in version 0.99.0:

- DaMiRseq package has been released!

[deepSNV](https://bioconductor.org/packages/deepSNV)
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

[DEGreport](https://bioconductor.org/packages/DEGreport)
---------

Version: 1.11.7
Text: 2017-04-08 Lorena Pantano <lorena.pantano@gmail.com> Fix: Add new
        contributor

Version: 1.11.6
Text: 2017-04-07 Lorena Pantano <lorena.pantano@gmail.com> Feature: Add
        function to plot genes in a wide format

Version: 1.11.4
Text: 02-17-2017 Lorena Pantano <lorena.pantano@gmail.com> Feature:
        Re-organize vignette. Feature: Ignore warnings when plotting
        Feature: Improve volcano plot

Version: 1.11.3
Text: 01-03-2017 Lorena Pantano <lorena.pantano@gmail.com> Fix: fix
        order of clusters figures that are not in the correct place in
        some cases with many groups.

Version: 1.11.2
Text: 12-09-2016 Lorena Pantano <lorena.pantano@gmail.com> Features:
        Add degCheckFactors functions to plot sizefactors used to
        normalize count data.

Version: 1.11.1
Text: 10-18-2016 Lorena Pantano <lorena.pantano@gmail.com> Fixes: print
        clusterProfiler output.

[derfinder](https://bioconductor.org/packages/derfinder)
---------

Changes in version 1.9.6:

BUG FIXES

- Fixed define_cluster() to match recent changes in BiocParallel and
  fixed an if clause in regionMatrix() that could lead to warnings in
  some situations.

Changes in version 1.9.3:

SIGNIFICANT USER-VISIBLE CHANGES

- regionMatrix() now has explicit arguments 'totalMapped' and
  'targetSize' so that users will almost always normalize by library
  size when using this function (if they see the help page) or in the
  steps prior to using regionMatrix().

Changes in version 1.9.2:

BUG FIXES

- Clarified the documentation of 'mc.cores' and 'mc.cores.load' in
  fullCoverage() thanks to feedback from Emily Burke
  https://github.com/emilyburke.

[DESeq2](https://bioconductor.org/packages/DESeq2)
------

Changes in version 1.16.0:

- DESeq() and nbinomWaldTest() the default setting will be
  betaPrior=FALSE, and the recommended pipeline will be to use
  lfcShrink() for producing shrunken LFC.

- Added a new function unmix(), for unmixing samples according to
  linear combination of pure components, e.g. "tissue deconvolution".

- Added a new size factor estimator, "poscounts", which evolved out of
  use cases in Paul McMurdie's phyloseq package.

- Ability to specify observation-specific weights, using
  assays(dds)&#91;&#91;"weights"&#93;&#93;. These weights are picked up by dispersion
  and NB GLM fitting functions.

Changes in version 1.15.40:

- Adding a new function unmix(), for unmixing samples according to pure
  components, e.g. "tissue deconvolution". The pure components are
  added on the gene expression scale (either normalized counts or
  TPMs), and the loss is calculated in a variance stabilized space.

Changes in version 1.15.39:

- Added a new size factor estimator, "poscounts", which evolved out of
  use cases in Paul McMurdie's phyloseq package.

Changes in version 1.15.36:

- Ability to specify observation-specific weights, using
  assays(dds)&#91;&#91;"weights"&#93;&#93;. These weights are picked up by dispersion
  and NB GLM fitting functions.

Changes in version 1.15.28:

- Remove some code that would "zero out" LFCs when both groups involved
  in a contrast had zero counts. This lead to inconsistency when
  similarly contrasts were performed by refactoring.

Changes in version 1.15.12:

- DESeq() and nbinomWaldTest() the default setting will be
  betaPrior=FALSE, and the recommended pipeline will be to use
  lfcShrink() for producing shrunken log2 fold changes for
  visualization and ranking. Explanation for the change is presented in
  the vignette section: "Methods changes since the 2014 DESeq2 paper"

Changes in version 1.15.9:

- Adding prototype function lfcShrink().

- Vignette conversion to Rmarkdown / HTML.

Changes in version 1.15.3:

- Removing betaPrior option for nbinomLRT, in an effort to clean up and
  reduce old un-used functionality.

[DiffBind](https://bioconductor.org/packages/DiffBind)
--------

Changes in version 2.4.0:

- Feature: add new plot - dba.plotVolcano

[diffHic](https://bioconductor.org/packages/diffHic)
-------

Changes in version 1.8.0:

- Streamlined filterDirect(), filterTrended(), and added tests for
  them. Also allowed specification of which assay to use for the data
  and reference objects.

- enrichedPairs() and neighborCounts() now return counts for
  neighbourhood regions, not just the enrichment values.

- filterPeaks() will compute (and optionally return) enrichment values
  from neighbourhood counts.

- normalizeCNV() and correctedContact() allow specification of which
  assay matrix to use from the SE objects.

- Refactored a great deal of the C++ code for improved clarity.

- Overhauled handling of DNase Hi-C data, so that pseudo-fragments are
  no longer necessary. Most functions now automatically recognise
  DNase-C data from an empty GRanges in param$fragments.  Deprecated
  segmentGenome() and prepPseudoPairs(), added the emptyGenome()
  function.

- Updated user's guide.

[DOSE](https://bioconductor.org/packages/DOSE)
----

Changes in version 3.1.3:

- output expected sample gene ID when input gene ID type not match
  <2017-03-27, Mon>

- dotplot for gseaResult <2016-11-23, Fri> +
  https://github.com/GuangchuangYu/DOSE/issues/20

Changes in version 3.1.2:

- in gseaplot, call grid.newpage only if dev.interactive() <2016-11-16,
  Wed>

- change minGSSize < geneSet_size & geneSet_size < maxGSSize to
  minGSSize <= geneSet_size & geneSet_size <= maxGSSize <2016-11-16,
  Wed>

- fixed show method issus of unknown setType in clusterProfiler::GSEA
  output <2016-11-15, Tue>

- throw more friendly error msg if fail to determine setType
  automatically in setReadable function <2016-11-15, Tue> +
  https://support.bioconductor.org/p/89445/#89479

- apply minGSSize and maxGSSize to fgsea <2016-11-14, Mon>

- change summary to as.data.frame in internal calls to prevent warning
  message <2016-11-14, Mon>

Changes in version 3.1.1:

- update startup message <2016-11-09, Wed>

- fixed parallel in Windows (not supported) <2016-10-24, Mon> +
  https://github.com/GuangchuangYu/DOSE/issues/16

- options(DOSE_workers = x) to set using x cores for GSEA analysis is
  removed. <2016-10-24, Mon> instead let MulticoreParam() to decide how
  many cores (can be set by `options(mc.cores = x)`).

[dupRadar](https://bioconductor.org/packages/dupRadar)
--------

Changes in version 1.5.3:

- Added support for Picard version >1.92, older versions no longer
  suported.

Changes in version 1.5.2:

- `CITATION` file was added.

[EBImage](https://bioconductor.org/packages/EBImage)
-------

Changes in version 4.18.0:

NEW FEATURES

- new arguments to 'display()' enabling control over layout and
  appearance of image grid in "raster" mode: 'nx' (number of frames in
  a row), 'drawGrid' (draw lines between frames), 'spacing' (separation
  between frames) and 'margin' (outer margin around the image)

- new function 'clahe()' for improving local contrast in images by
  performing Contrast Limited Adaptive Histogram Equalization

- re-introduced 'output.origin' argument to 'rotate()'

SIGNIFICANT USER-VISIBLE CHANGES

- object masks returned by `bwlabel()`, `propagate()`, and
  `watershed()`, as well as the result of thresh() are now of storage
  mode integer rather than double

- binary kernels constructed by 'makeBrush()' have storage mode integer
  (previously double)

- 'rmObjects()' and 'reenumerate()' now require input of storage mode
  integer

- 'untile()' and morphology operations preserve data storage mode

- modified boundary behaviour of 'thresh()' to reduce artifacts at
  image borders and to match the output of a corresponding call to
  'filter2()'

- added the ability for different boundary values for different frames
  in 'filter2()' linear mode (https://github.com/aoles/EBImage/pull/11)

- removed defunct '...GreyScale' family morphological functions

PERFORMANCE IMPROVEMENTS

- significantly improved performance of 'transpose()', 'getFrame()' and
  'getFrames()' by using C implementation

- numerous small improvements to execution time and memory consumption
  across the whole package, mostly by avoiding storage mode conversion
  and object duplication in C

BUG FIXES

- proper origin handling in 'resize()'

- import 'methods::slot' (fixes
  https://github.com/rstudio/blogdown/issues/17)

- fixed a bug in 'filter2()'
  (https://github.com/aoles/EBImage/issues/8)

- proper check of filter size in 'thresh()' and rectified behavior when
  filter dimensions are equal to image dimensions

- correct computation of 'selfComplementaryTopHat()'

- address PROTECT errors reported by Tomas Kalibera's 'maacheck' tool
  (https://stat.ethz.ch/pipermail/bioc-devel/2017-April/010771.html)

- fixed class retention in 'colorLabels()', 'colormap()', 'rgbImage()',
  'stackObjects()', 'tile() and untile()'

[edgeR](https://bioconductor.org/packages/edgeR)
-----

Changes in version 3.18.0:

- roast.DGEList(), mroast.DGEList(), fry.DGEList() and camera.DGEList()
  now have explicit arguments instead of passing arguments with ... to
  the default method.

- New function scaleOffset() to ensure scale of offsets are consistent
  with library sizes.

- Added decideTests() S3 methods for DGEExact and DGELRT objects. It
  now works for F-tests with multiple contrasts.

- Report log-fold changes for redundant contrasts in F-tests with
  multiple contrasts.

- Modified plotMD() S3 method for DGELRT and DGEExact objects. It now
  automatically uses decideTests() and highlights the DE genes on the
  MD plot.

- New argument 'plot' in plotMDS.DGEList().

- Removed S3 length methods for data objects.

- gini() now support NA values and avoids integer overflow.

[EGSEA](https://bioconductor.org/packages/EGSEA)
-----

Changes in version 1.3.2 (2017-03-16):

- Added: features to allow partitioning GO collections into GO domains

Changes in version 1.3.1 (2017-01-28):

- Fixed: bug in the row names of limmaTopTable. Thanks to Ali Jalali
  from MD Anderson for reporting it.

- Fixed: bug in plotSummaryHeatmap when there is only one single
  contrast.

- Added: "avg.logFC.Dir" to the EGSEA scores

- Improved: the plotSummaryHeatmap function to work with Direction
  scores

- Modified: EGSEA scores to be all small letters and updated
  egsea.sort() accordingly

- Added: median to combining p-values

- Added: plotBars function

- Fixed: bug in buildMSigDBIdx when no genes mapped to c5

- Improved: the wrapper I/O interfaces to become standard

- Improved: the implementation of GSVA by parallelizing the
  calculations on gene sets and calculating the gene set scores using
  the whole expression matrix

- Improved: the parallelization of several wrappers

- Added: ability to accept a design matrix with an intercept and a
  contrast vector of coefficient indexes.

- Added: a function to optimize the number of cores to be used for
  running EGSEA. It helps to avoid CPU overloading.

- Added: an information about the runnign time of the analysis.

- Added: a new way of report generation that completely depends on the
  EGSEAResults object. This allows users to re-generate their reports
  with different parameter values, e.g., display.top, sort.by,
  sum.plot.axis, sum.plot.cutoff.

- Added: summary heatmaps and bar plots to the report.

- Improved: the colour scheme of the summary heatmaps and bar plots.

- Fixed: bug in visualizations when log10(x) = Inf

- Added: fdr.cutoff to the calculation of Significance Score and
  Regulation Direction.

- Improved: the colour of summary heatmaps.

- Modified: buildMSigDBIdx to work with C5 collection of version 5.2

[EnrichedHeatmap](https://bioconductor.org/packages/EnrichedHeatmap)
---------------

Changes in version 1.5.1:

- anno_enriched(): can visualize positive signals and negative signals
  separatedly

- add rbind.normalizedMatrix function

[EnrichmentBrowser](https://bioconductor.org/packages/EnrichmentBrowser)
-----------------

Changes in version 2.4.1:

- Adding a min.cpm filter for RNA-seq data to de.ana

[ensembldb](https://bioconductor.org/packages/ensembldb)
---------

Changes in version 1.99.13:

USER VISIBLE CHANGES

- Most filter classes are now imported from the AnnotationFilter
  package.

- Parameter 'filter' supports now filter expression.

- Multiple filters can be combined with & and |.

- buildQuery is no longer exported.

Changes in version 1.99.11:

BUG FIXES

- ensDbFromGtf failed to fetch sequence length for some ensemblgenomes
  versions.

NEW FEATURES

- Retrieving also the taxonomy ID from the Ensembl databases and
  storing this information into the metadata table.

Changes in version 1.99.10:

BUG FIXES

- Fix problem on Windows systems failing to download files from Ensembl
  servers.

Changes in version 1.99.6:

BUG FIXES

- MySQL database name for useMySQL was not created as expected for
  GTF/GFF based EnsDbs.

Changes in version 1.99.5:

NEW FEATURES

- OnlyCodingTxFilter is now exported. This filter allows to query for
  protein coding genes.

Changes in version 1.99.3:

BUG FIXES

- Add two additional uniprot table columns to internal variable and fix
  failing unit test.

- Add two additional uniprot table columns to internal variable and fix
  failing unit test.

NEW FEATURES

- UniprotdbFilter and UniprotmappingtypeFilter.

USER VISIBLE CHANGES

- Fetching Uniprot database and the type of mapping method for Uniprot
  IDs to Ensembl protein IDs: database columns uniprot_db and
  uniprot_mapping_type.

Changes in version 1.99.2:

BUG FIXES

- Perl script is no longer failing if no chromosome info is available.

Changes in version 1.99.1:

BUG FIXES

- No protein table indices were created when inserting an EnsDb with
  protein data to MySQL.

Changes in version 1.99.0:

NEW FEATURES

- The perl script to create EnsDb databases fetches also protein
  annotations.

- Added functionality to extract protein annotations from the database
  (if available) ensuring backward compatibility.

- Add proteins vignette.

USER VISIBLE CHANGES

- Improved functionality to fetch sequence lengths for chromosomes from
  Ensembl or ensemblgenomes.

[ensemblVEP](https://bioconductor.org/packages/ensemblVEP)
----------

NOTE: As of Ensembl release 88 the name of the script was changed
      from variant_effect_predictor.pl to vep.

NEW FEATURES

    o add support for Ensembl release 85-88

MODIFICATIONS

    o document parseCSQToGRanges() behavior when no 'CSQ' data are found 

    o parseCSQToGRanges() returns mcols with names from CSQ field when 
      CSQ is present but empty 

    o add DBI perl module to SystemRequirements 

[epivizr](https://bioconductor.org/packages/epivizr)
-------

Changes in version 2.5:

- Add 'revisualize' method to add a new visualization using the same
  measurements as an existing visualization

- can save an 'EpivizApp' to disk as an 'rda' file and restart it using
  the 'restartEpiviz' function

- can use measurements from a remote epiviz UI server session to create
  visualizations from R. With this, remote epiviz UI sessions are now
  fully scriptable through R.

[epivizrData](https://bioconductor.org/packages/epivizrData)
-----------

Changes in version 999.999:

- This NEWS file is only a placeholder. The version 999.999 does not
  really exist. Please read the NEWS on Github: <URL:
  https://github.com/epiviz/epivizrData>

[epivizrServer](https://bioconductor.org/packages/epivizrServer)
-------------

Changes in version 999.999:

- This NEWS file is only a placeholder. The version 999.999 does not
  really exist. Please read the NEWS on Github: <URL:
  https://github.com/epiviz/epivizrServer>

[epivizrStandalone](https://bioconductor.org/packages/epivizrStandalone)
-----------------

Changes in version 999.999:

- This NEWS file is only a placeholder. The version 999.999 does not
  really exist. Please read the NEWS on Github: <URL:
  https://github.com/epiviz/epivizrStandalone>

[esetVis](https://bioconductor.org/packages/esetVis)
-------

Changes in version 1.1.2:

- fixes for ggplot 2.2.0: count for stat_bin_hex, center title

- fix use axis.text.y twice

Changes in version 1.1.1:

- add Rbuildignore

- test with rbokeh 0.5.0

- correct documentation to avoid warning R CMD check

[FamAgg](https://bioconductor.org/packages/FamAgg)
------

Changes in version 1.3.3:

NEW FEATURES

- Parameter id added to the findFounders method, which allows to find
  the founder couple for the pedigree of the specified individual.

Changes in version 1.3.2:

BUG FIXES

- FAData and pedigree<- ensure now that the IDs for individuals are
  unique, even across families.

Changes in version 1.3.1:

NEW FEATURES

- removeSingletons method for FAData objects.

SIGNIFICANT USER-VISIBLE CHANGES

- Additional argument family in buildPed.

[fCCAC](https://bioconductor.org/packages/fCCAC)
-----

Changes in version 1.1.1:

- Updated reference Madrigal (2016), Bioinformatics.

[fgsea](https://bioconductor.org/packages/fgsea)
-----

Changes in version 1.1.2:

- Fixed building issues

Changes in version 1.1.1:

- Fixed bug: not using Multicore on Windows

[flowPloidy](https://bioconductor.org/packages/flowPloidy)
----------

Changes in version 1.1.9 (2017-04-19):

User Visible Changes

- Users may now specify one or more internal standard sizes in
  functions FlowHist() and batchFlowHist(), via the argument
  'standards'. These values will be presented to the user in the
  browseFlowHist() viewer, and after being set by the user, the value
  will be used to calculate the GC size in pg in the function
  tabulateFlowHist.

- Users may now select which peak in the histogram to treat as the
  internal standard when calculating GC values.

- New vignettes added: "Getting Started", and "Histogram Tour".

- Old vignette "overview" removed.

- Internal help pages greatly expanded, including many internal
  functions. See ?flowPloidy for an overview

- Many minor bug fixes and GUI tweaks (for browseFlowHist).

Changes in version 1.1.3 (2016-11-25):

User Visible Changes

- Gating fully implemented in browseFlowHist! Major reorganization of
  the browseFlowHist layout.

Changes in version 1.1.2 (2016-11-22):

User Visible Changes

- Added support for processing files with two standards present. A new
  argument is available for functions that load FCS files
  (batchFlowHist, FlowHist etc.): `samples`. By default this is set to
  2, to account for a single unknown and an co-chopped standard. If you
  are using two co-chopped standards (or really anytime you have three
  distinct samples chopped together), set samples = 3. This can also be
  changed interactively in the browseFlowHist GUI.

- The layout of the browseFlowHist GUI has been re-arranged somewhat to
  accomodate the new features mentioned above.

- The linearity flag is over-ridden when no G2 peaks are present.
  Without a G2 peak, linearity can't be properly fit. This leads to
  linear gradients, because the linearity parameter is used in the
  S-phase components.

Changes in version 1.1.1 (2016-10-26):

Internal Changes

- Improved peak finding algorithm

- Reduced region searched for the starting bin for model fitting. Was
  originally 20, now set to 10. Need more data to establish best
  approach.

[FlowSOM](https://bioconductor.org/packages/FlowSOM)
-------

Changes in version 1.7.2:

NEW FEATURES

- Update to include the seed parameter in the vignette

Changes in version 1.7.1:

NEW FEATURES

- Update to make sure github and bioconductor contain the same
  functionality

[FunChIP](https://bioconductor.org/packages/FunChIP)
-------

Version: 1.0.1
Text:

[gCMAP](https://bioconductor.org/packages/gCMAP)
-----

Changes in version 1.19.4:

- CHANGE: As the bigmemory package is not available for Windows, only
  unix_type OS is supported.

Changes in version 1.19.1:

- BUGFIX: Updated NAMESPACE file to conform with R CMD check.

- BUGFIX: reactome2cmap function is available again.

- NEW: Added citation information.

[gCMAPWeb](https://bioconductor.org/packages/gCMAPWeb)
--------

Changes in version 1.15.1:

- NEW: Added citation information.

[gdsfmt](https://bioconductor.org/packages/gdsfmt)
------

Changes in version 1.12.0:

UTILITIES

- update liblzma to v5.2.3

- update lz4 to v1.7.5

- a new citation

Changes in version 1.10.1:

NEW FEATURES

- new data types (variable-length encoding of signed and unsigned
  integers)

[genbankr](https://bioconductor.org/packages/genbankr)
--------

Changes in version 1.3.2:

BUGFIXES

- Circularity (from LOCUS header information) is now applied to all
  sources in the file.

- Improved assertion related to id ordering in makeTxDbFromGenBank
  which passed in previous R version but was failing in recent ones

[geneClassifiers](https://bioconductor.org/packages/geneClassifiers)
---------------

Changes in version 1.0.0:

- The first publically available version of the geneClassifiers
  package. This packages currently contains a number of gene
  classifiers relating to survival in Multiple Myeloma.

[GeneNetworkBuilder](https://bioconductor.org/packages/GeneNetworkBuilder)
------------------

Changes in version 1.17.5:

- fix a bug in browseNetworkOutput

Changes in version 1.17.4:

- add cytoscape-exportbox

Changes in version 1.17.3:

- make the local copy of javascript library

Changes in version 1.17.2:

- add function exportNetwork.

- add cytoscape-searchbox

- update documents.

Changes in version 1.17.1:

- add function browseNetwork.

[geneplast](https://bioconductor.org/packages/geneplast)
---------

Changes in version 1.2.0:

- Improved documentation and vignette.

[GENESIS](https://bioconductor.org/packages/GENESIS)
-------

Changes in version 2.6.0:

- Major bug fix: assocTestSeq no longer drops some variants from
  aggregate tests in the case where the same variants are included in
  more than one aggregate unit.

- Added function for analysis of admixture mapping data.

[geNetClassifier](https://bioconductor.org/packages/geNetClassifier)
---------------

Changes in version 1.15.0:

- Update for compatibility with new R versions. The confusion matrix is
  now an object of class 'table'.

[genomation](https://bioconductor.org/packages/genomation)
----------

Changes in version 1.7.3:

IMPROVEMENTS AND BUG FIXES

- added OS check for tests involving BigWig files. They can not be read
  in windows OS.

Changes in version 1.7.2:

IMPROVEMENTS AND BUG FIXES

- Now we rely on Rsamtools idxStatsBam function to calculate rpm,
  previously it was a cpp function written by alexg9010.

Changes in version 1.7.1:

NEW FUNCTIONS AND FEATURES

- scoreMatrixBin() calculates coverage over windows that are not only
  GRanges, but also GRangesList. It's usefull for calculating
  transcript coverage of a set of exons.

- ScoreMatrix-like functions work with bigWig files and supplied
  weight.col and is.noCovNA=TRUE

IMPROVEMENTS AND BUG FIXES

- Added warning for rpm=TRUE and type='bigWig'

- type='auto' by default in ScoreMatrix-like functions

- narrowPeak() and broadPeak() are 0-based by default (#144 fixed)

- Fixed error in readGeneric when reading files with numeric
  chromosomes (#133 fixed)

- Show warning if windows fall off target

- Show error if windows have width 1

[GenomeInfoDb](https://bioconductor.org/packages/GenomeInfoDb)
------------

Changes in version 1.12.0:

NEW FEATURES

- Add function standardChromosomes()

- Seqlevels() setter now supports "fine" and "tidy" modes on
  GRangesList and GAlignmentsList objects

- Add assembly_accessions dataset

MODIFICATIONS

- Updated mapping table between UCSC and Ensembl to include recent
  builds

- Use https instead of http to fetch stuff from NCBI

- Replace 'force=TRUE' with 'pruning.mode="coarse"' in seqlevels()
  setter

- Add 'pruning.mode' argument to the keepSeqlevels(), dropSeqlevels(),
  and keepStandardChromosomes() functions. IMPORTANT NOTE: Like for the
  seqlevels() setter, the default pruning mode is "error", which means
  that now these functions fail when some of the seqlevels to drop from
  'x' are in use. The old behavior was to silently prune 'x' (doing
  "coarse" pruning)

- Update files in data directory

- Updated internal functions .lookup_refseq_assembly_accession() and
  fetch_assembly_report() for speed and efficiency

- move some files from GenomeInfoDb/data/ to GenomeInfoDbData
  annotation package

BUG FIXES

- fetch_assembly_summary() updated to work with recent changes to
  format of files assembly_summary_genbank.txt and
  assembly_summary_refseq.txt

[genomeIntervals](https://bioconductor.org/packages/genomeIntervals)
---------------

Changes in version 1.31.1:

- the signature of rank() has now ... and "last" as option in
  ties.method (for compatibility only -- not supported)

[GenomicFeatures](https://bioconductor.org/packages/GenomicFeatures)
---------------

Changes in version 1.28:

NEW FEATURES

- makeTxDbFromUCSC() supports new composite "NCBI RefSeq" track for
  hg38.

- Add 'metadata' argument to makeTxDbFromGFF().

- Add exonicParts() as an alternative to disjointExons(): -
  exonicParts() has a 'linked.to.single.gene.only' argument (FALSE by
  default) that is similar to the 'aggregateGenes' argument of
  disjointExons() but with opposite meaning. More precisely
  'exonicParts(txdb, linked.to.single.gene.only=TRUE)' returns the same
  exonic parts as 'disjointExons(txdb, aggregateGenes=FALSE)'.  -
  Unlike 'disjointExons(txdb, aggregateGenes=TRUE)', 'exonicParts(txdb,
  linked.to.single.gene.only=FALSE)' does NOT discard exon parts that
  are not linked to a gene.  - exonicParts() is almost twice more
  efficient than disjointExons().

- Add intronicParts(): similar to exonicParts() but returns intronic
  parts.

SIGNIFICANT USER-VISIBLE CHANGES

- Some work on distance,GenomicRanges,TxDb method: - pass ignore.strand
  to range() and distance() - return NA when 'id' cannot be collapsed
  into a single range or when 'id' is not found in 'y'

DEPRECATED AND DEFUNCT

- Argument 'force' of seqlevels() setters is deprecated in favor of new
  and more flexible 'pruning.mode' argument.

- Remove the 'vals' argument of the "transcripts", "exons", "cds", and
  "genes" methods for TxDb objects (was defunct in BioC 3.4).

BUG FIXES

- Fix bug in seqlevels() setter for TxDb objects reported here:
  https://support.bioconductor.org/p/90226/

[GenomicRanges](https://bioconductor.org/packages/GenomicRanges)
-------------

Changes in version 1.28.0:

NEW FEATURES

- Add coercion from ordinary list to GRangesList. Also the
  GRangesList() constructor function now accepts a list of GRanges as
  input (and just calls new coercion from list to GRangesList on it
  internally).

- seqlevels() setter now supports "fine" and "tidy" pruning modes on
  GRangesList objects (in addition to "coarse" mode, which is the
  default).

- "range" methods now have a 'with.revmap' argument (like "reduce" and
  "disjoin" methods).

- Add a bunch of range-oriented methods for GenomicRangesList objects.

SIGNIFICANT USER-LEVEL CHANGES

- Some changes/improvements to "precede" and "follow" methods for
  GenomicRanges objects motivated by discussion on support site:
  https://support.bioconductor.org/p/90664/

- Some changes/improvements to "rank" method for GenomicRanges objects:
  - now supports the same ties methods as base::rank() (was only
  supporting ties methods "first" and "min" until now) - default ties
  method now is "average", like base::rank() - now supports additional
  argument 'ignore.strand'.

DEPRECATED AND DEFUNCT

- Argument 'force' of seqinfo() and seqlevels() setters is deprecated
  in favor of new and more flexible 'pruning.mode' argument.

BUG FIXES

- Fix severe performance regression introduced in Bioconductor 3.3 in
  "intersect" and "setdiff" methods for GRangesList objects. Thanks to
  Jens Reeder <reeder.jens@gene.com> for catching and reporting this.

[GenomicScores](https://bioconductor.org/packages/GenomicScores)
-------------

Changes in version 0.99.0:

USER VISIBLE CHANGES

- Submission of the first version to the Bioconductor project. (start
  date: March 17, 2017)

[GenVisR](https://bioconductor.org/packages/GenVisR)
-------

Changes in version 1.4.8:

- fixed bug in lolliplot introduced by ensembl server changes

Changes in version 1.4.7:

- patch to fix possible conflicts between the rmvSilent=T and
  fileType="custom" parameter in waterfall()

Changes in version 1.4.6:

- fixed typo with new feature in cnFreq

Changes in version 1.4.5:

- fixed issue where plot would not render/take awhile for functions
  requiring a genome and using hg38

- Dramatically increased performance of cnFreq()

- cnFreq() no longer requires genomic segments to be identical across
  samples

- cnFreq() can now selectively plot chromosomes

- cnFreq() can now plot both frequency/proportion regardless of data
  input type

Changes in version 1.4.4:

- fixed bug where incorrect cell labels were overlayed with proper ones
  when removing silent mutations

Changes in version 1.4.3:

- fixed bug in cnFreq where cnFreq would produce an error (improperly)
  when looking for consistent windows

Changes in version 1.4.2:

- fixed bug caused by updates to gridExtra when using genCov

Changes in version 1.4.1:

- patch to fix alignments due to ggplot2 update

[ggtree](https://bioconductor.org/packages/ggtree)
------

Changes in version 1.7.11:

- remove layout.method parameter <2017-04-20, Thu> +
  https://github.com/GuangchuangYu/ggtree/issues/118#issuecomment-295130818
  + https://github.com/GuangchuangYu/ggtree/issues/125

Changes in version 1.7.10:

- add message for subview, inset, phylopic, theme_transparent and
  theme_inset <2017-03-23, Thu> + will be defunct in version >= 1.9.0 +
  user should use ggimage package to annotate tree with graphic object
  or image file

- update subview to support mainview produced by `ggplot() + layers`
  <2017-03-13, Mon>

Changes in version 1.7.9:

- fixed geom_range to support height_0.95_HPD <2017-03-03, Fri>

- fixed geom_tiplab(geom='label') <2017-03-02, Thu> +
  https://github.com/GuangchuangYu/ggtree/issues/115

Changes in version 1.7.8:

- get_taxa_name now sorted by taxa position and also support whole tree
  <2017-03-01, Wed>

- unrooted layout support branch.length="none", fixed #114 <2017-03-01,
  Wed>

- remove apeBootstrap and raxml object support as they were removed
  from treeio <2017-02-28, Tue>

Changes in version 1.7.7:

- supports parse="emoji" in geom_cladelabel, geom_text2, geom_label2,
  geom_tiplab, geom_tiplab2 <2017-02-16, Thu>

- aes(subset) now support logical vector contains NA <2017-02-16, Thu>

- add legend transparency to theme_transparent <2017-02-13, Mon> +
  <https://github.com/GuangchuangYu/ggtree/pull/112>

- update citation info <2017-01-20, Fri>

Changes in version 1.7.6:

- inset support reverse scale <2017-01-05, Thu> +
  https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!msg/bioc-ggtree/_JPfm71Z8nM/6gL93oxHFQAJ

Changes in version 1.7.5:

- disable labeling collapsed node as tip <2017-01-03, Tue> +
  https://groups.google.com/forum/#!topic/bioc-ggtree/nReqJatMvJQ

- fortify.phylo4d via converting phylo4d to treedata object
  <2016-12-28, Wed>

- improve viewClade function, use coord_cartesian instead of xlim
  <2016-12-28, Wed>

- remove codes that move to treeio and now ggtree depends treeio
  <2016-12-20, Tue>

Changes in version 1.7.4:

- is.ggtree function to test whether object is produced by ggtree
  <2016-12-06, Tue>

- now branch.length can set to feature available in phylo4d@data and
  yscale is supported for phylo4d object <2016-12-06, Tue>

- bug fixed of rm.singleton.newick, remove singleton parent instead of
  singleton <2016-12-01, Thu>

- reorder phylo to postorder before ladderrize <2016-11-28, Mon>

- allow yscale to use data stored in phylo4d object <2016-11-24, Thu> +
  https://github.com/GuangchuangYu/ggtree/issues/98

- groupOTU method now accept 'overlap = c("overwrite", "origin",
  "abandon")' parameter <2016-11-16, Wed> +
  https://groups.google.com/forum/#!topic/bioc-ggtree/Q4LnwoTf1DM

Changes in version 1.7.3:

- drop.tip method for NHX object <2016-11-11, Fri>

- update startup message <2016-11-09, Wed>

- reverse timescale x-axis <2016-11-07, Mon> +
  https://github.com/GuangchuangYu/ggtree/issues/87

Changes in version 1.7.2:

- make missing colors in gheatmap invisible (previously use 'white')
  <2016-11-03, Thu>

- xlim_expand for setting x axis limits of specific panel <2016-11-01,
  Tue> + xlim_tree is now a specific case of xlim_expand(xlim,
  panel='Tree')

- bug fixed of parsing tree text in beast file <2016-10-31, Mon> +
  https://github.com/GuangchuangYu/ggtree/issues/84

Changes in version 1.7.1:

- xlim_tree layer and test <2016-10-31, Mon> + set x axis limits for
  Tree panel for facet_plot

- update read.nhx <2016-10-30, Sun> + add tip numbers to @nhx_tags and
  add tests + https://github.com/GuangchuangYu/ggtree/pull/83 + store
  nhx_tags$node as numeric values <2016-10-31, Mon>

- facet_plot supports ggbio::geom_alignment <2016-10-26, Wed> +
  https://github.com/tengfei/ggbio/issues/83

- make tree stats available in facet_plot <2016-10-24, Mon>

[GISPA](https://bioconductor.org/packages/GISPA)
-----

Changes in version 0.99.0:

- First Submission to Bioconductor Bhakti Dwivedi and Jeanne Kowalski
  The Winship Cancer Institute, Emory University
  https://bbisr.winship.emory.edu/

[Glimma](https://bioconductor.org/packages/Glimma)
------

Changes in version 1.3.0:

- Added highlighting to bars in MDS plot.

- Added interaction with table when clicking on points in MD plot.

- Changed expressions to default to no transformation.

- Changed default colours.

- Changed style of table.

- Changed size of highlighted points.

[globalSeq](https://bioconductor.org/packages/globalSeq)
---------

Changes in version 1.4.0 (2016-04-25):

- Minor improvements

[GOexpress](https://bioconductor.org/packages/GOexpress)
---------

Changes in version 1.9.2:

BUG FIXES

- Remove function overlap_GO because the VennDiagram package has issues
  with condition tests on vectors with length greater than 1

Changes in version 1.9.1:

BUG FIXES

- Remove duplicated plot.title argument in ggplot call.

GENERAL UPDATES

- Minor code cleaning.

[GOpro](https://bioconductor.org/packages/GOpro)
-----

Version: 1.1.3
Text:

[GOSemSim](https://bioconductor.org/packages/GOSemSim)
--------

Changes in version 2.1.3:

- friendly error message for using IC method without IC computed
  <2017-02-17, Fri> +
  https://github.com/GuangchuangYu/GOSemSim/issues/11

- fixed https://github.com/GuangchuangYu/GOSemSim/issues/9 <2016-12-20,
  Tue>

Changes in version 2.1.2:

- use prettydoc for vignette <2016-11-30, Wed>

- remove using BiocStyle <2016-11-23, Wed>

Changes in version 2.1.1:

- update startup message <2016-11-09, Wed>

[gQTLstats](https://bioconductor.org/packages/gQTLstats)
---------

Version: 1.8
Category: NEW FEATURES
Text: AllAssoc() now reports a Z-score for HWE

Version: 1.8
Category: NEW FEATURES
Text: tqbrowser() facilitates interactive viewing of trans

Version: 1.8
Category: associations
Text:

Version: 1.8
Category: SIGNIFICANT USER-VISIBLE CHANGES
Text: reliance on GGtools has been eliminated

[graphite](https://bioconductor.org/packages/graphite)
--------

Changes in version 1.21.1 (2017-04-24):

- Updated all pathway data.

[GRridge](https://bioconductor.org/packages/GRridge)
-------

Changes in version 0.99.8 (2016-10-19):

- Added functions

- Added help pages

- Added vignette

[GSVA](https://bioconductor.org/packages/GSVA)
----

Changes in version 1.24:

BUG FIXES

- Bugfixes on the parallel execution of gsva() with bootstrap
  calculations.

[gtrellis](https://bioconductor.org/packages/gtrellis)
--------

Changes in version 1.7.1:

- size of axis and title are correctly calculated

- add_track supports raster image

[Gviz](https://bioconductor.org/packages/Gviz)
----

Changes in version 1.20.0:

BUG FIXES

- BiomartGeneRegionTracks can now deal with a featureMap list to
  provide alternative conditional mappings for different Biomarts.

[GWASTools](https://bioconductor.org/packages/GWASTools)
---------

Changes in version 1.21.1:

- Replace ZIP_RA with LZMA_RA for GDS compression.

- Default is no compression for genotypes.

[hiAnnotator](https://bioconductor.org/packages/hiAnnotator)
-----------

Changes in version 1.9.0:

- code spacing edits & GenomicRanges package update adjustments

[HIBAG](https://bioconductor.org/packages/HIBAG)
-----

Changes in version 1.11.1:

- change "hg20" to "hg38" according to the UCSC Genome Browser datasets
  and documentation

- add "DRB3" and "DRB4" to the HLA gene list

[hicrep](https://bioconductor.org/packages/hicrep)
------

Version: 0.99.1
Text: 1. added unit tests 2. adjusted spaces for the source code 3.
        added some missing part in documentation

Version: 0.99.0
Category: INITIAL RELEASE
Text:

[HilbertCurve](https://bioconductor.org/packages/HilbertCurve)
------------

Changes in version 1.5.2:

- `hc_map` supports to add labels under pixel mode

[hpar](https://bioconductor.org/packages/hpar)
----

Changes in version 1.17.2:

- Update to HPA version 16.1 (2017.01.31) <2017-02-14 Tue>

Changes in version 1.17.1:

- Using travis and codecov <2016-12-22 Thu>

- Migrate vignette to BiocStyle's html2 <2016-12-22 Thu>

Changes in version 1.17.0:

- Bioconductor devel 3.5

[HTSeqGenie](https://bioconductor.org/packages/HTSeqGenie)
----------

Changes in version 4.5.1:

- changed variant calling tests in accordance with changes made to
  bam_tally in gmapR

[HTSFilter](https://bioconductor.org/packages/HTSFilter)
---------

Changes in version 1.14.1:

- -- Fixed bug for missing parallel option when a CountDataSet is used
  with HTSFilter

[ideal](https://bioconductor.org/packages/ideal)
-----

Changes in version 0.99.0:

NEW FEATURES

- Ready for Bioc submission

- Completed the news

Changes in version 0.9.1:

NEW FEATURES

- Added Instructions fully from rendered version of the vignette to
  have available at runtime

- Added support for downloading all plots and tables

Changes in version 0.9.0:

NEW FEATURES

- Interactive tours are covering now all tabs, with extensive
  walkthroughs for the user

- Added all screenshots to vignette

Changes in version 0.6.2:

NEW FEATURES

- Interactive tours are now available, coded in external files

- Travis-CI is now supported

Changes in version 0.6.0:

NEW FEATURES

- Added MA plot with extra custom list to avoid manual selection of
  many genes

- MA plot function now automatically supports subset of gene to be
  extra plotted

- Added documentation with roxygen to all functions

- Heatmap functions for genes annotated to a GO term as signature

- Template report also provided

- Full draft of vignette now available, working towards bioc submission

- Added textual help to all sections, with collapsible element

- Added proof of principle to have interactive tours based on rintrojs

Changes in version 0.4.0:

NEW FEATURES

- Gene box info added, based on rentrez

- New look for MA plots and volcano plots

Changes in version 0.3.0:

NEW FEATURES

- Restructuring of the folders done, package can be correctly
  installed, loaded - namespace, description are set up

Changes in version 0.2.0:

NEW FEATURES

- Correct structure of the package

Changes in version 0.1.0:

NEW FEATURES

- Package created!

[ImpulseDE2](https://bioconductor.org/packages/ImpulseDE2)
----------

Changes in version 0.99.0 (2017-03-03):

- Initial release

[InPAS](https://bioconductor.org/packages/InPAS)
-----

Changes in version 1.7.5:

- add Julie as co-maintainer.

Changes in version 1.7.4:

BUG FIXES

- fix the bug if there is NA values for proximal site to be adjusted.

Changes in version 1.7.3:

BUG FIXES

- remove the code modified from voom. And directly use diffSplice.

Changes in version 1.7.2:

BUG FIXES

- Fix the memory trap if the dataset is huge when call CPsites.

Changes in version 1.7.1:

BUG FIXES

- Fix the bug for empty bedgraph.

[INSPEcT](https://bioconductor.org/packages/INSPEcT)
-------

Changes in version 1.5.5:

- updated the results on the simulated data in the vignette

Changes in version 1.5.4:

- modified the argument strandSpecific in makeRPKMs function, so that
  now the user can perform strand-specific read counting with this
  possible modes: 0 => unstranded 1 => stranded 2 => reversely stranded

Changes in version 1.5.3:

- modified internal functions inferKBetaFromIntegral,
  inferKBetaFromIntegralWithPre, inferKGammaFromIntegral in order to
  reduce the number of the missing values (NA) in the output

Changes in version 1.5.2:

- Updated the documentation

Changes in version 1.5.1:

- Fixed a bug that caused a mis-choiche of sigmoid or impulse function
  during modeling

[InteractionSet](https://bioconductor.org/packages/InteractionSet)
--------------

Changes in version 1.4.0:

- Deprecated anchors<- in favour of anchorIds<-, to avoid confusion
  about 'value' type.

- Added first(), second() functions for convenience.

- Updates to documentation, tests.

[IPO](https://bioconductor.org/packages/IPO)
---

Changes in version 1.1.2:

- vignette updated

- plot margins omitted

Changes in version 1.1.1:

- use package BiocParallel (via argument `BPPARAM`) instead of
  `nSlaves` to controll `xcms`-parallelization

- depends on xcms >= 1.50.0

- formatting `writeRScript` to output more beautifully

Changes in version 1.1.0:

- merge Bioconductor 1.0.0 release code with Github code.

[IRanges](https://bioconductor.org/packages/IRanges)
-------

Changes in version 2.10.0:

NEW FEATURES

- "range" methods now have a 'with.revmap' argument (like "reduce" and
  "disjoin" methods).

- Add coercion from list-like objects to IRangesList objects.

- Add "table" method for SimpleAtomicList objects.

- The "gaps" method for CompressedIRangesList objects now uses a chunk
  processing strategy if the input object has more than 10 million list
  elements. The hope is to reduce memory usage on very big input
  objects.

BUG FIXES

- Fix "setdiff" method for CompressedIRangesList for when all ranges
  are empty.

- Fix long standing bug in coercion from Ranges to PartitioningByEnd
  when the object to coerce has names.

DEPRECATED AND DEFUNCT

- Remove the RangedDataList and RDApplyParams classes, rdapply(), and
  the "split" and "reduce" methods for RangedData objects. All these
  things were defunct in BioC 3.4.

- Remove 'ignoreSelf' and 'ignoreRedundant' arguments (replaced by
  'drop.self' and 'drop.redundant') from findOverlaps,Vector,missing
  method (were defunct in BioC 3.4).

- Remove GappedRanges class (was defunct in BioC 3.4).

[isomiRs](https://bioconductor.org/packages/isomiRs)
-------

Changes in version 1.3.5:

FEATURES

- Add isomiRs naming to documentation

- Add design to the object to get better usability

- Remove non-template addition with C/G nucleotides by default
  (canonicalAdd)

- Remove sequences with mutations and more than one miRNA hit

Changes in version 1.3.4:

FIXES

- Fix removing false mutations from the raw files. Change sequences to
  correct the nucleotide at the specific position.

FEATURES

- Improve code to remove error sequencing from raw data

- Improve code to show the raw data with isoSelect

Changes in version 1.3.3:

FIXES

- Fix data with correct headers name

Changes in version 1.3.2:

OTHERS

- Preparing migration to new isomiRs naming using mirTOP naming system

Changes in version 1.3.1:

OTHERS

- Add option to IsomirDataSeqFromFiles to decide when to consider
  mutations as reals

[IVAS](https://bioconductor.org/packages/IVAS)
----

Changes in version 1.95.6:

- Replace the foreach function to the bplapply function in BiocParallel

Changes in version 1.95.5:

- Debugging and modify sQTLsFinder

Changes in version 1.9.0:

- Modify RatioFromFPKM, Splicingfinder, and sQTLsFinder

- Modify manual and tutorial

Changes in version 1.8.0:

- Defunct MsqtlFinder, calSignificant, and sqtlfinder.

- Adjust Chr names between GTF and SNP locus data.

- Change test of UTR region.

- sqtl finder reform

- Create new functions (RatioFromFPKM, Splicingfinder, and sQTLsFinder)

- Create ASclass object

[ldblock](https://bioconductor.org/packages/ldblock)
-------

Changes in version 1.5:

NEW FEATURES

- s3_1kg() generates TabixFile references to 1000 genomes VCF in AWS S3
  bucket

- ldByGene() obtains linkage information using snpStats ld() and erma
  genemodel() to retrieve focused information from VCF

[limma](https://bioconductor.org/packages/limma)
-----

Changes in version 3.32.0:

- New function cameraPR(), which implemented a pre-ranked version of
  camera().

- New function alias2SymbolUsingNCBI(), which converts gene aliases or
  synonyms into official gene symbols using an NCBI gene-info file.

- New function wsva() for weighted surrogate variable analysis.

- New function coolmap(). This is essentially a wrapper for the
  heatmap.2() function in the ggplots package, but with sensible
  default settings for genomic log-expression data.

- decideTests() is now an S3 generic function with a default method and
  a method for MArrayLM objects.  decideTests() now selects all null
  hypotheses as rejected if p.value=1.

- length() methods removed all limma data objects (objects of class
  EList, EListRaw, RGList, MAList or MArrayLM). length(x) will now
  return the number of list components in the object rather than the
  number of elements in the expression matrix.

- New argument 'style' for volcanoplot(). The default is now to use
  -log10(p-value) for the y-axis instead of the B-statistic.

- New argument 'xlab' for barcodeplot().

- New argument 'col' for plotSA().  plotSA() now longer plots a lowess
  curve trend, but if appropriate both high and low outlier variances
  are highlighted in a different color.

- Argument 'replace.weights' removed from voomWithQualityWeights().
  The function now always produces an EList, similar to voom(). The
  default behavior of the function is unchanged.

- barcodeplot() now ranks statistics from low to high, instead of from
  high to low, following the usual style of axes in R plots.  This
  means that left and right are now interchanged.

- plotSA() now plots quarter-root variances instead of log2(variances).

- Default for 'legend' argument of plotWithHighlights() changed from
  "topleft" to "topright".

- fitFDist() now estimates the scale by mean(x) when df2 is estimated
  to be Inf.  This will make the results from eBayes() less
  conservative than before when df.prior=Inf.

- plotSA() now indicates, by way of an open plotting symbol, any points
  that have low robust df.prior values.

- Clearer error message from fitFDistRobustly() when some variances are
  zero.

- C functions are now registered using R_registerRoutines.

- Bug fix for contrastAsCoef() when there is more than one contrast.
  Previously the coefficients for the transformed design matrix were
  correct only for the first contrast.

- Bug fix for kegga() when the universe is explicitly specified.

- Bug fix for fitFDistRobustly() when there is an extreme outlier.
  Previously floating point underflow for the outlier p-value could
  cause an error.

- Bug fix to mroast(), which was ignoring 'geneid' argument.

- Bug fix to printHead() for arrays with 1 column.

[Logolas](https://bioconductor.org/packages/Logolas)
-------

Changes in version 0.99.0:

- Release

[LymphoSeq](https://bioconductor.org/packages/LymphoSeq)
---------

Changes in version 1.3.0:

- Added function alignSeq

- Added funciton phyloSeq

- Added function exportFasta

- Added function differentialAbundance

- Added function clonalRelatedness

- Added function commonSeqsBar

[M3D](https://bioconductor.org/packages/M3D)
---

Changes in version 1.9.2-3:

- Removed parallel function due to incompatibility with Windows. I am
  investigation a cross platform solution in the meantime. Due to speed
  ups, the sequential ‘lite’ function should be fast enough for all
  practical needs.

Changes in version 1.9.1:

- Bug fix for readBedFiles. Thanks to Francesca Cairoli of the
  University of Trieste for the feedback.

[maftools](https://bioconductor.org/packages/maftools)
--------

Changes in version 1.2.0:

NEW FUNCTIONS

- mafSurvival - Performs survival analysis.

- tcgaCompare - Compares mutation load from given MAF against all 33
  TCGA cohorts.

- pancanComparision - Perform PacCancer analysis/comparision

- prepareMutSig - Prepares MAF file for MutSig analysis by fixing
  descrepencies in gene symbols.

SIGNIFICANT USER-LEVEL IMPROVEMENT

- plotmafSummary has argument titvRaw. You can set it to FALSE to plot
  TiTV fraction instead of raw counts.

NON SIGNIFICANT CHANGES

- Bug fixes and improvements.

[matter](https://bioconductor.org/packages/matter)
------

Changes in version 1.1.3:

NEW FEATURES

- Added class 'drle' for delta-run-length encoding vectors

- Added '+', '-', '*', '/', '^', 'exp', 'log', 'log2', and 'log10' as
  possible delayed operations to on-disk atoms

SIGNIFICANT USER-VISIBLE CHANGES

- Slots of 'atoms' class now use delta-run-length encoding

- Reduced metadata size by changing 'atoms' class to use groups rather
  then relying on a 'list' of 'atoms'

- The 'scale' method for 'matter_mat' now matches 'scale.default' more
  correctly when 'center = FALSE' and 'scale = TRUE'

Changes in version 1.1.2:

NEW FEATURES

- Added support for char, uchar, ushort, uint, and ulong datamodes

- Added support for raw (Rbyte) matter objects

SIGNIFICANT USER-VISIBLE CHANGES

- S4 methods for matrix-specific summary statistics are now only
  defined on matter_mat and its subclasses

BUG FIXES

- Dramatically improved speed of matrix multiplication

Changes in version 1.1.1 (2016-11-29):

NEW FEATURES

- Added 'crossprod' (t(x) %*% y) and 'tcrossprod' (x %*% t(y)) methods

- Added 'atomdata' accessor method, for which 'adata' is now an alias

BUG FIXES

- Added S3 versions of some S4 methods to fix scoping issues

- Removed Cardinal package from Suggests to avoid circular dependency

- Reduced memory consumption in bigglm-matter method

[MaxContrastProjection](https://bioconductor.org/packages/MaxContrastProjection)
---------------------

Changes in version 1.0.0:

- Uploaded package to Bioconductor

[metabomxtr](https://bioconductor.org/packages/metabomxtr)
----------

Changes in version 1.9.1:

- *modified package so that normalization function does not output data
  when any model parameters are not possible to estimate.  *added
  function metabPlot to plot pre vs. post normalization metabolite
  abundances.  *updated mixnorm vignette

[MetaboSignal](https://bioconductor.org/packages/MetaboSignal)
------------

Changes in version 1.5.2:

- MetaboSignal includes a new function: "MS_interactionType()". This
  function allows getting the interaction subtype between signaling
  nodes. The output matrix generated by this function can be used for
  "MetaboSignal_NetworkCytoscape()" and also for
  "MS_GetShortestpaths()".

- "MS_GetShortestpaths()" has been modified and now the output shortest
  path(s) can be represented as a network-table (i.e. 2-column matrix).

[metagenomeFeatures](https://bioconductor.org/packages/metagenomeFeatures)
------------------

Changes in version 1.6.0 (2017-04-14):

- Bioc release cleaned up documentation and debugging MgDb class
  definition with tree slot

Changes in version 1.5.1 (2017-03-20):

- editing documentation

[metavizr](https://bioconductor.org/packages/metavizr)
--------

Changes in version 1.3.20:

- Standalone mode introduced, a version of epiviz with reduced
  capabilities is now included as part of epivizr. The epiviz web app
  is run locally using 'httpuv's http server

- Add and remove seqinfo (e.g., chromosome info) to any epiviz session

Changes in version 1.3.11:

- Add NEWS file

- Update documentation on 'slideshow' function

Changes in version 1.3.10:

- Changed default on 'slideshow' to show all ranges

Changes in version 1.3.9:

- Added 'heatmapChart' convenience function

Changes in version 1.3.8:

- Fixed bug in 'startEpiviz' not sending 'seqName' parameter correctly

Changes in version 1.3.7:

- Fixed bug in 'EpivizBpData' that sent 'metadata' info in wrong format

Changes in version 1.3.6:

- Changes slots using lists in 'EpivizDeviceMgr' to environments to
  avoid crashing RStudio due to inspection of manager objects

Changes in version 1.3.5:

- Fails gracefully on daemonization request on Windows

- Deprecates the 'proxy' argument to 'startEpiviz'

Changes in version 1.3.4:

- Upgrading to Epiviz v2 webapp

[MetCirc](https://bioconductor.org/packages/MetCirc)
-------

Changes in version 1.1.5:

- extend unit tests for shinyApp, plottingFunctions and convert2MSP
  &#91;2017-04-03 Mon&#93;

Changes in version 1.1.4:

- change MSP class, create slots mz, rt, names, classes, information
  and adduct &#91;2017-01-28 Sat&#93;

- add tabPanels in shinyCircos (Main, Appearance) &#91;2017-01-28 Sat&#93;

- rearrange position of legend, implement option to show/not show l
  egend &#91;2017-01-28 Sat&#93;

- rescale plot when changing window size, allow for further
  scaling/descaling of the plot &#91;2017-01-28 Sat&#93;

- adjust convertMSP2MSP to new class MSP, create unit tests &#91;2017-01-29
  Sun&#93;

- include msp2msp matrix, a test data set for convertMSP2MSP
  &#91;2017-01-29 Sun&#93;

- set methods for names, classes, adduct and information &#91;2017-01-29
  Sun&#93;

- change the interactive shinyCircos such that the user can update the
  annotation data of an MSP object (name, class, information and adduct
  ion information) &#91;2017-01-29 Sun&#93;

Changes in version 1.1.3:

- use new email adress &#91;2016-12-05 Mon&#93;

- use option to calculate MSP-object from msp-file directly &#91;2016-12-05
  Mon&#93;

Changes in version 1.1.2:

- use absolute masses when calculating similarities in
  createSimilarityMatrix (bug fix) &#91;2016-11-17 Thu&#93;

- add option links in highlight, i.e. should links be plotted or not?
  &#91;2016-11-17 Thu&#93;

Changes in version 1.1.1:

- change slider input as such that one is able to select a lower and
  lower bound instead of only a lower bound &#91;2016-11-04 Fri&#93;

[MethylAid](https://bioconductor.org/packages/MethylAid)
---------

Changes in version 1.9:

- adapted to changes in minfi such as, read.meth.array

- now fully support for EPIC arrays

[methylKit](https://bioconductor.org/packages/methylKit)
---------

Changes in version 1.1.8:

IMPROVEMENTS AND BUG FIXES

- fix methSeg error when only one segment is returned from fastseg, add
  case handling for methSeg2bed

- check for user interruption in methCall to enable stop in execution

- changes to select() function: check for out-of-bound indices to
  prevent downstream errors

Changes in version 1.1.7:

IMPROVEMENTS AND BUG FIXES

- fix methCall segementation fault, added tests and test files from
  bismark

- fixed missing p.value at coercion of methylBase to GRanges

- changes to the pool() function: save.db=TRUE for methylBaseDB by
  default, differing lengths of given sample.ids and unique treatment
  lead to error, added tests

- fix osx related error when reading gzipped files with methRead

- change deprecated function names in test files

- fixed bug with dataSim() function updated the manual

- fix methSeg error when only one segment is returned from fastseg, add
  case handling for methSeg2bed

- check for user interruption in methCall to enable stop in execution

- changes to select() function: check for out-of-bound indices to
  prevent downstream errors

Changes in version 1.1.6:

IMPROVEMENTS AND BUG FIXES

- fixed a bug where tileMethylCounts() function did not work with small
  chromosomes/scaffolds with few bases covered.

Changes in version 1.1.5:

IMPROVEMENTS AND BUG FIXES

- fixes missing error messages during methRead :
  https://github.com/al2na/methylKit/pull/57

Changes in version 1.1.4:

IMPROVEMENTS AND BUG FIXES

- merging tabix files is fixed:
  https://github.com/al2na/methylKit/pull/56

Changes in version 1.1.3:

IMPROVEMENTS AND BUG FIXES

- Typos in the vignette are fixed, thanks to Marcin Kosinski

Changes in version 1.1.2:

IMPROVEMENTS AND BUG FIXES

- Bug fixes in SAM file reading process. If there were more than one
  header line there were problems in reading. Now this is fixed.
  https://github.com/al2na/methylKit/pull/51

Changes in version 1.1.1:

IMPROVEMENTS AND BUG FIXES

- Fisher's exact test now works as described in the manual. It is
  automatically applied in calculateDiffMeth() when there are only two
  groups with one replicate each.

- During logistic regression modeling, the samples without counts are
  removed from the model but the same filtering is not applied for
  covariates data.frame, which can cause errors if min.per.group
  argument is used.  Now this is fixed:
  https://github.com/al2na/methylKit/issues/50

[MIGSA](https://bioconductor.org/packages/MIGSA)
-----

Changes in version 0.99.0:

DOCUMENTATION

- `NEWS` file was added.

- First functional version

[minfi](https://bioconductor.org/packages/minfi)
-----

Changes in version 1.21:

- Moving RGChannelSet, MethylSet and RatioSet from building on eSet
  (from Biobase) to SummarizedExperiment (from SummarizedExperiment).
  Most important changes are that the constructor functions now uses
  the argument colData instead of pData; some of them have more
  arguments.  The updateObject methods have been extended to update to
  the new class backend.  While the pData, sampleNames, featureNames
  methods still work, we recommend (at least for package writers) to
  move to colData, colnames and rownames.

- Reverted the bugfix to preprocessQuantile mentioned under news for
  version 1.19. Our fix was wrong; the original code did not have a
  bug. Thanks to users who reported issues with the function (Frederic
  Fournier and David Martino).

- bugfix for getSnpBeta for subsetted (and combined) RGChannelSets
  (reported and diagnosed by Warren Cheung).

- Accessing the manifest or annotation now fails for an 'unknown'
  array.

- We now support gzipped IDAT files.

- Fixed a bug in read.metharray() which resulted in an error in some
  situations when running the function with argument force=TRUE to read
  IDAT files of different length.  Reported by Maria Calleja Cervantes
  <mcalleja@idibell.cat>.

[miRcomp](https://bioconductor.org/packages/miRcomp)
-------

Changes in version 1.5.1:

- Updated miRcomp Shiny app to final version prior to paper submission.

[miRNAmeConverter](https://bioconductor.org/packages/miRNAmeConverter)
----------------

Changes in version 1.3.1:

SIGNIFICANT USER-LEVEL CHANGES

- Updated citation

- Updated vignette

[MODA](https://bioconductor.org/packages/MODA)
----

Changes in version 1.1.2:

NEW FEATURES

- Multiple modules identification methods, avoid large modules

USER-LEVEL CHANGES

- Method name changes: WeightedModulePartitionHierarchical

Changes in version 1.0.1:

USER-LEVEL CHANGES

- R markdown vignette

BUG FIXES

- PartitionModularity

[monocle](https://bioconductor.org/packages/monocle)
-------

Changes in version 2.4.0:

- The default expressionFamily is now negbinomial.size, instead of
  Tobit. If you are using TPM or FPKM data, we urge you to convert it
  to relative transcript counts with relative2abs and use the negative
  binomial distribution in your CellDataSet objects.

- Revamped clusterCells functionality based on t-SNE and densityPeak

- New procedure for selected ordering genes called "dpFeature". See
  vignette for details.

[motifcounter](https://bioconductor.org/packages/motifcounter)
------------

Changes in version 0.99.0:

- Inital Bioconductor Submission

[motifStack](https://bioconductor.org/packages/motifStack)
----------

Changes in version 1.19.4:

NEW FEATURES

- update the css file for tooltip to function browseMotifs.

Changes in version 1.19.3:

NEW FEATURES

- add tooltip to function browseMotifs.

Changes in version 1.19.2:

NEW FEATURES

- add radialPhylog layout to function browseMotifs.

Changes in version 1.19.1:

NEW FEATURES

- add new function browseMotifs.

[msa](https://bioconductor.org/packages/msa)
---

Changes in version 1.7.2:

- fix for new clang 4 compiler on Mac OS

Changes in version 1.7.1:

- additional conversions implemented for msaConvert() function

- added a new method msaConsensusSequence() that extends the
  functionality provided by Biostring's consensusString() method

- added a new method msaConservationScore()

- print() method extended such that it now also allows for
  customization of the consensus sequence (via the new
  msaConsensusSequence() method)

- package now depends on Biostrings version >= 2.40.0 in order to make
  sure that consensusMatrix() also works correctly for masked
  alignments

- corresponding changes in documentation and vignette

Changes in version 1.7.0:

- new branch for Bioconductor 3.5 devel

[MSnbase](https://bioconductor.org/packages/MSnbase)
-------

Changes in version 2.1.18:

- suggest reshape2, as it's used in vignette <2017-04-18 Tue>

Changes in version 2.1.17:

- Update NEWS file <2017-04-11 Tue>

Changes in version 2.1.16:

- Remove timing test as it fails occasionally on the Bioconductor
  servers <2017-04-09 Sun>

Changes in version 2.1.15:

- Remove reshape2 dependency; see #201 <2017-04-06 Thu>

Changes in version 2.1.14:

- Internal rewrite and speedup of topN; Briefly multiple apply calls
  are avoided, `getTopIdx` and `subsetById` are replaced by `.topIdx`.
  See PR #199 for details. <2017-03-20 Mon>

- Fix mz calculation for terminal modifications and z > 1 in
  `calculateFragments`; closes #200 <2017-03-22 Wed>

- Fix errors and notes <2017-03-30 Thu>

Changes in version 2.1.13:

- Internal rewrite and speedup of plotNA <2017-02-26 Sun>

Changes in version 2.1.12:

- Import dist from stats <2017-02-25 Sat>

- Fix filing example <2017-02-25 Sat>

Changes in version 2.1.11:

- Fix breaks calculation for binning single (closes #191) and multiple
  (closes #190) spectra. The fix for single spectra (#191) could result
  in slightly different breaks on the upper end of the m/z values.
  <2017-02-10 Fri>

- New `aggvar` function, to assess aggregation variability <2017-02-11
  Sat>

Changes in version 2.1.10:

- New `diff.median` normalisation for `MSnSet`s. <2017-01-26 Thu>

- Fix combineFeatures message <2017-02-01 Wed>

Changes in version 2.1.9:

- When fully trimmed, an (empty) spectrum has peaksCount of 0L - see
  https://github.com/lgatto/MSnbase/issues/184 <2017-01-20 Fri>

- Add filterEmptySpectra,MSnExp method (see issue #181) <2017-01-20
  Fri>

- Add a section about notable on-disk and in-memory differences (was
  issue #165) <2017-01-20 Fri>

Changes in version 2.1.8:

- Remove order option altogether <2017-01-19 Thu> (superseeds setting
  default sorting using "auto" on R < 3.3 and "radix" otherwise
  <2017-01-03 Tue>)

Changes in version 2.1.7:

- Setting default sorting using "auto" on R < 3.3 and "radix" otherwise
  <2017-01-03 Tue>

- filterMz returns an empty spectrum when no data is within the mz
  range (see issue #181) <2017-01-16 Mon>

- Performance improvement: a new private .firstMsLevel will efficiently
  return the first MS level in an MSnExp and OnDiskMSnExp. See issue
  #183 for details/background <2017-01-18 Wed>

Changes in version 2.1.6:

- Migrate io and dev vignettes to BiocStyle's html_document2 style
  <2016-12-23 Fri>

- Update show method to display class.

- Migrated to NEWS.md <2016-12-23 Fri>

- Update DESCRIPTION (and README) to reflect wider usage of MSnbase
  (replaced MS-based proteomics by mass spectrometry and proteomics)
  <2016-12-23 Fri>

Changes in version 2.1.5:

- Fix (unexported) navMS example code <2016-12-14 Wed>

Changes in version 2.1.4:

- Jo added netCDF support <2016-11-30 Wed>

Changes in version 2.1.3:

- FeaturesOfInterest collections can now be assigned names - addresses
  issue #172 <2016-11-25 Fri>

Changes in version 2.1.2:

- Update readMSnSet2 to save filename <2016-11-09 Wed>

- Ensure that header information is read too if spectra data is loaded
  for onDiskMSnExp objects (see issue #170) <2016-11-24 Thu>

Changes in version 2.1.1:

- Fix typo in impute man page <2016-10-19 Wed>

- Cite Lazar 2016 in vignette imputation section <2016-10-28 Fri>

Changes in version 2.1.0:

- New version for Bioconductor devel

- New version for Bioconductor release version 3.4

[msPurity](https://bioconductor.org/packages/msPurity)
--------

Changes in version 1.1.1:

- Added pcalc functions to be used by user

- Added option to remove isotopes from calculation

[MWASTools](https://bioconductor.org/packages/MWASTools)
---------

Changes in version 1.0.0:

SIGNIFICANT USER-VISIBLE CHANGES

- Package introduced.

NEW FEATURES

- Package introduced.

[mzR](https://bioconductor.org/packages/mzR)
---

Changes in version 2.9.11:

- Restore -fpermissive flag on windows

Changes in version 2.9.10:

- Remove register keyword causing WARNING.

Changes in version 2.9.9:

- Remove C++ references to cout/cerr/abort in pwiz code (see issue #89)

Changes in version 2.9.8:

- Fix reading spectrum polarity from mzML using pwiz backend, closes
  #81

Changes in version 2.9.7:

- Fix compilation on macOS

Changes in version 2.9.6:

- Compile on macOS, but hdf5 path hard-coded

Changes in version 2.9.5:

- Add missing boost/config/platform/macos.hpp <2017-01-25 Wed>

Changes in version 2.9.4:

- New chromatogram accessors (for pwiz backend only) - see issue #73
  <2017-01-23 Mon>

Changes in version 2.9.3:

- bump to new Rcpp 0.12.8 version <2017-01-05 Thu>

Changes in version 2.9.2:

- cleanup CFLAGS and LIBS for libnetcdf

- add file missing for oaxaca (Apple clang 3.5svn / 600.0.57)

Changes in version 2.9.1:

- Delete RAMPAdapter pointer in pwiz backend (by jotsetung) <2016-11-20
  Sun>

- Use spectra in addition to peaks (see issue #15) <2016-12-09 Fri>

- New pwiz (commit 946d23d75dc70a7a4913d8e05e3d59b9255f278e)

Changes in version 2.9.0:

- Bioc devel 3.5

[NanoStringQCPro](https://bioconductor.org/packages/NanoStringQCPro)
---------------

Changes in version 1.7.1 (2017-04-10):

- Documentation improvements to address user-reported issues.

[netReg](https://bioconductor.org/packages/netReg)
------

Changes in version 1.0.0:

- First bioc release

Changes in version 0.99.0:

- Devel version 0.99.0

- `edgenet` penalization using CCD

- Implementation of model selection using Dlib

[omicade4](https://bioconductor.org/packages/omicade4)
--------

Changes in version 1.15.1:

- plotVar function updated

[OncoScore](https://bioconductor.org/packages/OncoScore)
---------

Changes in version 1.3.3 (2017-03-27):

- More bug fixes to the queries to PubMed.

- Implementing a function for combining genes.

- Support HTTPS.

[OncoSimulR](https://bioconductor.org/packages/OncoSimulR)
----------

Changes in version 2.6.0:

- Many additions to the vignette and documentation.

- LOD and POM (lines of descent, path of maximum, sensu Szendro et
  al.).

- Diversity of sampled genotypes.

- Genotyping error can be added in samplePop.

- fixation of a genotype/gene as stopping mechanism.

- rfitness: shifting by subtraction and mu of normal distribution.

- simOGraph: using proper transitive reduction.

- simOGraph can also output rT data frames.

- accessible genotypes now done in C++.

- Handling of trivial cases in genotFitness.

- Clarified McFarland parameterization.

- Better (and better explained) estimates of simulation error for McFL.

- AND of detectedSizeP and lastMaxDr.

- sampledGenotypes in user code.

- clonePhylog et al: deal with never any descendant.

- samplePop can handle failed simulations graciously.

- summary.oncosimulpop can handle failed simulations graciously.

- Citation shows Bioinformatics paper.

Changes in version 2.5.14 (2017-04-07):

- Fixed repeated entries in NEWS for BioC 3.5.

Changes in version 2.5.13 (2017-04-07):

- Updated NEWS for BioC 3.5.

Changes in version 2.5.12 (2017-02-18):

- rfitness: allow simple forcing of wt to 1, shifting by subtraction,
  and specifying mu of normal distribution.

- simOGraph: proper trm comparison.

- Citation now shows Bioinformatics reference.

Changes in version 2.5.11 (2017-01-27):

- Transitive reduction: must call transitive.closure first.

Changes in version 2.5.10 (2017-01-27):

- Transitive reduction: calling nem in simOGraph

Changes in version 2.5.9 (2017-01-09):

- Added code coverage comments to vignette.

Changes in version 2.5.8 (2016-12-17):

- Handle trivial cases in genotFitness.

Changes in version 2.5.7 (2016-12-15):

- Clarified McFarland parameterization.

Changes in version 2.5.6 (2016-12-14):

- Fixed a few typos in help files.

Changes in version 2.5.5 (2016-12-14):

- Vignette: miscell changes (typos, etc)

Changes in version 2.5.4 (2016-12-12):

- Vignette: miscell changes (order of examples, typos, etc)

Changes in version 2.5.3 (2016-12-12):

- Vignette uses pander in tables.

- Typos fixed and other enhancements in vignette.

Changes in version 2.5.2 (2016-12-10):

- Lots and lots of addition to vignette including benchmarks.

- Diversity of sampled genotypes.

- Genotyping error can be added in samplePop.

- LOD and POM (lines of descent, path of maximum, sensu Szendro et
  al.).

- simOGraph can also out rT data frames.

- Better (and better explained) estimates of simulation error for McFL.

Changes in version 2.5.1 (2016-11-12):

- AND of detectedSizeP and lastMaxDr.

- fixation as stopping mechanism.

- sampledGenotypes in user code.

- clonePhylog et al: deal with never any descendant.

- samplePop can handle failed simulations graciously.

- summary.oncosimulpop can handle failed simulations graciously.

- accessible genotypes now done in C++.

- OcurringDrivers should not be a factor.

- samplePop always returns gene names.

- to_Magellan is much faster with rfitness objects.

- Several improvements in vignette (English and additional
  explanations).

[Organism.dplyr](https://bioconductor.org/packages/Organism.dplyr)
--------------

Changes in version 1.0.0:

NEW FEATURES

- This package provides an interface to combined _Bioconductor_ org.*
  (identifier) and TxDb.* (genomic coordinate) annotation resources.
  The interface is implemented at several levels, including low-level
  'dplyr', org-like select(), and TxDb-like genes(), etc.

SIGNIFICANT USER-VISIBLE CHANGES

- Filters use strict CamelCase convention.

BUG FIXES

- *IdFilter and *RankFilter are numeric (integer), rather than
  character.

[OrganismDbi](https://bioconductor.org/packages/OrganismDbi)
-----------

Changes in version 1.18.0:

BUG FIXES

- avoid duplicate factor levels during compression of metadata for
  cdsBy and friends; previously introduced incorrectly empty metadata

[pathview](https://bioconductor.org/packages/pathview)
--------

Changes in version 1.15.1:

- fixed bug in node.map on single row/mapped data, and in mapped row
  numbers.

[Pbase](https://bioconductor.org/packages/Pbase)
-----

Changes in version 0.15.3:

- Working up new API (see issue #41) <2017-04-09 Sun>

- Temporarily remove some vignettes, until the API has stabilised.
  <2017-04-09 Sun>

Changes in version 0.15.2:

- New Proteins class implementation - see issue #38. <2016-12-10 Sat>

- Fixed many errors from current rewrite <2017-04-08 Sat>

Changes in version 0.15.1:

- Proteins,EnsDb method allowing to fetch a Proteins object from an
  EnsDb database.

- Fixes in the mapToGenome method: - Works also for negative strand
  encoded proteins (issue #29). - Supports a GRangesList object with
  arbitrary mcols. - Faster implementation (issue #30).

- mapToGenome,Proteins,EnsDb and pmapToGenome,Proteins,EnsDb methods
  that map peptide features to the genome using annotations fetched
  from an EnsDb.

- The seqnames of the Proteins object are used as names for the
  resulting GRangesList object from the mapToGenome and pmapToGenome
  methods (issue #34).

- Drop unique `seqnames` requirement; see #28, #32

- Create `names` as synonym for `seqnames`; close #32

[pbcmc](https://bioconductor.org/packages/pbcmc)
-----

Changes in version 1.3.2:

CODE

- *`test_subjectReport` update after subjectReport modification.

DEPENDENCIES

- *R (>= 3.4) was updated.

Changes in version 1.3.1:

CODE

- *`subjectReport` bug removed over axis.text.x duplicated parameter.

Changes in version 1.3.0:

VERSION

- Bump version after creating 3.5 devel branch

[pcaExplorer](https://bioconductor.org/packages/pcaExplorer)
-----------

Changes in version 2.2.0:

NEW FEATURES

- Added Demo data, loadable via demo button

BUG FIXES

- Plots work now without cutting out points when zooming in

OTHER NOTES

- Saved reactive values are now exported to dedicate environments
  (instead of assigning to global)

[PGA](https://bioconductor.org/packages/PGA)
---

Changes in version 1.5.1:

- Fix the bug in Outputaberrant2

[philr](https://bioconductor.org/packages/philr)
-----

Changes in version 1.1.0:

USER-VISIBLE CHANGES

- Inverse ILR ilrpInv is now implemented as well as the inverse clrp
  transform (clrpInv) and the inverse (shiftpInv) function.

- In order to untransformed more generally transformed PhILR data
  (e.g., with branch length weights), a philrInv function has been
  created, this is likely the most userfriendly way to invert any
  transformed data.

- Updated documentation

- Various Bugfixes

- Added updated citation information for package

[Pi](https://bioconductor.org/packages/Pi)
--

Changes in version 1.3.3:

NEW FEATURES

- Add a new function ('xPierROCR') for assessing the dTarget
  performance via ROC and Precision-Recall (PR) analysis

Changes in version 1.3.2:

NEW FEATURES

- Define S3 classes ('pNode', 'eTarget', 'dTarget', 'sTarget',
  'cTarget', etc)

[piano](https://bioconductor.org/packages/piano)
-----

Changes in version 1.14.5:

BUG FIXES

- Reset plot layout after networkPlot.

- Issue warning instead of error in writeFilesForKiwi when gene-level
  statistics are not p-values. This means that the GLS file will not be
  generated, but the GSS and GSC files are.

Changes in version 1.14.3:

BUG FIXES

- Required functions from the snow package, used by snowfall, are now
  properly loaded.

[Pigengene](https://bioconductor.org/packages/Pigengene)
---------

Changes in version 1.1.10 (2017-03-30):

Changes in existing functions

- Checking the pigengene input of module.heatmap().

- Issues in the balance() function (not exported) where Labels is a
  factor are resolved. Also, if all sampls have the same size,
  oversampling is automatically turned off.

- If Labels is a factor, it is now converted to a character vector in
  check. pigengene.input().

Changes in version 1.1.6 (2017-03-27):

Changes in existing functions

- The module.heatmap() function now has the doAddEigengene and
  scalePngs arguments.

- The compute.pigengene() function now reports also the size of modules
  in the pigengene_pvalue.csv output file.

[polyester](https://bioconductor.org/packages/polyester)
---------

Version: 1.99.3
Text: NB function now exported

Version: 1.99.3
Text: note that version 1.99.3 on GitHub was version 1.1.0 on
        Bioconductor.

Version: 1.99.2
Text: bug fix in fragment generation (last 2 bases of transcript were
        never sequenced)

Version: 1.99.1
Text:

[pqsfinder](https://bioconductor.org/packages/pqsfinder)
---------

Changes in version 1.3:

NEW FEATURES

- Added estimated time to complete (ETTC) to the progress line.

- New default values for scoring parameters as a result from training
  on G4-seq experimental data and co-testing on a set of known
  quadruplexes from literature.

- New PQS metadata reported: number of tetrads, bulges, mismatches and
  loop lengths. To get this data, use elementMetadata accessor
  function.

- New algorithm option: reporting of all overlapping PQS.

- Score distribution is reported. For each sequence position you get
  maximal score of PQS that was found overlapping the position. For
  details, see scoreDistribution function.

- Novel scoring parameter: exponent of loop length mean in the scoring
  equation to express non-linear dependency of the PQS propensity to
  the loop lengths.

BUG FIXES

- Minimal loop length is set back to 0 by default, but only one
  zero-length loop is allowed by the scoring system.

- Fixed unintentional cast of loop length mean factor and loop standard
  deviation factors from float to integer.

[pRoloc](https://bioconductor.org/packages/pRoloc)
------

Changes in version 1.15.9:

- update biomart attribute names to relect changes <2017-04-20 Thu>

Changes in version 1.15.8:

- New mrkConsProfiles function to calculate average/consensus marker
  profiles <2017-04-11 Tue>

Changes in version 1.15.7:

- Fix warnings and notes <2017-02-25 Sat>

- Add section about dimensionality methods reduction and t-SNE in the
  tutorial <2017-03-07 Tue>

- Fix error due to new uniprot attribute names <2017-04-06 Thu>

Changes in version 1.15.6:

- fix (unexported) remap function - see issue #92 <2016-12-15 Thu>

- plot2Ds now only adds segments when the featureNames are identical
  <2017-01-11 Wed>

- Import (rather than Suggest) hexbin <2017-01-18 Wed>

- Increase margin in QSep plotting (contributed by S. Gibb) <2017-02-07
  Tue>

Changes in version 1.15.5:

- Update plotDist to use sampleNames() to label x axis ticks and
  getUnknowncol() as default pcol (see issue #91) <2016-11-08 Tue>

- xlab and ylab are args in plotDist <2016-12-04 Sun>

Changes in version 1.15.4:

- Update human markers - see pRolocdata's issue 21
  https://github.com/lgatto/pRolocdata/issues/21 <2016-11-08 Tue>

Changes in version 1.15.3:

- Fix bug in plot2D to ignore fcol when using hexbin method <2016-11-04
  Fri>

Changes in version 1.15.2:

- Fix Arabidopsis parameters for biomaRt questions: 'TAIR locus ID' is
  now 'Stable gene ID'. Changed in several manual files and
  dunkley2006params. <2016-11-02 Wed>

Changes in version 1.15.1:

- new plot3D function <2016-10-27 Thu>

- update CITATION <2016-10-27 Thu>

- use predict:::predict.plsa, as it is not exported anymore <2016-11-01
  Tue>

[pRolocGUI](https://bioconductor.org/packages/pRolocGUI)
---------

Changes in version 1.9.4:

- fixed remap=FALSE bug in compare app <2017-01-12 Thu>

- Added mirrorX and mirrorY to the compare app <2017-01-12 Thu>

Changes in version 1.9.3:

- Update vignette use latest BiocStyle::html_document2() with floating
  table-f content <2016-12-30 Fri>

- Change to NEWS.md <2016-12-30 Fri>

- mirrorX and mirrorY are now ignored in pRolocVis - see issue #84
  <2017-01-11 Wed>

Changes in version 1.9.2:

- Remove accidental merging left-over <2016-12-13 Tue>

Changes in version 1.9.1:

- Remove accidental call to browser

Changes in version 1.9.0:

- Bioc devel 3.5

[Prostar](https://bioconductor.org/packages/Prostar)
-------

Changes in version 1.7.21:

NEW FEATURES

- Pre-release version

Changes in version 1.7.19:

NEW FEATURES

- Code restructured

Changes in version 1.7.17:

NEW FEATURES

- Interactive volcanoplot

Changes in version 1.7.13:

NEW FEATURES

- The normalization function has been modified to take into account the
  quantile normalization. Only works with DAPAR version >= 1.7.13. Back
  compatibilty with previous versions

Changes in version 1.7.11:

NEW FEATURES

- In the Desctiptive statistics panel, The NA values are colored in the
  table of quantitative data

Changes in version 1.7.9:

NEW FEATURES

- In the Desctiptive statistics panel, the variance distribution plot
  has been replaced by a CV distribution plot

Changes in version 1.7.7:

NEW FEATURES

- In the missing values imputation tool, two features have been
  implemented : - a serie of options for the imp4p method - a help text
  that describes the method selected by the user

Changes in version 1.7.5:

BUG FIXES

- The legend of the x-axis for the boxplots has been modified. With a
  lot of samples, there was a problem with the margins which were too
  large to display the plots. Now, the legend appears in one line
  instead of one line per type of information

Changes in version 1.7.3:

NEW FEATURES

- The agregation tool now deals with a dataset with missing values

[psichomics](https://bioconductor.org/packages/psichomics)
----------

Changes in version 1.1.10:

- Gene, protein and transcript information:

- Fix tooltip text presentation in transcript plot

- Fix JavaScript issues when zooming the transcript plot

- Fix error when plotting events associated with multiple genes

- Fix error when plotting single-exon transcripts

- Protein name, length and function are now presented when available

- Improved general presentation of the information

- Differential splicing analyses:

- Click and drag in the plot to zoom in and subsequently filter events
  shown in the table

- Decreased step of sliders

- Improve interface of previewed survival curves

- When clicking on a table link to navigate to differential splicing
  analyses of a single event, the appropriate analyses will now be
  automatically rendered with the respective options, as expected

- Settings (renamed to "Help"):

- Add links to tutorials and user feedback

- Add app information and acknowledgments

- Remove unused option for choosing cores (all performed operations are
  still single-core, given the difficulty of working with
  multiprocesses in Shiny)

- Improve dialogs regarding missing data and other minor interface
  elements

- Update documentation with volcano plot

Changes in version 1.1.9:

- Differential splicing analyses:

- Add volcano plot to represent events through selected attributes,
  such as p-values and descriptive statistics (e.g. median and
  variance) between groups of interest

- Transform values of the X and Y axis in the plot using log
  transformed, inverted and absolute values, for instance

- Highlight events in the plot based on values of the X and Y axis

- Table of differential analyses per alternative splicing event is
  filtered according to highlighted and selected events in the plot

- Gene, protein and transcript information:

- Transcript plot is now interactive and zoomable

- Protein are now rendered based on selected transcript alone

- Faster parsing of Uniprot's web API response

- Improve display of article information when data is missing

- Principal component analysis:

- Improve presentation of available options

- When clicking on previews of differential splicing and survival
  analyses, the appropriate analyses will now be automatically rendered
  with the respective options

- Fix buggy browser history when the user is directed to a different
  tab

- Consistently use Firebrowse and Firehose across the package

- Update documentation

Changes in version 1.0.8:

- Support GTEx data loading and analysis

- Fix clinical data dependency: - Fix error when trying to load a file
  containing alternative splicing quantification without first loading
  clinical data - Fix error where samples from junction quantification
  were matched to clinical information even if clinical data were not
  loaded - Inform user when clinical data is not loaded while trying to
  plot survival curves

- Improve data grouping: - Create sample groups like patient groups and
  perform set operations between any created groups - Create groups
  using patient and sample identifiers - Check number of patients and
  samples per group - Rename selected groups - Alert user when groups
  cannot be created due to missing data

- Differential splicing analysis: - Analyse all samples as one group

- Survival analysis: - Select any clinical attribute for
  starting/follow up and ending times

- Create table containing TCGA sample metadata when calculating or
  loading alternative splicing quantification

- Minor UI improvements

Changes in version 1.0.7:

- Survival analysis: - Fix error caused by some non-matched patients
  not being in the patient-sample matching matrix

Changes in version 1.0.6:

- Update tutorials with more relevant and complex examples

- Update minimum versions required of highcharter (0.5.0) and shiny
  (1.0.0): - Fix function usage as according to new version of
  highcharter - More options available when exporting plots (PNG, JPEG,
  SVG, XLS and CSV)

- Faster alternative splicing quantification

- Differential splicing analysis: - Fix major bug where samples could
  be placed in the wrong groups - Shorten speed of the calculation for
  the optimal PSI cut-off that minimises the survival difference - Fix
  not performing statistical tests for two selected sample types while
  analysing a single event with three or more sample types - Fix
  differential analysis on one splicing event not working when using
  `diffAnalyses()` function - Fix differential analysis not showing for
  individual events before navigating to the page where the analysis is
  performed for all events - Improve readability and information of
  statistical tests for single events

- Principal component analysis: - Shorten time taken to calculate
  principal components and to render the loadings plot - Fix loadings
  plot error when rendering some principal components

- Survival analysis: - Fix incorrect number of patients from the
  survival groups in the contextual information for the selected
  cut-off (below the slider) - Improve how alternative splicing
  quantification is assigned to patients based on their samples

- Protein annotation: - Warn user when trying to render proteins with
  no domains

Changes in version 1.0.5:

- Navigate history using the browser forward and back buttons

- Fix delay when displaying large data by removing columns containing
  missing values exclusively

- Principal component analysis: - Improve speed when calculating total
  contribution of each variable to the principal components

- Survival analysis: - Shorten calculation of optimal PSI that
  minimises the survival difference - Improve visual cues of optimal
  PSI cut-off and present p-value of selected PSI cut-off - Fix
  ambiguous error messages - Fix incorrect Cox model results for
  formula-based calculations - Fix null Cox models crashing the program

- Differential splicing analysis: - Select sample types for
  differential splicing analysis - Fix statistical tests not displaying
  for individual events after differentially analysing all events using
  the other statistical tests

Changes in version 1.0.4:

- Correctly load files and quantify alternative splicing for PRAD, OV
  and PAAD tumour types from The Cancer Genome Atlas (TCGA)

- Fix session disconnecting when exporting plots in Firefox

- Improve text and behaviour of fields to select datasets and AS events

- Fix author names and add contributor

Changes in version 1.0.3:

- Bug fixes regarding gene annotation: - Fix disabled gene selection
  when choosing a splicing event associated with a single gene after
  selecting an event related to multiple genes - Fix display of PubMed
  articles related to previously selected gene when selecting a
  single-gene-associated event after selecting an event related to
  multiple genes

- Bug fixes regarding groups: - Fix groups by rows not working - Fix
  group selection not working when only one group exists - Improve
  argument name of getGroupsFrom()

- Other minor improvements

Changes in version 1.0.2:

- Fix UTF-8 encoding in author list

Changes in version 1.0.1:

- Improve metadata (title, description, authors and vignette titles)

[PureCN](https://bioconductor.org/packages/PureCN)
------

Changes in version 1.6.0:

- Lots of improvements to command line scripts

- Improved somatic vs. germline status calling

- Better mapping bias estimation and correction

- Better integration into existing copy number pipelines

- Support for cell lines

- New GC-normalization for smaller gene panels

- Added sub-clonal SNV state (SOMATIC.M0)

- Polished plots, added new GC-normalization and volcano plots

- Better copy number normalization using multiple best normals

- Removed automatic curation, since the tuned likelihood model of
  runAbsoluteCN was hard to beat

- More control over homozygous deletions (significant portion of wrong
  maximum likelihood solutions had many homozygous deletions)

- Faster post.optimize=TRUE by not optimizing poor fits or unlikely
  solutions

- Automatic 50bp interval padding

- Tweaks to segmentationPSCBS

- seg.file can contain multiple samples

- Contamination rate estimation (experimental)

- Code cleanups (switch from inlinedocs to roxygen, from message/warn
  to futile.logger) API CHANGES

- runAbsoluteCN output from PureCN 1.2 cannot be analyzed with PureCN
  1.6 and needs to be re-run. We hope to avoid this in the future.

- Renamed functions: readCoverageGatk to readCoverageFile since future
  versions will likely support additional third-party tools.

- Deprecated functions: createSNPBlacklist, getDiploid,
  autoCurateResults

- Defunct functions: createExonWeightFile

- Changed defaults:

- min.normals 4 (from 10) in setMappingBiasVcf

- max.segments 300 (from 200) in runAbsoluteCN

- min.targeted.base 5 (from 4) in filterTargets

- max.homozygous.loss now a double(2) vector, with first element
  specifying the maximum fraction of genome deleted (default 5%) and
  the second value the maximum size of a homozygous loss (default
  10mb).

- prior somatic for variants in both dbSNP and COSMIC changed from 0.01
  and requiring 3 hits to 0.5 and requiring 4 hits.

- Other minor changes:

- Renamed some predictSomatic() output column names

- Removed "beta.model" from "SNV.posterior" slot since model is now an
  option

- Moved remove.off.target.snvs to filterVcfBasic

- Moved normalDB from filterTargets to runAbsoluteCN, since it is now
  used for more than target filtering

- Dropped BED file support in calculateGCContentByInterval Instead
  provide support for GRanges

- poolCoverage: w argument now used as provided, not normalized so that
  w[1] is 1

- Removed ... from runAbsoluteCN

- min.coverage removed from segmentation function, since this is now
  done by filterTargets

- Added centromeres to segmentation function

- Replaced contamination.cutoff with contamination.range in
  filterVcfBasic

- Removed verbose from most functions, since messages are now
  controlled with futile.logger

- Smoothing of log-ratios before segmentation now optionally done by
  runAbsoluteCN, not segmentation function

- setMappingBiasVcf now returns a list with elements bias (the old
  return value) and pon.count, the number of hits in the PON PLANNED
  FEATURES FOR 1.8

- Better sample summary statistics, like mutation burden, chromosomal
  instability

- Better performance in low purity samples

- Better performance in high purity samples with significant
  heterogeneity

- LOH database

- Switch to S4 data structures (maybe)

- Whole dataset visualizations (maybe)

- Better support for known, small deletions and amplifications (e.g.
  EGFRvIII, MYC)

- Support for GATK4

- Better runtime performance by ignoring unlikely solutions early

[pwOmics](https://bioconductor.org/packages/pwOmics)
-------

Changes in version 1.7.1:

- include signaling axes identification functions

- include phosphorylation information prefiltering in intersection
  analysis

- include direction of regulation prefiltering in intersection analysis

- include visualization of temporal correlations between phosphosite
  expression data and transcriptome data

[QDNAseq](https://bioconductor.org/packages/QDNAseq)
-------

Changes in version 1.12.0:

RELEASE

- Bioconductor 3.5

IMPROVEMENTS

- VCF and SEG file export have been implemented to allow use of
  downstream analysis tools such as Cartegenia (NGS) Bench.

- binReadCounts() now supports parallel computing

- calculateBlackListByRegions() has been implemented for convient bin
  overlap calculation of any set of regions.

[qpgraph](https://bioconductor.org/packages/qpgraph)
-------

Changes in version 2.10:

BUG FIXES

- Bugfix in the calculation of the p-values of qpPCC() when missing
  observations are present in the input data.

[QUBIC](https://bioconductor.org/packages/QUBIC)
-----

Changes in version 1.3.2:

- OpenMP with Macintosh operating systems is enabled

[R3CPET](https://bioconductor.org/packages/R3CPET)
------

Version: 1.8.0
Text: Updates: * Fixed some import issues * Fixed a bug in the
        visualizeCircos function * updated the documentation

[ramwas](https://bioconductor.org/packages/ramwas)
------

Changes in version 1.0.0:

NEW FEATURES

- Added joint methylation-genotype analysis

BUG FIXES

- package under active development

[RCyjs](https://bioconductor.org/packages/RCyjs)
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

[ReactomePA](https://bioconductor.org/packages/ReactomePA)
----------

Changes in version 1.19.1:

- update startup message <2016-11-09, Wed>

[recount](https://bioconductor.org/packages/recount)
-------

Changes in version 1.1.27:

NEW FEATURES

- Added the add_predictions() function which appends the predicted
  phenotypes to a RSE object downloaded with recount. The phenotypes
  were predicted by Shannon Ellis et al, 2017 (citation coming up
  soon!).

Changes in version 1.1.26:

SIGNIFICANT USER-VISIBLE CHANGES

- Changed the citation now that the recount2 paper has been published
  at http://www.nature.com/nbt/journal/v35/n4/full/nbt.3838.html.

Changes in version 1.1.25:

NEW FEATURES

- Added the function getRPKM() which can be used with
  RangedSummarizedExperiment objects from recount and from other
  sources.

Changes in version 1.1.24:

SIGNIFICANT USER-VISIBLE CHANGES

- recount_url now includes the URLs for the GTEx bigWig files.

Changes in version 1.1.19:

SIGNIFICANT USER-VISIBLE CHANGES

- coverage_matrix() now returns a RangedSummarizedExperiment object.
  This matches the behavior of recount.bwtool::coverage_matrix_bwtool()
  and is more consistent with the use of RSE objects in recount.

Changes in version 1.1.18:

BUG FIXES

- coverage_matrix()'s helper function .read_pheno() was failing for
  some projects.

Changes in version 1.1.16:

BUG FIXES

- Fixed a bug in the counts in coverage_matrix(). They were being
  incorrectly multiplied by 100.

Changes in version 1.1.14:

SIGNIFICANT USER-VISIBLE CHANGES

- Completed the change to Gencode v25 annotation for exon and gene
  counts.

Changes in version 1.1.13:

SIGNIFICANT USER-VISIBLE CHANGES

- We dropped TxDb.Hsapiens.UCSC.hg38.knownGene completely from recount
  and will be using Gencode v25 instead.

Changes in version 1.1.12:

BUG FIXES

- Updated snaptron_query() to comply with recent changes in Snaptron.

Changes in version 1.1.8:

SIGNIFICANT USER-VISIBLE CHANGES

- Updated the package so you can now access TCGA data. Now there's over
  8 terabytes of data available in the recount project!

Changes in version 1.1.6:

SIGNIFICANT USER-VISIBLE CHANGES

- snaptron_query() can now access GTEx and TCGA data.

Changes in version 1.1.5:

SIGNIFICANT USER-VISIBLE CHANGES

- Snaptron changed from stingray.cs.jhu.edu:8090 to snaptron.cs.jhu.edu
  so snaptron_query() has been changed accordingly.

Changes in version 1.1.2:

SIGNIFICANT USER-VISIBLE CHANGES

- The function reproduce_ranges() now has the 'db' argument. By default
  it's set to TxDb.Hsapiens.UCSC.hg38.knownGene to reproduce the actual
  information used in recount. But it can also be used with
  EnsDb.Hsapiens.v79 to use the ENSEMBL annotation. Then with
  coverage_matrix() you can get the counts for either an updated
  TxDb.Hsapiens.UCSC.hg38.knownGene or for EnsDb.Hsapiens.v79 at the
  exon and/or gene levels as shown in the vignette.

Changes in version 1.1.1:

SIGNIFICANT USER-VISIBLE CHANGES

- The vignette now describes how to download all the data, how to check
  exon-exon junctions by class, and how to use SciServer compute to
  access all the recount data (over 6 TB) via http://www.sciserver.org/

[RedeR](https://bioconductor.org/packages/RedeR)
-----

Changes in version 1.24.0:

- Implemented a new call-back method using R core infrastructure.

[regionReport](https://bioconductor.org/packages/regionReport)
------------

Changes in version 1.9.1:

SIGNIFICANT USER-VISIBLE CHANGES

- Changed the default style to BiocStyle::html_document2.

[ReportingTools](https://bioconductor.org/packages/ReportingTools)
--------------

Changes in version 2015-3-27:

- Updated email for Jessica L. Larson

Changes in version 2013-4-1:

- Minor bug fixes in NAMESPACE and DESCRIPTION and css

- Enhanced vignettes to include information on our website and
  publication

- Enhanced vignettes to clarify use of .modifyDF()

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

[rGREAT](https://bioconductor.org/packages/rGREAT)
------

Changes in version 1.7.1:

- submitGreatJob(): remove additional column by pintersect()

[rgsepd](https://bioconductor.org/packages/rgsepd)
------

Changes in version 1.7.1:

BUG FIXES

- sampleMeta may load as factors, and uses strings later. To prevent
  possible sample-swapping, need to make all inputs as
  strings/notfactors. We can force this internally.

[rhdf5](https://bioconductor.org/packages/rhdf5)
-----

Changes in version 2.20.0:

NEW FEATURES

- Indexing into spaces with more than .Machine$integer.max elements is
  supported using numeric (rather than integer) indexing; this provides
  exact indexing into spaces with about 51 bits of precision.

- Zero-length indexing is now supported (returning zero-length slabs).

BUG FIXES

- Using bit64conversion = "double" would always warn about loss of
  precision, but now only warns when precision is actually lost.

[RiboProfiling](https://bioconductor.org/packages/RiboProfiling)
-------------

Changes in version 1.5.3:

- Corrected a bug in utils.R, readStartCov1Aln function.
  normRange(listReadStartCov, fixedInterval), in case ixReverse null

Changes in version 1.5.2:

- Down-sampled the ctrlGAlignments.rda object from the data folder

- Changed conflicting sweave-knitr in the vignette

Changes in version 1.5.1:

- For the vignette the example bam files have been down-sampled and
  added to the extdata folder

- Changed the vignette accordingly

[RnBeads](https://bioconductor.org/packages/RnBeads)
-------

Changes in version 1.7.5:

- Several minor bugfixes and performance improvements

- added a vignette section on working with RnBSet objects

Changes in version 1.7.4:

- Several minor bugfixes and performance improvements

- Vignette installation instructions updated

- Reduce warnings in R CMD check

Changes in version 1.7.3:

- Age predictor (MethylAger) updates and documentation

- Support for external tools bedToBigBed and bedGraphToBigWig

- Minor bugfixes

Changes in version 1.7.2:

- Added genetic purity estimation based on SNP probes (option
  qc.snp.purity, microarrays only)

Changes in version 1.7.1:

- Added support for the ENmix.oob background subtraction method

- Several improvements in age prediction module

- Added conversion from minfi raw dataset to RnBeadRawSet

- Several minor bug fixes

[rols](https://bioconductor.org/packages/rols)
----

Changes in version 2.3.4:

- Ammend failing Ontologies unit test (related to changes made in
  version 2.3.3: the GO name reverted back to Gene Ontolgy) <2017-01-11
  Wed>

Changes in version 2.3.3:

- Add ctb to Authors@R <2016-12-28 Wed>

- use NEWS.md <2016-12-28 Wed>

- Ammend failing Ontologies unit test <2017-01-02 Mon>

Changes in version 2.3.2:

- Update test to reflect GO's new title <2016-12-21 Wed>

Changes in version 2.3.1:

- Fix failing unit test <2016-11-22 Tue>

Changes in version 2.3.0:

- Bioconductor devel 3.5

[ropls](https://bioconductor.org/packages/ropls)
-----

Changes in version 1.7.2:

INTERNAL MODIFICATION

- vignette now in pdf format

[rpx](https://bioconductor.org/packages/rpx)
---

Changes in version 1.11.2:

- Update unit test to reflect upstream changes <2017-01-25 Wed>

Changes in version 1.11.1:

- Migrate vignette to BiocStyle::html_document2 <2016-12-22 Thu>

- Use NEWS.md <2016-12-28 Wed>

Changes in version 1.11.0:

- Bioconductor devel 3.5

[Rqc](https://bioconductor.org/packages/Rqc)
---

Changes in version 1.10:

BUG FIXES

- The rqcQA function always returns list

[Rsamtools](https://bioconductor.org/packages/Rsamtools)
---------

Changes in version 1.27:

BUG FIXES

- qnameSuffixStart<-(), qnamePrefixEnd<-() accept 'NA' (bug report from
  Peter Hickey).

- scanBam() accepts a single tag mixing 'Z' and 'A' format. See
  https://support.bioconductor.org/p/94553/

[Rsubread](https://bioconductor.org/packages/Rsubread)
--------

Changes in version 1.26.0:

NEW FEATURES

- Gene annotation can be provided to align() and subjunc() to improve
  exon junction detection.

- Improve sanity checking for input and output data for align(),
  subjunc() and featureCounts().

- Resolve inconsistency between runs for align() and subjunc() when
  more than one CPU thread is used.

[RTCA](https://bioconductor.org/packages/RTCA)
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

[RTN](https://bioconductor.org/packages/RTN)
---

Changes in version 1.14.0:

- Concluded the VSE/EVSE pipeline, including examples.

- Improved computational performance and RTN workflows.

- In order to improve stability and portability, all dependencies have
  been extensively revised and, when available, replaced with more
  stable options.

- In order to simplify documentation and usability, some pipelines have
  been revised and should be distributed as separated packages, or 'on
  demand', focused on the main workflow, as for example the new
  'RTNduals' package.

[RTNduals](https://bioconductor.org/packages/RTNduals)
--------

Changes in version 1.0.0:

- 1st Bioconductor release of RTNduals &#91;2017-03-01&#93;.

[S4Vectors](https://bioconductor.org/packages/S4Vectors)
---------

Changes in version 0.14.0:

NEW FEATURES

- Add Linteger vectors: similar to ordinary integer vectors (int values
  at the C level) but store "large integers" i.e. long long int values
  at the C level. These are 64-bit on Intel platforms vs 32-bit for int
  values. See ?Linteger for more information. This is in preparation
  for supporting long Vector derivatives (planned for BioC 3.6).

- Default "rank" method for Vector objects now supports the same ties
  method as base::rank() (was only supporting ties methods "first" and
  "min" until now).

- Support x&#91;&#91;i,j&#93;&#93; on DataFrame objects.

- Add "transform" methods for DataTable and Vector objects.

SIGNIFICANT USER-VISIBLE CHANGES

- Rename union classes characterORNULL, vectorORfactor,
  DataTableORNULL, and expressionORfunction -> character_OR_NULL,
  vector_OR_factor, DataTable_OR_NULL, and expression_OR_function,
  respectively.

- Remove default "xtfrm" method for Vector objects. Not needed and
  introduced infinite recursion when calling order(), sort() or rank()
  on Vector objects that don't have specific order/sort/rank methods.

DEPRECATED AND DEFUNCT

- Remove compare() (was defunct in BioC 3.4).

- Remove elementLengths() (was defunct in BioC 3.4).

BUG FIXES

- Make showAsCell() robust to nested lists.

- Fix bug where subsetting a List object 'x' by a list-like subscript
  was not always propagating 'mcols(x)'.

[sampleClassifier](https://bioconductor.org/packages/sampleClassifier)
----------------

Changes in version 0.99.0:

- First submission to BioConductor

[scater](https://bioconductor.org/packages/scater)
------

Changes in version 1.3.49:

- plotRLE() function to make relative log expression plots to assess
  and compare normalizations

- Refactored newSCESet() with defined hierarchy of data types

- read10XResults() to read in results from 10x Chromium CellRanger
  output

- Refined QC metrics

- Bug fixes, efficiency improvements and more tests

[scone](https://bioconductor.org/packages/scone)
-----

Changes in version 1.0.0:

- Faster SVD with rARPACK.

- New zero mode options in scone.

- Newly exported get functions.

- Updated scone scaling function defaults.

- Error handling and documentation updates.

- Bug fixes to sample filtering functions.

Changes in version 0.99.0 (2016-11-14):

- Major release: many changes due to porting to S4.

- Widespread updates to documentation, including examples.

- Compatibility with Bioconductor format, including biocViews.

- Updated cell-cycle genes.

- Default BPPARAM value, now passed as argument to scone.

- Removed zinb function.

- Test different parallel back ends.

- Added class SconeExperiment based on SummarizedExperiment.

- Added constructor sconeExperiment to create SconeExperiment objects.

- Added many helper methods to retrieve content of slots.

- Wrapper get_normalized() to extract/compute single normalization from
  scone object.

- Wrapper get_design() to extract the design matrix associated with a
  given normalization.

- Method select_methods() to get a smaller SconeExperiment object
  containing only the requested normalization schemes.

- biplot_interactive() now works with SconeExperiment objects.

- Dropped date from DESCRIPTION as better practice is to add date to
  NEWS file.

- Added wrapper for scran normalization, removed FQP.

Changes in version 0.0.8 (2016-09-26):

- Added sconeReport() function for shiny browser of SCONE results.

- scone() outputs are now sorted according to mean *rank* of scores
  rather than mean scores. Although this is a relative measure it
  accounts for differing variability of metrics.

- RLE_IQR now quantifies the variability of the IQR RLE across samples,
  rather than the mean.

- Bug Fix: Previously imputation could be applied to the wrong subset
  of params rows when user passed params arguments. Imputation
  functions are now indexed by name to avoid this error.

- Bug Fix: Negative infinity expected likelihood is temporarily
  permitted in the ziber loop.

- "ezbake" script and scone_easybake() function added for pipelined
  SCONE commands.

- Revised documentation and removed old scripts.

- Various other bug fixes.

Changes in version 0.0.6 (2016-07-22):

- Added option for restoring zeroes after scaling step.

- New argument format for imputation via impute_args.

- Simplified ziber fnr estimation function - requires control genes.

- Fixed bug when using plot functionality of filtering functions.

- "Conditional" pam replaced with "Stratified" pam, clustering each
  bio-cross-batch condition separately, rather than simply each bio
  condition.

- Simple FNR for filtering is now based on medians of natural log
  expression in "expressing" cells / robust to convergence issues.

- Removed all sufficient thresholds for metric sample filter.

- Added option to write normalized matrices to HDF5 file.

- Added wrapper function get_normalized() to retrieve normalized data.

- New biplot_interactive function to explore the results.

Changes in version 0.0.5:

- Modified biplot to handle general coloring schemes.

- Limit number of WV and UV factors to eval_pcs in computing WV and UV
  scores.

- Updated dependencies.

- Added error-handling to sample filter.

- Removed var preserved measure due to length of running time.

Changes in version 0.0.4:

- Fixed a few bugs and documentation mismatches.

- Removed stability evaluation (redundant with sil width and slow).

- Removed clusterExperiment dependency.

- Removed RUV correlation score. UV correlation now takes ruv control
  genes as default.

- Added RLE measures to scone evaluation.

- Added FQT_FN to implement careful ties handling by FQ.

- Better handling of plots in sample filter functions.

- Mean score rather than median rank to evaluate normalizations.

- Default value for imputation.

- Minor optimizations to evaluation functions.

Changes in version 0.0.3:

- Fixed various bugs.

- Added “compactness” measure for stability evaluation.

- Compute RUV factors only when needed.

- Fixed Github issues #11, #12, #13, #14, #21, #28.

- zinb now works for non-integer whole numbers.

- Updated tests.

- Added documentation for datasets.

- Added biplot_colored function.

[scran](https://bioconductor.org/packages/scran)
-----

Changes in version 1.4.0:

- Switched default BPPARAM to SerialParam() in all functions.

- Added run argument to selectorPlot(). Bug fix to avoid adding an
  empty list.

- Added exploreData() function for visualization of scRNA-seq data.

- Minor bug fix to DM() when extrapolation is required.

- Added check for centred size factors in trendVar(), decomposeVar()
  methods. Refactored trendVar() to include automatic start point
  estimation, location rescaling and df2 estimation.

- Moved spike-in specification to the scater package.

- Deprecated isSpike<- to avoid confusion over input/output types.

- Generalized sandbag(), cyclone() to work for other classification
  problems.

- Added test="f" option in testVar() to account for additional scatter.

- Added per.gene=FALSE option in correlatePairs(), expanded accepted
  value types for subset.row. Fixed an integer overflow in
  correlatePairs(). Also added information on whether the permutation
  p-value reaches its lower bound.

- Added the combineVar() function to combine results from separate
  decomposeVar() calls.

- Added protection against all-zero rows in technicalCV2().

- Added the improvedCV2() function as a more stable alternative to
  technicalCV2().

- Added the denoisePCA() function to remove technical noise via
  selection of early principal components.

- Removed warning requiring at least twice the max size in
  computeSumFactors(). Elaborated on the circumstances surrounding
  negative size factors. Increased the default number of window sizes
  to be examined. Refactored C++ code for increased speed.

- Allowed quickCluster() to return a matrix of ranks for use in other
  clustering methods. Added method="igraph" option to perform
  graph-based clustering for large numbers of cells.

- Added the findMarkers() function to automatically identify potential
  markers for cell clusters.

- Added the overlapExprs() function to compute the overlap in
  expression distributions between groups.

- Added the buildSNNGraph() function to build a SNN graph for cells
  from their expression profiles.

- Added the correctMNN() function to perform batch correction based on
  mutual nearest neighbors.

- Streamlined examples when mocking up data sets.

[semisup](https://bioconductor.org/packages/semisup)
-------

Changes in version 1.0.0 (2016-04-25):

- Added functions

- Added help pages

- Added vignette

[SeqArray](https://bioconductor.org/packages/SeqArray)
--------

Changes in version 1.16.0:

- a new argument 'intersect' in `seqSetFilter()` and
  `seqSetFilterChrom()`

- a new function `seqSetFilterCond()`

- `seqVCF2GDS()` allows arbitrary numbers of different alleles if REF
  and ALT in VCF are missing

- optimize internal indexing for FORMAT annotations to  avoid reloading
  the indexing from the GDS file

- a new CITATION file

- 'LZMA_RA' is the default compression method in `seqBED2GDS()` and
  `seqSNP2GDS()`

- `seqVCF_Header()` correctly calculates ploidy with missing genotypes

Changes in version 1.15.0-1.15.6:

- the version number was bumped for the Bioconductor develop version
  3.4

Changes in version 1.14.1:

- The default compression setting in `seqVCF2GDS()` and `seqMerge()` is
  changed from "ZIP_RA" to "LZMA_RA"

- `seqVCF2GDS()`: variable-length encoding method is used to store
  integers in the FORMAT field of VCF files to reduce the file size and
  compression time

[seqTools](https://bioconductor.org/packages/seqTools)
--------

Changes in version 1.9.1:

NEW FEATURES

- (none)

SIGNIFICANT USER-VISIBLE CHANGES

- (none)

BUG FIXES

- Added entry in NAMESPACE

[SGSeq](https://bioconductor.org/packages/SGSeq)
-----

Changes in version 1.10.0:

- Bug fixes and documentation improvements

[signeR](https://bioconductor.org/packages/signeR)
------

Changes in version 1.1.7:

- Reorder method was renamed to Reorder_signatures.

- New methods Reorder_samples and Reorder_mutations was added.

Changes in version 1.1.2:

- Pre-calculate plot data during object construction

[SIMAT](https://bioconductor.org/packages/SIMAT)
-----

Changes in version 1.7.1:

- bug related to ggplot package fixed &#91;2016-12-07 Wed&#93;

[SNPRelate](https://bioconductor.org/packages/SNPRelate)
---------

Changes in version 1.10.0:

- new functions `snpgdsAdmixPlot()` and `snpgdsAdmixTable()`

- `snpgdsPCASNPLoading()` and `snpgdsPCASampLoading()` support the
  eigen results of `snpgdsEIGMIX()` allowing projecting new samples to
  the existing coordinate

- `snpgdsFst()` provides W&C84 mean Fst together with weighted Fst

- a new argument 'outgds' in `snpgdsPCACorr()` allows exporting
  correlations to a gds file

- a friendly warning is given when openning a SeqArray file with
  `snpgdsOpen()`

- a new option "Corr" in `snpgdsGRM()` for scaled GRM

[specL](https://bioconductor.org/packages/specL)
-----

Changes in version 1.9.15 (2017-04-21):

- added prozor (>= 0.2.2) to the Suggest list.

- added more specific R package version numbers in DESCRIPTION file.

- in plot.specLSet (normalized RT versus RT) use pch=16 and color with
  parameter alpha=0.1.

- fixed issue #22 by including the iRTs in the ionlibrary; LIB <-
  genSwathIonLib(data=peptideStd, data.fit=peptideStd.redundant);
  LIB@input.parameter$iRTpeptides.

- fixed issue #19.

- removed par command in specLset plot function.

- added vignettes/report.Rmd file, see also <URL:
  http://bioconductor.org/packages/devel/bioc/vignettes/specL/inst/doc/report.html>.

[splatter](https://bioconductor.org/packages/splatter)
--------

Changes in version 0.99.16 (2017-04-23):

- Splatter is a package for the simple simulation of single-cell
  RNA-seq data, including:

- Multiple simulation models

- Parameter estimation from real data

- Functions for comparing simulations and real datasets

- Simulation of complex groups and differentiation paths

Changes in version 0.99.0 (2016-12-05):

- Package prepared for Bioconductor submission.

[spliceSites](https://bioconductor.org/packages/spliceSites)
-----------

Changes in version 1.23.5:

- Corrected two errors (addMaxEnt function and addGenomeData function)
  Using the old version, errors in calculation of mxe_ps3, mxe_ms5 and
  mxe_ms3 have occured. Also erroneus wgis values were calculated. BUG
  FIXES

- (none)

[statTarget](https://bioconductor.org/packages/statTarget)
----------

Version: 1.5.7
Category: The result of Fold Change will be calculated accroding to the
        raw data in
Text:

Version: 1.5.6
Category: NEW FEATURES
Text: transX and transCode was added to generate statTarget inputs from
        Mass Spectrometry Data softwares, like XCMS.

[SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment)
--------------------

Changes in version 1.6.0:

NEW FEATURES

- Add saveHDF5SummarizedExperiment() and loadHDF5SummarizedExperiment()
  for saving/loading HDF5-based SummarizedExperiment objects to/from
  disk.

DEPRECATED AND DEFUNCT

- Remove SummarizedExperiment0 class (was introduced to ease transition
  from old SummarizedExperiment class defined in GenomicRanges to new
  RangedSummarizedExperiment class defined in SummarizedExperiment
  package).

[SWATH2stats](https://bioconductor.org/packages/SWATH2stats)
-----------

Changes in version 1.5.8:

NEW FEATURES

- filter_on_max_peptides: add removeDecoyProteins and
  unifyProteinGroupLabels to function

- filter_on_min_peptides: add removeDecoyProteins and
  unifyProteinGroupLabels to function

BUG FIXES

- filter_on_max_peptides: unique selected peptides to prevent
  duplication of rows

Changes in version 1.5.7:

NEW FEATURES

- featurealigner2msstats_withRT.py: reads now csv and tab delimited
  file and checks for right data structure

Changes in version 1.5.6:

NEW FEATURES

- plot.fdr_cube: add option to select mscore levels to plot FDR
  estimation.

- assess_fdr_byrun: add option to select mscore levels to plot FDR
  estimation.

Changes in version 1.5.5:

NEW FEATURES

- sample_annotation: added "fixed" option to grep function to increase
  speed.

BUG FIXES

- convert4mapDIA: diagnostic output in case there were several values
  was not displaying correct rows.

Changes in version 1.5.4:

BUG FIXES

- removeDecoyProteins and unifyProteinGroupLabels: fixed bug that
  affected protein groups of 10-19 proteins.

Changes in version 1.5.3:

NEW FEATURES

- plot_variation: make function to work also if comparison contains
  more than 3 elements

- plot_variation_vs_total: make function to work also if comparison
  contains more than 3 elements

- sample_annotation: introduce fail-safe if input data is in the
  data.table format.

Changes in version 1.5.2:

NEW FEATURES

- unifyProteinGroupLabels: unifies different ProteinGroupLabels

- removeDecoyProteins, rmDecoyProt: Removes decoy protein labels from
  protein Group label

Changes in version 1.5.1:

NEW FEATURES

- SWATH2stats in BioC 3.5 development release

Changes in version 1.4.1:

NEW FEATURES

- SWATH2stats in BioC 3.4 release

[swfdr](https://bioconductor.org/packages/swfdr)
-----

Changes in version 0.99.0:

- First version of swfdr package submitted to Bioconductor

[synapter](https://bioconductor.org/packages/synapter)
--------

Version: 1.99.2
Text: update show,MasterPeptides to check of there's a fragmentlibrary
        slot before trying to access it <2017-04-19 Wed>

Version: 1.99.1
Text: update NEWS file <2017-04-11 Tue>

Version: 1.99.0
Category: NEW FEATURES
Text:

Version: 1.99.0
Category: ION MOBILITY/GRID SEARCH
Text: Replace 2D grid search (retention time, *m/z*) of synapter1 by 3D
        grid search (retention time, *m/z*, ion mobility); set argument
        `imdiff = Inf` to get the original 2D grid search; closes #33.

Version: 1.99.0
Category: ION MOBILITY/GRID SEARCH
Text: Add `{set,get}ImDiff` methods.

Version: 1.99.0
Category: ION MOBILITY/GRID SEARCH
Text: `getGrid` returns an array instead of a matrix (because of the
        new 3D grid search) &#91;2014-05-16 Fri&#93;.

Version: 1.99.0
Category: ION MOBILITY/GRID SEARCH
Text: `plotFeatures(..., what = "all")` gains a new argument:
        "ionmobilty" to plot *m/z* vs ionmobility as well. &#91;2014-05-16
        Fri&#93;

Version: 1.99.0
Category: ION MOBILITY/GRID SEARCH
Text: `plotGrid` gets a new argument "maindim" to decided which of the
        three dimension should be used. &#91;2014-05-16 Fri&#93;

Version: 1.99.0
Category: ION MOBILITY/GRID SEARCH
Text: Add `filterNonUniqueIdentMatches` to remove matches of multiple
        identification data to a single quantitation entry (see #111
        for details) &#91;2016-02-22 Mon&#93;.

Version: 1.99.0
Category: FRAGMENT MATCHING
Text: Load identification fragments (`final_fragments.csv`) and
        quantitation spectra (`Spectrum.xml`) via `Synapter`
        constructor.

Version: 1.99.0
Category: FRAGMENT MATCHING
Text: New functions: `fragmentMatchingPerformance`,
        `filterUniqueMatches`, `filterNonUniqueMatches`,
        `filterFragments`, `plotCumulativeNumberOfFragments`,
        `plotFragmentMatchingPerformance`, `getIdentificationFragments`
        and `getQuantitationSpectra`.

Version: 1.99.0
Category: FRAGMENT MATCHING
Text: Integrate a fragment library into *master* objects; closes #63
        and #74.

Version: 1.99.0
Category: MISC
Text: Allow to use an RDS instead of a fasta file as 'Unique Peptides
        Database', adds `createUniquePeptideDbRds`; closes #55
        &#91;2014-04-29 Tue&#93;.

Version: 1.99.0
Category: MISC
Text: Introduce `IisL` argument to `dbUniquePeptideSet` which treats
        I/L as same aminoacid if `IisL == TRUE` (default: `IisL =
        FALSE`); closes #60 &#91;2014-04-30 Wed&#93;.

Version: 1.99.0
Category: MISC
Text: Add `rescueEMRTs` functions; replaces the argument `mergedEMRTs`
        in `findEMRTs`; closes #93 &#91;2015-07-26 Sun&#93;.

Version: 1.99.0
Category: MISC
Text: Add `synergise2` which combines the integrates the new 3D grid
        search, the fragment matching; and uses slightly different
        default arguments than `synergise1`; closes #119 &#91;2016-10-25
        Di&#93;.

Version: 1.99.0
Category: MISC
Text: Load isotopic distributions from Pep3D data and also export them
        to MSnSet, to allow the correction of detector saturation;
        closes #39 &#91;2015-03-29 Sun&#93;.

Version: 1.99.0
Category: MISC
Text: Add `synapterPlgsAgreement` to find agreement between *synapter*
        and *PLGS*; closes #73.

Version: 1.99.0
Category: MISC
Text: Introduce `modelIntensity` to correct systematic intensity shifts
        (similar to `modelRt`); closes #116.

Version: 1.99.0
Category: IMPROVEMENTS
Text: Extract the ion that was used for identification (`isFid == 1`)
        from the Pep3D file instead of the first instance &#91;2014-05-13
        Tue&#93;.

Version: 1.99.0
Category: IMPROVEMENTS
Text: Add `updateObject` and `validObject` method &#91;2014-11-16 Sun&#93;.

Version: 1.99.0
Category: IMPROVEMENTS
Text: Rename `QuantPep3DData$Function` column into
        `QuantPep3DData$matchedEMRTs`; closes #67 &#91;2015-07-26 Sun&#93;.

Version: 1.99.0
Category: IMPROVEMENTS
Text: Use just unique peptides in master creation (see #107)
        &#91;2016-01-23 Sat&#93;.

Version: 1.99.0
Category: IMPROVEMENTS
Text: New `rmarkdown` based reports for `synergise1` (synonym to
        `synergise`) and `synergise2`.

Version: 1.99.0
Category: BUGFIXES
Text: Use new loess model in master creation (now based on m-estimator
        instead of least squares, identical to retention time model in
        classical synergise workflow; see #107 for details) &#91;2016-01-23
        Sat&#93;

Version: 1.99.0
Category: BUGFIXES
Text: Fix retention time model calculation in `plotFeatures(...,
        what="some")` &#91;2014-04-28 Mon&#93;.

Version: 1.99.0
Category: INTERNAL CHANGES
Text: Add `testthat` to Suggests &#91;2014-04-25 Fri&#93;.

Version: 1.99.0
Category: INTERNAL CHANGES
Text: Add recommended biocView &#91;2014-06-05 Thu&#93;.

Version: 1.99.0
Category: INTERNAL CHANGES
Text: Replace `any(is.na(...)` by `anyNA(...)`; *synapter* depends on
        `R >= 3.1.0` now &#91;2014-11-01 Sat&#93;.

Version: 1.99.0
Category: INTERNAL CHANGES
Text: Add `ClassVersion` field to `Synapter` class &#91;2014-11-21 Fri&#93;.

Version: 1.99.0
Category: INTERNAL CHANGES
Text: Add `Versioned` class as parent class to `MasterPeptides` and
        `MasterFdrResults` &#91;2014-11-22 Sat&#93;.

Version: 1.99.0
Category: INTERNAL CHANGES
Text: Adapt `synergise` to new grid search (closes #81) &#91;2016-10-16
        So&#93;.

Version: 1.99.0
Category: INTERNAL CHANGES
Text: Replace `hwriter` by `rmarkdown` report in `synergise`; closes
        #120. &#91;2016-10-17 Mon&#93;

Version: 1.99.0
Category: REMOVED FUNCTIONS/ARGUMENTS
Text: Remove `synapterGUI`.

Version: 1.99.0
Category: REMOVED FUNCTIONS/ARGUMENTS
Text: Remove unused internal functions: `filterCommonSeq`,
        `filterKeepUniqueSeq`, `filterKeepUniqueProt` &#91;2014-11-27 Thu&#93;.

Version: 1.99.0
Category: REMOVED FUNCTIONS/ARGUMENTS
Text: Remove "mergedEMRTs" argument from `findEMRTs`. Now `rescueEMRTs`
        has to be called manually at the end of the processing; close
        #93 &#91;2015-07-26 Sun&#93;

Version: 1.99.0
Category: REMOVED FUNCTIONS/ARGUMENTS
Text: Remove "light" version of `writeMergedPeptides` and
        `writeMachtedPeptides` (now always the full `data.frame` is
        saved; see #95) &#91;2016-10-16 Sun&#93;

Version: 1.99.0
Category: REMOVED FUNCTIONS/ARGUMENTS
Text: Update `synapterTiny` and `synapterTinyData` &#91;2016-10-16 So&#93;

[TarSeqQC](https://bioconductor.org/packages/TarSeqQC)
--------

Changes in version 1.5.1:

CODE

- Fix error of duplicated parameter in plotRegion function

- Changes in TargetExperiment contstructor to avoid errors related to
  unmapped reads in the alignment BAM files

- The checkBedFasta function was added to perform a control of the Bed
  and Fasta file consistency

- Changes in plotRegion methods to allow filter noise SNPs

- Changes in plotGeneAttrPerFeat method to incorporate the exploration
  of overlapped amplicons

VIGNETTE

- Troubleshoot explaining the way to exclude unmapped reads with older
  TarSeqQC versions

[TCGAbiolinksGUI](https://bioconductor.org/packages/TCGAbiolinksGUI)
---------------

Version: 0.99.0
Text: FIRST VERSION. Demo:
        https://tcgabiolinksgui.shinyapps.io/TCGAbiolinksGUI/

[TFBSTools](https://bioconductor.org/packages/TFBSTools)
---------

Changes in version 3.5:

NEW FEATURES

- Add function to parse MEME output.

Changes in version 3.4:

BUG FIXES

- Fix a bug in PFMSimilarity.

- Fix an error when there are multiple classes for motif matrx.

NEW FEATURES

- readJASPARMatrix function to add the JASPAR format file.

Changes in version 3.3:

BUG FIXES

- Adapt the runMEME to work with meme 4.10.x version.

- Fix the scientific notation in run_MEME

- Better error handling of MEME wrappe

[tigre](https://bioconductor.org/packages/tigre)
-----

Changes in version 1.30.0:

BUG FIXES

- RSQLite 1.1 compatibility update

[tofsims](https://bioconductor.org/packages/tofsims)
-------

Version: 1.3.1
Category: SIGNIFICANT USER-VISIBLE CHANGES
Text: Import of Ulvac-Phi Raw data from WinCadence V1.18.0.22 is now
        supported

Version: 1.3.1
Category: INTERNALS
Text:

Version: 1.3.1
Category: BUGFIXES
Text:

Version: 099.1
Category: SIGNIFICANT USER-VISIBLE CHANGES
Text: changed function behvaiour in the whole package from call-by-ref
        to call-by value. Adjusted accordingly all examples and the
        vignette.

Version: 099.1
Category: INTERNALS
Text: depends now on ProtGenerics from which it uses 'mz'

Version: 099.1
Category: INTERNALS
Text: exchanged various print() with message()

Version: 099.1
Category: BUGFIXES
Text:

[TPP](https://bioconductor.org/packages/TPP)
---

Changes in version 3.2.2:

- Bug fix in the TPP-TR analysis: - It is now possible to leave the
  'otherRequirements' slot empty in the object containing the
  normalization criteria.

Changes in version 3.2.0:

- Changes to the Excel output of TPP-TR experiments: - Plot paths in
  excel output are stored in separate columns for spline and melting
  curve fits - Reported p-values for Tm based analysis are now
  increased in value after removing a bug in their calculation. The
  thresholds applied for determining significant hits have been updated
  justed accordingly.

- Improvements to the 2D-TPP analysis: - Split rows with multiple
  identifiers (e.g. separated by '|') into separate rows (was already
  present in TR- and CCR analysis) - Sort result table rows by
  temperature for each ID

Changes in version 3.1.3:

- Bug fix: ensure correct handling of drug concentrations when they are
  imported in scientific notation (xx.xxE-x)

- Deprecate functions and arguments, primarily those used in the first
  version of TPP-2D analysis workflow.

- Bug fixes: - sort x and y values before DR curve fitting because the
  functions for initial parameter estimation are dependent on the
  correct ordering - supress useless VD logfiles when creating VENN
  diagrams

[trackViewer](https://bioconductor.org/packages/trackViewer)
-----------

Changes in version 1.11.8:

- support break x-axis for viewTracks.

Changes in version 1.11.7:

- improve the ideogramPlot function.

Changes in version 1.11.6:

- improve the ideogramPlot function.

Changes in version 1.11.5:

- add ideogramPlot function.

Changes in version 1.11.4:

- remove the blank bottom space for lolliplot.

Changes in version 1.11.3:

- fix a bug when tracks has multiple max values.

Changes in version 1.11.2:

- add jitter label only for lolliplot.

Changes in version 1.11.1:

- adjust the x postion of lolliplot.

[treeio](https://bioconductor.org/packages/treeio)
------

Changes in version 0.99.11:

- bug fixed in get.fields method for paml_rst <2017-03-20, Mon>

- fixed raxml2nwk for using treedata as output of read.raxml
  <2017-03-17, Fri>

- taxa_rename function <2017-03-15, Wed>

- phyPML method moved from ggtree <2017-03-06, Mon>

Changes in version 0.99.10:

- remove raxml class, now read.raxml output treedata object
  <2017-02-28, Tue>

- bug fixed of read.beast <2017-02-27, Mon>

Changes in version 0.99.9:

- read.newick for parsing node.label as support values <2017-01-03,
  Tue>

- read.beast support MrBayes output <2016-12-30, Fri>

- export as.phylo.ggtree <2016-12-30, Fri>

Changes in version 0.99.8:

- as.treedata.ggtree <2016-12-30, Fri>

- as.treedata.phylo4 & as.treedata.phylo4d <2016-12-28, Wed>

Changes in version 0.99.7:

- groupOTU, groupClade, gzoom methods from ggtree <2016-12-21, Wed>

Changes in version 0.99.6:

- add unit test of NHX (move from ggtree) <2016-12-14, Wed>

Changes in version 0.99.3:

- fixed BiocCheck by adding examples <2016-12-07, Wed>

Changes in version 0.99.1:

- fixed link in DESCRIPTION <2016-12-06, Tue>

Changes in version 0.99.0:

- add vignette <2016-12-06, Tue>

- move parser functions from ggtree <2016-12-06, Tue>

Changes in version 0.0.1:

- read.nhx from ggtree <2016-12-06, Tue>

- as.phylo.treedata to access phylo from treedata object <2016-12-06,
  Tue>

- as.treedata.phylo to convert phylo to tree data object <2016-12-06,
  Tue>

- treedata class definition <2016-12-06, Tue>

[TRONCO](https://bioconductor.org/packages/TRONCO)
------

Changes in version 2.7.7:

- RNA Seq validation

- Random restart on Hill Climbing added to CAPRI algorithm

- Minor fixes to algorithms and error model

Changes in version 2.7.3:

- Assignment to .GlobalEnv removed

Changes in version 2.6.1:

- Last stable releade with bugfix

[TVTB](https://bioconductor.org/packages/TVTB)
----

Changes in version 1.1.12 (2017-04-10):

Major changes

- Update GRangesFilter following update of the condition argument to
  type.

Changes in version 1.1.11 (2017-04-08):

Major changes

- Update NAMESPACE imports following move of GRangesFilter and
  GenenameFilter from AnnotationFilter to ensembldb.

Changes in version 1.1.10 (2016-11-18):

Bug fix

- Disambiguated the variable metric that was used for two different
  things in the same method and produced incorrect DataTrack name in
  plotInfo.

Major changes

- Vignette uses alternate allele frequency to demonstrate pairsInfo and
  plotInfo methods.

Changes in version 1.1.9 (2016-11-18):

Major changes

- New argument zero.rm in method plotInfo to hide data points with
  value of zero on the Y axis. Intended to reduce overplotting of
  variants absent from phenotype levels, and instead emphasise variants
  of low frequency.

Changes in version 1.1.8 (2016-11-14):

New features

- New method pairsInfo to visualise a matrix of pairwise plots that
  displays a metric calculated in levels of a given phenotype, and
  stored in columns of the info slot of a VCF object.

- ggpairs method imported from the GGally package.

Bug fix

- Updated pattern used to detect INFO columns for a given metric, to be
  more specific.

Major changes

- Internal method .findInfoMetricColumns moved to utils.R, as it is now
  used in two different user-visible methods.

- Updated Introduction vignette to better present usage of the
  addFrequencies method, better present _VCF filter rules_, and
  introduce the new method pairsInfo.

Changes in version 1.1.7 (2016-11-11):

Major changes

- Plotting type(s) cab be selected for the DataTrack of plotInfo.

Minor changes

- Update to README.

- _AppVeyor_ caches R library.

Changes in version 1.1.6 (2016-11-10):

New features

- plotInfo method to visualise a metric calculated in levels of a given
  phenotype, and stored in columns of the info slot.

- Methods imported from Gviz package.

- More methods imported from ensembldb package.

Major changes

- Added new VCF file and associated preprocessing script in extdata/
  for gene ADH1B.

- Introduction vignette updated to change activated filter rules.

- Introduction vignette updated to introduce the plotInfo method.

Minor changes

- Moved content of the table of motivations to implement _VCF filter
  rules_ to a CSV file in misc/.

- BED file and VCF files for gene ADH1B.

- Updated Shell script to preprocess VCF files with the VEP script.

- Ignore .svn/.

Changes in version 1.1.5 (2016-11-08):

Bug fix

- Updated reference to renamed object in _Shiny_ app.

- When no phenotypes are supplied, set phenoData slot to a DataFrame
  with rownames set to colnames(vcf) and 0 columns, instead of the
  default behaviour of the VariantAnnotation package which is to create
  a column named Sample filled with seq_along.

- NEWS file closing brackets.

Major changes

- autodetectGTimport setting available in tSVE method.

- New checkbox in _Shiny_ app to update selected genotypes after
  importing variants and autodetection of genotypes present in the
  data.

Minor changes

- Do not ignore *.Rproj files.

- Removed commented lines in _AppVeyor_ YAML file.

- Removed files in misc/.

- Display list of error messages in a new session panel of the _Shiny_
  app.

Changes in version 1.1.4 (2016-11-04):

Bug fix

- Bumped version number to try and update the Bioconductor GitHub
  mirror to display the latest code instead of version 0.1.7.

Minor changes

- Added TVTB.Rproj to tracked files.

- Deleted deprecated and misc files in inst/.

- Deleted commented lines from AppVeyor settings.

Changes in version 1.1.3 (2016-11-03):

New features

- The autodetectGenotypes method creates or updates the genotypes
  defined in the TVTBparam that is stored in the metadata slot of a VCF
  object.

- The argument autodetectGT of the readVcf method may be used to call
  the new autodetectGenotypes method immediately after a VCF object is
  initialised from the parsed VCF file.

Major changes

- vepInPhenoLevel returns a GRanges instead of a data.frame; the key
  advantage is that ranges may have non-unique names.

- Genotypes objects can now be initialised without specifying ref, het,
  and alt genotype vectors (with a warning).  A default Genotypes
  object is created with ref, het, and alt slots set to NA_character_.
  The new autodetectGenotypes method may be used to populate those
  slots after variants are imported (see _New features_ section).

- TVTBparam objects can now be initialised without supplying a
  Genotypes object (with a warning).  A default Genotypes object is
  created (see above).

- Constructors for classes Genotypes and TVTBparam are now high-level
  methods, *not* S4 methods methods anymore.

- Default settings of the _Shiny_ app are stored as an environment that
  can be overriden by arguments of the tSVE method.

- _Shiny_ app stores more objects in reactiveValues.

- _Shiny_ app stores more error messages in reactiveValues to better
  deal with optional inputs and better help users to resovle sources of
  errors.

Minor changes

- The show method throws warning messages for TVTBparam and Genotypes
  objects that have not fully defined all genotypes.

- Better layout of badges in README.

- Non-reactive settings of the _Shiny_ app stored in hidden objects.

- Helper methodS getEdb, tryParsePheno, tryParseBed, tryParseVcfHeader,
  tryParseMultipleVcf, and tryParseSingleVcf removed and integrated
  into the server side of the _Shiny_ app.

- Massive cleaning of messages in the global.R file of the _Shiny_ app.

- GRanges, Genotypes, and Phenotypes panels removed from Session panel
  of the _Shiny_ app.

- Table reporting status of BiocParallel configurations of the _Shiny_
  app on various system stored as an RDS file.

- _Shiny_ app displays a warning at the top of the screen if the
  genotypes are not fully defined.

- Tab width of _Shiny_ files set to 2.

- Branches tracked by Travis CI.

- Added a couple of files in inst/badexamples folder.

- Added YAML file for _AppVeyor_.

- Added pander in Suggests section of DESCRIPTION, to render vignette
  tables.

Changes in version 1.1.2 (2016-10-21):

Minor changes

- README indicates status on BioC-release, BioC-devel, and Travis CI.

Changes in version 1.1.1 (2016-10-21):

Minor changes

- Updates to README: weblinks, installation, unit tests.

- Branches tracked by Travis CI.

- Coverage: exclude AllClasses.R, tSVE.R.

- Four-space indents in DESCRIPTION.

[tximport](https://bioconductor.org/packages/tximport)
--------

Changes in version 1.3.8:

- Support for inferential replicates written by Rob Patro! Works for
  Salmon, Sailfish and kallisto. See details in ?tximport.

Changes in version 1.3.6:

- Now, 'countsFromAbundance' not ignored when txOut=TRUE.

Changes in version 1.3.4:

- Support for kallisto HDF5 files thanks to Andrew Parker Morgan and
  Ryan C Thompson

- Removing 'reader' argument, leaving only 'importer' argument. In
  addition, read_tsv will be used by default if readr package is
  installed.

- Messages from the importing function are captured to avoid screen
  clutter.

[variancePartition](https://bioconductor.org/packages/variancePartition)
-----------------

Changes in version 1.5.7:

- include splines in foreach .packages

Changes in version 1.5.6:

- compatibility with tximport v1.3.5

Changes in version 1.5.5:

- Decrease computing time of effective sample size with ESS() by
  additional ~10x with sparse solver

- fix margins for plotPercentBars()

- Fix bug for getVarianceComponents() when correlated continous
  variables are included

- compatibility with ggplot2 2.2.0

- center plot titles

- fix order of bars in plotPercentBars()

- legend background to transparent

- set text to be black

- include lme4 in foreach .packages

- change residuals color to not be transparent

- add CITATION information

- plotCorrMatrix now shows dendrogram by default

- Estimate run time for fitExtractVarPartModel() / fitVarPartModel()

- improve warnings for plotPercentBar()

- improve warnings for plotCorrStructure()

- define ylab for plotVarPart() - add as.matrix.varPartResults()
  (hidden) - define isVaryingCoefficientModel() (hidden)

[VariantAnnotation](https://bioconductor.org/packages/VariantAnnotation)
-----------------

Changes in version 1.22.0:

NEW FEATURES

- add import() wrapper for VCF files

- add support for Number='R' in vcf parsing

- add indexVcf() and methods for character,VcfFile,VcfFileList

MODIFICATIONS

- throw message() instead of warning() when non-nucleotide variations
  are set to NA

- replace 'force=TRUE' with 'pruning.mode="coarse"' in seqlevels()
  setter

- add 'pruning.mode' argument to keepSeqlevels() in man page example

- idempotent VcfFile()

- add 'idType' arg to IntergenicVariants() constructor

- modify locateVariants man page example to work around issue that
  distance,GRanges,TxDb does not support gene ranges on multiple
  chromosomes

- modify VcfFile() constructor to detect index file if not specified

- order vignettes from intro to advanced; use BiocStyle::latex2()

- remove unused SNPlocs.Hsapiens.dbSNP.20110815 from the Suggests field

- follow rename change in S4Vectors from vector_OR_factor to
  vector_or_factor

- pass classDef to .RsamtoolsFileList; VariantAnnotation may not be on
  the search path

BUG FIXES

- fix expansion of 'A' fields when there are multiple columns

[VariantFiltering](https://bioconductor.org/packages/VariantFiltering)
----------------

Changes in version 1.12:

USER VISIBLE CHANGES

- The VariantFilteringParam constructor is restricted to one
  (multisample) VCF file.

- mafByOverlaps() returns now a GRanges object with minor allele
  frequency values in the metadata columns.

[vsn](https://bioconductor.org/packages/vsn)
---

Changes in version 3.43.10:

- Vignette now is in HTML, using BiocStyle::html_document2

[wateRmelon](https://bioconductor.org/packages/wateRmelon)
----------

Version: 1.19.1
Text: and methods will return 'MethylSet' instead of beta matrix.

[wiggleplotr](https://bioconductor.org/packages/wiggleplotr)
-----------

Changes in version 1.0.0:

- First submission to Bioconductor.

[xcms](https://bioconductor.org/packages/xcms)
----

Changes in version 1.51.11:

NEW FEATURES

- Parameter "filled" for featureValues (issue #157).

- Parameters "rt" and "mz" in chromPeaks method allowing to extract
  chromatographic peaks from the specified ranges (issue #156).

BUG FIXES

- Fixed possible memory problem in obiwarp (issue #159).

- Update getPeaks to use non-deprecated API (issue #163).

Changes in version 1.51.10:

NEW FEATURES

- filterRt for Chromatogram class (issue #142).

- adjustRtimePeakGroups function (issue #147).

- adjustRtime,XCMSnExp,PeakGroupsParam and do_adjustRtime_peakGroups
  support use of pre-defined matrix to perform alignment (issue #153).

- plotAdjustedRtime to visualize alignment results (issue #141).

USER VISIBLE CHANGES

- featureDefinitions and featureValues return DataFrame and matrix with
  rownames corresponding to arbitrary feature IDs (issue #148).

- New peakGroupsMatrix slot for PeakGroupsParam class (issue #153).

BUG FIXES

- Issue #146: ensure adjusted retention times returned by the
  peakGroups method to be in the same order than the raw retention
  times.

Changes in version 1.51.9:

NEW FEATURES

- fillChromPeaks, dropFilledChromPeaks methods and FillChromPeaksParam
  class.

- featureValues method.

USER VISIBLE CHANGES

- Extended new_functionality vignette.

- Change default backend for reading mzML files to pwiz.

BUG FIXES

- Issue #135: fix peak signal integration for centWave.

- Issue #139: problem with expand.mz and expand.rt in fillPeaks.chrom.

- Issue #137: Error in findChromPeaks if no peaks are found.

Changes in version 1.51.8:

NEW FEATURES

- Add Chromatogram class and extractChromatograms method.

BUG FIXES

- Issue #118: failing unit test on Windows build machine.

- Issue #133: error with c() and xcmsSet without peaks.

- Issue #134: xcmsSet constructor endless loop.

Changes in version 1.51.7:

USER VISIBLE CHANGES

- Major renaming of methods and classes to follow the naming
  convention:

- chromatographic peak (chromPeak): the peaks identified in rt
  dimension.

- feature: mz-rt feature, being the grouped chromatographic peaks
  within and across samples.

BUG FIXES

- Issue #127: failing unit test on Windows build machine.

Changes in version 1.51.6:

NEW FEATURES

- groupFeatures and adjustRtime methods for XCMSnExp objects.

- New Param classes for groupFeatures and adjustRtime analysis methods:
  FeatureDensityParam, MzClustParam, NearestFeaturesParam,
  FeatureGroupsParam and ObiwarpParam.

BUG FIXES

- Issue #124 (filterRt,XCMSnExp returned empty object).

Changes in version 1.51.5:

NEW FEATURES

- MsFeatureData and XCMSnExp objects.

- features, features<-, adjustedRtime, adjustedRtime<-, featureGroups,
  featureGroups<-, hasAlignedFeatures, hasAdjustedRtime and
  hasDetectedFeatures methods.

- dropFeatures, dropFeatureGroups and dropAdjustedRtime methods.

- filterMz, filterRt, filterFile etc implemented.

- mz, intensity and rtime methods for XCMSnExp allowing to return
  values grouped by sample.

BUG FIXES

- Issue #99 (rtrange outside of retention time range in
  getEIC,xcmsSet).

- Issue #101 (xcmsRaw function returns NULL if mslevel = 1 is
  specified).

- Issue #102 (centWave returns empty matrix if scales not OK). Thanks
  to J. Stanstrup.

- Issue #91 (warning instead of error if no peaks in ROI). Thanks to J.
  Stanstrup.

Changes in version 1.51.4:

BUG FIXES

- added deepCopy to avoid corrupting the original object, thanks to J.
  Stanstrup, closes #93

Changes in version 1.51.3:

NEW FEATURES

- binYonX binning function.

- imputeLinInterpol function providing linear interpolation of missing
  values.

- breaks_on_binSize and breaks_on_nBins functions to calculate breaks
  defining bins.

- New vignette "new_functionality.Rmd" describing new and modified
  functionality in xcms.

- Add do_detectFeatures_matchedFilter function.

- Add do_detectFeatures_centWave function.

- Add do_detectFeatures_centWaveWithPredIsoROIs function and unit test.

- Implement a new data import function.

- Add do_detectFeatures_MSW function and unit test.

- Argument stopOnError in xcmsSet function that allows to perform
  feature detection on all files without stopping on errors.

- Method showError for xcmsSet objects that list all errors during
  feature detection (if stopOnError = FALSE in the xcmsSet function).

- [ method to subset xcmsRaw objects by scans.

- profMat method to extract/create the profile matrix from/for an
  xcmsRaw.

- Add new detectFeatures methods for MSnExp and OnDiskMSnExp objects
  from the MSnbase package.

- Add new CentWaveParam, MatchedFilterParam, MassifquantParam, MSWParam
  and CentWavePredIsoParam parameter class to perform method dispatch
  in the detectFeatures method.

- retcor.obiwarp uses the new binning methods for profile matrix
  generation.

- scanrange,xcmsRaw reports always a scanrange of 1 and
  length(object@scantime).

- scanrange,xcmsSet reports the scanrange eventually specified by the
  user in the xcmsSet function.

- Fixed bug in rawMat (issue #58).

- Fix issue #60: findPeaks.massifquant always returns a xcmsPeaks
  object.

Changes in version 1.51.2:

USER VISIBLE CHANGES

- As suggested by Jan Stanstrup, do not error if a centWave ROI
  contains no data, closes #90

Changes in version 1.51.1:

BUG FIXES

- Fix incorrrect indexing getEIC function reported by Will Edmands,
  closes #92

[xps](https://bioconductor.org/packages/xps)
---

Changes in version 3.3:

VERSION xps-1.35.4

- eliminate dependency on ROOTSYS

- update configure.in file

- update Makefile.arch to add ROOTSYS an update -rpath

VERSION xps-1.35.3

- eliminate dependency on DYLD_LIBRARY_PATH and LD_LIBRARY_PATH

- update configure.in file

- update Makefile.arch to use -rpath with ld

VERSION xps-1.35.2

- update Makefile for MacOS Sierra

- update README file

VERSION xps-1.35.1

- update import in NAMESPACE

[yamss](https://bioconductor.org/packages/yamss)
-----

Changes in version 1.1:

- Chunk matrix multiplications in density estimation for faster run
  times.

- Change vignette format from PDF to HTML.

- Fix sessionInfo format in vignette and triggering of data.table
  warnings with nomatch.

- Update citation.

- Fix accent in citation.

Deprecated and Defunct Packages
===============================

One software package (betr) was removed from this release (after being
deprecated in BioC 3.4).

Nine software packages (AtlasRDF, coRNAi, saps, MeSHSim, GENE.E, mmnet,
CopyNumber450k, GEOsearch, pdmclass) are deprecated in this release and
will be removed in BioC 3.6.

Two experimental data packages (encoDnaseI, ggtut) were removed from this
release (after being deprecated in BioC 3.4).

One experimental data package (CopyNumber450kData) is deprecated in this
release and will be removed in BioC 3.6.

