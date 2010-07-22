# ![](/images/icons/blog.gif)BioC 2.6 Released #

**Bioconductor 2.6, consisting of 389 packages and designed to work with R
version 2.11, was released on April 23, 2010.**

We are pleased to announce the release of Bioconductor 2.6. This release
includes 37 new software packages, and many changes to existing packages.
Bioconductor 2.6 consists of 389 software packages and is compatible with
the recently released R 2.11.0.

Bioconductor 2.6 is supported on Linux, 32-bit Windows, Mac OS X 10.5
(Leopard), and Mac OS X 10.6 (Snow Leopard). Bioconductor 2.6 also has
experimental 64-bit Windows builds for most of its packages.

## Contents ##

* Getting Started with Bioconductor 2.6
* New Software Packages
* Additional Software Package Changes

## Getting Started with Bioconductor 2.6 ##

To install Bioconductor 2.6

1. Install R 2.11.  Bioconductor 2.6 has been designed expressly for this
   version of R.
1. Follow the [installation](/install/) instructions.

## New Software Packages ##

There are 37 new packages in this release of Bioconductor.

New sequence analysis tools address infrastructure (`GenomicRanges`,
`Rsamtools`, `girafe`); ChIP-seq (`BayesPeak`, `CSAR`, `PICS`); digital gene
expression and RNA-seq (`DESeq`, `goseq`, `segmentSeq`); and motif discovery
(`MotIV`, `rGADEM`).

Microarray analysis includes new packages for pre-process and
technology-specific assays (`affyILM`, `frma`, `frmaTools`, `BeadDataPackR`,
`MassArray`); analysis of specific experimental protocols (`charm`, `genoCN`,
`iChip`, `methVisual`); and novel statistical methods (`ConsensusClusterPlus`,
`ExpressionView`, `eisa`, `GSRI`, `PROMISE`, `tigre`).

Flow cytometry packages include `SamSPECTRAL`, `flowMeans`, `flowTrans`, and
`iFlow`.

Annotation and integrative analysis are facilitated by new packages interfacing
with GEO (`GEOsubmission`), the Sequence Read Archive (`SRAdb`), and tabulation
of genome sequence project data (`genomes`); the GSRI package to estimate
differentially expressed genes in a gene set; PCA and CCA dependency modeling
(`pint`); and updated access to exon array annotations (`xmapcore`).

### Packages in detail ###

1. [affyILM](/packages/2.6/bioc/html/affyILM.html) -
   Linear Model of background subtraction and the Langmuir isotherm
1. [BayesPeak](/packages/2.6/bioc/html/BayesPeak.html) -
   Bayesian Analysis of ChIP-seq Data
1. [BeadDataPackR](/packages/2.6/bioc/html/BeadDataPackR.html) -
   Compression of Illumina BeadArray data
1. [charm](/packages/2.6/bioc/html/charm.html) -
   Analysis of DNA methylation data from CHARM microarrays
1. [ConsensusClusterPlus](/packages/2.6/bioc/html/ConsensusClusterPlus.html) -
   Algorithm for determining cluster count and membership by stability evidence
   in unsupervised analysis
1. [CSAR](/packages/2.6/bioc/html/CSAR.html) -
   Statistical tools for the analysis of ChIP-seq data
1. [DESeq](/packages/2.6/bioc/html/DESeq.html) -
   Digital gene expresion analysis based on the negative binomial distribution
1. [eisa](/packages/2.6/bioc/html/eisa.html) -
   Expression data analysis via the Iterative Signature Algorithm
1. [ExpressionView](/packages/2.6/bioc/html/ExpressionView.html) -
   Visualize biclusters identified in gene expression data
1. [flowMeans](/packages/2.6/bioc/html/flowMeans.html) -
   Non-parametric Flow Cytometry Data Gating
1. [flowTrans](/packages/2.6/bioc/html/flowTrans.html) -
   Parameter Optimization for Flow Cytometry Data Transformation
1. [frma](/packages/2.6/bioc/html/frma.html) -
   Frozen RMA and Barcode
1. [frmaTools](/packages/2.6/bioc/html/frmaTools.html) -
   Frozen RMA Tools
1. [genoCN](/packages/2.6/bioc/html/genoCN.html) -
   genotyping and copy number study tools
1. [genomes](/packages/2.6/bioc/html/genomes.html) -
   Genome sequencing project metadata
1. [GenomicRanges](/packages/2.6/bioc/html/GenomicRanges.html) -
   Representation and manipulation of genomic intervals
1. [GEOsubmission](/packages/2.6/bioc/html/GEOsubmission.html) -
   Prepares microarray data for submission to GEO
1. [girafe](/packages/2.6/bioc/html/girafe.html) -
   Genome Intervals and Read Alignments for Functional Exploration
1. [goseq](/packages/2.6/bioc/html/goseq.html) -
   Gene Ontology analyser for RNA-seq and other length biased data
1. [GSRI](/packages/2.6/bioc/html/GSRI.html) -
   Gene Set Regulation Index
1. [hyperdraw](/packages/2.6/bioc/html/hyperdraw.html) -
   Visualizing Hypergaphs
1. [iChip](/packages/2.6/bioc/html/iChip.html) -
   Bayesian Modeling of ChIP-chip Data Through Hidden Ising Models
1. [iFlow](/packages/2.6/bioc/html/iFlow.html) -
   GUI based visualization for flow cytometry
1. [keggorthology (replaces keggortho)](/packages/2.6/bioc/html/keggorthology.html) -
   graph support for KO, KEGG Orthology
1. [MassArray](/packages/2.6/bioc/html/MassArray.html) -
   Analytical Tools for MassArray Data
1. [methVisual](/packages/2.6/bioc/html/methVisual.html) -
   Methods for visualization and statistics on DNA methylation data
1. [MotIV](/packages/2.6/bioc/html/MotIV.html) -
   Motif Identification and Validation
1. [PICS](/packages/2.6/bioc/html/PICS.html) -
   Probabilistic inference of ChIP-seq
1. [pint](/packages/2.6/bioc/html/pint.html) -
   Pairwise INTegration of functional genomics data
1. [PROMISE](/packages/2.6/bioc/html/PROMISE.html) -
   PRojection Onto the Most Interesting Statistical Evidence
1. [rGADEM](/packages/2.6/bioc/html/rGADEM.html) -
   De novo motif discovery
1. [Rsamtools](/packages/2.6/bioc/html/Rsamtools.html) -
   Import aligned BAM file format sequences into R / Bioconductor
1. [SamSPECTRAL](/packages/2.6/bioc/html/SamSPECTRAL.html) -
   Identifies cell population in flow cytometry data
1. [segmentSeq](/packages/2.6/bioc/html/segmentSeq.html) -
   Takes high-throughput sequencing data and uses it to define segments of the
   genome to which a high density of reads align
1. [SRAdb](/packages/2.6/bioc/html/SRAdb.html) -
   A compilation of metadata from NCBI SRA and tools
1. [tigre](/packages/2.6/bioc/html/tigre.html) -
   Transcription factor Inference through Gaussian process Reconstruction of
   Expression
1. [xmapcore](/packages/2.6/bioc/html/xmapcore.html) -
   Core access to the xmap database (installed separately)

## Additional Software Package Changes ##

`keggorth` has been renamed `keggorthology`.
