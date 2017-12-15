# ![](/images/icons/help.gif)Bioconductor Workflows

Bioconductor provides software to help analyze diverse high-throughput
genomic data. Common workflows include:

<h2 id="basic">Basic</h2>

* __[Sequence Analysis](sequencing/)__  
  Import fasta, fastq, BAM, gff, bed, wig, and other sequence formats.
  Trim, transform, align, and manipulate sequences. Perform quality
  assessment, ChIP-seq, differential expression, RNA-seq, and other
  workflows.  Access the Sequence Read Archive.

* __[Oligonucleotide Arrays](arrays/)__  
  Import Affymetrix, Illumina, Nimblegen, Agilent, and other
  platforms.  Perform quality assessment, normalization, differential
  expression, clustering, classification, gene set enrichment,
  genetical genomics and other workflows for expression, exon, copy
  number, SNP, methylation and other assays.  Access GEO,
  ArrayExpress, Biomart, UCSC, and other community resources.

* __[Changing genomic coordinate systems with rtracklayer::liftOver](/help/workflows/liftOver/)__  
  The liftOver facilities developed in conjunction with the UCSC
  browser track infrastructure are available for transforming
  data in GRanges formats.  This is illustrated here with
  an image of the NHGRI GWAS catalog that is, as of Oct. 31 2014,
  distributed with coordinates defined by NCBI build hg38.


<h2 id="annotation">Annotations</h2>

* __[Annotation Resources](annotation/Annotation_Resources/)__  
  Introduction to using gene, pathway, gene ontology, homology annotations
  and the AnnotationHub. Access GO, KEGG, NCBI, Biomart, UCSC, vendor,
  and other sources.

* __[Annotating Genomic Ranges](annotation/Annotating_Genomic_Ranges/)__  
  Represent common sequence data types (e.g., from BAM, gff, bed, and
  wig files) as genomic ranges for simple and advanced range-based
  queries.

* __[Annotating Genomic Variants](variants/)__  
  Read and write VCF files. Identify structural location of variants
  and compute amino acid coding changes for non-synonymous
  variants. Use SIFT and PolyPhen database packages to predict
  consequence of amino acid coding changes.


<h2 id="rnaseq">RNA Sequencing</h2>

* __[RNA-Seq workflow: gene-level exploratory analysis and differential expression](/help/workflows/rnaseqGene/)__  
  This lab will walk you through an end-to-end RNA-Seq differential
  expression workflow, using DESeq2 along with other Bioconductor
  packages. We will start from the FASTQ files, show how these were
  aligned to the reference genome, prepare gene expression values
  as a count matrix by counting the sequenced fragments, perform
  exploratory data analysis (EDA), perform differential gene
  expression analysis with DESeq2, and visually explore the results.

* __[RNA-seq analysis is easy as 1-2-3](/help/workflows/RNAseq123/)__  
  This workflow demonstrates how to analyse RNA-sequencing data using the edgeR,
  limma and Glimma packages. The edgeR package is first used to import,
  organise, filter and normalise the data, followed by the limma package with
  its voom method, linear modelling and empirical Bayes moderation to assess
  differential expression and perform gene set testing. This pipeline is further
  enhanced by the Glimma package which enables interactive exploration of the
  results so that individual samples and genes can be examined by the user.

* __[Gene-level RNA-seq differential expression and pathway analysis](/help/workflows/RnaSeqGeneEdgeRQL/)__  
  Gene-level RNA-seq differential expression and pathway analysis using
  Rsubread and the edgeR quasi-likelihood pipeline

* __[recountWorkflow](/help/workflows/recountWorkflow/)__  
  The recount2 resource is composed of over 70,000 uniformly processed human
  RNA-seq samples spanning TCGA and SRA, including GTEx. The processed data can
  be accessed via the recount2 website and the recount Bioconductor
  package. This workflow explains in detail how to use the recount package and
  how to integrate it with other Bioconductor packages for several analyses that
  can be carried out with the recount2 resource. In particular, we describe how
  the coverage count matrices were computed in recount2 as well as different
  ways of obtaining public metadata, which can facilitate downstream
  analyses. Step-by-step directions show how to do a gene level differential
  expression analysis, visualize base-level genome coverage data, and perform an
  analyses at multiple feature levels. This workflow thus provides further
  information to understand the data in recount2 and a compendium of R code to
  use the data.

* __[EGSEA123](/help/workflows/EGSEA123/)__  
  R package that supports the F1000Research workflow article `Easy and efficient
  ensemble gene set testing with EGSEA', Alhamdoosh et al. (2017).


<h2 id="singlecell">Single-cell Workflows</h2>

* __Low-level analyses of single-cell RNA-sequencing data__  
  [Introduction](/help/workflows/simpleSingleCell/intro/) | [Part 1](/help/workflows/simpleSingleCell/part1/) | [Part 2](/help/workflows/simpleSingleCell/part2/) | [Part 3](/help/workflows/simpleSingleCell/part3/)  
  This workflow implements a low-level analysis pipeline for scRNA-seq
  data using scran, scater and other Bioconductor packages. It describes
  how to perform quality control on the libraries, normalization of
  cell-specific biases, basic data exploration and cell cycle phase
  identification. Procedures to detect highly variable genes,
  significantly correlated genes and subpopulation-specific marker genes
  are also shown. These analyses are demonstrated on a range of publicly
  available scRNA-seq data sets.

* __[CyTOF workflow: differential discovery in high-throughput high-dimensional
  cytometry datasets](/help/workflows/cytofWorkflow/)__  
  High dimensional mass and flow cytometry (HDCyto) experiments have become a
  method of choice for high throughput interrogation and characterization of
  cell populations. Here, we present an R-based pipeline for differential
  analyses of HDCyto data, largely based on Bioconductor packages. We
  computationally define cell populations using FlowSOM clustering, and
  facilitate an optional but reproducible strategy for manual merging of
  algorithm-generated clusters. Our workflow offers different analysis paths,
  including association of cell type abundance with a phenotype or changes in
  signaling markers within specific subpopulations, or differential analyses of
  aggregated signals. Importantly, the differential analyses we show are based
  on regression frameworks where the HDCyto data is the response; thus, we are
  able to model arbitrary experimental designs, such as those with batch
  effects, paired designs and so on. In particular, we apply generalized linear
  mixed models to analyses of cell population abundance or
  cell-population-specific analyses of signaling markers, allowing
  overdispersion in cell count or aggregated signals across samples to be
  appropriately modeled. To support the formal statistical analyses, we
  encourage exploratory data analysis at every step, including quality control
  (e.g. multi-dimensional scaling plots), reporting of clustering results
  (dimensionality reduction, heatmaps with dendrograms) and differential
  analyses (e.g. plots of aggregated signals).


<h2 id="variants">Genomic Variants</h2>

* __[Variant Calling](/help/course-materials/2014/BioC2014/Lawrence_Tutorial.pdf)__  
  This presentation illustrates a typical variant calling workflow starting
  with FASTQ data working through alignment, filtering, tallying, and calling.
  QC issues such as alignment coverage, mappability, and problematic
  homopolymers are explored. Called variants are exported as a vcf file and
  compared against published genotypes for concordance.  Final variants are
  annotated with with coding consequence and disease association.

* __[Nucleotide Tallies](/help/course-materials/2014/CSAMA2014/3_Wednesday/labs/Tutorial.pdf)__  
  Managing sequence data of large cohorts for population level analysis has
  become increasingly difficult with current file formats such as BAM, VCF,
  BCF, GTF, etc. Many studies work exclusively on the level of preprocessed
  variant calls stored in VCF/MAF file simply because there is no way to look
  at the data with reasonable resource usage. This tutorial presents an HDF5
  alternative for storing variant tallies from BAM files. This intermediate
  file format stores nucleotide tallies rather than alignments and provides
  efficient random access to cohort-level data. Once created, the tally files
  can be easily manipulated and used to create custom reports and plots.


<h2 id="advanced">Domain Specific</h2>

* __[High Throughput Assays](/help/workflows/highthroughputassays/)__  
  Import, transform, edit, analyze and visualize flow cytometric, mass
  spec, HTqPCR, cell-based, and other assays.

* __[Mass spectrometry and proteomics](/help/workflows/proteomics/)__  
  This lab demonstrates how to access data from proteomics data
  repositories, how to parse various mass spectrometry data formats, how
  to identify MS2 spectra and analyse the search results, how to use the
  high-level infrastructure for raw mass spectrometry and quantitative
  proteomics experiments and quantitative data processing and analysis.

* __[Transcription Factor Binding](/help/workflows/generegulation/)__  
  Finding Candidate Binding Sites for Known Transcription Factors via
  Sequence Matching.

* __[Cloud-enabled cis-eQTL search and annotation](/help/workflows/eQTL/)__  
  Bioconductor can be used to perform detailed analyses of
  relationships between DNA variants and mRNA abundance.  Genotype
  (potentially imputed) and expression data are organized in packages
  prior to analysis, using very concise representations.  SNP and
  probe filters can be specified at run time. Transcriptome-wide
  testing can be carried out using multiple levels of concurrency
  (chromosomes to nodes, genes to cores is a common approach).
  Default outputs of the cloud-oriented interface ciseqByCluster
  include FDR for all SNP-gene pairs in cis, along with locus-specific
  annotations of genetic and genomic contexts.

* __[Differential Binding from ChIP-seq data](/help/workflows/chipseqDB/)__  
  This workflow describes an analysis pipeline for de novo detection of
  differential binding (DB) from ChIP-seq data, from read alignment to
  interpretation of putative DB regions. It will be based on the use of sliding
  windows in the csaw package, with statistical modelling performed using
  methods in the edgeR package. Analyses will be demonstrated on real histone
  mark and transcription factor ChIP-seq data.

* __[Gene Expression Normalization Workflow](/help/workflows/ExpressionNormalizationWorkflow/)__  
  This workflow elucidates a customizable strategy to identify the effects of
  technical and confounding factors on gene expression data and normalize it
  while preserving the underlying biological features of interest. The example
  analysis demonstrated here explores how certain technical covariates 
  influence the interpretation of the impact of Coronary Artery Disease on
  peripheral blood gene expression.

* __[Methylation Array Analysis](/help/workflows/methylationArrayAnalysis/)__  
  Methylation in the human genome is known to be associated with development and
  disease. The Illumina Infinium methylation arrays are by far the most common
  way to interrogate methylation across the human genome. This Bioconductor
  workflow uses multiple packages for the analysis of methylation array
  data. Specifically, we demonstrate the steps involved in a typical
  differential methylation analysis pipeline including: quality control,
  filtering, normalization, data exploration and statistical testing for
  probe-wise differential methylation. We further outline other analyses such as
  differential methylation of regions, differential variability analysis,
  estimating cell type composition and gene ontology testing. Finally, we
  provide some examples of how to visualise methylation array data.

* __[TCGA Workflow: Analyze cancer genomics and epigenomics data](/help/workflows/TCGAWorkflow/)__  
  Biotechnological advances in sequencing have led to an explosion of publicly
  available data via large international consortia such as The Cancer Genome
  Atlas (TCGA), The Encyclopedia of DNA Elements (ENCODE), and The NIH Roadmap
  Epigenomics Mapping Consortium (Roadmap). These projects have provided
  unprecedented opportunities to interrogate the epigenome of cultured cancer
  cell lines as well as normal and tumor tissues with high genomic
  resolution. The Bioconductor project offers more than 1,000 open-source
  software and statistical packages to analyze high-throughput genomic
  data. However, most packages are designed for specific data types
  (e.g. expression, epigenetics, genomics) and there is no one comprehensive
  tool that provides a complete integrative analysis of the resources and data
  provided by all three public projects. A need to create an integration of
  these different analyses was recently proposed. In this workflow, we provide a
  series of biologically focused integrative analyses of different molecular
  data. We describe how to download, process and prepare TCGA data and by
  harnessing several key Bioconductor packages, we describe how to extract
  biologically meaningful genomic and epigenomic data. Using Roadmap and ENCODE
  data, we provide a work plan to identify biologically relevant functional
  epigenomic elements associated with cancer. To illustrate our workflow, we
  analyzed two types of brain tumors: low-grade glioma (LGG) versus high-grade
  glioma (glioblastoma multiform or GBM). 


<h2 id="Contribute">Contribute a Workflow</h2>

See the [HOWTO Creating Workflow Vignettes](/developers/how-to/workflows/)__  
for information on contributing your own workflow.
