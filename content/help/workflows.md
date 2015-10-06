# ![](/images/icons/help.gif)Bioconductor Workflows

Bioconductor provides software to help analyze diverse high-throughput
genomic data. Common workflows include:

<h2 id="basic">Basic Workflows</h2>

* [Sequence Analysis](sequencing/)
  Import fasta, fastq, BAM, gff, bed, wig, and other sequence formats.
  Trim, transform, align, and manipulate sequences. Perform quality
  assessment, ChIP-seq, differential expression, RNA-seq, and other
  workflows.  Access the Sequence Read Archive.

* [Oligonucleotide Arrays](arrays/)
  Import Affymetrix, Illumina, Nimblegen, Agilent, and other
  platforms.  Perform quality assessment, normalization, differential
  expression, clustering, classification, gene set enrichment,
  genetical genomics and other workflows for expression, exon, copy
  number, SNP, methylation and other assays.  Access GEO,
  ArrayExpress, Biomart, UCSC, and other community resources.

* [Annotation Resources](annotation/Annotation_Resources/)
  Introduction to using gene, pathway, gene ontology, homology annotations
  and the AnnotationHub. Access GO, KEGG, NCBI, Biomart, UCSC, vendor,
  and other sources.

* [Annotating Genomic Ranges](annotation/Annotating_Genomic_Ranges/)
  Represent common sequence data types (e.g., from BAM, gff, bed, and
  wig files) as genomic ranges for simple and advanced range-based
  queries.

* [Annotating Genomic Variants](variants/)
  Read and write VCF files. Identify structural location of variants
  and compute amino acid coding changes for non-synonymous
  variants. Use SIFT and PolyPhen database packages to predict
  consequence of amino acid coding changes.

* [Changing genomic coordinate systems with rtracklayer::liftOver](/help/workflows/liftOver/)
  The liftOver facilities developed in conjunction with the UCSC
  browser track infrastructure are available for transforming
  data in GRanges formats.  This is illustrated here with
  an image of the NHGRI GWAS catalog that is, as of Oct. 31 2014,
  distributed with coordinates defined by NCBI build hg38.

<h2 id="advanced">Advanced Workflows</h2>

* [High Throughput Assays](/help/workflows/highthroughputassays/)
  Import, transform, edit, analyze and visualize flow cytometric, mass
  spec, HTqPCR, cell-based, and other assays.

* [RNA-Seq workflow: gene-level exploratory analysis and differential expression](/help/workflows/rnaseqGene/)
  This lab will walk you through an end-to-end RNA-Seq differential 
  expression workflow, using DESeq2 along with other Bioconductor 
  packages. We will start from the FASTQ files, show how these were 
  aligned to the reference genome, prepare gene expression values 
  as a count matrix by counting the sequenced fragments, perform
  exploratory data analysis (EDA), perform differential gene 
  expression analysis with DESeq2, and visually explore the results.

* [Mass spectrometry and proteomics](/help/workflows/proteomics/)
  This lab demonstrates how to access data from proteomics data
  repositories, how to parse various mass spectrometry data formats, how
  to identify MS2 spectra and analyse the search results, how to use the
  high-level infrastructure for raw mass spectrometry and quantitative
  proteomics experiments and quantitative data processing and analysis.

* [Transcription Factor Binding](/help/workflows/generegulation/)
  Finding Candidate Binding Sites for Known Transcription Factors via
  Sequence Matching.

* [Cloud-enabled cis-eQTL search and annotation](/help/workflows/eQTL/)
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

* [Differential Binding from ChIP-seq data](/help/workflows/chipseqDB/)
  This workflow describes an analysis pipeline for de novo detection of
  differential binding (DB) from ChIP-seq data, from read alignment to
  interpretation of putative DB regions. It will be based on the use of sliding
  windows in the csaw package, with statistical modelling performed using
  methods in the edgeR package. Analyses will be demonstrated on real histone
  mark and transcription factor ChIP-seq data.

* [Variant Calling](/help/course-materials/2014/BioC2014/Lawrence_Tutorial.pdf)
  This presentation illustrates a typical variant calling workflow starting
  with FASTQ data working through alignment, filtering, tallying, and
  calling. QC issues are discussed such as alignment coverage and mappability,
  and problematic homopolymers. Called variants are exported as a vcf file and 
  compared against published genotypes for concordance and functionally 
  annotation with genomic content, coding consequence and disease association.

* [Nucleotide Tallys](/help/course-materials/2014/CSAMA2014/3_Wednesday/labs/Tutorial.pdf)
  Managing sequence data of large cohorts for population level analysis has
  become increasingly difficult with current file formats such as BAM, VCF,
  BCF, GTF, etc. Many studies work exclusively on the level of preprocessed
  variant calls stored in VCF/MAF file simply because there is no way to look
  at the data with reasonable resource usage. This tutorial presents an HDF5
  alternative for storing variant tallies from BAM files. This intermediate
  file format stores nucleotide tallies rather than alignments and provides
  efficient random access to cohort-level data. Once created, the tally files
  can be easily manipulated and used to create custom reports and plots.

<h2 id="Contribute">Contribute a Workflow</h2>

See the [HOWTO Creating Workflow Vignettes](/developers/how-to/workflows/) 
for information on contributing your own workflow.
