# ![](/images/icons/help.gif)Bioconductor Workflows

Bioconductor provides software to help analyze diverse high-throughput
genomic data. Common workflows include:

* [Sequence Analysis](high-throughput-sequencing/)  
  Import fasta, fastq, BAM, gff, bed, wig, and other sequence formats.
  Trim, transform, align, and manipulate sequences. Perform quality
  assessment, ChIP-seq, differential expression, RNA-seq, and other
  workflows.  Access the Sequence Read Archive.

* [RNAseq Differential Expression](/packages/release/data/experiment/html/parathyroidSE.html)  
  Use the <em>parathyroidSE</em> ExperimentData package and vignette
  to learn how to count reads and perform other common operations
  required for differential expression analysis.

* [Oligonucleotide Arrays](arrays/)  
  Import Affymetrix, Illumina, Nimblegen, Agilent, and other
  platforms.  Perform quality assessment, normalization, differential
  expression, clustering, classification, gene set enrichment,
  genetical genomics and other workflows for expression, exon, copy
  number, SNP, methylation and other assays.  Access GEO,
  ArrayExpress, Biomart, UCSC, and other community resources.

* [Variants](variants/)  
  Read and write VCF files. Identify structural location of variants
  and compute amino acid coding changes for non-synonymous
  variants. Use SIFT and PolyPhen database packages to predict
  consequence of amino acid coding changes.

* [Accessing Annotation Data](annotation/annotation/)  
  Use microarray probe, gene, pathway, gene ontology, homology and
  other annotations.  Access GO, KEGG, NCBI, Biomart, UCSC, vendor,
  and other sources.

* [Annotating Ranges](annotation/AnnotatingRanges)  
  Represent common sequence data types (e.g., from BAM, gff, bed, and
  wig files) as genomic ranges for simple and advanced range-based
  queries.

* [High Throughput Assays](/help/workflows/highthroughputassays/)  
  Import, transform, edit, analyze and visualize flow cytometric, mass
  spec, HTqPCR, cell-based, and other assays.

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

* [RNA-Seq workflow: gene-level exploratory analysis and differential expression](/help/workflows/rnaseqGene/)
  This lab will walk you through an end-to-end RNA-Seq differential expression workflow, using DESeq2 along with other Bioconductor packages. We will start from the FASTQ files, show how these were aligned to the reference genome, prepare gene expression values as a count matrix by counting the sequenced fragments, perform exploratory data analysis (EDA), perform differential gene expression analysis with DESeq2, and visually explore the results.

* [Changing genomic coordinate systems with rtracklayer::liftOver](/help/workflows/liftOver/)
  The liftOver facilities developed in conjunction with the UCSC
  browser track infrastructure are available for transforming
  data in GRanges formats.  This is illustrated here with
  an image of the NHGRI GWAS catalog that is, as of Oct. 31 2014,
  distributed with coordinates defined by NCBI build hg38.

* [Mass spectrometry and proteomics](/help/workflows/proteomics/)  
  This lab demonstrates how to accessing data from proteomics data
  repositories, how to parse various mass spectrometry data formats, how
  to identify MS2 spectra and analyse the search results, how to use the
  high-level infrastructure for raw mass spectrometry and quantitative
  proteomics experiments and quantitative data processing and analysis.

See the [HOWTO](/developers/how-to/workflows/) for information on contributing
your own workflow.

