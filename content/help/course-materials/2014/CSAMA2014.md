## Monday, June 22

Morning talks

- [pdf](1_Monday/lectures/140623-morgan-introduction.pdf) Introduction
  to R and Bioconductor
- [pdf](1_Monday/lectures/HTS_intro__Anders.pdf) Basics of
  high-throughput sequencing technologies and short read aligners
- [pdf](1_Monday/lectures/140623-gentleman-statistics-intro.pdf)
  Elements of statistics 1: t-test and linear model
- [pdf](1_Monday/lectures/140623-brixen-huber-multtestindepfilt.pdf)
  Elements of statistics 2: multiple testing, false discovery rates,
  independent filtering

Afternoon labs

- [zip](1_Monday/labs/R-Basics.zip) R introduction/refresher: data
  types, reading and writing files and spreadsheets, plotting,
  programming, functions and packages.
- [pdf](1_Monday/labs/vizlab14.pdf) [R](1_Monday/labs/vizlab14.R)
  Exploratory data analysis and visualization
  ([pdf](1_Monday/labs/solutionsToVizlab.pdf) solutions)
- [html](1_Monday/labs/IntermediateR.html)
  [R](1_Monday/labs/IntermediateR.R) Intermediate R 1: accessing
  resources - packages, classes, methods, and efficient code. Download
  [IntermediateR1_1.0.0.tar.gz](1_Monday/labs/IntermediateR1_1.0.0.tar.gz)
  and install as:

      source("http://bioconductor.org/biocLite.R")
      biocLite(c("IRanges", "GenomicRanges", "microbenchmark"))
      install.packages("IntermediateR1_1.0.0.tar.gz", repos=NULL, type="source")
    
- [pdf](1_Monday/labs/ScalableComputing.pdf)
  [R](1_Monday/labs/ScalableComputing.pdf) Intermediate R 2: scalable
  / performant computing. (Large files needed for some of this lab are
  NOT available for download). Download
  [CSAMA2014ScalableComputingLab_0.0.1.tar.gz](1_Monday/labs/CSAMA2014ScalableComputingLab_0.0.1.tar.gz)
  and install as:

      source("http://bioconductor.org/biocLite.R")
      biocLite(c("IRanges", "GenomicRanges", "Rsamtools", "ShortRead", 
          "rtracklayer", "GenomicAlignments", "GEOquery", "microbenchmark",
          "BiocParallel", "ggbio", "Biobase", "GenomicFiles"))
      install.packages("CSAMA2014ScalableComputingLab_0.0.1.tar.gz",
          repos=NULL, type="source")

## Tuesday

Morning talks

- [pdf](2_Tuesday/lectures/DESeq2-Anders.pdf) RNA-Seq 1: differential
  expression analysis - GLMs and testing
- RNA-Seq 2: shrinkage, empirical Bayes, FC estimation
- [pdf](2_Tuesday/lectures/Visualization_in_Statistical_Genomics-Carey.pdf)
  Visualisation
- [pdf](2_Tuesday/lectures/Ranges_Sequences_and_Alignments-Lawrence.pdf)
  Computing with genomic ranges, sequences and alignments

Afternoon labs

- [pdf](2_Tuesday/labs/RNA-Seq-Analysis-Lab.pdf)
  [R](2_Tuesday/labs/RNA-Seq-Analysis-Lab.R)
  ([Rnw](2_Tuesday/labs/RNA-Seq-Analysis-Lab.Rnw),
  [bib](2_Tuesday/labs/RNA-Seq-Analysis-Lab.bib)) A complete RNA-Seq
  differential expression workflow
  [DESeq2_result_table.RData](2_Tuesday/labs/DESeq2_result_table.RData)
  [Homo_sapiens.GRCh37.75.subset.gtf.gz](2_Tuesday/labs/Homo_sapiens.GRCh37.75.subset.gtf.gz)

## Wednesday

Morning talks

- [pdf](3_Wednesday/lectures/VariantCallingLecture.pdf) DNA-Seq 1:
  Variant calling
- [pdf](3_Wednesday/lectures/CSAMA2014.variant.visualisation.and.qc.pdf)
  DNA-Seq 2: visualisation and quality assessment of variant calls
- [pdf](3_Wednesday/lectures/GSEA2014.pdf) Gene set enrichment analysis
- [pdf](3_Wednesday/lectures/Annotation_slides.pdf)
  [R](3_Wednesday/lectures/Annotation_slides_script.R) Working with
  gene and genome annotations

Afternoon labs

- [pdf](3_Wednesday/labs/Tutorial.pdf)
  [R](3_Wednesday/labs/Tutorial.R) Variant tallies, visualisation,
  HDF5 [ExampleData.zip](3_Wednesday/labs/ExampleData.zip)
  [NRAS.tally.hfs5](3_Wednesday/labs/NRAS.tally.hfs5)

## Thursday

Morning talks

- [pdf](4_Thursday/lectures/isoforms_dexseq_talk.pdf) RNA-Seq 3:
  alternative exon usage
- [html](4_Thursday/lectures/thursClust.html) Elements of statistics
  3: Classification and clustering - basic concepts
- [pdf](4_Thursday/lectures/140626-brixen-machinelearn-huber.pdf)
  Elements of statistics 4: regularisation & kernels
- [pdf](4_Thursday/lectures/ChIPSeq_slides.pdf)
  [R](4_Thursday/lectures/ChIPSeq_slides_script.R) ChIP-Seq

Afternoon labs

- [pdf](4_Thursday/labs/Annotation_lab.pdf)
  [R](4_Thursday/labs/Annotation_lab.R) Working with the Ranges
  infrastructure: annotating and understanding regions

## Friday

Morning talks

- [pdf](5_Friday/lectures/ExpDesign_Anders.pdf) Elements of statistics
  5: experimental design
- [pdf](5_Friday/lectures/vjcEQTLfri.pdf) eQTL / molecular-QTL
  analyses
- [pdf](5_Friday/lectures/proteomics_slides.pdf) Proteomics
- Emerging topic -- [pdf](5_Friday/lectures/ImageData.pdf) image
  analysis

Afternoon labs

- Reporting your analysis - authoring knitr/Rmarkdown, ReportingTools,
  shiny. [rauthoring_0.2.2.tar.gz](5_Friday/labs/rauthoring_0.2.2.tar.gz)
