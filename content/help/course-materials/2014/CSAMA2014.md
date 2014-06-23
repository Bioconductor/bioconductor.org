## Monday, June 22

Morning talks

- [Introduction to R and Bioconductor](Introduction.pdf)
- Basics of high-throughput sequencing technologies and short read
  aligners
- Elements of statistics 1: t-test and linear model
- [Elements of statistics 2](140623-brixen-huber-multtestindepfilt.pdf):
  multiple testing, false discovery rates, independent filtering

Afternoon labs

- R introduction/refresher: data types, reading and writing files and
  spreadsheets, plotting, programming, functions and packages.
- Exploratory data analysis and visualization
- [Intermediate R 1](IntermediateR.html) ([R](IntermediateR.R)):
  accessing resources - packages, classes, methods, and efficient
  code. Download
  [IntermediateR1_1.0.0.tar.gz](IntermediateR1_1.0.0.tar.gz) and
  install as:

      source("http://bioconductor.org/biocLite.R")
      biocLite(c("IRanges", "GenomicRanges", "microbenchmark"))
      install.packages(IntermediateR1_1.0.0.tar.gz", repos=NULL, type="source")
    
- [Intermediate R 2](CSAMA2014ScalableComputingLab.pdf): scalable /
  performant computing. (Large files needed for some of this lab are
  NOT available for download). Download
  [CSAMA2014ScalableComputingLab_0.0.1.tar.gz](CSAMA2014ScalableComputingLab_0.0.1.tar.gz)  and install as:

      source("http://bioconductor.org/biocLite.R")
      biocLite(c("IRanges", "GenomicRanges", "Rsamtools", "ShortRead", 
          "rtracklayer", "GenomicAlignments", "GEOquery", "microbenchmark",
          "BiocParallel", "ggbio", "Biobase", "GenomicFiles"))
      install.packages(CSAMA2014ScalableComputingLab_0.0.1.tar.gz",
          repos=NULL, type="source")

## Tuesday

Morning talks

- RNA-Seq 1: differential expression analysis - GLMs and testing
- RNA-Seq 2: shrinkage, empirical Bayes, FC estimation
- Visualisation
- Computing with genomic ranges, sequences and alignments

Afternoon labs

- A complete RNA-Seq differential expression workflow

## Wednesday

Morning talks

- DNA-Seq 1: Variant calling
- DNA-Seq 2: visualisation and quality assessment of variant calls
- Gene set enrichment analysis
- Leveraging annotation and data integration: working with gene and
  genome annotations

## Thursday

Morning talks

- RNA-Seq 3: alternative exon usage
- Elements of statistics 3: Classification and clustering - basic
  concepts
- Elements of statistics 4: regularisation & kernels
- ChIP-Seq

Afternoon labs

- Working with the Ranges infrastructure: annotating and understanding
  regions
- Variant tallies, visualisation, HDF5

## Friday

Morning talks

- Elements of statistics 5: experimental design
- eQTL / molecular-QTL analyses
- Proteomics
- Emerging topic

Afternoon labs

- Reporting your analysis - authoring knitr/Rmarkdown, ReportingTools, shiny
