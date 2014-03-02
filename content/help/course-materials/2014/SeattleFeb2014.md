Download and install the package (containing all material) for use
with R-3.1.0 / Bioconductor 2.14.

* [BiocIntro_0.0.4.tar.gz](BiocIntro_0.0.4.tar.gz)

Install the course package with

    source("http://bioconductor.org/biocLite.R")
    dependencies <- c("Biostrings", "ShortRead", "ggplot2")
    biocLite(dependencies)
    install.packages("BiocIntro_0.0.3.tar.gz", repos=NULL)

Optionally, install suggested packages (used in exercises, etc) with

    source("http://bioconductor.org/biocLite.R")
    suggested <- c("BiocStyle", "knitr", "AnnotationHub",
        "BSgenome.Hsapiens.UCSC.hg19", "BiocParallel", "Biostrings",
        "GenomicAlignments", "GenomicFeatures", "GenomicRanges",
        "Gviz", "IRanges", "PSICQUIC", "RNAseqData.HNRNPC.bam.chr14",
        "TxDb.Hsapiens.UCSC.hg19.knownGene", "VariantAnnotation",
        "biomaRt", "knitr", "org.Hs.eg.db", "parallel", "rtracklayer")
    biocLite(suggested)

Explore the material through the following documents:

Introduction

* [pdf](Introduction.pdf), [Rnw](Introduction.Rnw) Abstract and tentative schedule 
* [pdf](Introduction_slides.pdf), [Rnw](Introduction_slides.Rnw) Slides 

Working with R

* [pdf](R_slides.pdf), [R](R_slides.R), [Rnw](R_slides.Rnw) Slides
* [pdf](R.pdf), [R](R.R), [Rnw](R.Rnw) Exercise

Sequencing work flows

* [pdf](Sequencing_workflows_slides.pdf), [Rnw](Sequencing_workflows_slides.Rnw) Slides 

Bioconductor for Sequence Analysis

* [pdf](Bioconductor_slides.pdf), [R](Bioconductor_slides.R), [Rnw](Bioconductor_slides.Rnw) Bioconductor - Slides 
* [pdf](Ranges_slides.pdf), [R](Ranges_slides.R), [Rnw](Ranges_slides.Rnw) Genomic Ranges - Slides 
* [pdf](Bioconductor_sequences.pdf), [R](Bioconductor_sequences.R), [Rnw](Bioconductor_sequences.Rnw) Working with DNA Sequences
* [pdf](Bioconductor.pdf), [R](Bioconductor.R), [Rnw](Bioconductor.Rnw) Working with FASTQ, BAM, and VCF files
* [pdf](Ranges.pdf), [R](Ranges.R), [Rnw](Ranges.Rnw) Working with Genomic Ranges

RNA-Seq

* [pdf](RNASeq.pdf), [R](RNASeq.R), [Rnw](RNASeq.Rnw) RNASeq Analysis 
* Exercises: see the [RNA-seq Differential Expression Lab](/help/course-materials/2013/EMBOBGI/)

Annotation and visualization

* [pdf](Annotation_slides.pdf), [R](Annotation_slides.R), [Rnw](Annotation_slides.Rnw) Slides 
* [pdf](Annotation.pdf), [R](Annotation.R), [Rnw](Annotation.Rnw) Working with Annotations
* [pdf](Visualization.pdf), [R](Visualization.R), [Rnw](Visualization.Rnw) Visualization of Genomic Data

