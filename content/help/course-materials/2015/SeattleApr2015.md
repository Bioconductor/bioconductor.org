## Use R / Bioconductor for Sequence Analysis
==========================================

Fred Hutchinson Cancer Research Center, Seattle, WA<br />
6-7 April, 2015

Contact: Martin Morgan
([mtmorgan@fredhutch.org](mailto:mtmorgan@fredhutch.org))

This **INTERMEDIATE** course is designed for individuals comfortable
using _R_, and with some familiarity with _Bioconductor_. It consists
of approximately equal parts lecture and practical sessions addressing
use of _Bioconductor_ software for analysis and comprehension of
high-throughput sequence and related data. Specific topics include use
of central Bioconductor classes (e.g., _GRanges_,
_SummarizedExperiment_), RNASeq gene differential expression, ChIP-seq
and methylation work flows, approaches to management and integrative
analysis of diverse high-throughput data types, and strategies for
working with large data. Participants are required to bring a laptop
with wireless internet access and a modern version of the Chrome or
Safari web browser.


Schedule (tentative)
--------------------

Day 1 (9:00 - 12:30; 1:30 - 5:00)

- [A. Introduction](vignettes/A_Introduction.Rmd). _Bioconductor_ and
  sequencing work flows
- [B. Genomic Ranges](vignettes/B_GenomicRanges.Rmd). Working with Genomic
  Ranges and other _Bioconductor_ data structures (e.g., in the
  [GenomicRanges](http://bioconductor.org/packages/devel/bioc/html/GenomicRanges.html).
  package).  
- [C. Differential Gene Expression](vignettes/C_DifferentialExpression.Rmd). RNA-Seq
  known gene differential expression with
  [DESeq2](http://bioconductor.org/packages/devel/bioc/html/DESeq2.html)
  and
  [edgeR](http://bioconductor.org/packages/devel/bioc/html/edgeR.html).
- [D. Machine Learning](vignettes/D_MachineLearning.Rmd).
- [E. Gene Set Enrichment](vignettes/E_GeneSetEnrichment.Rmd).

Day 2 (9:00 - 12:30; 1:30 - 5:00)

- [F. ChIP-seq](vignettes/F_ChIPSeq.Rmd) ChIP-seq with
  [csaw](http://bioconductor.org/packages/devel/bioc/html/csaw.html)
- [G. Methylation](vignettes/G_Methylation.Rmd) and regulatory work flows with
  [minfi](http://bioconductor.org/packages/devel/bioc/html/minfi.html).
- [H. Integrative Data Analysis](vignettes/H_IntegrativeAnalysis.Rmd) -- emerging
  approaches
- [I. Large Data](vignettes/I_LargeData.Rmd) -- efficient, parallel, and cloud
  programming with
  [BiocParallel](http://bioconductor.org/packages/devel/bioc/html/BiocParallel.html),
  [GenomicFiles](http://bioconductor.org/packages/devel/bioc/html/GenomicFiles.html),
  and other resources.
