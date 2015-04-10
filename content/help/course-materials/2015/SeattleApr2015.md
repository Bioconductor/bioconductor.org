[Course Package](UseBioconductor_0.1.0.tar.gz)


<h4 id="ami">To launch an Amazon Machine Image (AMI) for this course:</h4>

* [Create an Amazon Web Services (AWS) Account](https://aws.amazon.com/) if you
  don't already have one.
* Start the instance <%= ami_url("ami-705b6718") %>; See the [documentation for this](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/launching-instance.html). Make sure your [security group](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/using-network-security.html) has
port 80 accessible.
* Paste the Public DNS name of your instance into a web browser.
* Log in to RStudio with username *ubuntu* and password *bioc* .
* Be sure and [terminate your instance](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/terminating-instances.html) when you are done using it, in order to avoid excessive charges.



## Use R / Bioconductor for Sequence Analysis

Fred Hutchinson Cancer Research Center, Seattle, WA<br />
6-7 April, 2015

Schedule
--------

Day 1 (9:00 - 12:30; 1:30 - 5:00)

- [A. Introduction](A_Introduction.html). _Bioconductor_ and
  sequencing work flows
- [B. Genomic Ranges](B_GenomicRanges.html). Working with Genomic
  Ranges and other _Bioconductor_ data structures (e.g., in the
  [GenomicRanges](http://bioconductor.org/packages/devel/bioc/html/GenomicRanges.html).
  package).  
- [C. Differential Gene Expression](C_DifferentialExpression.html). RNA-Seq
  known gene differential expression with
  [DESeq2](http://bioconductor.org/packages/devel/bioc/html/DESeq2.html)
  and
  [edgeR](http://bioconductor.org/packages/devel/bioc/html/edgeR.html).
- [D. Machine Learning](D_MachineLearning.html).
- [E. Gene Set Enrichment](E_GeneSetEnrichment.html).

Day 2 (9:00 - 12:30; 1:30 - 5:00)

- [F. ChIP-seq](F_ChIPSeq.html) ChIP-seq with
  [csaw](http://bioconductor.org/packages/devel/bioc/html/csaw.html)
- [G. Methylation](G_Methylation.html) and regulatory work flows with
  [minfi](http://bioconductor.org/packages/devel/bioc/html/minfi.html).
- [H. Integrative Data Analysis](H_IntegrativeAnalysis.html) -- emerging
  approaches
- [I. Large Data](I_LargeData.html) -- efficient, parallel, and cloud
  programming with
  [BiocParallel](http://bioconductor.org/packages/devel/bioc/html/BiocParallel.html),
  [GenomicFiles](http://bioconductor.org/packages/devel/bioc/html/GenomicFiles.html),
  and other resources.
