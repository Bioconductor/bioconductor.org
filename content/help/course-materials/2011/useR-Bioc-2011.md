Welcome!

* To prepare for this workshop, please install R version
  [2.13.1](http://cran.r-project.org).

* Exercises require additional packages. Install these by starting an R
  session and entering (e.g., cut-and-paste from below)

      source("http://bioconductor.org/biocLite.R")
      biocLite(c("GenomicFeatures", "ShortRead", "edgeR",
                 "org.Dm.eg.db", "BSgenome.Dmelanogaster.UCSC.dm3"))

* The workshop uses moderately sized data; your laptop (Mac, Windows,
  or Linux) should have 2 Gb or more of memory and ample disk space.

* If using Linux or MacOS, download the [Source
  package](useR2011_0.1.0.tar.gz) and install with

      install.packages("useR2011_0.1.0.tar.gz", repos=NULL,
          type="source")

  On windows, download the [Windows package](useR2011_0.1.0.zip) and
  install with

      install.packages("useR2011_0.1.0.tar.gz", repos=NULL)

  In either case, the first argument should be the path to the file
  that you downloaded!

* Once installed, load the library and read the extensive vignette

      library(useR2011)
      browseVignettes("useR2011")

  the package contains an extensive vignette (PDF) and R scripts for
  each chapter.

If you have questions, please do not hesitate to contact me, mtmorgan
at fhcrc dot org.
