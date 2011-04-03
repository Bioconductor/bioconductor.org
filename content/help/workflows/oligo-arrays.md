![](/images/icons/help.gif)Using Bioconductor for Microarray Analysis
=====================================================================

Bioconductor has advanced facilities for analysis of microarray
platforms including Affymetrix, Illumina, Nimblegen, Agilent, and
other one- and two-color technologies.

Bioconductor includes extensive support for analysis of expression
arrays, and well-developed support for exon, copy number, SNP,
methylation, and other assays.

Major workflows in Bioconductor include pre-processing, quality
assessment, differential expression, clustering and classification,
gene set enrichment analysis, and genetical genomics.

Bioconductor offers extensive interfaces to community resources,
including GEO, ArrayExpress, Biomart, genome browsers, GO, KEGG, and
diverse annotation sources.

* [Sample Workflow](#sample-workflow)  
* [Installation and Use](#install-and-use)
* [Exploring Package Content](#exploring-package-content)
* [Pre-Processing Resources](#pre-processing-resources)

<h2 id="sample-workflow"> Sample Workflow</h2>

The following psuedo-code illustrates a typical R / Bioconductor
session. It uses RMA from the `affy` package to pre-process Affymetrix
arrays, and the `limma` package for assessing differential expression.

    ## Load packages
    > library(affy)   # Affymetrix pre-processing
    > library(limma)  # two-color pre-processing; differential
                      # expression
                    
    ## import "phenotype" data, describing the experimental design
    > phenoData <- read.AnnotatedDataFrame("sample-description.csv")
    
    ## RMA normalization
    > eset <- justRMA("/celfile-directory", phenoData=phenoData)
    
    ## differential expression
    > design <-                   # describe model to be fit
          model.matrix(~ Disease, pData(eset))
    > fit <- lmFit(eset, design)  # fit each probeset to model
    > efit <- eBayes(fit)        # empirical Bayes adjustment
    > topTable(efit, coef=2)      # table of differentially expressed probesets
    
A top table resulting from a more complete analysis, described in
Chapter 7 of [Bioconductor Case
Studies](/help/publications/books/bioconductor-case-studies/), 
is shown below. The table enumerates Affymetrix probes, the log-fold
difference between two experimental groups, the average expression
across all samples, the t-statistic describing differential
expression, the unadjusted and adjusted (controlling for false
discovery rate, in this case) significance of the difference, and
log-odds ratio. These results can be used in further analysis and
annotation.

          ID logFC AveExpr    t  P.Value adj.P.Val     B
    636_g_at  1.10    9.20 9.03 4.88e-14  1.23e-10 21.29
    39730_at  1.15    9.00 8.59 3.88e-13  4.89e-10 19.34
     1635_at  1.20    7.90 7.34 1.23e-10  1.03e-07 13.91
     1674_at  1.43    5.00 7.05 4.55e-10  2.87e-07 12.67
    40504_at  1.18    4.24 6.66 2.57e-09  1.30e-06 11.03
    40202_at  1.78    8.62 6.39 8.62e-09  3.63e-06  9.89
    37015_at  1.03    4.33 6.24 1.66e-08  6.00e-06  9.27
    32434_at  1.68    4.47 5.97 5.38e-08  1.70e-05  8.16
    37027_at  1.35    8.44 5.81 1.10e-07  3.08e-05  7.49
    37403_at  1.12    5.09 5.48 4.27e-07  1.08e-04  6.21
   

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="install-and-use">Installation and Use</h2>

Follow [installation instructions](/install/) to start using these
packages.  The `affy` and `limma` packages are part of the core
Bioconductor packages, and are installed automatically with

    > source("http://bioconductor.org/biocLite.R")
    > biocLite()

To install additional packages, such as the annotations associated
with the Affymetrix Human Genome U95A 2.0, use

    > source("http://bioconductor.org/biocLite.R")
    > biocLite("hgu95av2.db")

Package installation is required only once per R installation. View a
/packagesfull list of [available packages](/packages/release/bioc/).

To use the `affy` and `limma` packages, evaluate the commands

    > library("affy")
    > library("limma")

These commands are required once in each R session.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="exploring-package-content">Exploring Package Content</h2>

Packages have extensive help pages, and include vignettes highlighting
common use cases. The help pages and vignettes are available from
within R. After loading a package, use syntax like

    > help(package="limma")
    > ?topTable

to obtain an overview of help on the `limma` package, and the
`topTable` function, and

    > browseVignettes(package="limma")

to view vignettes (providing a more comprehensive introduction to
package functionality) in the `limma` package. Use

    > help.start()

to open a web page containing comprehensive help resources.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="pre-processing-resources">Pre-Processing Resources</h2>

The following provide a brief overview of packages useful for
pre-processing. More comprehensive workflows can be found in
documentation (available from
[package descriptions](/packages/release/bioc/))
and in Bioconductor [Books and monographs](/help/publications/).

### Affymetrix 3'-biased Arrays ###

[affy](/packages/release/bioc/html/affy.html),
[gcrma](/packages/release/bioc/html/gcrma.html),
[affyPLM](/packages/release/bioc/html/affyPLM.html)

* Require cdf package, probe package and annotation package
* All these packages are available from Bioconductor via `biocLite()`

[xps](/packages/release/bioc/html/xps.html)

* Requires installation of [ROOT](http://root.cern.ch/)
* Uses data files from Affymetrix (.CDF, .PGF, .CLF, .CSV) directly

### Affymetrix Exon ST Arrays ###

[oligo](/packages/release/bioc/html/oligo.html)

* Requires a pdInfoPackage built using
  [pdInfoBuilder](/packages/release/bioc/html/pdInfoBuilder.html)
* This package collates cdf, probe, annotation data together
* These packages are available from Bioconductor via `biocLite()`
* Most cases will require a 64-bit computer running Linux and &gt;= 8Gb RAM

[exonmap](/packages/release/bioc/html/exonmap.html)

* Requires installation of MySQL and Ensembl core database tables
* Requires specially modified
  [cdf](http://xmap.picr.man.ac.uk/download/) and affy package
* Requires a 64-bit computer running Linux and &gt;= 8 Gb RAM

[xps](/packages/release/bioc/html/xps.html)

* Requires installation of [ROOT](http://root.cern.ch/)
* Uses data files from Affymetrix (.CDF, .PGF, .CLF, .CSV) directly
* Will run on conventional desktop computers

### Affymetrix Gene ST Arrays ###

[oligo](/packages/release/bioc/html/oligo.html)

* Requires a pdInfoPackage built using
  [pdInfoBuilder](/packages/release/bioc/html/pdInfoBuilder.html)
* This package collates cdf, probe, annotation data together
* These packages are available from Bioconductor via `biocLite()`

[xps](/packages/release/bioc/html/xps.html)

* Requires installation of [ROOT](http://root.cern.ch/)
* Uses data files from Affymetrix (.CDF, .PGF, .CLF, .CSV) directly

### Affymetrix SNP Arrays ###

[oligo](/packages/release/bioc/html/oligo.html)

* Requires a pdInfoPackage built using
  [pdInfoBuilder](/packages/release/bioc/html/pdInfoBuilder.html)
* This package collates cdf, probe, annotation and HapMap data
* These packages are available from Bioconductor via `biocLite()`
* Not yet capable of processing CNV regions in SNP5.0 and SNP6.0 </ul>

### Affymetrix Tiling Arrays ###

[oligo](/packages/release/bioc/html/oligo.html)

* Requires a pdInfoPackage built using
  [pdInfoBuilder](/packages/release/bioc/html/pdInfoBuilder.html)
* This package collates data from bpmap and cif files

### Nimblegen Arrays ###

[oligo](/packages/release/bioc/html/oligo.html)

* Requires a `pdInfoPackage` built using
  [pdInfoBuilder](/packages/release/bioc/html/pdInfoBuilder.html)

### Illumina Expression Microarrays ###

[lumi](/packages/release/bioc/html/lumi.html)

* Requires lumi-specific mapping and annotation packages (e.g.,
  `lumiHumanAll.db` and `lumiHumanIDMapping`)

[beadarray](/packages/release/bioc/html/beadarray.html)

* Requires beadarray-specific mapping and annotation packages (e.g.,
  `illuminaHumanv1BeadID.db` and `illuminaHumanV1.db`)

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>
