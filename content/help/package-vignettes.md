# Package Vignettes

Each Bioconductor package contains at least one vignette, a document
that provides a task-oriented description of package
functionality. Vignettes contain executable examples and are intended
to be used interactively.

You can access the compiled version of a vignette for any installed package
from inside R using the `browseVignettes` function. For example, to
view the vignettes in the Biostrings package enter the following at an
R prompt:

```r
browseVignettes(package = "Biostrings")
```

This will open a web browser with links to the compiled vignette, as well as
a plain-text R file containing the code used in the vignette.

The vignette files, both the compiled and the sources documents,
are located in the `doc` directory of an installed package
(inst/doc for an uninstalled package tarball). You can discover the
location of an installed package as follows:

```r
system.file(package = "Biostrings")
```

Older vignettes are mostly written in Sweave format with extension Rnw while
newer vignttes mostly take advantage of Rmarkdown with extension Rmd. You can
extract all of the R code from a vignette source file from a package. The result
will be a `.R` file for each vignette source file written to your current
working directory. You can use these R files to cut and paste into an R
session so that you can work the examples in the vignette yourself.

The steps to extract R code depend if the source code is a Rnw Sweave format or
Rmd Rmarkdown format. These commands would be run at an R prompt.

For Sweave we use <code>tools::Stangle</code>:

```r
library("tools")
vigSrc = list.files(pattern = "Rnw$",
                    system.file("doc", package = "Biostrings"),
                    full.names = TRUE)
vigSrc
for (v in vigSrc) Stangle(v)
```

For Rmarkdown we use <code>knitr::purl</code>:

```r
library("knitr")
vigSrc = list.files(pattern = "Rmd$",
                    system.file("doc", package = "BiocFileCache"),
                    full.names = TRUE)
vigSrc
for (v in vigSrc) purl(v)
```

You can also view or download the vignette for a package that isn't installed
on your system by visiting the package's description page linked from
the [list of Bioconductor packages](/packages/release/bioc/)

_Term of use for the vignettes_: You are welcome to use these
materials for instructional purposes. However, you may not include
these in separately published works (articles, books, websites). When
using all or parts of the Bioconductor course materials (slides,
vignettes, scripts), we would appreciate it if you would cite the
authors and refer your audience to the Bioconductor website.

Please also remember to check package documention help files. You can access
help files for a package within R with <code>?</code>.

```r
library(Biostrings)
?DNAString
?reverseComplement
```