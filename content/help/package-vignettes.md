# ![](/images/icons/help.gif)Package Vignettes #

Each Bioconductor package contains at least one vignette, a document
that provides a task-oriented description of package
functionality. Vignettes contain executable examples and are intended
to be used interactively.

You can access the PDF version of a vignette for any installed package
from inside R using the `browseVignettes` function.  For example, to
view the vignettes in the Biostrings package enter the following at an
R prompt:

    browseVignettes(package = "Biostrings")

This will open a web browser with links to the vignette PDF as well as
a plain-text R file containing the code used in the vignette.
    
The vignette files, both the PDF and the Rnw sources document, are
located in the `doc` directory of an installed package
(inst/doc for an uninstalled package tarball). You can discover the
location of an installed package as follows:

    system.file(package = "Biostrings")

You can extract all of the R code from a vignette source file (these
files are in Sweave format) as shown below. The result will be a `.R`
file for each vignette source file written to your current working
directory. You can use these R files to cut and paste into an R
session so that you can work the examples in the vignette yourself.

    library("tools")
    vigSrc = list.files(pattern = "Rnw$",
                        system.file("doc", package = "Biostrings"),
                        full.names = TRUE)
    vigSrc
    for (v in vigSrc) Stangle(v)

You can also download the vignette for a package that isn't installed
on your system by visiting the package's description page linked from
the [list of Bioconductor packages](/packages/release/bioc/)

*Term of use for the vignettes*: You are welcome to use these
materials for instructional purposes. However, you may not include
these in separately published works (articles, books, websites). When
using all or parts of the Bioconductor course materials (slides,
vignettes, scripts), we would appreciate it if you would cite the
authors and refer your audience to the Bioconductor website.
