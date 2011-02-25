* [Install Bioconductor Packages](#install-bioconductor-packages)  
* [Find Bioconductor Packages](#find-bioconductor-packages)  
* [Update Bioconductor Packages](#update-bioconductor-packages)
* [Install R](#install-R)  

<h2 id="install-bioconductor-packages">Install Bioconductor Packages</h2>

Use the `biocLite.R` script to install a selection of core
Bioconductor packages, or to install individual packages. To install
core packages, type the following in an R command window:

    source("http://bioconductor.org/biocLite.R")
    biocLite()

This installs the following packages and their dependencies: `affy`,
`affydata`, `affyPLM`, `affyQCReport`, `annaffy`, `annotate`,
`Biobase`, `biomaRt`, `Biostrings`, `DynDoc`, `gcrma`, `genefilter`,
`geneplotter`, `GenomicRanges`, `hgu95av2.db`, `limma`, `marray`,
`multtest`, `vsn`, and `xtable`.  After downloading and installing
these packages, the script prints "Installation complete" and "TRUE".

The `biocLite.R` script has arguments that change its default behavior:

    pkgs
        Character vector of Bioconductor packages to install.
    destdir
        File system directory for downloaded packages.
    lib
        R library where packages are installed.

[biocViews](/help/bioc-views/release/BiocViews.html) and the list of
[software packages](/help/bioc-views/release/bioc/) are ways to
discover Bioconductor packages.  To install a particular package,
e.g., EBImage, use

    source("http://bioconductor.org/biocLite.R")
    biocLite("EBImage")

Install "pkg1" and "pkg2" with

    biocLite(c("pkg1", "pkg2"))

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>


<h2 id="find-bioconductor-packages">Find Bioconductor Packages</h2>

Visit the [Workflows](/help/workflows/) help page and
[BiocViews](/help/bioc-views/<%=config[:release_version]%>/BiocViews.html)
taxonomy to discover available software.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="update-bioconductor-packages">Update Installed Bioconductor Packages</h2>

Bioconductor packages, especially those in the development branch, are
updated fairly regularly. To identify packages requiring update, start
a new session of R and enter

    source("http://bioconductor.org/biocLite.R")
    old.packages(repos=biocinstallRepos())

To update all installed packages that are out of date, start a new
session of R and enter

    source("http://bioconductor.org/biocLite.R")
    update.packages(repos=biocinstallRepos(), ask=FALSE, checkBuilt=TRUE)

Read the help page for `update.packages` for additional details.
Recompiling installed Bioconductor packages

Rarely, underlying changes in the operating system require ALL
installed packages to be recompiled for source (C or Fortran)
compatibility. One way to address this might be to start a new R
session and enter

    source("http://bioconductor.org/biocLite.R")
    pkgs <- rownames(installed.packages())
    biocLite(pkgs)

As this will reinstall all currently installed packages, it likely
involves a significant amount of network bandwidth and compilation
time. All packages are implicitly updated, and the cumulative effect
might introduce wrinkles that disrupt your workflow.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="install-R">Install R</h2>

1. Download the most recent version of [R][].  The [R FAQs][] and the [R
Installation and Administration Manual][1] contain detailed instructions
for installing R on various platforms (Linux, OS X, and Windows being
the main ones).

[R]: http://www.r-project.org/
[R FAQs]: http://cran.r-project.org/faqs.html
[1]: http://cran.r-project.org/doc/manuals/R-admin.html

2. Start the R program; on Windows and OS X, this will usually mean
   double-clicking on the R application, on UNIX-like systems, type
   "R" at a shell prompt.
   
3. As a first step with R, start the R help browser by typing
   `help.start()` in the R command window. For help on any
   function, e.g. the "mean" function, type `? mean`.

