![](/images/icons/install.gif)Install
===========================================================

* [Install R](#install-R)  
* [Find Bioconductor Packages](#find-bioconductor-packages)  
* [Install Bioconductor Packages](#install-bioconductor-packages)  
* [Update Bioconductor Packages](#update-bioconductor-packages)


<h2 id="install-R"">Install R</h2>

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


<h2 id="find-bioconductor-packages">Find Bioconductor Packages</h2>

Find a package

<h2 id="install-bioconductor-packages">Install Base Bioconductor Packages</h2>


Install BioConductor packages using the biocLite.R installation
script. In an R command window, type the following:

<p class="code_box">
    source("http://bioconductor.org/biocLite.R")<br/>
    biocLite()
</p>

This installs the following packages: affy, affydata, affyPLM,
annaffy, annotate, Biobase, Biostrings, DynDoc, gcrma, genefilter,
geneplotter, hgu95av2.db, limma, marray, matchprobes, multtest, ROC,
vsn, xtable, affyQCReport. After downloading and installing these
packages, the script prints "Installation complete" and TRUE.

The biocLite script has arguments that change the default behavior:

<p class="code_box">
    <b>pkgs</b><br/>
    &nbsp&nbsp&nbsp Character vector of BioConductor packages to install.<br/>
    <b>destdir</b><br/>
    &nbsp&nbsp&nbsp File system directory for downloaded packages.<br/>
    <b>lib</b><br/>
    &nbsp&nbsp&nbsp R library where packages are installed.<br/>
</p>

<h2>Install additional Bioconductor packages</h2>

There are many Bioconductor and R packages in addition to those in the
default installation of biocLite. A catalog of the Bioconductor
packages is available at BiocViews. To install a new package, e.g.,
EBImage, use

<p class="code_box">
    source("http://bioconductor.org/biocLite.R")<br/>
    biocLite("EBImage")
</p>

Install "pkg1" and "pkg2" with

<p class="code_box">
    biocLite(c("pkg1", "pkg2"))
</p>


<h2 id="update-bioconductor-packages">Update installed Bioconductor packages</h2>

Bioconductor packages, especially those in the development branch, are
updated fairly regularly. To identify packages requiring update, start
a new session of R and enter

<p class="code_box">
    source("http://bioconductor.org/biocLite.R")<br/>
    old.packages(repos=biocinstallRepos())
</p>

To update all installed packages that are out of date, start a new
session of R and enter

<p class="code_box">
    source("http://bioconductor.org/biocLite.R")<br/>
    update.packages(repos=biocinstallRepos(), ask=FALSE)
</p>

Read the help page for update.packages for additional details.
Recompiling installed Bioconductor packages

Rarely, underlying changes in the operating system require ALL
installed packages to be recompiled for source (C or Fortran)
compatibility. One way to address this might be to start a new R
session and enter

<p class="code_box">
    source("http://bioconductor.org/biocLite.R")<br/>
    pkgs <- rownames(installed.packages())<br/>
    biocLite(pkgs)
</p>

As this will reinstall all currently installed packages, it likely
involves a significant amount of network bandwidth and compilation
time. All packages are implicitly updated, and the cumulative effect
might introduce wrinkles that disrupt your work flow.

