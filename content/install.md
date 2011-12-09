
<ul class="inline_list">
    <li><a href="#bioc-version">Getting The Latest Version of Bioconductor</a></li>
	<li><a href="#install-bioconductor-packages">Install Packages</a></li> 
	<li><a href="#find-bioconductor-packages">Find Packages</a></li> 
	<li><a href="#update-bioconductor-packages">Update Packages</a></li> 
	<li><a href="#install-R">Install R</a></li> 
</ul>


<h2 id="bioc-version">Getting The Latest Version of Bioconductor</h2>

If you have installed the latest release of R, you will automatically
get packages from the latest version of Bioconductor by following the steps
below. The current release version of R is
<%= config[:r_version_associated_with_release]%>, and the currently
released Bioconductor version is <%= config[:release_version] %>.

<h2 id="install-bioconductor-packages">Install Bioconductor Packages</h2>

Use the `biocLite.R` script to install Bioconductor packages.  To
install a particular package, e.g., limma, type the following in an R
command window:

    source("http://bioconductor.org/biocLite.R")
    biocLite("limma")

Install several packages, e.g., "GenomicFeatures" and "AnnotationDbi", with

    biocLite(c("GenomicFeatures", "AnnotationDbi"))

To install a selection of core Bioconductor packages, use

    biocLite()

Packages and their dependencies installed by this usage are: 
`Biobase`, `IRanges`, and `AnnotationDbi`.


The BiocInstaller package (installed by the`biocLite.R` script)
has arguments that change its default behavior:

    pkgs
        Character vector of Bioconductor packages to install.
    destdir
        File system directory for downloaded packages.
    lib
        R library where packages are installed.


You can type `?biocLite` for further installation help.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="find-bioconductor-packages">Find Bioconductor Packages</h2>

Visit the [Workflows](/help/workflows/) page,
/packages[BiocViews](/packages/<%=config[:release_version]%>/BiocViews.html)
taxonomy, and [software package list](/packages/release/bioc/)
to discover available packages.

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
    biocLite(character(), ask=FALSE)

Read the help page for `?biocLite` for additional details.

<h3>Recompiling installed Bioconductor packages</h3>

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

