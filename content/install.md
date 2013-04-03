
<ul class="inline_list">
    <li><a href="#bioc-version">Getting The Latest Version of Bioconductor</a></li>
	<li><a href="#install-bioconductor-packages">Install Packages</a></li> 
	<li><a href="#find-bioconductor-packages">Find Packages</a></li> 
	<li><a href="#update-bioconductor-packages">Update Packages</a></li> 
	<li><a href="#troubleshoot-bioconductor-packages">Troubleshoot Package Installations</a></li> 
	<li><a href="#install-R">Install&nbsp;R</a></li> 
</ul>


<h2 id="bioc-version">Getting The Latest Version of Bioconductor</h2>

After installing the latest release of R, get packages from the latest
version of Bioconductor by following the steps below. The current
release version of R is <%=
config[:r_version_associated_with_release]%>, and the currently
released Bioconductor version is <%= config[:release_version] %>.

<h2 id="install-bioconductor-packages">Install Bioconductor Packages</h2>

Use the `biocLite.R` script to install Bioconductor packages. To
install core packages, type the following in an R command window:

    source("http://bioconductor.org/biocLite.R")
    biocLite()

Install specific packages, e.g., "GenomicFeatures" and "AnnotationDbi", with

    biocLite(c("GenomicFeatures", "AnnotationDbi"))

To update all your installed packages, use:

    biocLite(character())

The `biocLite()` function (in the BiocInstaller package installed by
the`biocLite.R` script) has arguments that change its default
behavior; type `?biocLite` for further help.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="find-bioconductor-packages">Find Bioconductor Packages</h2>

Visit the [Workflows](/help/workflows/) page
and [software package list](/packages/release/BiocViews.html#___Software)
to discover available packages.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="update-bioconductor-packages">Update Installed Bioconductor Packages</h2>

Bioconductor packages, especially those in the development branch, are
updated fairly regularly. To identify packages requiring update, start
a new session of R and enter

    source("http://bioconductor.org/biocLite.R")
	biocLite()              ## R version 3.0 or later

Use the argument `ask=FALSE` to update old packages without being
prompted.  For older versions of `R`, use the command
`biocLite(NULL)`.  Read the help page for `?biocLite` for additional
details.

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

<h2 id="troubleshoot-bioconductor-packages">Troubleshoot Package Installations</h2>

Use the commands

    library(BiocInstaller)
	biocValid()             ## R version 3.0 or later

to flag packages that are either out-of-date or too new for your
version of Bioconductor. The output suggests ways to solve identified
problems.

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

