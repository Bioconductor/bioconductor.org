
<ul class="inline_list">
    <li><a href="#bioc-version">Using Bioconductor</a></li>
	<li><a href="#install-R">Install&nbsp;R</a></li>
	<li><a href="#install-bioconductor-packages">Install Packages</a></li> 
	<li><a href="#find-bioconductor-packages">Find Packages</a></li> 
	<li><a href="#update-bioconductor-packages">Update Packages</a></li> 
	<li><a href="#troubleshoot-bioconductor-packages">Troubleshoot Package Installations</a></li> 
</ul>


<h2 id="bioc-version">Using <em>Bioconductor</em></h2>

The current release of Bioconductor is version
<%=config[:release_version] %>; it works with R version
<%=config[:r_version_associated_with_release]%>. Users of older R and
Bioconductor versions must update their installation to take advantage
of new features.

[Install](#install-R) the latest release of R, then get the latest version of
Bioconductor by starting R and entering the commands

    source("http://bioconductor.org/biocLite.R")
    biocLite()

Details, including instructions to
[install additional packages](#install-bioconductor-packages") and
[update packages](#update-bioconductor-packages"), are provided below.
A [devel](/developers/how-to/useDevel/) version of Bioconductor is
available.

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

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="install-bioconductor-packages">Install <em>Bioconductor</em> Packages</h2>

Use the `biocLite.R` script to install Bioconductor packages. To
install core packages, type the following in an R command window:

    source("http://bioconductor.org/biocLite.R")
    biocLite()

Install specific packages, e.g., "GenomicFeatures" and "AnnotationDbi", with

    biocLite(c("GenomicFeatures", "AnnotationDbi"))

The `biocLite()` function (in the BiocInstaller package installed by
the `biocLite.R` script) has arguments that change its default
behavior; type `?biocLite` for further help.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="find-bioconductor-packages">Find <em>Bioconductor</em> Packages</h2>

Visit the [Workflows](/help/workflows/) page
and [software package list](/packages/release/BiocViews.html#___Software)
to discover available packages.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="update-bioconductor-packages">Update Installed <em>Bioconductor</em> Packages</h2>

Bioconductor packages, especially those in the development branch, are
updated fairly regularly. To identify packages requiring update within
your version of Bioconductor, start a new session of R and enter

    source("http://bioconductor.org/biocLite.R")
	biocLite()                  ## R version 3.0 or later

Use the argument `ask=FALSE` to update old packages without being
prompted.  For older versions of `R`, use the command
`biocLite(NULL)`.  Read the help page for `?biocLite` for additional
details.

<h3>Upgrading installed <em>Bioconductor</em> packages</h3>

Some versions of R support more than one version of Bioconductor. To
use the latest version of Bioconductor for your version of R, enter

    source("http://bioconductor.org/biocLite.R")
	biocLite("BiocUpgrade")     ## R version 2.15 or later

Read the help page for `?BiocUpgrade` for additional details. Remember
that more recent versions of Bioconductor may be available if your
version of R is out-of-date.

<h3>Recompiling installed <em>Bioconductor</em> packages</h3>

Rarely, underlying changes in the operating system require ALL
installed packages to be recompiled for source (C or Fortran)
compatibility. One way to address this might be to start a new R
session and enter

    source("http://bioconductor.org/biocLite.R")
    pkgs <- rownames(installed.packages())
    biocLite(pkgs, type="source")

As this will reinstall all currently installed packages, it likely
involves a significant amount of network bandwidth and compilation
time. All packages are implicitly updated, and the cumulative effect
might introduce wrinkles that disrupt your work flow. It also requires
that you have the necessary compilers installed.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="troubleshoot-bioconductor-packages">Troubleshoot Package Installations</h2>

Use the commands

    library(BiocInstaller)
	biocValid()             ## R version 3.0 or later

to flag packages that are either out-of-date or too new for your
version of Bioconductor. The output suggests ways to solve identified
problems, and the help page `?biocValid` lists arguments influencing
the behavior of the function.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

