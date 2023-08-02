
<ul class="inline_list">
    <li><a href="#bioc-version">Using Bioconductor</a></li>
	<li><a href="#install-R">Install&nbsp;R</a></li>
	<li><a href="#install-bioconductor-packages">Install&nbsp;Packages</a></li>
	<li><a href="#find-bioconductor-packages">Find&nbsp;Packages</a></li>
	<li><a href="#update-bioconductor-packages">Update&nbsp;Packages</a></li>
	<li><a href="#troubleshoot-bioconductor-packages">Troubleshoot&nbsp;Package&nbsp;Installations</a></li>
	<li><a href="#why-biocmanagerinstall">Why&nbsp;BiocManager::install()?</a></li>
	<li><a href="#preconfigured">Pre-configured&nbsp;Bioconductor</a></li>
	<li><a href="#Legacy">Legacy&nbsp;and&nbsp;Older&nbsp;R&nbsp;Versions</a></li>
</ul>


<h2 id="bioc-version">Using <em>Bioconductor</em></h2>

The current release of _Bioconductor_ is version
<%=config[:release_version] %>; it works with _R_ version
<%=config[:r_version_associated_with_release]%>. Users of older R and
_Bioconductor_ must update their installation to take advantage
of new features and to access packages that have been added to
_Bioconductor_ since the last release.

The development version of _Bioconductor_ is version
<%=config[:devel_version] %>; it works with _R_ version
<%=config[:r_version_associated_with_devel]%>. More recent 'devel'
versions of _R_ (if available) will be supported during the next
_Bioconductor_ release cycle.

[Install](#install-R) the latest release of R, then get the latest version of
_Bioconductor_ by starting R and entering the commands

    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install(version = "<%=config[:release_version] %>")

It may be possible to change the _Bioconductor_ version of an existing
installation; see the '[Changing version][]' section of the BiocManager
vignette.

[Changing version]: https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html

Details, including instructions to
[install additional packages](#install-bioconductor-packages) and to
[update](#update-bioconductor-packages),
[find](#find-bioconductor-packages), and
[troubleshoot](#troubleshoot-bioconductor-packages) are provided
below.  A [devel](/developers/how-to/useDevel/) version of
_Bioconductor_ is available. There are good
[reasons for using `BiocManager::install()`](#why-biocmanagerinstall) for
managing _Bioconductor_ resources.

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

To install core packages, type the following in an R command window:

    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install()

Install specific packages, e.g., "GenomicFeatures" and "AnnotationDbi", with

    BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))

The `install()` function (in the BiocManager package) has arguments that change
its default behavior; type `?install` for further help.

For a more detailed explanation on using BiocManager and its advanced usage,
such as version switching, please refer to the
[BiocManager vignette](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html).

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="find-bioconductor-packages">Find <em>Bioconductor</em> Packages</h2>

Visit the [software package list](/packages/release/BiocViews.html#___Software)
to discover available packages.

To search through available packages programmatically, use the following:

    BiocManager::available()

For example, using a "^org" search pattern will show all of the available
organism annotation packages.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="update-bioconductor-packages">Update Installed <em>Bioconductor</em> Packages</h2>

_Bioconductor_ packages, especially those in the development branch, are
updated fairly regularly. To identify packages requiring update within
your version of _Bioconductor_, start a new session of R and enter

    BiocManager::install()

Use the argument `ask=FALSE` to update old packages without being
prompted.  Read the help page for `?install` for additional details.

<h3 id="upgrade-bioconductor-packages">Upgrading installed <em>Bioconductor</em> packages</h3>

Due to the development cycle, all versions of R will eventually support more
than one version of _Bioconductor_. To use the latest version of _Bioconductor_
for your version of R, enter

    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install()

Remember that more recent versions of _Bioconductor_ may be available if your
version of R is out-of-date. BiocManager will notify you when your version
of R is out-of-date.

For more details on <em>Bioconductor</em> approaches to versioning, see
the <a
href="https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html#managing-multiple-versions">advanced section</a>
in the vignette and version numbering in the <a
href="/developers/how-to/version-numbering/">developer reference</a> section.

<h3>Recompiling installed <em>Bioconductor</em> packages</h3>

Rarely, underlying changes in the operating system require ALL
installed packages to be recompiled for source (C or Fortran)
compatibility. One way to address this might be to start a new R
session and enter

    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    pkgs <- rownames(installed.packages())
    BiocManager::install(pkgs, type = "source", checkBuilt = TRUE)

As this will reinstall all currently installed packages, it likely
involves a significant amount of network bandwidth and compilation
time. All packages are implicitly updated, and the cumulative effect
might introduce wrinkles that disrupt your work flow. It also requires
that you have the necessary compilers installed.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="troubleshoot-bioconductor-packages">Troubleshoot Package Installations</h2>

Use the commands

    BiocManager::valid()     ## R version 3.5 or later

to flag packages that are either out-of-date or too new for your
version of _Bioconductor_. The output suggests ways to solve identified
problems, and the help page `?valid` lists arguments influencing
the behavior of the function.

<h3 id="troubleshoot-biocmanager">Troubleshoot BiocManager</h3>

One likely reason for BiocManager not working on your system could
be that your version of R is too old for `BiocManager`. In order
avoid this issue, please ensure that you have the latest version of R
installed in your system. BiocManager supports R versions from 3.5.0
and above.


<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="why-biocmanagerinstall">Why use BiocManager::install()?</h2>

`BiocManager::install()` is the recommended way to install _Bioconductor_
packages. There are several reasons for preferring this to the
'standard' way in which R pacakges are installed via
`install.packages()`.

_Bioconductor_ has a repository and release schedule that differs from R
(_Bioconductor_ has a 'devel' branch to which new packages and updates
are introduced, and a stable 'release' branch emitted once every 6
months to which bug fixes but not new features are introduced).

A consequence of the mismatch between R and _Bioconductor_ release
schedules is that the _Bioconductor_ version identified by
`install.packages()` is sometimes not the most recent 'release'
available. For instance, an R minor version may be introduced some
months before the next Bioc release. After the Bioc release the users
of the R minor version will be pointed to an out-of-date version of
_Bioconductor_.

A consequence of the distinct 'devel' branch is that
`install.packages()` sometimes points only to the 'release'
repository, whereas _Bioconductor_ developers and users wanting
leading-edge features wish to access the _Bioconductor_ 'devel'
repository. For instance, the _Bioconductor_ 3.0 release is available
for R.3.1.x, so _Bioconductor_ developers and leading-edge users need to
be able to install the devel version of _Bioconductor_ packages into the
same version (though perhaps different instance or at least library
location) of R that supports version 2.14 of _Bioconductor_.

An indirect consequence of _Bioconductor_'s structured release is that
packages generally have more extensive dependencies with one another,
both explicitly via the usual package mechanisms and implicitly
because the repository, release structure, and _Bioconductor_ community
interactions favor re-use of data representations and analysis
concepts across packages. There is thus a higher premium on knowing
that packages are from the same release, and that all packages are
current within the release.

The BiocManager package serves as the primary way to ensure that
the appropriate _Bioconductor_ installation is used with respect
to the version of R in use regardless of the R and _Bioconductor_
release cycles.

    > library(BiocManager)
    Bioconductor version 3.9 (BiocManager 1.30.4), ?BiocManager::install
    for help

The `install()` function is provided by BiocManager. This is a
wrapper around `install.packages`, but with the repository chosen
according to the version of _Bioconductor_ in use, rather than to the
version relevant at the time of the release of R.

`install()` also nudges users to remain current within a release, by
default checking for out-of-date packages and asking if the user would
like to update

    > BiocManager::install()
    Bioconductor version 3.9 (BiocManager 1.30.4), R 3.6.0 Patched
    (2019-05-02 r76454)
    Update old packages: 'BBmisc', 'genefilter', 'GenomicAlignments',
      'GenomicRanges', 'IRanges', 'MASS', 'reshape2', 'Rgraphviz',
      'RJSONIO', 'rtracklayer'
    Update all/some/none? [a/s/n]:

The BiocManager package provides facilities for switching to the
'devel' version of _Bioconductor_

    > BiocManager::install(version = "devel")
    Upgrade 89 packages to Bioconductor version '3.10'? [y/n]: y
    Installing package(s) 'BiocVersion'
    trying URL 'https://bioconductor.org/packages/3.10/bioc/src/contrib/BiocVersion_3.10.0.tar.gz'
    Content type 'application/x-gzip' length 987 bytes
    ==================================================
    downloaded 987 bytes

    * installing *source* package 'BiocVersion' ...
    ** help
    *** installing help indices
    ** building package indices
    ** testing if installed package can be loaded
    * DONE (BiocVersion)

    ...
    Bioconductor version 3.10 (BiocManager 1.30.4), ?BiocManager::install for
    help

(at some points in the R / _Bioconductor_ release cycle use of 'devel'
requires use of a different version of R itself, in which case the
attempt to install devel fails with an appropriate message).

The BiocManager package also provides `valid()` to test that the
installed packages are not a hodgepodge from different _Bioconductor_
releases (the 'too new' packages have been installed from source
rather than a repository; regular users would seldom have these).

    > BiocManager::valid()

    * sessionInfo()

    R version 3.6.0 Patched (2019-05-02 r76454)
    Platform: x86_64-pc-linux-gnu (64-bit)
    ...

    Bioconductor version '3.9'

      * 2 packages out-of-date
      * 0 package too new
    ...
    create a valid installation with

      BiocManager::install(c(
        "GenomicFeatures", "AnnotationDbi"
      ), update = TRUE, ask = FALSE)

    more details: BiocManager::valid()$too_new, BiocManager::valid()$out_of_date

    Warning message:
    2 packages out-of-date; 0 packages too new

For users who spend a lot of time in _Bioconductor_, the features
outlined above become increasingly important and `install()` is much
preferred to `install.packages()`.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="preconfigured">Pre-configured <em>Bioconductor</em></h2>

_Bioconductor_ is also available as [Docker images](/help/docker/) or available
for use in the [AnVIL](https://anvil.bioconductor.org/)

<h2 id="Legacy">Legacy and Older R Versions</h2>

It is always recommended to update to the most current version of <em>R</em> and
<em>Bioconductor</em>.  If this is not possible and `R < 3.5.0` , please use the following
for installing <em>Bioconductor</em> packages

To install core packages, type the following in an R command window:

    source("https://bioconductor.org/biocLite.R")

Install specific packages, e.g., "GenomicFeatures" and "AnnotationDbi", with

    BiocInstaller::biocLite(c("GenomicFeatures", "AnnotationDbi"))


<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>
