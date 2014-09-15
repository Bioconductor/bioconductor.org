# ![](/images/icons/magnifier.gif)FAQ

This is the FAQ, currently under construction.

* [Package Installation](#install-packages-faq)
* [Package Use](#use-packages-faq)
* [Annotations](#annotation-faq)
* [Developing Packages](#developer-faq)
* [Citation](#citation-faq)

<!--

   Add these above when first FAQ added

* [Package-Specific FAQs](#package-specific-faq)

   Numbering scheme is meant to accommodate future changes in FAQ
   without requiring existing FAQs to be renumbered. FAQs w/in the
   first section are 100, 200, etc. New sections bisect current
   sections, e.g., will be 50, 20, 10 if sections were to be
   added always immediately after the original first
   section.

-->

<h2 id="install-packages-faq">Package Installation</h2>

<h3 class="faq" id="30">biocLite warns that a package is not available</h3>

Most Bioconductor packages are available for Windows, Mac OS, and
Linux operating systems. A few packages are not available on one or
more platforms. This usually occurs because the package relies on
additional software that is not available for the operating
system. For instance, a user trying to install `GeneRfold`
encountered this message:

    > biocLite("GeneRfold")
    Using R version 2.11.1, biocinstall version 2.6.8.
    Installing Bioconductor version 2.6 packages:
    [1] "GeneRfold"
    Please wait...
    Warning message:
    In getDependencies(pkgs, dependencies, available, lib) :
      package 'GeneRfold' is not available

Visiting the list of [package home pages][home-pages]
shows that the package was not available on Windows, the platform on
which the user was trying to install the package. If the package
Description does not indicate why the package is not available, please
feel free to ask on the [Bioconductor support site](/help/support/).

It's useful to check that the package name is spelt correctly, with
correct capitalization!

<h3 class="faq" id="70">Package XXX fails to install</h3>

A common reason for a package to fail to install is that `R` or
`Bioconductor` software dependencies are not satisfied, as shown here
from a user trying to install the `affyPLM` package:

    ...
    ** inst
    ** preparing package for lazy loading
    Error: package 'affy' required by 'affyPLM' could not be found
    Execution halted
    ERROR: lazy loading failed for package 'affyPLM'

Be sure to use `biocLite` to [install packages][2] that are
appropriate for your system and version of `R`.  Be sure that your
installed packages are up-to-date by following [update packages][1].

Less commonly, packages may install but then fail to load, as here
with the `Rsamtools` package:

    Error in dyn.load(file, DLLpath = DLLpath, ...) :
    unable to load shared library 
    '/usr/local/lib64/R/library/Rsamtools/libs/Rsamtools.so':
      /usr/local/lib64/R/library/Rsamtools/libs/Rsamtools.so: undefined symbol: ecrc32

This is likely a system configuration issue, e.g., a Windows `PATH`
environment variable is not set correctly, or the Linux `ldconfig`
program or `LD_LIBRARY_PATH` environment variable is incorrect.

Packages may also fail to install because third party software is not
available. This often happens during the `configure` part of the
package installation process, as illustrated here with the `XML`
package:

    * Installing *source* package 'XML' ...
    ...
    checking for xml2-config... no
    Cannot find xml2-config
    ERROR: configuration failed for package 'XML'

These types of errors can sometimes be easily solved (installing
necessary libraries or other software, perhaps referenced on the
[package home page][home-pages]). It will often be necessary to
understand your system more thoroughly than you'd like, perhaps with
the assistance of the Bioconductor [support site](/help/support/).

<h3 class="faq" id="80">Installing packages from the development repository</h3>

To install a package from the development version of the Bioconductor repository,
install and use the development version of R. <code>biocLite()</code> will automatically
download packages from the correct repository. Similarly, if using older versions
of R, <code>biocLite()</code> will install the appropriate version of packages.

<h2 id="use-packages-faq">Package Use</h2>

<h3 class="faq" id="50">How can I find information about using a package?</h3>

There are three main steps to using a package. (1) Identify an
appropriate package. Do this using [biocViews][/packages/release/BiocViews.html] to browse
available software. (2) Explore overall package functionality and work
flows. Do this by reading the package vignettes, listed on the page
describing the package and available from biocViews. For instance,
locate [IRanges vignettes][iranges-landing-page]. (3) Find help on
particular functions, e.g.,

    library(IRanges)
    help(package="IRanges")  ## overview
    ?findOverlaps            ## specific function

For a more exploratory interface to the help system, try

    help.start()

If you are new to `R`, then it will help to enter into the process
knowing that some basic R skills are assumed by the vignettes and help
pages; spend some time learning `R` itself.

<!--

<h2 id="package-specific-faq">Package-specific questions</h2>

-->


<h2 id="annotation-faq">Annotations</h2>

<h3 class="fq" id="50">Different sources (e.g., annotation packages,
biomaRt, manufacturer, UCSC, GO) disagree. Why?</h3>

Different sources take different approaches to managing
annotations. The annotation packages in Bioconductor are based on
downloads obtained shortly before each Bioconductor release, and so
can lag by six months (at the end of the release cycle) compared to
on-line resources. The advantage of this approach is that the
annotations do not change unexpectedly during development of an
analysis, while the disadvantage is that the resource is not quite
up-to-date with current understanding. To find out information about
data sources used for each annotation package, try a command analogous
to

    library(org.Hs.eg.db)
    org.Hs.eg_dbInfo()

Bioconductor packages can help further investigate discrepancies,
e.g., AnnotationDbi, rtracklayer, Biostrings, and BSgenome (e.g.,
BSgenome.Hsapiens.UCSC.hg17).

<h2 id="developer-faq">Developing Packages</h2>

<h3 class="faq" id="50">What packages belong in the Depends:,
    Imports:, or Suggests: fields?</h3>

Two relevant mailing list posts
([a](https://stat.ethz.ch/pipermail/r-devel/2008-December/051602.html),
[b](https://stat.ethz.ch/pipermail/bioc-devel/2010-September/002310.html))
address this. Generally, packages whose code you use in your own
package should where ever possible be Import:'ed. Packages required
for vignettes are often Suggest:'ed. Depends: is appropriate for
packages that cannot be Import:'ed (e.g., because they do not have a
NAMESPACE) or for packages that provide essential functionality needed
by the user of your package, e.g., your functions always return
`GRanges` objects, so the user will necessarily need `GenomicRanges`
on their search path.

<h2 id="citation-faq">Citation</h2>

<h3 class="faq">How do I cite Bioconductor?</h3>

Many packages provide a citation available from within R. Try

    citation('DESeq')

for any installed package. To cite the project as a whole, use

    citation('Biobase')

[1]: /install/index.html#update-bioconductor-packages
[2]: /install/index.html#install-bioconductor-packages
[mailing-list]: /help/mailing-list/
[home-pages]: /packages/release/bioc/
[bioc-views]: /packages/release/BiocViews.html
[iranges-landing-page]: /packages/release/bioc/html/IRanges.html
