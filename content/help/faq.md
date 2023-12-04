# ![](/images/icons/magnifier.gif) FAQ

* [Package Installation](#install-packages-faq)
  - [Package Not Available](#30)
  - [Fails To Install](#70)
  - [Install From Devel](#80)
* [Package Use](#use-packages-faq)
  - [Find A Package or Package Information](#50)
  - [ERROR: No Package Called](#51)
* [Annotations](#annotation-faq)
  - [Annotations Disagree](#disagree)
* [Developing Packages](#developer-faq)
  - [Guidelines To Follow](#guidelines)
  - [How To Declare Dependencies](#dependencies)
  - [Troubleshooting Build Report](#buildreport)
* [Citation](#citation-faq)
  - [How To Cite Bioconductor](#cite)
* [Connecting](#connecting)
  - [Where to Ask Help](#help)
  - [Connect With Community](#social)

<h2 id="install-packages-faq">Package Installation</h2>

<h3 class="faq" id="30">BiocManager::install() warns that a package is not available</h3>

_Bioconductor_ has a 'release' and a 'devel' branch. The 'release'
branch is intended for most users, and is available using the
**current version of _R_**. New packages are first added to the
'devel' branch, and then become available after the next
_Bioconductor_ release, typically in April and October.  See
[Why use BiocManager::install()?][1] for more information about our release
policy, and [these instructions][2] if you wish to use the devel
branch.

Most Bioconductor packages are available for Windows, Mac OS, and
Linux operating systems. A few packages are not available on one or
more platforms. This usually occurs because the package relies on
additional software that is not available for the operating
system. For instance, a user trying to install `Rcollectl`
encountered this message:

```r
    > BiocManager::install("Rcollectl")
    Bioconductor version 3.18 (BiocManager 1.30.22), R 4.3.2 Patched (2023-11-07
      r85494)
    Installing package(s) 'Rcollectl'
    Warning message:
    package 'Rcollectl' is not available (for R version 4.3.2 Patched)
```

Visiting the list of [package home pages][home-pages]
shows that the package was not available on Windows or Mac OS, the platform on
which the user was trying to install the package. If the package
Description does not indicate why the package is not available, please
feel free to ask on the [Bioconductor support
site][support]. You can also see if any platform is
not supported on the [package build report
pages](https://bioconductor.org/checkResults/release/bioc-LATEST/long-report.html)


It's useful to check that the package name is spelt correctly, with
correct capitalization!

<h3 class="faq" id="70">Package XXX fails to install</h3>

A common reason for a package to fail to install is that `R` or
`Bioconductor` software dependencies are not satisfied, as shown here
from a user trying to install the `affyPLM` package:

```r
    ...
    ** inst
    ** preparing package for lazy loading
    Error: package 'affy' required by 'affyPLM' could not be found
    Execution halted
    ERROR: lazy loading failed for package 'affyPLM'
```

Be sure to use `BiocManager::install` to [install packages][4] that are
appropriate for your system and version of `R`.  Be sure that your
installed packages are up-to-date by following [update packages][3].

Less commonly, packages may install but then fail to load, as here
with the `Rsamtools` package:

```r
    Error in dyn.load(file, DLLpath = DLLpath, ...) :
    unable to load shared library
    '/usr/local/lib64/R/library/Rsamtools/libs/Rsamtools.so':
      /usr/local/lib64/R/library/Rsamtools/libs/Rsamtools.so: undefined symbol: ecrc32
```

This is likely a system configuration issue, e.g., a Windows `PATH`
environment variable is not set correctly, or the Linux `ldconfig`
program or `LD_LIBRARY_PATH` environment variable is incorrect.

Packages may also fail to install because third party software is not
available. This often happens during the `configure` part of the
package installation process, as illustrated here with the `XML`
package:

```r
    * Installing *source* package 'XML' ...
    ...
    checking for xml2-config... no
    Cannot find xml2-config
    ERROR: configuration failed for package 'XML'
```

These types of errors can sometimes be easily solved (installing
necessary libraries or other software, perhaps referenced on the
[package home page][home-pages]). It will often be necessary to
understand your system more thoroughly than you'd like, perhaps with
the assistance of the Bioconductor [support site][support].

<h3 class="faq" id="80">Installing packages from the development repository</h3>

To install a package from the development version of the Bioconductor repository,
install and use the development version of R and
<code>BiocManager::install(version='devel')</code>. See [Using 'Devel' Version
of Bioconductor](http://contributions.bioconductor.org/use-devel.html)


<h2 id="use-packages-faq">Package Use</h2>

<h3 class="faq" id="50">How can I find information about using a package?</h3>

There are three main steps to using a package. (1) Identify an
appropriate package. Do this using [biocViews](/packages/release/BiocViews.html) to browse
available software. (2) Explore overall package functionality and work
flows. Do this by reading the package vignettes, listed on the page
describing the package and available from biocViews. For instance,
locate [IRanges vignettes][iranges-landing-page]. (3) Find help on
particular functions, e.g.,

```r
    library(IRanges)
    help(package="IRanges")  ## overview
    ?findOverlaps            ## specific function
```

For a more exploratory interface to the help system, try

```r
    help.start()
```

If you are new to `R`, then it will help to enter into the process
knowing that some basic R skills are assumed by the vignettes and help
pages; spend some time learning `R` itself.

<h3 class="faq" id="51">Library ERROR: No Package Called xxx</h3>

A package must be [installed](4) before it can be loaded with
<code>library()</code>. Also all package dependencies should be installed. If a
package or package dependency is not installed you may often see something
similar to the following

```r
	library(BiocFileCache)
	Error in library(BiocFileCache) : there is no package called 'BiocFileCache'

```

You must first install the package with BiocManager

```r
	BiocManager::install("BiocFileCache")
	library(BiocFileCache)
```


It's useful to check that the package name is spelt correctly, with
correct capitalization! Remember R is case sensitive! If ERROR persists, reach
out on the [Bioconductor support site][support]


<h2 id="annotation-faq">Annotations</h2>

<h3 class="fq" id="disagree">Different sources (e.g., annotation packages,
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

<h3 class="faq" id="guidelines">Bioconductor Package Development Guidelines</h3>

See the [guidelines](https://contributions.bioconductor.org/) for Bioconductor package development, maintenance, and peer
review. 

<h3 class="faq" id="dependencies">What packages belong in the Depends:,
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

<h3 class="faq" id="buildreport">Troubleshooting the Build Report</h3>

Package maintainers are expected to monitor their packages. If a package starts
failing on the [Bioconductor build system](/checkResults/), it is a maintainers
responsibility to investigate and push any changes to the Bioconductor package
repository to ensure a clean build report. Some tips on reproducing ERRORs and
commonly known issues can be found on the [Troubleshooting Build
Report](http://contributions.bioconductor.org/troubleshooting-build-report.html)
section of the contributions guidelines.

<h2 id="citation-faq">Citation</h2>

<h3 id = "cite" class="faq">How do I cite Bioconductor?</h3>

Many packages provide a citation available from within R. Try the following on
any installed package:

```r
    citation('DESeq')
```

To cite the project as a whole, use

```r
    citation('Biobase')
```

<h2 id="connecting">Connecting With Bioconductor Community</h2>

<h3 id = "help" class="faq">Where To Ask Questions</h3>

Most questions are appropriate to ask on the [Bioconductor support
site][support]. Package development questions are best asked on the [bioc-devel
mailing list](https://stat.ethz.ch/mailman/listinfo/bioc-devel). The
[Bioconductor slack](https://slack.bioconductor.org/) also has specialized
channels for areas of research that may be helpful.

<h3 id = "social" class="faq">Connect With Bioconductor Community</h3>

See [Support Forums](/help/support/) for a listing of ways to connect with the
community.


[1]: /install/index.html#why-biocmanagerinstall
[2]: http://contributions.bioconductor.org/use-devel.html
[3]: /install/index.html#update-bioconductor-packages
[4]: /install/index.html#install-bioconductor-packages
[mailing-list]: /help/mailing-list/
[home-pages]: /packages/release/bioc/
[bioc-views]: /packages/release/BiocViews.html
[iranges-landing-page]: /packages/release/bioc/html/IRanges.html
[support]: https://support.bioconductor.org/