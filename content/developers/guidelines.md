![](/images/icons/magnifier.gif)Package Guidelines
==================================================

* [Introduction](#introduction)
* [Correctness, Space and Time](#correctness)
* [Package Name](#name)
* [License](#license)
* [Package Content](#content)
* [Package Dependencies](#dependencies)
* [S4 Classes and Methods](#classes)
* [Vectorized Calculations](#vectorized)
* [End-User Messages](#messages)
* [The Sweave Vignette](#vignettes)
* [Citations](#citations)
* [Version Numbering](#versions)
* [C or Fortran code](#c-code)
* [Duplication of Packages in CRAN and Bioconductor](#duplications)
* [Package Maintainer Responsibilities](#responsibilities)

<h2 id="introduction">Introduction</h2>

The Bioconductor project strives to promote high-quality, well documented,
and interoperable software. These guidelines help to achieve this object;
they are not meant to put undue burden on package authors, and authors
having difficultly satisfying guidelines should seek advice on the
[bioc-devel](/help/mailing-list/) mailing list.

Package maintainers are urged to follow these guidelines as closely as
possible when developing Bioconductor packages.

<h2 id="correctness">Correctness, Space and Time</h2>

Bioconductor packages must pass `R CMD build` (or `R CMD INSTALL --build`)
and pass `R CMD check` with no errors and no warnings using a recent R-devel.
Authors should also try to address all notes that arise during build or check.

Do not use filenames that differ only in case, as not all file systems are
case sensitive.

The source package resulting from running `R CMD build` should occupy less than
2MB on disk. The package should require less than 5 minutes to run R CMD check.
This includes the time required to build the Sweave vignette.

<h2 id="name">Package Name</h2>

Choose a descriptive name. An easy way to check whether your name is already
in use is to check that the following command fails

    source("http://bioconductor.org/biocLite.R")
    biocLite("mypackage")

<h2 id="license">License</h2>

The "License:" field in the DESCRIPTION file should preferably refer to a
standard license (see [opensource.org](http://wiki.fhcrc.org/bioc/opensource.org)
or [wikipedia](http://en.wikipedia.org/wiki/Comparison_of_free_software_licences))
using one of R's standard specifications. Be specific about any version that
applies (e.g., GPL-2). Core Bioconductor packages are typically licensed under
Artistic-2.0. To specify a non-standard license, include a file named LICENSE
in your package (containing the full terms of your license) and use the string
"file LICENSE" (without the double quotes) in the "License:" field of your
DESCRIPTION file.

<h2 id="content">Package Content</h2>

Packages must

* Contain a Sweave-style
  [vignette](http://cran.fhcrc.org/doc/manuals/R-exts.html#Writing-package-vignettes)
  that demonstrates how to use the package to accomplish a task (more on this
  below).
* Include examples in all
  [man pages](http://cran.fhcrc.org/doc/manuals/R-exts.html#Rd-format).
* Specify one or more
  [biocViews categories](http://wiki.fhcrc.org/bioc/biocViews_categories).
* Contain a
  [NAMESPACE](http://cran.fhcrc.org/doc/manuals/R-exts.html#Package-name-spaces)
  file to define the functions, classes, and methods that are imported into the
  name space, and exported for users.
* Contain (literature) references to the methods used as well as to other
  similar or related packages.
* Make use of appropriate existing packages (e.g., biomaRt, AnnotationDbi,
  Biostrings) and classes (e.g., ExpressionSet, AnnotatedDataFrame, RangedData,
  RLE, DNAStringSet), and avoid duplication of functionality available in other
  Bioconductor packages.
* Document data structures used and, if different from data structures used by
  similar packages, explain why a different data structure was used.
* Contain only code that can be redistributed according to the package license.
  In particular, packages may not include any code from
  [Numerical Recipes](http://www.nr.com/).
* Not contain unnecessary files such as .DS_Store, .project, .svn, cache file,
  log files, etc.

<h2 id="dependencies">Package Dependencies</h2>

Reuse, rather than re-implement or duplicate, well-tested functionality from
other packages. Specify package dependencies in the DESCRIPTION file, listed
as follows

* **Imports**: is for packages that provide functions, methods, or classes
  that are used inside your package name space. Most dependencies are listed
  here.
* **Depends**: is appropriate when the package whose functionality you are
  using does not have a name space. In this case, use fully qualified variables
  (pkg::variable). Depends: is also appropriate when a package is used in the
  example section of a man page. It is very unusual for a package to list more
  than three packages as 'Depends:'.
* **Suggests**: is appropriate for packages used in your vignette.

Packages should specify the R version on which they depend. This is usually the
current development version.

<h2 id="classes">S4 Classes and Methods</h2>

We recommend the following structure/layout:

* All class definitions in R/AllClasses.R
* All generic function definitions in R/AllGenerics.R
* Methods are defined in a file named by the generic function. For example, all
  `show` methods would go in R/show-methods.R. This is not written in stone,
  but tends to provide a useful organization. Sometimes a collection of methods
  that provide the interface to a class are best put in a SomeClass-accessors.R
  file.

A Collates: field in the DESCRIPTION file may be necessary to order class and
method definitions appropriately during package installation.

<h2 id="vectorized">Vectorized Calculations</h2>

Many R operations are performed on the whole object, not just the elements of
the object (e.g., sum(x), not x[1] + x[2] + ...). In particular, relatively few
situations require an explicit for loop.

<h2 id="messages">End-User Messages</h2>

* `message()` communicates diagnostic messages (e.g., progress during lengthy
  computations) during code evaluation.
* `warning()` communicates unusual situations handled by your code.
* `stop()` indicates an error condition.
* `cat()` or `print()` are used only when displaying an object to the user,
  e.g., in a `show` method.

<h2 id="vignettes">The Sweave Vignette</h2>

A vignette demonstrates how to accomplish non-trivial tasks embodying the core
functionality of your package. A Sweave vignette is an .Rnw file that contains
LaTeX and chunks of R code. The R code chunk starts with a line <<>>=, and ends
with @. Each chunk is evaluated during `R CMD build`, prior to LaTeX
compilation. Refer to
[Writing package vignettes](http://cran.fhcrc.org/doc/manuals/R-exts.html#Writing-package-vignettes)
for technical details.

A vignette provides reproducibility: the vignette produces the same results as
copying the corresponding commands into an R session. It is therefore essential
that the vignette embed R code between <<>>= and @; short-cuts (e.g., using a
LaTeX verbatim environment, or using the Sweave eval=FALSE flag) undermine the
benefit of vignettes.

All packages are expected to have at least one Sweave vignette.

<h2 id="citations">Citations</h2>

Appropriate citations must be included in help pages (e.g., in the see also
section) and vignettes; this aspect of documentation is no different from any
scientific endeavor. The file inst/CITATION can be used to specify how a
package is to be cited.

<h2 id="versions">Version Numbering</h2>

All Bioconductor packages use an x.y.z version scheme. The following rules
apply:

* x is usually 0 for packages that have not yet been released.
* y is even for packages in release, and odd for packages in devel.
* z is incremented whenever committing changes to a package.

For more details, see
[Version Numbering Standards](http://wiki.fhcrc.org/bioc/Version_Numbering_Standards)

<h2 id="c-code">C or Fortran code</h2>

If the package contains C or Fortran code, it should adhere to the standards
and methods described in the
[System and foreign language interfaces](http://cran.r-project.org/doc/manuals/R-exts.html#System-and-foreign-language-interfaces)
section of the Writing R Extensions manual. In particular:

* Use internal R functions, e.g., R_alloc and random number generators, over
  system supplied ones.
* Use C function registration (See the
  [Registering native routines](http://cran.fhcrc.org/doc/manuals/R-exts.html#Registering-native-routines))
* Use R_CheckUserInterrupt in C level loops when there is a chance that they
  may not terminate for certain parameter settings or when their run time
  exceeds 10 seconds with typical parameter settings, and the method is
  intended for interactive use.
* Make judicious use of Makevars and Makefile. These are often not required at
  all (See the
  [Configure and cleanup](http://cran.fhcrc.org/doc/manuals/R-exts.html#Configure-and-cleanup)).

Use of external libraries whose functionality is redundant with libraries
already supported is strongly discouraged. In cases where the external library
is complex the author may need to supply pre-built binary versions for some
platforms.

<h2 id="duplications">Duplication of Packages in CRAN and Bioconductor</h2>

Authors are strongly discouraged from placing their package into both CRAN and
Bioconductor. This avoids burdening the author with extra work and confusing
the user.

<h2 id="responsibilities">Package Maintainer Responsibilities</h2>

Acceptance of packages into Bioconductor brings with it ongoing responsibility
for package maintenance. These responsibilities include:

* Subscription to the bioc-devel mailing list.
* Response to bug reports and questions from users regarding your package, as
  posted on the bioconductor mailing list.
* Package Maintenance through software release cycles, including prompt updates
  to software and documentation necessitated by underlying changes in R.
