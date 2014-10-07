![](/images/icons/magnifier.gif)Package Guidelines
==================================================
<a name="top"></a>

* [Introduction](#introduction)
* [Types of Packages](#package-types)
* [Version of Bioconductor and R](#version)
* [Correctness, Space and Time](#correctness)
* [Package Name](#name)
* [License](#license)
* [Package Content](#content)
* [Package Dependencies](#dependencies)
* [S4 Classes and Methods](#classes)
* [Vectorized Calculations](#vectorized)
* [End-User Messages](#messages)
* [Graphics Device](#graphical)
* [Vignette(s)](#vignettes)
* [Citations](#citations)
* [Version Numbering](#versions)
* [C/C++ or Fortran code](#c-code)
* [Unit tests](#unitTests)
* [Videos](#videos)
* [Duplication of Packages in CRAN and Bioconductor](#duplications)
* [Package Author and Maintainer Responsibilities](#responsibilities)

<h2 id="introduction">Introduction</h2>

The Bioconductor project promotes high-quality, well documented, and
interoperable software. These guidelines help to achieve this objective;
they are not meant to put undue burden on package authors, and authors
having difficultly satisfying guidelines should seek advice on the
[bioc-devel](/help/mailing-list/) mailing list.

Package maintainers are urged to follow these guidelines as closely as
possible when developing Bioconductor packages.

General instructions for producing packages can be found in the
[Writing R Extensions](http://cran.r-project.org/doc/manuals/R-exts.html)
manual, available from within R (`RShowDoc("R-exts")`) or on the [R web
site](http://cran.fhcrc.org/manuals.html).

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>


<h2 id="package-types">Types of Packages</h2>

Most packages contributed by users are [software][software-pkgs]
packages that perform analytic calculations. Users also contribute
[annotation][annotation-pkgs] and [experiment data][exptdata-pkgs]
packages. Annotation packages are database-like packages that provide
information linking identifiers (e.g., Entrez gene names or Affymetrix
probe ids) to other information (e.g., chromosomal location, Gene
Ontology category). Experiment data packages provide data sets that
are used, often by software packages, to illustrate particular
analyses. An excellent practice is to develop a software package, and
to provide or use an existing experiment data package to give a
comprehensive illustration of the methods in the software package. The
guidelines below apply to all packages, but annotation and experiment
data packages are not required to conform to the space limitations of
software packages. Developers wishing to contribute annotation or
experiment data packages should seek [additional support][support]
associated with package submission.

[software-pkgs]: /packages/release/bioc/
[annotation-pkgs]: /packages/release/data/annotation/
[exptdata-pkgs]: /packages/release/data/experiment/
[support]: /developers/package-submission/#support

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>


<h2 id="version">Version of Bioconductor and R</h2>

Package developers should always use the 
[devel version of Bioconductor](/packages/devel/BiocViews.html#___Software) when developing and testing packages to be contributed.

Depending on the R release cycle, using Bioconductor devel may or may
not involve also using the devel version of R. See the how-to on
[using devel version of Bioconductor](/developers/how-to/useDevel/)
for up-to-date information.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>


<h2 id="correctness">Correctness, Space and Time</h2>

Bioconductor packages must pass `R CMD build` (or 
<code>R CMD INSTALL &#45;&#45;build</code>)
and pass `R CMD check` with no errors and no warnings using a recent R-devel.
Authors should also try to address all notes that arise during build or check.

Do not use filenames that differ only in case, as not all file systems are
case sensitive.

The source package resulting from running `R CMD build` should occupy 
less than 4MB on disk. The package should require less than 5 minutes to run
<code>R CMD check &#45;&#45;no&#45;build&#45;vignettes</code>. 
Using the <code>&#45;&#45;no&#45;build&#45;vignettes</code>
option ensures that the vignette is built only once.

Vignette and man page examples should not use more than 2GB of memory
since R cannot allocate more than this on 32-bit Windows.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="name">Package Name</h2>

Choose a descriptive name. An easy way to check whether your name is already
in use is to check that the following command fails

    source("http://bioconductor.org/biocLite.R")
    biocLite("MyPackage")

Avoid names that are easily confused with existing package names, or
that imply a temporal (e.g., `ExistingPackage2`) or qualitative (e.g.,
`ExistingPackagePlus`) relationship.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

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

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="content">Package Content</h2>

Packages must

* Contain a 
  [vignette](http://cran.fhcrc.org/doc/manuals/R-exts.html#Writing-package-vignettes)
  that demonstrates how to use the package to accomplish a task (more on this
  below).
* Include examples in all
  [man pages](http://cran.fhcrc.org/doc/manuals/R-exts.html#Rd-format).
* Specify one or more
  [biocViews categories](/packages/release/BiocViews.html#___Software).
* Contain a
  [NAMESPACE](http://cran.fhcrc.org/doc/manuals/R-exts.html#Package-namespaces)
  file to define the functions, classes, and methods that are imported into the
  name space, and exported for users.
* Contain (literature) references to the methods used as well as to other
  similar or related packages.
* Make use of appropriate existing packages (e.g., biomaRt, AnnotationDbi,
  Biostrings) and classes (e.g., ExpressionSet, AnnotatedDataFrame, GRanges,
  Rle, DNAStringSet), and avoid duplication of functionality available in other
  Bioconductor packages. Note that, starting with BioC 2.12, the use of
  RangedData or RangedDataList objects (those classes are defined in the
  IRanges package) is discouraged so new contributed packages should use
  GRanges/GRangesList objects instead (those classes are defined in the
  GenomicRanges package).
* Document data structures used and, if different from data structures used by
  similar packages, explain why a different data structure was used.
* Contain only code that can be redistributed according to the package license.
  Be aware of the licensing agreements for packages you are depending on in your package. 
  Not all packages are open source even if they are publicly available.
  In particular, packages may not include any code from
  [Numerical Recipes](http://www.nr.com/).
* Not contain unnecessary files such as .DS_Store, .project, .svn, cache file,
  log files, etc.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="dependencies">Package Dependencies</h2>

Packages you depend on must be available via Bioconductor or CRAN;
users and the automated build system has no way to install packages
from any other source.

Reuse, rather than re-implement or duplicate, well-tested functionality from
other packages. Specify package dependencies in the DESCRIPTION file, listed
as follows

* **Imports**: is for packages that provide functions, methods, or classes
  that are used inside your package name space. Most dependencies are listed
  here.
* **Depends**: is appropriate when a package is used in the
  example section of a man page. It is very unusual for a package to list more
  than three packages as 'Depends:'.
* **Suggests**: is appropriate for packages used in your vignette.
* **Enhances**: is for packages such as Rmpi or parallel that enhance
    the performance of the current package, but are not strictly
    needed for its functionality.
* **SystemRequirements**: is for listing any external software which 
   is required, but not automatically installed by the normal package
   installation process. If the installation process is non-trivial,
   a top-level README file should be included to document the process.

A package may rarely offer optional functionality, e.g., visualization
with `rgl` when that package is available. Authors then list the
package in the **Suggests** field, and use `requireNamespace()` (or
`loadNamespace()`) to condition code execution. Functions from the
loaded namespace should be accessed using `::` notation, e.g.,

    x <- sort(rnorm(1000))
    y <- rnorm(1000)
    z <- rnorm(1000) + atan2(x,y)
    if (requireNamespace("rgl", quietly=TRUE)) {
        rgl::plot3d(x, y, z, col=rainbow(1000))
    } else {
        ## code when "rgl" is not available
    }

This approach does not alter the user `search()` path, and ensures
that the necessary function (`plot3d()`, from the `rgl` package) is
used.  Such conditional code increases complexity of the package and
frustrates users who do not understand why behavior differs between
installations, so is often best avoided.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="classes">S4 Classes and Methods</h2>

Re-use existing S4 classes and generics where possible. This
encourages interoperability and simplifies your own package
development. If your data requires a new representation or function,
carefully design an S4 class or generic so that other package
developers with similar needs will be able to re-use your hard work,
and so that users of related packages will be able to seamlessly use
your data structures. Do not hesitate to ask on the Bioc-devel mailing
list for advice.

Implement a constructor (typically a simple function) if the user is
supposed to be able to create an instance of your class. Write short
accessors (functions or methods) if the user needs to extract from or
assign to slots in the class. Constructors and accessors help separate
the interface seen by the user from the implementation details
relevant to the developer.

The following layout is sometimes used to organize classes and
methods; other approaches are possible and acceptable.

* All class definitions in R/AllClasses.R
* All generic function definitions in R/AllGenerics.R
* Methods are defined in a file named by the generic function. For example, all
  `show` methods would go in R/show-methods.R.

A Collates: field in the DESCRIPTION file may be necessary to order class and
method definitions appropriately during package installation.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="vectorized">Vectorized Calculations</h2>

Many R operations are performed on the whole object, not just the elements of
the object (e.g., sum(x), not x[1] + x[2] + ...). In particular, relatively few
situations require an explicit for loop.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="messages">End-User Messages</h2>

* `message()` communicates diagnostic messages (e.g., progress during lengthy
  computations) during code evaluation.
* `warning()` communicates unusual situations handled by your code.
* `stop()` indicates an error condition.
* `cat()` or `print()` are used only when displaying an object to the user,
  e.g., in a `show` method.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="graphical">Graphics Device</h2>

Use `dev.new()` to start a graphics device if necessary. Avoid using `x11()`
or `X11()` for it can only be called on machines that have access to an X
server. 

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="vignettes">Vignette(s)</h2>

A vignette demonstrates how to accomplish non-trivial tasks embodying
the core functionality of your package. There are two common types of
vignettes. A _Sweave_ vignette is an .Rnw file that contains LaTeX and
chunks of R code. The R code chunk starts with a line <<>>=, and ends
with @. Each chunk is evaluated during `R CMD build`, prior to LaTeX
compilation to a PDF document. An _R markdown_ vignette is similar to
a Sweave vignette, but uses
[markdown](http://daringfireball.net/projects/markdown/) instead of
LaTeX for structuring text sections and resulting in HTML output. The
[knitr](http://yihui.name/knitr/) package can process most Sweave and
all R markdown vignettes, producing pleasing output. Refer to
[Writing package vignettes](http://cran.fhcrc.org/doc/manuals/R-exts.html#Writing-package-vignettes)
for technical details. See the
[BiocStyle](/packages/devel/bioc/html/BiocStyle.html) package for a
convenient way to use common macros and a standard style.

A vignette provides reproducibility: the vignette produces the same
results as copying the corresponding commands into an R session. It is
therefore essential that the vignette embed R code between <<>>= and
@; short-cuts (e.g., using a LaTeX verbatim environment, or using the
Sweave eval=FALSE flag, or equivalent tricks in markdown) undermine
the benefit of vignettes.

All packages are expected to have at least one vignette.
Vignettes go in the `vignettes` directory of the package.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="citations">Citations</h2>

Appropriate citations must be included in help pages (e.g., in the see also
section) and vignettes; this aspect of documentation is no different from any
scientific endeavor. The file inst/CITATION can be used to specify how a
package is to be cited.

Whether or not a CITATION file is present, an automatically-generated
citation will appear on the package landing page on the
Bioconductor web site. For optimal formatting of author names
(if a CITATION file is not present),
specify the package author and maintainer using the Authors@R
field as described in
[Writing R Extensions](http://cran.r-project.org/doc/manuals/r-release/R-exts.html#The-DESCRIPTION-file).


<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="versions">Version Numbering</h2>

All Bioconductor packages use an x.y.z version scheme. The following rules
apply:

* x is usually 0 for packages that have not yet been released.
* y is even for packages in release, and odd for packages in devel.
* z is incremented whenever committing changes to a package.

When first submitted to Bioconductor, a package usually has version 0.99.0.
For more details, see
[Version Numbering](../version-numbering/)

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="c-code">C or Fortran code</h2>

If the package contains C or Fortran code, it should adhere to the standards
and methods described in the
[System and foreign language interfaces](http://cran.r-project.org/doc/manuals/R-exts.html#System-and-foreign-language-interfaces)
section of the Writing R Extensions manual. In particular:

* Use internal R functions, e.g., R_alloc and random number
  generators, over system supplied ones.
* Use C function registration (See the [Registering native
  routines](http://cran.fhcrc.org/doc/manuals/R-exts.html#Registering-native-routines))
* Use R_CheckUserInterrupt in C level loops when there is a chance
  that they may not terminate for certain parameter settings or when
  their run time exceeds 10 seconds with typical parameter settings,
  and the method is intended for interactive use.
* Make judicious use of Makevars and Makefile within a package. These
  are often not required at all (See the [Configure and
  cleanup](http://cran.fhcrc.org/doc/manuals/R-exts.html#Configure-and-cleanup)).
<a name="development-compiler-settings"></a>
* During package development, enable all warnings and disable
  optimizations. If you plan to [use a
  debugger](/developers/how-to/c-debugging/), tell the compiler to
  include debugging symbols. The easiest way to enforce these is to
  create a user-level Makevars file user's home directory in a
  sub-directory called '.R'). See examples below for flags for common
  toolchains. Consult the Writing R Extensions Manual for [details
  about Makevars
  files](http://cran.r-project.org/doc/manuals/R-exts.html#Using-Makevars).

  - Example for gcc/g++:

        CFLAGS=-Wall -Wextra -pedantic -O0 -ggdb
        CXXFLAGS=-Wall -Wextra -pedantic -O0 -ggdb
        FFLAGS=-Wall -Wextra -pedantic -O0 -ggdb

  - Example for clang/clang++:

        CFLAGS=-Weverything -O0 -g
        CXXFLAGS=-Weverything -O0 -g
        FFLAGS=-Wall -Wextra -pedantic -O0 -g

<h4 id="third-party-code">Third-party code</h4>

Use of external libraries whose functionality is redundant with
libraries already supported is strongly discouraged. In cases where
the external library is complex the author may need to supply
pre-built binary versions for some platforms.

By including third-party code a package maintainer assumes
responsibility for maintenance of that code. Part of the maintenance
responsibility includes keeping the code up to date as bug fixes and
updates are released for the mainline third-party project.

For guidance on including code from some specific third-party sources,
see the [external code sources
section](/developers/how-to/mavericks-howto/#external-code-sources) of
the C++ Best Practices guide.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="unitTests">Unit Tests</h2>

Unit tests are highly recommended.  We find them indispensable for 
both package development and maintenance.  Examples and explanations are provided 
[here](http://www.bioconductor.org/developers/unitTesting-guidelines).


<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="videos">Videos</h2>

You can submit an instructional video along with 
your package.
In the DESCRIPTION file of your package, add a "Video:" line
which contains the link to your video. We will then feature your 
video on our
[Bioconductor YouTube Channel](https://www.youtube.com/user/bioconductor).

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>



<h2 id="duplications">Duplication of Packages in CRAN and Bioconductor</h2>

Authors are strongly discouraged from placing their package into both CRAN and
Bioconductor. This avoids burdening the author with extra work and confusing
the user.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="responsibilities">Package Author and Maintainer Responsibilities</h2>

Acceptance of packages into Bioconductor brings with it ongoing
responsibility for package maintenance. These responsibilities
include:

* Subscription to the [bioc-devel](/help/mailing-list/) mailing list.
* Response to bug reports and questions from users regarding your
  package, as posted on the Bioconductor mailing list or directly to
  developers. Add a `BugReports:` field to the DESCRIPTION file if
  reports should be directed to a particular web page rather than the
  package maintainer.
* Package maintenance through software release cycles, including
  prompt updates to software and documentation necessitated by
  underlying changes in R.

All authors mentioned in the package DESCRIPTION file are entitled to
modify package source code. Changes to package authorship require
consent of all authors.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>
