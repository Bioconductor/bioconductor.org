# Package Guidelines

_version 1.0 - 2021-04-08_

This page gives details concerning guiding principles and formatting
required for _Bioconductor_ packages. See also [Package Submission][]
for an overview of the submission process and what is expected as a
_Bioconductor_ maintainer.

[Package Submission]: http://bioconductor.org/developers/package-submission/

<a name="top"></a>

- [Introduction](#intro)
- [General Package Development](#devel)

  * [Versions of R and _Bioconductor_](#version)
  * [Correctness, Space, and Time](#correctness)
  * [R CMD check environment](#checkingenv)

- [DESCRIPTION](#description)
- [NAMESPACE](#namespace)
- [NEWS](#news)
- [CITATION](#citation)
- [Including Data](#data)
- [Package Documentation](#documentation)

  * [Vignettes](#vignettes)
  * ['man' pages](#manpages)

- [Unit Tests](#unittest)
- [R Code](#rcode)
- [C/Fortran Code](#ccode)
- [.gitignore](#gitignore)
- [Conclusion](#closing)

<a name="intro"></a>

## Introduction

The _Bioconductor_ project promotes high-quality, well documented, and
interoperable software. These guidelines help to achieve this
objective; they are not meant to put undue burden on package authors,
and authors having difficultly satisfying guidelines should seek
advice on the [bioc-devel](/help/mailing-list/) mailing list.

Package maintainers are urged to follow these guidelines as closely as
possible when developing _Bioconductor_ packages.

General instructions for producing packages can be found in the
[Writing R Extensions][] manual, available from within R
(`RShowDoc("R-exts")`) or on the [R web site][].

[Writing R Extensions]: http://cran.r-project.org/doc/manuals/R-exts.html
[R web site]: http://cran.fhcrc.org/manuals.html

Remember these are the minimum requirements for package acceptance and
the package will still be subject to other guidelines below and a
formal technical review by trained _Bioconductor_ package reviewer.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="devel"></a>

## 1. General Package Development

<a name="version"></a>

### 1.1 Version of _Bioconductor_ and R

Package developers should always use the
[devel version of _Bioconductor_][1] when developing and testing
packages to be contributed.

Depending on the R release cycle, using _Bioconductor_ devel may or may
not involve also using the devel version of R. See the how-to on
[using devel version of _Bioconductor_][2] for up-to-date information.

[1]: /packages/devel/BiocViews.html#___Software
[2]: /developers/how-to/useDevel/

<a name="correctness"></a>

### 1.2 Correctness, Space and Time

1.2.1 _Bioconductor_ packages must minimally pass `R CMD build` (or `R CMD
   INSTALL --build`) and pass `R CMD check` with no errors and no
   warnings using a recent R-devel. Authors should also try to
   address all notes that arise during build or check.[^1]

1.2.2 Packages must also pass `BiocCheckGitClone()` and `BiocCheck()` with no
   errors and no warnings. The [BiocCheck][] package is a set of tests that
   encompass _Bioconductor_ Best Practices. Every effort should be made
   to address any notes that arise during this build or check.[^1]

1.2.3 Do not use filenames that differ only in case, as not all file
   systems are case sensitive.

1.2.4 The source package resulting from running `R CMD build` should
   occupy less than 5MB on disk.

1.2.5 The package should require less than 10 minutes to run `R CMD check
   --no-build-vignettes`. Using the `--no-build-vignettes` option
   ensures that the vignette is built only once.[^2]

1.2.6 Vignette and man page examples should not use more than 3GB of
   memory since R cannot allocate more than this on 32-bit Windows.

1.2.7 For software packages, individual files must be <= 5MB. This
   restriction exists even after the package is accepted and added to
   the `_Bioconductor_` repository.

1.2.8 The raw package directory should not contain unnecessary files,
   system files, or hidden files such as .DS_Store, .project, .git,
   cache file, log files, .Rproj, .so, etc.. These files may be
   present in your local directory but should not be commited to git
   (see [.gitignore](#gitignore)).

[BiocCheck]: https://bioconductor.org/packages/devel/BiocCheck

[^1]: The _Bioconductor_ team member assigned to review the package during the submission process will expect all ERROR, WARNINGS, and NOTES to be addressed. If there are any remaining, a justification of why they are not corrected will be expected.

[^2]: This is true for Software Packages. Experiment Data, Annotation, and Workflow packages are allowed additional space and check time.

<a name="checkingenv"></a>

### 1.3 R CMD check environment

It is possible to activate or deactivate a number of options in `R CMD build`
and `R CMD check`. Options can be set as individual environment variables or
they can be [listed in a file](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Checking-and-building-packages).
Descriptions of all the different options available can be found [here](https://cran.r-project.org/doc/manuals/r-devel/R-ints.html#Tools).
_Bioconductor_ has chosen to customize some of these options for incoming
submission during `R CMD check`. The file of utilized flags can be downloaded
from
[Github](https://github.com/Bioconductor/packagebuilder/blob/master/check.Renviron). The
file can either be place in a default directory as directed
[here](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Checking-and-building-packages)
or can be set through environment variable `R_CHECK_ENVIRON` with a command
similar to

```
export R_CHECK_ENVIRON = <path to downloaded file>
```

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="description"></a>

### 1.4 README file

If a README file is included in the package or the github repository and it
provides also installation instructions, those instruction should also include
_Bioconductor_ installation instructions.

## 2 DESCRIPTION file

The DESCRIPTION file must be properly formatted. The following section
will review some important notes regarding DESCRIPTION fields and
associated files.

2.1 "Package:" field: This is the name of the package. This should
   match the github repository name and is case sensitive. A package
   name should be descriptive and not already exist as a current
   package (case-insensitive) in [_Bioconductor_][BioC] or
   [CRAN][cran]. Avoid names that are easily confused with existing
   package names, or that imply a temporal (e.g., `ExistingPackage2`)
   or qualitative (e.g., `ExistingPackagePlus`) relationship. An easy
   way to check whether your name is already in use is to check that
   the following command fails

    ```
    if (!requireNamespace("BiocManager"))
	install.packages("BiocManager")
    BiocManager::install("MyPackage")
    ```

2.2 "Title:" field: Have a brief but descriptive Title

2.3 "Version:" field: All _Bioconductor_ packages use an x.y.z version
   scheme. See [Version Numbering][vernum] for specifics to how the
   release and devel _Bioconductor_ versioning proceeds. When first
   submitted to _Bioconductor_, a package should have pre-release
   version 0.99.0. The following rules apply:

* x is usually 0 for packages that have not yet been released.
* y is even for packages in release, and odd for packages in
  devel. Generally, do not bump this number especially in
  pre-release.
* z is incremented whenever committing changes to a package.

2.4 "Description:" field: The description should be a relatively short
   but detailed overview of what the package functionality entails. It
   should be at least three complete sentences.

2.5 "Authors@R:" field: The `Authors@R` field should be used. A maintainer
   designation (`cre` for Authors@R ) is required with an actively maintained
   email. This email will be used for contact regarding any issues that arise
   with your package in the future. For persons with an ORCID identifier (See 
   [ORCiD](https://orcid.org/) for more information) provide the identifier via
   an element named "ORCID" in the comment argument of person(). Example: 
   `person("Lori", "Shepherd", email=Lori.Shepherd@roswellpark.org,
   role=c("cre", aut"), comment = c(ORCID = "0000-0002-5910-4010"))`.

   Only one person should be listed as the `Maintainer` to ensure a single point
   of contact. This person by default will have commit access to the git
   repository on git.bioconductor.org. Commit access can be given to other
   developers by request on the `bioc-devel` mailing list.  Another option is
   to add collaborators to the github repository. This approach enables
   development by many but restricts push access to git.bioconductor.org.

2.6 "License:" field: should preferably refer to a standard license
   (see [wikipedia][wikiLic]) using one of R's standard
   specifications. Be specific about any version that applies (e.g.,
   GPL-2). Licenses restricting use, e.g., to academic or non-profit
   researchers, are not suitable for _Bioconductor_. Core
   _Bioconductor_ packages are typically licensed under
   Artistic-2.0. To specify a non-standard license, include a file
   named LICENSE in your package (containing the full terms of your
   license) and use the string "file LICENSE" (without the double
   quotes) in this "License:" field. The package should contain only
   code that can be redistributed according to the package license. Be
   aware of the licensing agreements for packages you are depending on
   in your package.  Not all packages are open source even if they are
   publicly available.

2.7 "LazyData:" field: For packages that include data, we recommend not
   including `LazyData: TRUE`. This rarely proves to be a good
   thing. In our experience it only slows down the loading of packages
   with large data.

2.8 "Depends/Imports/Suggests/Enhances:" fields:

* All packages must be available via [_Bioconductor_][BioC] or
  [CRAN][cran]; the use of `Remotes` is not supported hence dependencies
  only available on e.g. github are not allowed.
* Reuse, rather than re-implement or duplicate, well-tested
  functionality from other packages. Make use of appropriate
  existing packages (e.g., biomaRt, AnnotationDbi, Biostrings) and
  classes (e.g., SummarizedExperiment, GRanges, Rle,
  DNAStringSet), and avoid duplication of functionality available
  in other _Bioconductor_ packages. See
  [Common _Bioconductor_ Methods and Classes][CommonMethods].
  _Bioconductor_ Reviewers are very strict on this point! New
  packages should be interoperable with existing _Bioconductor_
  classes and not reimplement functionality especially with
  regards to importing/reading data.
* A package can be listed only once between
  Depends/Imports/Suggests/Enhances. Determine placement of
  package based on the following guidelines:

  + **Imports:** is for packages that provide functions, methods,
    or classes that are used inside your package name space. Most
    packages are listed here.
  + **Depends:** is for packages that provide essential
    functionality for users of your package, e.g., the
    `GenomicRanges` package is listed in the Depends: field of
    `GenomicAlignments`.  It is unusual for more than three
    packages to be listed as 'Depends:'.
  + **Suggests:** is for packages used in vignettes, examples,
    and in conditional code. Commonly, annotation and experiment
    packages (e.g., `TxDb*`) used in vignette and example code are
    included in this field thus avoiding users a costly download.
    In the case where an external one-off function is required for
    package code, external package availability can be checked via
    `if (!requireNamespace('extraPKG')) stop(...)`.
  + **Enhances:** is for packages such as `Rmpi` or `parallel`
    that enhance the performance of your package, but are not
    strictly needed for its functionality.

* It is seldom necessary to specify _R_ or specific versions as
  dependencies, since the _Bioconductor_ release strategy and
  standard installation instructions guarantee these
  constraints. Repositories mirrored outside _Bioconductor_ should
  include branches for each _Bioconductor_ release, and may find it
  useful to fully specify versions to enforce constraints otherwise
  guaranteed by _Bioconductor_ installation practices.

2.9 "SystemRequirements:" field: This field is for listing any external
   software which is required, but not automatically installed by the
   normal package installation process. If the installation process is
   non-trivial, a top-level README file should be included to document
   the process.

2.10 "biocViews:" field: REQUIRED! Specify at least two leaf node
    [biocViews][BioC]. Multiple leaf terms are encouraged but
    terms must come from the same trunk or package type (i.e.,
    Software, AnnotationData, ExperimentData, or Workflow).

2.11 "BugReports:" field: It is encouraged to include the relevant
    links to Github for reporting Issues.

2.12 "URL:" field: This field directs users to source code
    repositories, additional help resources, etc; details are provided
    in "Writing _R_ Extensions", `RShowDoc("R-exts")`.

2.13 "Video:" field: This field displays links to instructional videos.

2.14 "Collates:" field: This may be necessary to order class and method
    definitions appropriately during package installation.

2.15 "BiocType" field: This is required if submitting a `Docker` or
    `Workflow`. Otherwise this field could optionally define the type of
    Bioconductor package `Software`, `ExperimentData`, `Annotation`.

[BioC]: https://bioconductor.org/packages/release/BiocViews.html
[cran]: https://cran.r-project.org/web/packages/available_packages_by_name.html
[wikiLic]: http://en.wikipedia.org/wiki/Comparison_of_free_software_licences
[vernum]: https://bioconductor.org/developers/how-to/version-numbering/
[CommonMethods]: https://bioconductor.org/developers/how-to/commonMethodsAndClasses/

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="namespace"></a>

## 3 NAMESPACE

A [Namespace][namespace] file defines the functions, classes, and
methods that are imported into the name space, and exported for
users. _Bioconductor_ reviewers will be looking for:

3.1 Exported functions should use camel case or underscoring and not
   include "." indicate S3 dispatch.

3.2 Generally `importFrom()` is encouraged over importing an entire
   package, however if there are many functions from a single package,
   `import()` is okay.

3.3 Exporting all functions with `exportPattern("^[[:alpha:]]+")` is strongly discouraged.

[namespace]: http://cran.fhcrc.org/doc/manuals/R-exts.html#Package-namespace

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="news"></a>

## 4 NEWS

A NEWS file should be included to keep track of changes to the code
from one version to the next. It can be a top level file or in the
inst/ directory. Only one NEWS file should exist. The following are acceptable
formats and locations:

:--|:--------------|:-----------------------------
1. |./inst/NEWS.Rd |  latex
2. |./inst/NEWS    |  formatted text see ?news
3. |./inst/NEWS.md |  mardown
4. |./NEWS.md      |  markdown
5. |./NEWS         |  formatted text see ?news


Specifics on formatting can be found on the help page
for `?news`.  _Bioconductor_ uses the NEWS
file to create the semi-annual release announcement. It must include list
elements and **cannot** be a plain text file. An example format:

```
Changes in version 0.99.0 (2018-05-15)
+ Submitted to Bioconductor

Changes in version 1.1.1 (2018-06-15)
+ Fixed bug. Begin indexing from 1 instead of 2
+ Made the following significant changes
  o added a subsetting method
  o added a new field to database
```

After you install your package, the following can be run to see if the NEWS is
properly formatted:

```
utils::news(package="<name of your package>")
```

The output should look similar to the following

```
Changes in version 1.1.1 (2018-06-15):

    o   Fixed bug. Begin indexing from 1 instead of 2

    o   Made the following significant changes
	o added a subsetting method
	o added a new field to database

Changes in version 0.99.0 (2018-05-15):

    o   Submitted to Bioconductor

```

If you get something like the following there are formatting ERRORS that need to
be corrected:

```
Version: 0.99.0
Date: 2018-05-15
Text: Submitted to Bioconductor

Version: 1.1.1
Date: 2018-06-15
Text: Fixed bug. Begin indexing from 1 instead of 2

Version: 1.1.1
Date: 2018-06-15
Text: Made the following significant changes o added a subsetting
	method o added a new field to database
```

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="citation"></a>

## 4 CITATION

Appropriate citations must be included in help pages (e.g., in the see
also section) and vignettes; this aspect of documentation is no
different from any scientific endeavor. The file `inst/CITATION` can
be used to specify how a package is to be cited. If this option is utilized, 
a maintainer can check proper formatting of the CITATION file by running
`readCitationFile("inst/CITATION")`; This must run without ERROR for the 
CITATION to be accurately displayed on the package landing pages. 

Whether or not a CITATION file is present, an automatically-generated
citation will appear on the package landing page on the _Bioconductor_
web site. For optimal formatting of author names (if a CITATION file
is not present), specify the package author and maintainer using the
Authors@R field as described in [Writing R Extensions][].

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="data"></a>

## 5 Including Data

An excellent practice is to develop a software package, and to provide
or use an existing [experiment data package][expDpkg],
[annotation data][annDpkg] or data in the [ExperimentHub][exphub] or
[AnnotationHub][annhub] to give a comprehensive illustration of the
methods in the software package.

If existing data is not available or applicable, or a new smaller
dataset is needed for examples in the package, data can be included
either as a separate data package (for larger amounts of data) or
within the package (for smaller datasets).

### 5.1 Additional Experiment Data Package

Experimental data packages contain data specific to a particular
analysis or experiment. They often accompany a software package for
use in the examples and vignettes and in general are not updated
regularly. If you need a general subset of data for workflows or
examples first check the AnnotationHub resource for available files
(e.g., BAM, FASTA, BigWig, etc.). _Bioconductor_ strongly encourages creating an
experiment data package that utilizes ExperimentHub or AnnotationHub
(See [Creating an Experiment Hub Package][createHubExp] or
[Creating an Annotation Hub Package][createHubAnn]) but a traditional
package that encapsulates the data is also okay. See the Package
Submission package for submitting related packages.

### 5.2 Adding Data to Existing Package

_Bioconductor_ strongly encourages the use of existing datasets but if
not available data can be included directly in the package for use in
the examples found in man pages, vignettes, and tests of your
package. This is a good reference by Hadley Wickham
[concerning data][wickhamData]. As mentioned _Bioconductor_ however
does not encourage using `LazyData: True` despite its recommendataion
in this article. Some key points are summarized below.

### 5.3 Exported Data and the `data/` Directory

Data in `data/` is exported to the user and readily available. It is
made available in an R session through the use of `data()`.  It will
require documentation concerning its creation and source
information. It is most often a `.RData` file created with `save()`
but other types are acceptible as well, see `?data()`. Please remember
to compress the data.

### 5.4 Raw Data and the `inst/extdata/` Directory

It is often desirable to show a workflow which involves parsing or
loading of raw files. _Bioconductor_ recommends finding existing raw
data already provided in another package or the hubs, however if this
is not applicable, raw data files should be included in the
`inst/extdata`. Files of these type are often accessed utilizing
`system.file()`. _Bioconductor_ requires documentation on these files in
an `inst/script/` directory.

### 5.5 Internal Data

Rarely, a package may require parsed data that is used internal but
should not be exported to the user. An `R/sysdata.rda` is often the
best place to include this type of data.


[expDpkg]: https://bioconductor.org/packages/release/BiocViews.html#___ExperimentData
[exphub]: https://bioconductor.org/packages/release/bioc/html/ExperimentHub.html
[annDpkg]: https://bioconductor.org/packages/release/BiocViews.html#___AnnotationData
[annhub]: https://bioconductor.org/packages/release/bioc/html/AnnotationHub.html
[createHubExp]: https://bioconductor.org/packages/release/bioc/vignettes/ExperimentHub/inst/doc/CreateAnExperimentHubPackage.html
[createHubAnn]: https://bioconductor.org/packages/release/bioc/vignettes/AnnotationHub/inst/doc/CreateAnAnnotationPackage.html
[wickhamData]: http://r-pkgs.had.co.nz/data.html

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="documentation"></a>

## 6 Package Documentation

Package documentation is important for users to understand how to work
with your code. _Bioconductor_ requires a [vignette][vig] with
executable code that demonstrates how to use the package to accomplish
a task, [man pages][man] for all exported functions with runnable
examples, well documented data structures especially if not a
[pre-exiting class][preclass], and well documented datasets for data
in `data` and in `inst/extdata`. References to the methods used as
well as to other simlar or related project/packages is also
expected. If data structures differ from similar packages,
_Bioconductor_ reviewers will expect some justification as to why. Keep
in mind it is always possible to extend existing classes.

<a name="vignettes"></a>

### 6.1 Vignettes

A vignette demonstrates how to accomplish non-trivial tasks embodying
the core functionality of your package. There are two common types of
vignettes. A _Sweave_ vignette is an .Rnw file that contains LaTeX and
chunks of R code. The R code chunk starts with a line <<>>=, and ends
with @. Each chunk is evaluated during `R CMD build`, prior to LaTeX
compilation to a PDF document. An _R markdown_ vignette is similar to
a Sweave vignette, but uses [markdown][] instead of LaTeX for
structuring text sections and resulting in HTML output. The [knitr][]
package can process most Sweave and all R markdown vignettes,
producing pleasing output. Refer to [Writing package vignettes][] for
technical details. See the [BiocStyle][] package for a convenient way
to use common macros and a standard style.

[markdown]: http://daringfireball.net/projects/markdown/
[knitr]: http://yihui.name/knitr/
[Writing package vignettes]: http://cran.fhcrc.org/doc/manuals/R-exts.html#Writing-package-vignettes
[BiocStyle]: /packages/devel/bioc/html/BiocStyle.html

A vignette provides reproducibility: the vignette produces the same
results as copying the corresponding commands into an R session. It is
therefore **essential** that the vignette embed executed R
code. short-cuts (e.g., using a LaTeX verbatim environment, or using
the Sweave eval=FALSE flag, or equivalent tricks in markdown)
undermine the benefit of vignettes and are **generally not allowed**;
exceptions can be made with proper justification and are at the
_Bioconductor_ Reviewers discretion.

All packages are required to have at least one vignette.  Vignettes go
in the `vignettes` directory of the package. Vignettes are often used
as stand-alone documents, so best practices are to include an
informative _title_, the primary _author_ of the vignette, the _last
modified date_ of the vignette, and a link to the package landing
page. We encourage the use of BiocSytle for formatting.

Some best coding practices for Biocondcutor vigenttes are as follow:

6.1.1 Add an "Introduction" section that serves as an abstract to
   introduce the objective, models, unique functions, key points, etc
   that distinguish the package from other packages of similar type.

6.1.2 Add an "Installation" section that show to users how to download
   and load the package from _Bioconductor_.

6.1.3 If appropriate, we strongly encourage a table of contents

6.1.4 Non-trival executable code is a must!!! Static vignettes are not
   acceptable.

6.1.5 Include a section with the `SessionInfo()`

6.1.6 Only the vignette file (.Rnw or .Rmd) and any necessary static
   images should be in the vignette directory. No intermediate files
   should be present.

6.1.7 Remember to include any relavent references to methods.

<a name="manpages"></a>

### 6.2 'man' Pages

See Writing R Extensions section on [man pages][man] for detailed instruction or
format information for documenting a package, functions, classes, and data
sets. All help pages should be comprehensive.

6.2.1 All exported functions and classes need will have a man page. Man pages
   describing new classes must be very detailed on the structure and what
   type of information is stored.

6.2.2 _Bioconductor_ encourages having a package man page with an overview of
   the package and links to the main functions.

6.2.3 Data man pages must include source information and data structure 
   information.  
   
6.2.4 All man pages should have an runnable examples. `donttest and dontrun` 
   are discouraged and generally not allowed; exceptions can be made with 
   proper justification and are at the _Bioconductor_ Reviewers discretion. If 
   this option is used it will also be preferrable to use `donttest` instead 
   of `dontrun`; `donttest` requires valid R code while `dontrun` does not. 

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

### 6.3 inst/script/

The scripts in this directory can vary. Most importantly if data was
included in the `inst/extdata/`, a related script must be present in
this directory documenting very clearly how the data was generated. It
should include source urls and any important information regarding
filtering or processing. It can be executable code, sudo code, or a
text description. A user should be able to download and be able to
roughly reproduce the file or object that is present as data.

[vig]: http://cran.fhcrc.org/doc/manuals/R-exts.html#Writing-package-vignettes
[man]: http://cran.fhcrc.org/doc/manuals/R-exts.html#Rd-format
[preclass]: http://bioconductor.org/developers/how-to/commonMethodsAndClasses/

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="unittest"></a>

## 7 Unit Tests

Unit tests are highly recommended.  We find them indispensable for
both package development and maintenance. Two of the main frameworks
for testing are `RUnit` and `testthat`. Examples and explanations are
provided [here][]. There is also the opportunity to create a full testing suite
more in depth than traditional testing guidelines but this will require the use
of [long tests][]. If a package developer is considering the use of long tests
we highly recommend reaching out on the bioc-devel mailing list to ensure proper
use and justification.

[here]: http://bioconductor.org/developers/unitTesting-guidelines
[long tests]: /developers/how-to/long-tests

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="rcode"></a>

## 8 R Code and Best Practices

Everyone has their own coding style and formats. There are however
some best practice guidelines that _Bioconductor_ will look for (see
[coding style][style]). There are also some other key points:

8.1 Only contain code that can be distributed under the license
   specified (see also 2.6).

8.2 Many common coding and sytax issues are flagged in `R CMD check`,
   and `BiocCheck()`. (see the `R CMD check` [cheatsheet][] and
   [BiocCheck][BiocCheck-vignette] vignette. Some of
   the more promenient offenders:

+ Use `vapply()` instead of `sapply()` and use the various apply
  functions instead of for loops.
+ Use `seq_len()` or `seq_along()` instead of `1:...`
+ Use TRUE/FALSE instead of T/F
+ Use of numeric indices (rather than robust named indices)
+ Avoid `class()==` and `class()!=` instead use `is()`
+ Use `system2()` instead of `system`
+ Do not use `set.seed` in any internal R code.
+ No `browser()` calls should be in code
+ Avoid the use of `<<-`.
+ Avoid use of direct slot access with `@` or `slot()`. Accessor
  methods should be created and utilized
+ Downloading of external data from unsanctioned providers such
  as GitHub, Dropbox, etc. We recommend the use of `*Hub` packages.

8.3 Some additional formatting and syntax guidelines

+ Use `<-` instead of `=` for assignment
+ Function names should be camel case or utilize the underscore
  `_` and not have a dot `.` which indicates S3 dispatch.
+ Use `dev.new()` to start a graphics drive if necessary. Avoid
  using `x11()` or `X11()` for it can only be called on machines
  that have access to an X server.
+ Use `message`/`warning`/`error` and not `cat` method (except for
  customized `show` methods). `paste0` should generally not be used in these
  methods except for collapsing multiple values from a variable.

8.4 Avoid re-implementing functionality or classes (see also 2.8). Make use of
   appropriate existing packages (e.g., biomaRt, AnnotationDbi,
   Biostrings, GenomicRanges) and classes (e.g., `SummarizedExperiment`,
   `AnnotatedDataFrame`, `GRanges`, `DNAStringSet`) to avoid duplication of
   functionality available in other _Bioconductor_ packages. See also
   [Common _Bioconductor_ Methods and Classes][CommonMethods]. This
   encourages interoperability and simplifies your own package
   development. If new representation is needed, see the
   [Essential S4 interface][RECEssentialS4] section of
   [Robust and Efficient Code][]. In general, Bioconductor will insist on
   interoperability with [Common Classes][CommonMethods] for acceptance.

8.5 Avoid large chunks of repeated code. If code is being repeated this
   is generally a good indication a helper function could be
   implemented.

8.6 Excessively long functions should also be avoided. Write small
   functions.  It's best if each function has only one job that needs
   to do.  And it's also best if it does that job in as few lines of
   code as possible.  If you find yourself writing great big functions
   that wrap on for more than a screen then you should probably take a
   moment to split it up into smaller helper functions.  Smaller
   functions are easier to read, debug and to reuse.

8.7 Argument names to functions should be descriptive and well
   documented. Arguments should generally have default values. Check
   arguments against a validity check.

8.8 Vectorize! Many R operations are performed on the whole object, not
   just the elements of the object (e.g., `sum(x)`, not `x[1] + x[2] +
   x[2] + ...`).  In particular, relatively few situations require an
   explicit `for` loop. See the [Vectorize][RECVectorize] section of
   [Robust and Efficient Code][] for additional detail.

8.9 Follow guiding principles on [Querying Web Resources][webre] if applicable

8.10 For parallel implementation please use
    [BiocParallel](/packages/devel/BiocParallel). See also the
    [Parallel Recommendations][RECPara] section of
    [Robust and Efficient Code][]. A minimal number of cores (1 or 2) should be
    set as a default.

[cheatsheet]: http://r-pkgs.had.co.nz/check.html
[BiocCheck-vignette]:  https://bioconductor.org/packages/devel/bioc/vignettes/BiocCheck/inst/doc/BiocCheck.html

[style]: http://bioconductor.org/developers/how-to/coding-style/
[webre]: http://bioconductor.org/developers/how-to/web-query/
[RECVectorize]: /developers/how-to/efficient-code/#vectorize
[Robust and Efficient Code]: /developers/how-to/efficient-code
[RECEssentialS4]: /developers/how-to/efficient-code/#essential-s4-interface
[RECPara]: /developers/how-to/efficient-code/#parallel-recommendations

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="ccode"></a>

## 9 C or Fortran code

If the package contains C or Fortran code, it should adhere to the
standards and methods described in the
[System and foreign language interfaces][interface] section of the Writing R
Extensions manual. In particular:

9.1 Use internal R functions, e.g., `R_alloc` and random number
   generators, over system supplied ones.
   
9.2 Use C function registration (See the [Registering native routines][]).

9.3 Use R_CheckUserInterrupt in C level loops when there is a chance
   that they may not terminate for certain parameter settings or when
   their run time exceeds 10 seconds with typical parameter settings,
   and the method is intended for interactive use.

9.4 Make judicious use of Makevars and Makefile within a package. These
   are often not required at all (See the
   [Configure and cleanup][candc]).

9.5 During package development, enable all warnings and disable
   optimizations. If you plan to [use a debugger][debug], tell the
   compiler to include debugging symbols. The easiest way to enforce
   these is to create a user-level Makevars file user's home directory
   in a sub-directory called '.R'). See examples below for flags for
   common toolchains. Consult the Writing R Extensions Manual for
   [details about Makevars files][makevars].

  - Example for gcc/g++:
```
    CFLAGS=-Wall -Wextra -pedantic -O0 -ggdb
	CXXFLAGS=-Wall -Wextra -pedantic -O0 -ggdb
	FFLAGS=-Wall -Wextra -pedantic -O0 -ggdb
```
  - Example for clang/clang++:
```
	CFLAGS=-Weverything -O0 -g
	CXXFLAGS=-Weverything -O0 -g
	FFLAGS=-Wall -Wextra -pedantic -O0 -g
```
[Registering native routines]: http://cran.fhcrc.org/doc/manuals/R-exts.html#Registering-native-routines

<a name="third-party-code"></a>

## 10 Third-party code

Use of external libraries whose functionality is redundant with
libraries already supported is strongly discouraged. In cases where
the external library is complex the author may need to supply
pre-built binary versions for some platforms.

By including third-party code a package maintainer assumes
responsibility for maintenance of that code. Part of the maintenance
responsibility includes keeping the code up to date as bug fixes and
updates are released for the mainline third-party project.

For guidance on including code from some specific third-party sources,
see the [external code sources section][external] of the C++ Best
Practices guide.

[debug]: /developers/how-to/c-debugging/
[candc]: http://cran.fhcrc.org/doc/manuals/R-exts.html#Configure-and-cleanup
[interface]: http://cran.r-project.org/doc/manuals/R-exts.html#System-and-foreign-language-interfaces
[makevars]: http://cran.r-project.org/doc/manuals/R-exts.html#Using-Makevars
[external]: /developers/how-to/mavericks-howto/#external-code-sources

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="gitignore"></a>

## 11 The .gitignore File

_Bioconductor_ requires a git repository for submission.  There are
certain system files that should not be git tracked and are
unacceptable to include. These files can remain on a local system but
should be excluded from the git repository which is possible by
including a `.gitignore` file.

The following are files that are checked by _Bioconductor_ and flagged
as unacceptable:

```
hidden_file_ext = (
    ".renviron", ".rprofile", ".rproj", ".rproj.user", ".rhistory",
    ".rapp.history", ".o", ".sl", ".so", ".dylib", ".a", ".dll",
    ".def", ".ds_store", "unsrturl.bst", ".log", ".aux", ".backups",
    ".cproject", ".directory", ".dropbox", ".exrc", ".gdb.history",
    ".gitattributes", ".gitmodules", ".hgtags", ".project", ".seed",
    ".settings", ".tm_properties"
)
```

<a name="closing"></a>

## Conclusion

The following exercise
[How to Build _Bioconductor_ Package with RStudio][] may also be
helpful.

Remember that every _Bioconductor_ package goes through a formal review
process and may still receive technical feedback from the assigned
_Bioconductor_ reviewer. An overview of the submission process may
be found [here][submission] and a package may be submitted to the
[new package tracker][tracker].

[How to Build _Bioconductor_ Package with RStudio]: http://bioconductor.org/developers/how-to/buildingPackagesForBioc/
[submission]: http://bioconductor.org/developers/package-submission/
[tracker]: https://github.com/Bioconductor/Contributions
