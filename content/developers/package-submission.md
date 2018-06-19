Package Submission
==================================================

This page gives an overview of the submission process along with key principles
to follow. See also [Package
Guidelines](http://bioconductor.org/developers/package-guidelines/) for package
specific guidelines and requirement and the [Bioconductor new package submission
tracker](https://github.com/Bioconductor/Contributions). 

<a name="top"></a>

[**Introduction**](#intro)

[**Types of Packages**](#type)

[**Author/Maintainer Expectations**](#author)

[**Submission**](#submission)  

[**Review Process**](#whattoexpect)

[**Following Acceptance**](#afteraccept)

[**Additonal Support**](#support)

<a name="intro"></a>

Introduction
------------

Bioconductor Packages should

* Address areas of high-throughput genomic analysis where
  Bioconductor already makes significant contributions, e.g.,
  sequencing, expression and other microarrays, flow cytometry, mass
  spectrometry, image analysis; see
  [biocViews](http://bioconductor.org/packages/devel/BiocViews.html#___Software).
* Interoperate with other Bioconductor packages, re-using common data
  structures ([S4 classes and methods][]) and existing infrastructure
  (e.g., `rtracklayer::import()` for input of common genomic files).
* Adopt software best practices that enable reproducible research and
  use, such as full documentation and vignettes (including fully
  evaluated code) as well as commitment to long-term user support
  through the Bioconductor [support site](https://support.bioconductor.org).
* Not exist on CRAN. A package can only be submitted to one or the other. 
* Comply with [Package
  Guidelines](http://bioconductor.org/developers/package-guidelines/)
* Your package cannot depend on any package (or version of a package) that is
  not (yet) available on CRAN or Bioconductor.


[S4 classes and methods]: http://bioconductor.org/developers/how-to/commonMethodsAndClasses/

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>






<a name="type"></a>

Types of Packages
------------

Bioconductor packages are broadly defined by three main package types: Software,
Experiment Data, and Annotation. Most packages contributed by users are
[software][software-pkgs] packages that perform analytic calculations. Users
also contribute [annotation][annotation-pkgs] and [experiment
data][exptdata-pkgs] packages.

Annotation packages are database-like packages that provide
information linking identifiers (e.g., Entrez gene names or Affymetrix
probe ids) to other information (e.g., chromosomal location, Gene
Ontology category). It is also encouraged to utilize AnnotationHub for
storage and access to large raw data files and their conversion to standard
R formats. Instructions for adding data to AnnotationHub and designing a
annotaiton package to use AnnotationHub can be found here: [Creating
AnnotationHub Packages][annoHowTo].

Experiment data packages provide data sets that are used, often by software
packages, to illustrate particular analyses. These packages contain curated
data from an experiment, teaching course or publication and in most cases
contain a single data set. It is also encouraged to utilize ExperimentHub for
storage and access to larger data files. ExperimentHub is also particularly
useful for hosting collections of related data sets. Instructions for adding
data to ExperimentHub and designing an experiment data package to use
ExperimentHub can be found here: [Creating ExperimentHub Packages][expHowTo].

See [Package Guidelines][] for details on package format and syntax.

[software-pkgs]: /packages/release/bioc/
[annotation-pkgs]: /packages/release/data/annotation/
[exptdata-pkgs]: /packages/release/data/experiment/
[annoHowTo]: /packages/devel/bioc/vignettes/AnnotationHub/inst/doc/CreateAnAnnotationPackage.html
[expHowTo]: /packages/devel/bioc/vignettes/ExperimentHub/inst/doc/CreateAnExperimentHubPackage.html
[Package Guidelines]: http://bioconductor.org/developers/package-guidelines/

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>






<a name="author"></a>

Author/Maintainer Expectations
------------

Acceptance of packages into Bioconductor brings with it ongoing
responsibility for package maintenance. These responsibilities
include:

* Be familiar with the ‘devel’ and ‘release’ branch concepts used in
  the project.  New packages and features are added to the ‘devel’
  branch. The current devel branch becomes the next release, with a
  release in April and October. Once your package has been accepted,
  it will initially be in the ‘devel’ branch. Most users are expected
  to use the release branch, so will not immediately have access to
  your package.
* Realize Bioconductor, unlike CRAN, maintains all package source code
  under git version control. This means that you make  changes to your package
  using [git][11].
* Package maintenance through software release cycles, including
  prompt updates to software and documentation necessitated by
  underlying changes in R.
* Subscription to the [bioc-devel](/help/mailing-list/) mailing list.
* [Registration](https://support.bioconductor.org/accounts/signup/)
  on the [support site](https://support.bioconductor.org/).
* Response to bug reports and questions from users regarding your
  package, as posted on the [Bioconductor support site](https://support.bioconductor.org/)
  or directly to developers. Add a `BugReports:` field to the DESCRIPTION file if
  reports should be directed to a particular web page rather than the
  package maintainer. You should register on the
  [support site](https://support.bioconductor.org/)
  and edit your profile, changing the "Watched Tags" field to
  include all packages you maintain, so you will be notified
  when anyone posts a question about your package.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>






<a name="submission"></a>

Submission
------------


* Submit by opening a new issue at the Bioconductor
  [Contributions](https://github.com/Bioconductor/Contributions/issues/new) repository.
* Read the [Contribution Guidelines](https://github.com/Bioconductor/Contributions/blob/master/CONTRIBUTING.md) for full instructions.
* Assuming your package is in a
  [GitHub Repository](https://help.github.com/articles/create-a-repo/), under a default 'master' branch, 
  add the link to your repository to the issue you are submitting.


### Experiment Data Packages ###
Experimental data packages contain data specific to a particular
analysis or experiment. They often accompany a software package for use
in the examples and vignettes and in general are not updated regularly.
If you need a general subset of data for workflows or examples first check the
AnnotationHub resource for available files (e.g., BAM, FASTA, BigWig, etc.).

If you have an associated data package for your software package, please do
*NOT* create a separate issue in the our tracker repository for that. Instead, please add the
data package repository to the same issue as the software package.
The process for doing this is documented
[here](https://github.com/Bioconductor/Contributions/blob/master/CONTRIBUTING.md#submitting-related-packages).

### Annotation Packages ###

Annotation packages contain lightly or non-curated data from a public
source and are updated with each Bioconductor release (every 6 months).
They are a source of general annotation for one or many organisms and
are not specific to a particular experiment.  When possible, they
should support the select() interface from AnnotationDbi.

Annotation packages should *NOT* be posted to the tracker repository.
Instead send an email to <packages@bioconductor.org> with a description
of the proposed annotation package and futher instructions of where to
send the package will be provided.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>







<a name="whattoexpect"></a>

Review Process
------------

* A new package is initially labeled as '1. awaiting moderation'.
  A _Bioconductor_ team member will take a very brief look at your
  package, to ensure that it is intended for _Bioconductor_. Appropriate
  packages will be re-labelled '2. review in progress' and a reviewer will be
  automatically assigned. Your assigned reviewer will address your concerns and
  help you through the review process. The entire review process typically takes
  between 2 and 5 weeks.

* The package will be submitted to the _Bioconductor_ build
  system. The system will check out your package from GitHub. It will
  then run `R CMD build` to create a 'tarball' of your source code,
  vignettes, and man pages. It will run `R CMD check` on the tarball,
  to ensure that the package conforms to standard _R_ programming best
  practices. Finally, the build system will run `R CMD BiocCheck` to
  ensure that the package conforms to _Bioconductor_ [BiocCheck][4]
  standards. The system will perform these steps using the
  ['devel' version](https://bioconductor.org/developers/how-to/useDevel/)
  of _Bioconductor_, on three platforms (Linux, Mac OS X, and
  Windows).  After these steps are complete, a link to a build report
  will be appended to the new package issue. Avoid surprises by
  running these checks on your own computer, under the 'devel' version
  of _Bioconductor_, before submitting your package.

* If the build report indicates problems, modify your package and
  commit changes to the default branch of your GitHub repository.  If
  there are problems that you do not understand, seek help on the
  [bioc-devel][9] mailing list.

* To trigger a new build, include a version bump in your commit, e.g.,
  from `Version: 0.99.0` to `Version: 0.99.1`. Pre-release versions utilize the
  `0.99.z` format. When accepted and released, your package's version number
  will be automatically incremented to 1.0.0.

* Once your package builds and checks without errors or (avoidable)
  warnings, a _Bioconductor_ team member will provide a technical
  review of your package.  Other _Bioconductor_ developers and users
  with domain expertise are encouraged to provide additional community
  commentary.  Reviewers will add comments to the issue you created.

* Respond to the issues raised by the reviewers. You _must_ respond to
  the primary reviewer, and are strongly encouraged to consider
  community commentary. Typically your response will involve code
  modifications; commit these to the default branch of your GitHub
  repository to trigger subsequent builds. When you have addressed all
  concerns, add a comment to the issue created in step 2 to explain
  your response.

* The reviewer will assess your responses, perhaps suggesting further
  modifications or clarification. The reviewer will then accept your
  package for inclusion in _Bioconductor_, or decline it. The label
  '2. review in progress' will be replaced by '3a. accepted' or
  '3b. declined'.

* If your package is accepted, it will be added to _Bioconductor_'s
  Git repository and to the nightly 'devel'
  builds. All packages in the 'devel' branch of the repository are
  'released' to the user community once every six months, in
  approximately April and October.

* Once the review process is complete, the issue you created will be
  closed. All updates to your package will be through the
  [Bioconductor Git Server][11].

[4]: https://bioconductor.org/packages/devel/bioc/html/BiocCheck.html
[9]: https://stat.ethz.ch/mailman/listinfo/bioc-devel
[11]: http://bioconductor.org/developers/how-to/git/

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>






<a name="afteraccept"></a>

Following Acceptance
------------

Following acceptance of a package:

* Packages accepted on the tracker repository are added to the 'devel' branch of
  the Bioconductor GIT repository, with the current version number of
  the accepted package.
* Packages are then built by the Bioconductor nightly build
  process. If the build is successful, the package has its own
  'landing page' created, and the package is made available to users
  of the 'devel' branch of Bioconductor via `biocLite()`.
* Developers may continue to make changes to their package, but now do
  so to the version in the [Bioconductor git server][11].
* Developers should bump the `z` portion of their version number every
  time they commit changes to their package, following the
  [Version numbering](/developers/how-to/version-numbering/) guidelines. If
  developers don't bump the version, the changes made to their package
  *do not propagate* to the Bioconductor web site and package
  repository.


<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>






<a name="support"></a>

Additional Support
------------

We are eager to enhance the quality and interoperability of
Bioconductor software and will provide additional support when
requested by package developers. Example areas of assistance include
use of appropriate S4 structures, specific guidance on efficient
implementation, guidance on code structure, and critical assessment of
package documentation and structure.  Use the
[bioc-devel](/help/mailing-list/) mailing list or email
<packages@bioconductor.org> to obtain additional support.

* Support Email: <packages@bioconductor.org>


<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>
