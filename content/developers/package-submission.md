# Package Submission

This page gives an overview of the submission process along with key
principles to follow. See also [Package Guidelines][guidelines] for
package specific guidelines and requirement and the
[_Bioconductor_ new package submission tracker][tracker].

[guidelines]: https://contributions.bioconductor.org/
[tracker]: https://github.com/Bioconductor/Contributions

<a name="top"></a>

- [Introduction](#intro)
- [Types of Packages](#type)
- [Package Naming Policy](#naming)
- [Author/Maintainer Expectations](#author)
- [Submission](#submission)
- [Experiment data package](#experPackage)
- [Annotation package](#annPackage)
- [Review Process](#whattoexpect)
- [Following Acceptance](#afteraccept)
- [Additional Support](#support)

<a name="intro"></a>

## Introduction

To submit a package to _Bioconductor_ the package should:

* Address areas of high-throughput genomic analysis, e.g., sequencing,
  expression and other microarrays, flow cytometry, mass spectrometry,
  image analysis; see [biocViews][].
* Interoperate with other _Bioconductor_ packages by _re-using common data
  structures_ (see [S4 classes and methods][]) and existing infrastructure
  (e.g., `rtracklayer::import()` for input of common genomic files).
* Adopt software best practices that enable reproducible research and
  use, such as full documentation and vignettes (including fully
  evaluated code) as well as commitment to long-term user support
  through the _Bioconductor_ [support site][support].
* Not exist on CRAN. A package can only be submitted to one or the other.
* Comply with [Package Guidelines][guidelines].
* Your package cannot depend on any package (or version of a package)
  that is not (yet) available on CRAN or _Bioconductor_.

[biocViews]: /packages/devel/BiocViews.html
[S4 classes and methods]: /developers/how-to/commonMethodsAndClasses/

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="type"></a>

## Types of Packages

_Bioconductor_ packages are broadly defined by three main package types:
[**Software**][software-pkgs], [**Experiment Data**][exptdata-pkgs], and
[**Annotation**][annotation-pkgs].

* Most packages contributed by users are software packages. Software packages
  provide implementation of algorithms (e.g. statistical analysis), access to
  resources (e.g. biomart, or NCBI) or visualizations (e.g. volcano plots,
  pathways plots). Instructions for creating Software packages can be found
  here: [Package guidelines][guidelines].

* [Annotation packages](#annPackage) are database-like packages that provide
  information linking identifiers (e.g., Entrez gene names or Affymetrix
  probe ids) to other information (e.g., chromosomal location, Gene
  Ontology category). It is also encouraged to utilize AnnotationHub for
  storage and access to large raw data files and their conversion to
  standard R formats. Instructions for adding data to AnnotationHub and
  designing a annotation package to use AnnotationHub can be found here:
  [Creating A Hub Packages][HubHowTo].

* [Experiment data packages](#experPackage) provide data sets that are used,
  often by software packages, to illustrate particular analyses. These packages
  contain curated data from an experiment, teaching course or publication and
  in most cases contain a single data set. It is also encouraged to utilize
  ExperimentHub for storage and access to larger data files. ExperimentHub is
  also particularly useful for hosting collections of related data sets.
  Instructions for adding data to ExperimentHub and designing an experiment data
  package to use ExperimentHub can be found here:
  [Creating A Hub Packages][HubHowTo].

See [Package Guidelines][guidelines] for details on package format and syntax.

[software-pkgs]: /packages/release/bioc/
[annotation-pkgs]: /packages/release/data/annotation/
[exptdata-pkgs]: /packages/release/data/experiment/
[HubHowTo]: /packages/devel/bioc/vignettes/HubPub/inst/doc/CreateAHubPackage.html

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>


<a name="naming"></a>

## Package Naming Policy

Package naming: i) Ownership of package name. Bioconductor follows [CRAN's
policy][] in requiring that contributors give the right to use the package name to Bioconductor at time of submission, so that the Bioconductor team can orphan the package and allow another maintainer to take it over in the event that the package contributor discontinues package maintenance. See Bioconductor's package [end-of-life policy][] for more details. ii) Uniqueness of package name. Packages should be named in a way that does not conflict (irrespective of case) with any current or past BIOCONDUCTOR package, nor any current CRAN package.

[CRAN's policy]: https://cran.r-project.org/web/packages/policies.html
[end-of-life policy]: http://bioconductor.org/developers/package-end-of-life/

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>


<a name="author"></a>

## Author / Maintainer Expectations

Acceptance of packages into _Bioconductor_ brings with it ongoing
package maintenance responsibilities. Package authors are expected to:

* Follow Bioconductor guidelines

  These include standard [guidelines](guidelines), [version numbering][versioning],
  coding style, code performance requirements, memory usage, using existing
  data classes, and the other requirements described below.

* Follow the release cycle of Bioconductor

  There are two releases each year, around April and October. The
  [release schedule][release-schedule] will indicate the timetables and
  deadlines for each release. A release cycle typically produces two
  versions of packages, ‘devel’ and ‘release’. It is important to be familiar
  with these branch concepts. Once your package has been accepted, it will
  initially be in the ‘devel’ branch.  The current devel branch becomes the
  next release. Most users are expected to use the release branch, so they will
  not immediately have access to your package until the next release.
  Bug fixes can be fixed in both branches, while new features should only
  be added to the ‘devel’ branch.

* Maintain the package using version control

  Realize that _Bioconductor_, unlike CRAN, maintains all package source code
  under git version control. This means that you make changes to your
  package using [git][11]. If your package is accepted, you will receive
  instructions with typical git operations (see the after
  [acceptance section](#afteraccept)).  Package maintenance through software
  release cycles, including prompt updates to software and documentation,
  is needed due to possible underlying changes in R and/or other package
  dependencies.

* Subscribe to the [bioc-devel](/help/mailing-list/) mailing list

  The Bioconductor team communicates with developers through this list.
  It is also a good channel to communicate changes to other developers.
  Addressing Bioconductor team requests in a timely manner guarantees that
  your package remains available through Bioconductor.

* [Register][support-register] on and use the [support site][support]

  The support site is the official support channel for users. Users and even
  developers may ask questions regarding your package on this platform.
  Be sure to include all the packages that you maintain in the "Watched Tags"
  section of your support site profile. This will notify you of any questions
  posted regarding your package(s).
  It is important to promptly respond to bug reports and questions either on
  the [_Bioconductor_ support site][support] post or directly to developers.
  Some maintainers prefer to indicate a `BugReports:` field in their package's
  DESCRIPTION file. This field indicates a particular web page for submitting
  bug reports and questions.

[release-schedule]: /developers/release-schedule/
[support-register]: https://support.bioconductor.org/accounts/signup/
[support]: https://support.bioconductor.org/

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="submission"></a>

## Submission

* Read and follow the full [Contributor Guidelines][guidelines] section and
  make yourself familiar with the
  [Bioconductor Developer](/developers/) pages.

* Submit by opening a new issue in the _Bioconductor_
  [Contributions][issues] repository, following the [guidelines][tracker] of
  the `README.md` file.
  Assuming that your package is in a [GitHub Repository][git-repo-create] and
  under the default branch, add the link to your repository to
  the issue you are opening.  You cannot specify any alternative branches; the
  default branch is utilized. The default branch must contain only package
  code. Any files or directories for other applications (Github Actions,
  devtool, etc) should be in a different branch.

[issues]: https://github.com/Bioconductor/Contributions/issues/new
[git-repo-create]: https://help.github.com/articles/create-a-repo/

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="experPackage"></a>

## Experiment Data Packages

Experimental data packages contain data specific to a particular
analysis or experiment. They often accompany a software package for
use in the examples and vignettes and in general are not updated
regularly.  If you need a general subset of data for workflows or
examples first check the AnnotationHub resource for available files
(e.g., BAM, FASTA, BigWig, etc.) or ExperimentHub for available processed
example data set already included in _Bioconductor_. If no current files or data
sets are appropriate consider an associated Experiment Data Package that
utilizes [ExperimentHub][].

If you have an associated data package for your software package,
please do *NOT* create a separate issue in the our tracker repository
for that. Instead, please add the data package repository to the same
issue as the software package.  The process for doing this is
documented [here][].

[ExperimentHub]: /packages/devel/bioc/vignettes/HubPub/inst/doc/CreateAHubPackage.html
[here]: https://github.com/Bioconductor/Contributions/blob/master/CONTRIBUTING.md#submitting-related-packages

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="annPackage"></a>

## Annotation Packages

Annotation packages contain lightly or non-curated data from a public
source and are updated with each _Bioconductor_ release (every 6
months).  They are a source of general annotation for one or many
organisms and are not specific to a particular experiment.  When
possible, they should support the `select()` interface from
AnnotationDbi.

Annotation packages should *NOT* be posted to the tracker repository.
Instead send an email to <packages@bioconductor.org> with a
description of the proposed annotation package and futher instructions
of where to send the package will be provided. Whenever possible Annotation
Packages should use the [AnnotationHub][] for managing files.

[AnnotationHub]: /packages/devel/bioc/vignettes/HubPub/inst/doc/CreateAHubPackage.html

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
  system (BBS). The system will check out your package from GitHub. It will
  then run `R CMD build` to create a 'tarball' of your source code,
  vignettes, and man pages. It will run `R CMD check` on the tarball,
  to ensure that the package conforms to standard _R_ programming best
  practices. _Bioconductor_ has chosen to utilize a custom `R CMD check`
  environment; See [R CMD check environment][] for more details. Finally, the
  build system will run `BiocCheckGitClone()` and `BiocCheck()` to ensure that
  the package conforms to _Bioconductor_ [BiocCheck][4] standards. The system
  will perform these steps using the ['devel'
  version](/developers/how-to/useDevel/) of _Bioconductor_, on three platforms
  (Linux, Mac OS X, and Windows).  After these steps are complete, a link to a
  build report will be appended to the new package issue. Avoid surprises by
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
  [_Bioconductor_ Git Server][11].

[4]: /packages/BiocCheck/
[9]: https://stat.ethz.ch/mailman/listinfo/bioc-devel
[11]: /developers/how-to/git/
[R CMD check environment]: https://bioconductor.org/developers/package-guidelines/#checkingenv

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="afteraccept"></a>

## Following Acceptance

Following acceptance of a package:

* Packages accepted on the tracker repository are added to the 'devel'
  branch of the _Bioconductor_ GIT repository, with the current version
  number of the accepted package.
* Packages are then built by the _Bioconductor_ nightly build
  process. If the build is successful, the package has its own
  'landing page' created, and the package is made available to users
  of the 'devel' branch of _Bioconductor_ via `BiocManager::install()`.
* Changes to their package (if any), should be done to version in the 
  [_Bioconductor_ git server][11].
* Developers should bump the `z` portion of their version number every
  time they commit changes to their package, following the
  [Version numbering][versioning] guidelines. If developers don't 
  bump the version, the changes made to their package *do not propagate* 
  to the _Bioconductor_ web site and package repository.
  
  [versioning]:/developers/how-to/version-numbering/

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="support"></a>

## Additional Support

We are eager to enhance the quality and interoperability of
_Bioconductor_ software and will provide additional support when
requested by package developers. Example areas of assistance include
use of appropriate S4 structures, specific guidance on efficient
implementation, guidance on code structure, and critical assessment of
package documentation and structure.  Use the
[bioc-devel](/help/mailing-list/) mailing list or email
<packages@bioconductor.org> to obtain additional support.

* Support Email: <packages@bioconductor.org>

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>
