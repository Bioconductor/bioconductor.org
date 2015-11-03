![](/images/icons/magnifier.gif)Package Submission
==================================================

* [Introduction](#introduction)
* [Checklist](#checklist)
* [Submission](#submission)
* [Review Process](#review)
* [Additional Support](#support)

<h2 id="introduction">Introduction</h2>

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

[S4 classes and methods]: http://bioconductor.org/developers/package-guidelines/#classes

Bioconductor Authors should

* Be familiar with the ‘devel’ and ‘release’ branch concepts used in
  the project.  New packages and features are added to the ‘devel’
  branch. The current devel branch becomes the next release, with a
  release in April and October. Once your package has been accepted,
  it will initially be in the ‘devel’ branch. Most users are expected
  to use the release branch, so will not immediately have access to
  your package.
* Realize Bioconductor, unlike CRAN, maintains all package source code
  under version control (‘SVN’; ‘git’ is also possible). This means
  that, once your package is accepted, you will make additional
  changes to your package using SVN or git rather than submitting a
  new tarball.
* Be committed to maintaining your package across multiple release
  cycles.

Other avenues for distributing your package include

* [CRAN](http://www.r-project.org/) (for packages only tangentially
  related to areas of Bioconductor emphasis) and repositories such as
* [r-forge](https://r-forge.r-project.org/) or
* [GitHub](https://github.com) (for packages in early stages of
  development).

Many Bioconductor packages import or depend on CRAN
packages. CRAN packages importing or depending on many Bioconductor
packages can be problematic, because of the different approaches to
repository structure and release schedules.

<h2 id="checklist">Checklist</h2>

Packages must satisfy the following checklist:

* Pass `R CMD build`, `R CMD check`, and `R CMD BiocCheck` (see the
  [R CMD check](http://r-pkgs.had.co.nz/check.html) cheatsheet and the
  [BiocCheck](/packages/devel/bioc/html/BiocCheck.html) package) on all
  supported platforms (Windows, Macintosh, Linux) with no errors or warnings,
  using an appropiate version of R. To work out which version that is, see
  [useDevel](how-to/useDevel).
* The result of `R CMD build` must be less than 4MB;
* `R CMD check` must complete within 5 minutes.
* Contain a DESCRIPTION file with valid contact information, an informative
  title and description, correct license specification, appropriate biocViews
  terms, valid version number.
* Set Version: 0.99.0 in the DESCRIPTION.  Subsequent versions created
  during the review process will be numbered 0.99.1, 0.99.2, etc.
  When released, your package's version number will be automatically
  incremented to 1.0.0.
* Contain a NAMESPACE that imports all symbols used in the package, and
  exports just those symbols the package author identifies as appropriate.
  Use of a NAMESPACE implies that appropriate packages are mentioned in the
  Imports: field of the DESCRIPTION file.
* Contain a vignette that illustrates the major uses of the
  package. The vignette must be *evaluated* during package installation; a
  static vignette is not acceptable.
* Contain comprehensive help pages. This includes accurate description
  of function parameter and return values, and meaningful examples.
* Make use of appropriate existing packages (e.g., biomaRt, AnnotationDbi,
  Biostrings) and classes (e.g., ExpressionSet, AnnotatedDataFrame,
  RangedData, Rle, DNAStringSet) to avoid duplication of functionality
  available in other Bioconductor packages.
* Contain no extraneous files (e.g., '.DS_Store', '.project', '.svn', etc.),
  files with invalid names (e.g., differing only in case), or code that
  cannot be distributed under the license specified by the author.
* Packages should have a descriptive name that is not already in use.
  See if it is by running <code>biocLite("myPackageName")</code>. You
  cannot have a package name that is case-insensitively equal to
  an existing package name in CRAN or Bioconductor.
* Follow the [Package Guidelines](/package-guidelines) for details on
  appropriate use.
* Include an `inst/NEWS` file for providing users with information on package
  updates.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>


<h2 id="submission">Submission</h2>
* Submit by opening a new issue at https://tracker.bioconductor.org
* You will need to be registered in order to use the tracker.
* Once you've registered, log on and create a new issue.
* On the `New Issue Editing` page:
  * attach the package tarball (the .tar.gz file from R CMD build)
  * click Submit
  * (all other fields are automatically filled in, or filled in by us later)

*Note - your DESCRIPTION file will be automatically attached to the first
message, as a convenient summary of your package.*

### Experiment Data Packages ###
Experimental data packages contain data specific to a particular
analysis or experiment. They often accompany a software package for use
in the examples and vignettes and in general are not updated regularly.
If you need a general subset of data for workflows or examples first check the
AnnotationHub resource for available files (e.g., BAM, FASTA, BigWig, etc.).

If you have an associated data package for your software package, please do
*NOT* create a separate issue in the tracker for that. Instead, please add the
data package tarball to the same issue as the software package.

### Annotation Packages ###

Annotation packages contain lightly or non-curated data from a public
source and are updated with each Bioconductor release (every 6 months).
They are a source of general annotation for one or many organisms and
are not specific to a particular experiment.  When possible, they
should support the select() interface from AnnotationDbi.

Annotation packages should *NOT* be uploaded to the tracker. Instead send
an email to <packages@bioconductor.org> with a description of the proposed
annotation package and futher instructions of where to send the package will
be provided.

<h2 id="review">Review Process</h2>

After you submit a tarball a comment will be posted on the issue with the
result of `R CMD build`, `R CMD check` and `R CMD BiocCheck` on all four
platforms. Please address all the Warnings from `R CMD check` and all
'Required' and 'Recommended' issues from `R CMD BiocCheck`. Your issue is
assigned a reviewer, who will addresses your concerns and help you through the
review process. Reviewers are assigned daily, if after submission a reviewer is
not assigned within 3 working days please contact <packages@bioconductor.org>.
The entire review process typically takes between 2 and 5 weeks.

A typical review works as follows.

* The package developer submits first version of the package (0.99.0).
* Build system returns check results.
* The package developer fixes any issues found, runs `R CMD build`,
  `R CMD check` and `R CMD BiocCheck` on their local machine, and
  uploads the new version (0.99.1).
* A reviewer is assigned to the package.
* A detailed package review is returned to the developer within a few weeks.
* The package developer updates their package incorporating the
  reviewer comments, runs `R CMD build`, `R CMD check` and
  `R CMD BiocCheck` on their local machine, and uploads the new version
  (0.99.2).
* The process is repeated, with appropriate version bumps, until the
  package is accepted to Bioconductor.

Following acceptance of a package:

* Packages accepted on the tracker are added to the 'devel' branch of
  the Bioconductor SVN repository, with the current version number of
  the accepted package.
* Packages are then built by the Bioconductor nightly build
  process. If the build is successful, the package has its own
  'landing page' created, and the package is made available to users
  of the 'devel' branch of Bioconductor via `biocLite()`.
* Developers may continue to make changes to their package, but now do
  so to the version in the subversion repository.
* Developers should bump the `z` portion of their version number every
  time they commit changes to their package, following the
  [Version numbering](/developers/how-to/version-numbering/) guidelines. If
  developers don't bump the version, the changes made to their package
  *do not propagate* to the Bioconductor web site and package
  repository.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="support">Additional Support</h2>

We are eager to enhance the quality and interoperability of
Bioconductor software and will provide additional support when
requested by package developers. Example areas of assistance include
use of appropriate S4 structures, specific guidance on efficient
implementation, guidance on code structure, and critical assessment of
package documentation and structure.  Use the
[bioc-devel](/help/mailing-list/) mailing list or email
<packages@bioconductor.org> to obtain additional support.

* [Webinar on Package Submission](https://www.youtube.com/watch?v=QfqaK_BHebU)
* Support Email: <packages@bioconductor.org>

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>
