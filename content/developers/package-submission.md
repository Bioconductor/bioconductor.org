![](/images/icons/magnifier.gif)Package Submission
==================================================

* [Introduction](#introduction)
* [Checklist](#checklist)
* [Additional Support](#support)
* [Contact Information](#contact-info)

<h2 id="introduction">Introduction</h2>

Users are welcome to contribute new packages to Bioconductor. Packages
must meet a [checklist](#checklist) of standards of functionality,
documentation, and interoperability. See
[Contact information](#contact-info) to obtain instructions for
submitting new packages.

Authors receive initial feedback in 1 to 3 weeks. Packages will be checked
for adherence to the Bioconductor package guidelines and the checklist below
by a member of Bioconductor team. Package developers should consult the full
Bioconductor [Package Guidelines](/developers/package-guidelines/). Submission
implies commitment to package maintenance across multiple release cycles.

<h2 id="checklist">Checklist</h2>

Packages must satisfy the following checklist:

* Pass `R CMD build` and `R CMD check` on all supported platforms (Windows,
  Macintosh, Linux) with no errors or warnings, using a recent R-devel.
  The result of `R CMD build` must be less than 4MB; `R CMD check` must
  complete within 5 minutes.
* Contain a DESCRIPTION file with valid contact information, an informative
  title and description, correct license specification, appropriate biocViews
  terms, valid version number.
* Set Version to 0.99.0 in the DESCRIPTION.  Subsequent versions created 
  during the review process will be numbered 0.99.1, 0.99.2, etc.  When 
  released, your package's version number will be automatically incremented 
  to 1.0.0.
* Contain a NAMESPACE that imports all symbols used in the package, and
  exports just those symbols the package author identifies as appropriate.
  Use of a NAMESPACE implies that appropriate packages are mentioned in the
  Imports: field of the DESCRIPTION file.
* Contain a Sweave-style vignette that illustrates the major uses of the
  package. The vignette must be evaluated during package installation; a
  static vignette (e.g., PDF) is not acceptable.
* Contain complete help pages. This includes accurate description of function
  parameter and return values, and meaningful examples.
* Make use of appropriate existing packages (e.g., biomaRt, AnnotationDbi,
  Biostrings) and classes (e.g., ExpressionSet, AnnotatedDataFrame,
  RangedData, Rle, DNAStringSet) to avoid duplication of functionality
  available in other Bioconductor packages.
* Contain no extraneous files (e.g., '.DS_Store', '.project', '.svn', etc.),
  files with invalid names (e.g., differing only in case), or code that
  cannot be distributed under the license specified by the author.
* Packages must not already be available on CRAN.
* Packages should have a descriptive name that is not already in use.
  See if it is by running <code>biocLite("myPackageName")</code>. You
  cannot have a package name that is case-insensitively equal to
  an existing package name in CRAN or Bioconductor.

Packages should also conform to the following:

* Use existing S4 classes and generic functions; see the
  [Package Guidelines](/developers/package-guidelines) for details on
  appropriate use.
* Include an inst/NEWS file for providing users with information on package
  updates.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="support">Additional Support</h2>

We are eager to enhance the quality and interoperability of Bioconductor
software and will provide additional support when requested by package
developers. Example areas of assistance include use of appropriate S4
structures, specific guidance on efficient implementation, guidance on code
structure, and critical assessment of package documentation and structure.
Use the [bioc-devel](/help/mailing-list/) mailing list or the [Contact
information](#contact-information) below to obtain additional support.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="contact-info">Contact Information</h2>

To submit a package or obtain additional support, contact Marc
Carlson (email: mcarlson NEAR fhcrc POINT org)


<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>
