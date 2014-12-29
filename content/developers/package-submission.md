![](/images/icons/magnifier.gif)Package Submission
==================================================

* [Introduction](#introduction)
* [Checklist](#checklist)
* [From submission to acceptance](#acceptance)
* [Additional Support](#support)
* [Contact Information](#contact-info)

<h2 id="introduction">Introduction</h2>

Consider contributing your completed R package to Bioconductor if the
package

* Addresses areas of high-throughput genomic analysis where
  Bioconductor already makes significant contributions, e.g.,
  sequencing, expression and other microarrays, flow cytometry, mass
  spectrometry, image analysis; see
  [biocViews](http://www.bioconductor.org/packages/devel/BiocViews.html#___Software).
* Interoperates with other Bioconductor packages, re-using common data
  structures (see
  [S4 classes and methods](/developers/how-to/S4-classes/)) and
  existing infrastructure (e.g., `rtracklayer::import()` for input of
  common genomic files).
* Adopts software best practices that enable reproducible research and
  use, such as full documentation and vignettes (including fully
  evaluated code) as well as commitment to long-term user support
  through the Bioconductor
  [support site](https://support.bioconductor.org).

Other avenues for distributing your package include
[CRAN](http://www.r-project.org/) (for packages only tangentially
related to areas of Bioconductor emphasis) and repositories such as
[r-forge](https://r-forge.r-project.org/) or
[github](https://github.com) (for packages in early stages of
development).  Many Bioconductor packages import or depend on CRAN
packages. CRAN packages importing or depending on many Bioconductor
packages can be problematic, because of the different approaches to
repository structure and release schedules.

Packages submitted to Bioconductor must meet a [checklist](#checklist)
of standards of functionality, documentation, and
interoperability. See [Contact information](#contact-info) to obtain
instructions for submitting new packages.

Authors receive initial feedback in 1 to 3 weeks. Packages will be checked
for adherence to the Bioconductor package guidelines and the checklist below
by a member of Bioconductor team. Package developers should consult the full
Bioconductor [Package Guidelines](/developers/package-guidelines/). Submission
implies commitment to package maintenance across multiple release cycles.

<h2 id="checklist">Checklist</h2>

Packages must satisfy the following checklist:

* Pass `R CMD build`, `R CMD check`, and `R CMD BiocCheck` (see the
  [BiocCheck](/packages/devel/bioc/html/BiocCheck.html) package) on
  all supported platforms (Windows, Macintosh, Linux) with no errors
  or warnings, using a recent R-devel.  The result of `R CMD build`
  must be less than 4MB; `R CMD check` must complete within 5 minutes.
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

<h2 id="acceptance">From submission to acceptance</h2>

* Package developers submit their first version of a package with the
  version number at 0.99.0
* A detailed package review is returned to the developer
* The package developer updates their package, increments the version
  number to 0.99.1, builds and checks the package on the local
  machine, and re-submitts the package to Bioconductor.
* The process is repeated, with version bumps on each iteration, until
  the package is accepted to Bioconductor.

Example  

* A package was submitted as 'DemoPackage_0.99.0.tar.gz'
* During package review it was recommended that biocViews terms in the
  DESCRIPTION file, a NEWS file and a vignette be added.
* The package developer made all the recommended changes, bumped the
  version number to 0.99.1 in the DESCRIPTION file, built and checked
  the package on their local machine.
* The developer uploaded the tar ball 'DemoPackage_0.99.1.tar.gz' to
  the tracker.

Following acceptance of a package:

* Packages accepted on the tracker are added to the devel branch of
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
  **do not propagate** to the Bioconductor web site and package
  repository.

Example  

* The accepted package tarball is added to the devel branch by the
  Bioconductor group members. The package is built on all supported
  platforms during the nightly build. If the build is successful, the
  package has its own 'landing page' created. Users of the devel
  branch can now download the "0.99.1" version of "DemoPackage" using
  `biocLite("DemoPackage")`.
* The package developer wants to add a new function to his
  package. The following steps are recommended.

  * Check out the package from the svn repository (using the link,
    username and password) emailed to the Developer.
  * Add the functionality to the package.
  * Increment the version number from 0.99.1 to 0.99.2.
  * Build and check the package.
  * Commit changes to the svn repository.
  
* The following day, after the build report is generated - the users
  can now access the "0.99.2" version from the Bioconductor
  repository/website (via `biocLite("DemoPackage")`).

If the developer had made an error and not bumped the version from
0.99.1 to 0.99.2, then the package would be built by the nightly build
process, but the new version would not propogate to the public
repository. `biocLite("DemoPackage")` would return the old version
0.99.1 of the package.

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

To submit a package or obtain additional support, contact Sonali Arora
(email: bioconductorseattle NEAR gmail POINT com)


<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>
