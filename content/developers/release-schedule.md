# ![](/images/icons/magnifier.gif)Bioconductor 2.9 Release Schedule


## Release objectives

The next release of Bioconductor will be BioC 2.9 and is scheduled for
Tuesday, November 1st, 2011.

The following schedule outlines major goals: 

* October 10: All packages conform to Bioconductor
  [guidelines][guidelines]...

* October 17: ... and pass R CMD build / R CMD check without
  errors.

* October 24: All packages pass R CMD build / R CMD check without
  warnings.

[guidelines]: /developers/package-guidelines


## Tentative schedule

### October 3

R:

* R grand feature freeze (alpha).

Bioconductor:

* Packages for deprecation identified.

* All Bioconductor package developers should start using R-devel alpha
  for testing their packages.

### October 10

Bioconductor:

* Deadline for new package submissions.

* Release candidate of annotation packages built and posted to devel
  annotation data repository.

* All packages conform to Bioconductor [guidelines][guidelines].

### October 17

R:

* R FEATURE FREEZE (beta).

Bioconductor:

* Grand Feature Freeze (alpha):  no new packages added to BioC
  release roster.

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without error.

### October 20

Bioconductor:

* Bioconductor Feature Freeze (beta): no API changes to BioC.

### October 24

R:

* R CODE FREEZE (release candidate).

Bioconductor:

* Deadline for packages passing R CMD build and check without warning.
  Some warnings will be accepted, clarification on the bioc-devel mailing
  list.

### October 26

Bioconductor:

* Bioconductor release candidate 1.  Package maintainers should limit
  changes to "show-stopper" bugs and documentation improvements.

* Annotation data packages finalized.

### October 28

R:

* R PRERELEASE.

### October 31

R:

* R 2.14.0 RELEASE.

Bioconductor:

* Creation of the BioC 2.9 release branch. Development can resume on
  trunk, but changes will not be part of the release.

* Build final release repositories.

* Test install scripts, GUI installation.

### November 1st

* Bioconductor release.


## External Resources

For more information on the release details for R, visit the [R
Developer Page](http://developer.r-project.org).
