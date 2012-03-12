# ![](/images/icons/magnifier.gif)Bioconductor 2.10 Release Schedule


## Release objectives

The next release of Bioconductor will be BioC 2.10 and is scheduled for
Monday, April 2nd, 2012.

The following schedule outlines major goals: 

* March 9: All packages conform to Bioconductor
  [guidelines][guidelines]...

* March 16: ... and pass R CMD build / R CMD check without
  errors.

* March 23: All packages pass R CMD build / R CMD check without
  warnings.

[guidelines]: /developers/package-guidelines


## Tentative schedule

### March 2

R:

* R GRAND-FEATURE FREEZE (alpha).

Bioconductor:

* Packages for deprecation identified.

* All Bioconductor package developers should start using R-devel alpha
  for testing their packages.

### March 9

Bioconductor:

* Deadline for new package submissions.

* Release candidate of annotation packages built and posted to devel
  annotation data repository.

* All packages conform to Bioconductor [guidelines][guidelines].

### March 16

R:

* R FEATURE FREEZE (beta).

Bioconductor:

* Grand Feature Freeze (alpha):  no new packages added to BioC
  release roster.

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without error.

### March 20

Bioconductor:

* Bioconductor Feature Freeze (beta): no API changes to BioC.

### March 23

R:

* R CODE FREEZE (release candidate).

Bioconductor:

* Deadline for packages passing R CMD build and check without warning.
  Some warnings will be accepted, clarification on the bioc-devel mailing
  list.
  
* Stop release (BioC 2.9) builds. Commits to this branch will be disabled.

### March 26

Bioconductor:

* Bioconductor release candidate 1.  Package maintainers should limit
  changes to "show-stopper" bugs and documentation improvements.

* Annotation data packages finalized.

### March 27

R:

* R PRERELEASE.

### March 30

R:

* R 2.15.0 RELEASE.

Bioconductor:

* Creation of the BioC 2.10 release branch. Development can resume on
  trunk, but changes will not be part of the release.

* Build final release repositories.

* Test install scripts, GUI installation.

### April 2

* Bioconductor release.


## External Resources

For more information on the release details for R, visit the [R
Developer Page](http://developer.r-project.org).
