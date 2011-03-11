# ![](/images/icons/magnifier.gif)Bioconductor 2.8 Release Schedule

## Release objectives

The next release of R is not yet officially scheduled, but is likely in
around mid-April, 2011. An emphasis in this release cycle is to bringing
existing packages into conformance with Bioconductor [guidelines][guidelines]
and marking as 'deprecated' packages that are no longer useful or
whose functionality is better implemented elsewhere.

The following schedule outlines major goals: 

* March 4: identify packages for deprecation.

* March 11: All packages conform to Bioconductor
  [guidelines][guidelines]

* March 18: and pass R CMD build / R CMD check without
  errors.

* March 25: All packages pass R CMD build / R CMD check without
  warnings.

[guidelines]: /developers/package-guidelines

## Tentative schedule

### March 16

R:

* R grand feature freeze (alpha).

Bioconductor:

* Packages for deprecation identified.

* All Bioconductor package developers should start using R-devel alpha
  for testing their packages.

### March 23

Bioconductor:

* Deadline for new package submissions.

* Release candidate of annotation packages built and posted to devel
  annotation data repository.

* All packages conform to Bioconductor [guidelines][guidelines].

### March 30

R:

* R FEATURE FREEZE (beta).

Bioconductor:

* Grand Feature Freeze (alpha):  no new packages added to BioC
  release roster.

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without error.

### April 3

Bioconductor:

* Bioconductor Feature Freeze (beta): no API changes to BioC.

### April 6

R:

* R CODE FREEZE (release candidate).

Bioconductor:

* Deadline for packages passing R CMD build and check without warning.
  Some warnings will be accepted, clarification on the bioc-devel mailing
  list.

### April 10

R:

* R PRERELEASE.

Bioconductor:

* Bioconductor release candidate 1.  Package maintainers should limit
  changes to "show-stopper" bugs and documentation improvements.

* Annotation data packages finalized.

### April 13

R:

* R 2.13.0 RELEASE.

Bioconductor:

* Creation of the BioC 2.8 release branch. Development can resume on
  trunk, but changes will not be part of the release.

* Build final release repositories.

* Test install scripts, GUI installation.


### April 14

* Bioconductor release.

## External Resources

For more information on the release details for R, visit the [R
Developer Page](http://developer.r-project.org).
