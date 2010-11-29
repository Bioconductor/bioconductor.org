# ![](/images/icons/magnifier.gif)Bioconductor 2.7 Release Schedule

## Release objectives

The next release of R is not yet scheduled, but is likely in early
April, 2011. An emphasis in this release cycle is to bringing existing
packages into conformance with Bioconductor [guidelines][guidelines]
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

## Release template

The following template will be filled in as the release approaches;
dates are relative to the announced R release date:
 
### -28 days

R:

* R grand feature freeze (alpha).

Bioconductor:

* Packages for deprecation identified.

* All Bioconductor package developers should start using R-devel alpha
  for testing their packages.

### -21 days

Bioconductor:

* Deadline for new package submissions.

* Release candidate of annotation packages built and posted to devel
  annotation data repository.

* All packages conform to Bioconductor [guidelines][guidelines].

### -14 days

R:

* R FEATURE FREEZE (beta).

Bioconductor:

* Grand Feature Freeze (alpha):  no new packages added to BioC
  release roster.

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without error.

### -10 days

Bioconductor:

* Bioconductor Feature Freeze (beta): no API changes to BioC.

### -7 days

R:

* R CODE FREEZE (release candidate).

Bioconductor:

* Deadline for packages passing R CMD build and check without warning.
  Some warnings will be accepted, clarification on the bioc-devel mailing
  list.

### -3 days

R:

* R PRERELEASE.

Bioconductor:

* Bioconductor release candidate 1.  Package maintainers should limit
  changes to "show-stopper" bugs and documentation improvements.

* Annotation data packages finalized.

### R Release

R:

* R 2.12.0 RELEASE.

Bioconductor:

* Creation of the BioC 2.7 release branch. Development can resume on
  trunk, but changes will not be part of the release.

* Build final release repositories.

* Test install scripts, GUI installation.


### +1 day

* Bioconductor release.

## External Resources

For more information on the release details for R, visit the [R
Developer Page](http://developer.r-project.org).
