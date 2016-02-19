# ![](/images/icons/magnifier.gif)Bioconductor 3.3 Release Schedule

This release will use R-3.3.0 (Supposedly Educational).

## Tentative schedule

### March 16, 2016

* Packages for deprecation identified.

### March 17

* (**R**): Grand Feature Freeze (3.3.0 alpha)

### March 23

* All packages conform to Bioconductor [guidelines][guidelines].

* Deadline for first-pass .db0 packages to be available for developers.

[guidelines]: /developers/package-guidelines


### March 25

* Deadline for new package submissions.

* Deadline for release candidate .db0 packages to be available
  for developers.


### March 30

* Release candidate of annotation packages built and posted to devel
  annotation data repository.

### March 31

* (**R**): Feature Freeze (3.3.0 beta)

### April 1

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without error.

### April 6

* Grand Feature Freeze (alpha):  no new packages added to BioC
  release roster.

### April 7

* (**R**): Code Freeze (3.3.0 RC)

### April 8

* Bioconductor Feature Freeze (beta): no API changes to BioC.

* Deadline for packages passing R CMD build and check without warning.
   Some warnings will be accepted, clarification on the bioc-devel mailing
   list.

* Stop release (BioC 3.2) builds. Commits to this branch will be disabled.

* Deadline for annotation package contributors to upload updated packages.


### April 11

* Bioconductor release candidate.  Package maintainers should limit
   changes to "show-stopper" bugs and documentation improvements.

* Annotation data packages finalized.

* Package NEWS files updated. Latest NEWS will be collated and included
  in release announcement.

* (**R**): Prerelease.  

### April 14

* Creation of the BioC 3.3 release branch. Development can resume on
   trunk, but changes will not be part of the release.

* Build final release repositories.

* Test install scripts, GUI installation.

* (**R**): Release

### April 15

* Bioconductor release.
