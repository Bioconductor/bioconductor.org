# ![](/images/icons/magnifier.gif)Bioconductor 3.1 Release Schedule

This release will use R-3.2.0 (Full of Ingredients). 

## Tentative schedule

### March 19, 2015

* (**R**): Grand Feature Freeze (3.2.0 alpha)

### March 25


* All packages conform to Bioconductor [guidelines][guidelines].

* Deadline for first-pass .db0 packages to be available for developers.

[guidelines]: /developers/package-guidelines


### March 27

* Deadline for new package submissions.

* Deadline for release candidate .db0 packages to be available
  for developers.


### March 31

* Release candidate of annotation packages built and posted to devel
  annotation data repository.


### April 2

* (**R**): Feature Freeze (3.2.0 beta)

### April 3

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without error.


### April 7

* Grand Feature Freeze (alpha):  no new packages added to BioC
  release roster.


### April 9

* (**R**): Code Freeze (3.2.0 RC)

* Bioconductor Feature Freeze (beta): no API changes to BioC.

* Deadline for packages passing R CMD build and check without warning.
   Some warnings will be accepted, clarification on the bioc-devel mailing
   list.

* Stop release (BioC 3.0) builds. Commits to this branch will be disabled.

* Deadline for annotation package contributors to upload updated packages.


### April 13

* (**R**): Prerelease

* Bioconductor release candidate.  Package maintainers should limit
   changes to "show-stopper" bugs and documentation improvements.

* Annotation data packages finalized.

* Package NEWS files updated. Latest NEWS will be collated and included
  in release announcement.

### April 16

* (**R**): Release (3.2.0)

* Creation of the BioC 3.1 release branch. Development can resume on
   trunk, but changes will not be part of the release.


### April 17

* Bioconductor release.
