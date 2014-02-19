# ![](/images/icons/magnifier.gif)Bioconductor 2.14 Release Schedule

This release will use R-3.1.0 (Spring Dance). 

## Tentative schedule

### March 11, 2014

* Packages for deprecation identified.

### March 18

* All packages conform to Bioconductor [guidelines][guidelines].

* Deadline for .db0 packages to be available for developers.

[guidelines]: /developers/package-guidelines


### March 20

* Deadline for new package submissions.


### March 24

* Release candidate of annotation packages built and posted to devel
  annotation data repository.


### March 26

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without error.
   
### April 1

* Grand Feature Freeze (alpha):  no new packages added to BioC
  release roster.

### April 3

* Bioconductor Feature Freeze (beta): no API changes to BioC.

* Deadline for packages passing R CMD build and check without warning.
   Some warnings will be accepted, clarification on the bioc-devel mailing
   list.

* Stop release (BioC 2.13) builds. Commits to this branch will be disabled.

* Deadline for annotation package contributors to upload updated packages.


### April 7

* Bioconductor release candidate.  Package maintainers should limit
   changes to "show-stopper" bugs and documentation improvements.

* Annotation data packages finalized.

* Package NEWS files updated. Latest NEWS will be collated and included
  in release announcement.

### April 11

* Creation of the BioC 2.14 release branch. Development can resume on
   trunk, but changes will not be part of the release.

* Build final release repositories.

* Test install scripts, GUI installation.

### April 14

* Bioconductor release.
