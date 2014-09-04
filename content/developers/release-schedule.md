# ![](/images/icons/magnifier.gif)Bioconductor 3.0 Release Schedule

This release will use R-3.1.1 (Sock it to Me). 

## Tentative schedule

### September 16, 2014

* Packages for deprecation identified.

### September 23

* All packages conform to Bioconductor [guidelines][guidelines].

* Deadline for first-pass .db0 packages to be available for developers.

[guidelines]: /developers/package-guidelines


### September 25

* Deadline for new package submissions.

* Deadline for release candidate .db0 packages to be available
  for developers.


### September 29

* Release candidate of annotation packages built and posted to devel
  annotation data repository.


### October 2

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without error.

### October 6

* Grand Feature Freeze (alpha):  no new packages added to BioC
  release roster.

### October 8

* Bioconductor Feature Freeze (beta): no API changes to BioC.

* Deadline for packages passing R CMD build and check without warning.
   Some warnings will be accepted, clarification on the bioc-devel mailing
   list.

* Stop release (BioC 2.14) builds. Commits to this branch will be disabled.

* Deadline for annotation package contributors to upload updated packages.


### October 10

* Bioconductor release candidate.  Package maintainers should limit
   changes to "show-stopper" bugs and documentation improvements.

* Annotation data packages finalized.

* Package NEWS files updated. Latest NEWS will be collated and included
  in release announcement.

### October 13

* Creation of the BioC 3.0 release branch. Development can resume on
   trunk, but changes will not be part of the release.

* Build final release repositories.

* Test install scripts, GUI installation.

### October 14

* Bioconductor release.
