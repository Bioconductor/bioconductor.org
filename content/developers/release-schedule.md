# ![](/images/icons/magnifier.gif)Bioconductor 2.13 Release Schedule

## Tentative schedule

### September 17, 2013

* Packages for deprecation identified.

### September 24

* All packages conform to Bioconductor [guidelines][guidelines].

* Deadline for .db0 packages to be available for developers.

[guidelines]: /developers/package-guidelines


### September 26

* Deadline for new package submissions.


### September 30

* Release candidate of annotation packages built and posted to devel
  annotation data repository.


### October 2

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without error.
   
### October 7

* Grand Feature Freeze (alpha):  no new packages added to BioC
  release roster.

### October 9

* Bioconductor Feature Freeze (beta): no API changes to BioC.

* Deadline for packages passing R CMD build and check without warning.
   Some warnings will be accepted, clarification on the bioc-devel mailing
   list.

* Stop release (BioC 2.12) builds. Commits to this branch will be disabled.

* Deadline for annotation package contributors to upload updated packages.


### October 11

* Bioconductor release candidate.  Package maintainers should limit
   changes to "show-stopper" bugs and documentation improvements.

* Annotation data packages finalized.

* Package NEWS files updated. Latest NEWS will be collated and included
  in release announcement.

### October 14

* Creation of the BioC 2.13 release branch. Development can resume on
   trunk, but changes will not be part of the release.

* Build final release repositories.

* Test install scripts, GUI installation.

### October 15

* Bioconductor release.
