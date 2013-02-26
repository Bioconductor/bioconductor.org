# ![](/images/icons/magnifier.gif)Bioconductor 2.12 Release Schedule

## Tentative schedule

### March 4, 2013

* Packages for deprecation identified.

### March 11

* All packages conform to Bioconductor [guidelines][guidelines].

* Deadline for .db0 packages to be available for developers.

[guidelines]: /developers/package-guidelines


### March 13

* Deadline for new package submissions.


### March 15

* Release candidate of annotation packages built and posted to devel
  annotation data repository.


### March 18

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without error.
   
### March 22

* Grand Feature Freeze (alpha):  no new packages added to BioC
  release roster.

### March 25

* Bioconductor Feature Freeze (beta): no API changes to BioC.

* Deadline for packages passing R CMD build and check without warning.
   Some warnings will be accepted, clarification on the bioc-devel mailing
   list.

* Stop release (BioC 2.11) builds. Commits to this branch will be disabled.

* Deadline for annotation package contributors to upload updated packages.


### March 28

* Bioconductor release candidate.  Package maintainers should limit
   changes to "show-stopper" bugs and documentation improvements.

* Annotation data packages finalized.

* Package NEWS files updated. Latest NEWS will be collated and included
  in release announcement.

### April 3

* Creation of the BioC 2.12 release branch. Development can resume on
   trunk, but changes will not be part of the release.

* Build final release repositories.

* Test install scripts, GUI installation.

### April 4

* Bioconductor release.
