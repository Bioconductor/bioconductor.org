# ![](/images/icons/magnifier.gif)Bioconductor 2.10 Release Schedule


## Release objectives

The next release of Bioconductor will be BioC 2.11 and is scheduled for
Tuesday, October 2nd, 2012.

The following schedule outlines major goals:

* September 14: All packages conform to Bioconductor
   [guidelines][guidelines]...

* September 21: ... and pass R CMD build / R CMD check without
   errors.

* September 28: All packages pass R CMD build / R CMD check without
   warnings.

[guidelines]: /developers/package-guidelines


## Tentative schedule

### September 3

* Packages for deprecation identified.

### September 10

* Deadline for new package submissions.

* Release candidate of annotation packages built and posted to devel
   annotation data repository.

* All packages conform to Bioconductor [guidelines][guidelines].

* Deadline for .db0 packages to be available for developers.

### September 17

* Grand Feature Freeze (alpha):  no new packages added to BioC
   release roster.

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
   without error.
   
### September 21

* Bioconductor Feature Freeze (beta): no API changes to BioC.

### September 24

* Deadline for packages passing R CMD build and check without warning.
   Some warnings will be accepted, clarification on the bioc-devel mailing
   list.

* Stop release (BioC 2.10) builds. Commits to this branch will be disabled.

* Deadline for contributed annotation package authors to get working 
packages

### September 27

* Bioconductor release candidate 1.  Package maintainers should limit
   changes to "show-stopper" bugs and documentation improvements.

* Annotation data packages finalized.

### October 1

* Creation of the BioC 2.11 release branch. Development can resume on
   trunk, but changes will not be part of the release.

* Build final release repositories.

* Test install scripts, GUI installation.

### October 2

* Bioconductor release.
