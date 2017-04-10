# ![](/images/icons/magnifier.gif)Bioconductor 3.5 Release Schedule

This release will use R-3.4.0 ("You Stupid Darkness").

## Tentative schedule

### Friday March 24 

* Packages for deprecation identified and announced.

* Start building/checking workflows in R devel.

### Friday March 31

* Deadline for new package submissions.

* .db0 packages available for developers.

### Friday April 7

* Maintainers update package NEWS files. Latest NEWS will be collated 
  and included in release announcement.

* Annotations Deadline:
  - Core packages built and posted to devel annotation data repository.
  - Contributed packages posted to devel annotation data repository.

* Bioconductor Feature Freeze: 
  - No API changes to BioC 3.5.
  - No new packages added to BioC 3.5 roster.

* Stop current release (BioC 3.4) builds.  Commits to this branch will be
  disabled.

### Friday April 14

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without errors or warnings.  Some warnings will be accepted, clarification 
  on the bioc-devel mailing list.

* Deadline for workflows to build successfully.

* Bioconductor release candidate.  Package maintainers should limit
  changes to "show-stopper" bugs and documentation improvements.

* Packages 'deprecated' in Bioconductor 3.4 are marked as 'defunct' and 
  removed from the nightly builds for Bioconductor 3.5.
  See [End of Life](/developers/package-end-of-life) for details.

### Monday April 24 

* Creation of the BioC 3.5 release branch.  Development can resume on
  trunk, but changes there will not be part of the release.

* Build final release repositories.

* Test install scripts, GUI installation.

### Tuesday April 25 

* Bioconductor release.


## Post-release

### Wednesday April 26 

* Build AMIs for release and devel

* Build Dockers for release and devel

* Update Chef recipes
