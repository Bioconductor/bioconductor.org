# Bioconductor 3.6 Release Schedule - DRAFT

This release will use R-3.4.2 ("Short Summer").

### Friday Sept 29 

* Announce deprecated / defunct packages for BioC 3.6.

* Start building/checking workflows in R devel.

### Friday Oct 6 

* Deadline for new package submissions.

* .db0 packages available for developers.

### Friday Oct 13

* Maintainers update package NEWS files. Latest NEWS will be collated 
  and included in release announcement.

* TxDb, OrgDb packages built and posted to devel annotation data repository.

### Tuesday Oct 17

* Contributed annotation packages posted to devel annotation data repository.

* Bioconductor Feature Freeze: 
  - No API changes to BioC 3.6.
  - No new packages added to BioC 3.6 roster.

* Stop current release (BioC 3.5) builds.  Commits to this branch will be
  disabled.

* Start building BioC 3.7 new devel

### Tuesday October 24 

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without errors or warnings.  Some warnings will be accepted, clarification 
  on the bioc-devel mailing list.

* Deadline for workflows to build successfully.

* Identify packages to be 'deprecated' in the new devel, Bioconductor 3.7. This 
  includes packages in 3.6 with no Landing Page (have not built this 
  devel cycle) or those with errors and unresponsive maintainers. 

* Packages 'deprecated' in Bioconductor 3.5 are marked as 'defunct' and 
  removed from the nightly builds for Bioconductor 3.6.
  See [End of Life](/developers/package-end-of-life) for details.

* Deadline for NEWS files to be updated.

* Bioconductor release candidate.  Package maintainers should limit
  changes to "show-stopper" bugs and documentation improvements.

### Monday October 30 

* Creation of the BioC 3.6 release branch.  Development can resume on
  trunk, but changes there will not be part of the release.

* Build final release repositories.


### Tuesday October 31 

* Bioconductor release.


## Post-release

### Wednesday November 1 

* Build AMIs for release and devel

* Build Dockers for release and devel

* Update Chef recipes
