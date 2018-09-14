# ** DRAFT ** 
# Bioconductor 3.8 Release Schedule

This release will use R-3.5.1 ("Feather Spray").

### Monday September 24 

* Announce new package submission deadline of October 5.

### Monday October 1 

* Announce deprecated packages for BioC 3.8. They'll be removed from BioC 3.8.
  See [End of Life](/developers/package-end-of-life) for details.

### Friday October 5 

* Deadline for new package submissions.

### Monday October 8 

* .db0 packages available for developers.

* Update package NEWS files. Latest NEWS will be collated 
  and included in release announcement.

### Monday October 15 

* Stop building BioC 3.7, current release. Commits to this branch will be
  disabled.

* Start building BioC 3.8, new devel.

### Wednesday October 24 

* No API changes to BioC 3.8.

* Deadline to add new packages to the BiocC 3.8 manifest.

* TxDb, OrgDb packages built and posted to devel annotation data repository.

* Contributed annotation packages posted to devel annotation data repository.

### Friday October 26 

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without errors or warnings. This includes software, data experiment
  and workflow packages. Some warnings will be accepted, clarification 
  on the bioc-devel mailing list.

* Deadline for NEWS files to be updated.

* Bioconductor release candidate.  Package maintainers should limit
  changes to "show-stopper" bugs and documentation improvements.

### Monday October 29

* Last day to commit changes to the Bioc 3.8 devel branch. The branch will be
  frozen prior to creating the release branch on Tuesday, October 30.
  Committing last minute changes could break your package in both release and
  devel! Be sure to run 'R CMD build' and 'R CMD check' locally before
  committing any changes.

### Tuesday October 30 

* Creation of the BioC 3.8 release branch. Development can resume on
  trunk, but changes there will not be part of the release.

### Wednesday October 31 

* Bioconductor release.


## Post-release

### Thursday November 1 - Wednesday May 2 

* Build AMIs for release and devel

* Build Dockers for release and devel

* Packages marked as deprecated in BioC 3.8 are now removed from the
  BioC 3.8 nightly builds.

* Identify packages to be deprecated in the new devel, BioC 3.8.
  This includes packages with errors and unresponsive maintainers.
