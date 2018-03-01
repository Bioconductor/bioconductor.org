# **DRAFT**  
# Bioconductor 3.7 Release Schedule

This release will use R-3.5.0 ("Unsuffered Consequences").

### Monday March 26 

* Announce new package deadline of April 9.

### Monday April 2 

* Announce deprecated / defunct packages for BioC 3.7.

### Monday April 9

* Deadline for new package submissions.

* .db0 packages available for developers.

* Update package NEWS files. Latest NEWS will be collated 
  and included in release announcement.

### Monday April 16

* TxDb, OrgDb packages built and posted to devel annotation data repository.

* Contributed annotation packages posted to devel annotation data repository.

* Bioconductor Feature Freeze: 
  - No API changes to BioC 3.7.
  - No new packages added to BioC 3.7 roster.

* Stop building BioC 3.6, current release.  Commits to this branch will be
  disabled.

* Start building BioC 3.8, new devel.

### Wednesday April 25 

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without errors or warnings. This includes software, data experiment
  and workflow packages. Some warnings will be accepted, clarification 
  on the bioc-devel mailing list.

* Identify packages to be 'deprecated' in the new devel, Bioconductor 3.8.
  This includes packages in 3.7 with no Landing Page (have not built this devel
  cycle) or those with errors and unresponsive maintainers. 

* Packages 'deprecated' in Bioconductor 3.6 are marked as 'defunct' and 
  removed from the nightly builds for Bioconductor 3.7.
  See [End of Life](/developers/package-end-of-life) for details.

* Deadline for NEWS files to be updated.

* Bioconductor release candidate.  Package maintainers should limit
  changes to "show-stopper" bugs and documentation improvements.

### Sunday April 29

* Last day to commit changes to the Bioc 3.7 devel branch. The branch will be
  frozen prior to creating the release branch on Monday, April 30.  Committing
  last minute changes could break your package in both release and devel! Be
  sure to run 'R CMD build' and 'R CMD check' locally before committing any
  changes.

### Monday April 30 

* Creation of the BioC 3.7 release branch. Development can resume on
  trunk, but changes there will not be part of the release.

### Tuesday May 1 

* Bioconductor release.


## Post-release

### Wednesday May 2 

* Build AMIs for release and devel

* Build Dockers for release and devel
