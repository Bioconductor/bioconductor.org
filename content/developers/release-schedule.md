# Bioconductor 3.9 Release Schedule 

This release will use R-3.6.0 ("Planting of a Tree").
The official release date is schedule for Tuesday April 30th. 
The following highlights important deadlines for the release.

### Friday March 22

* Announce new package submission deadline of April 5.

### Monday April 1

* Announce deprecated packages for Bioc 3.9. They'll be removed from Bioc 3.10.
  See [End of Life](/developers/package-end-of-life) for details.

### Friday April 5

* Deadline for new package submissions.

### Monday April 8

* .db0 packages available for developers.

### Monday April 15

* Stop building Bioc 3.8, current release. Commits to this branch will be
  disabled.

* Start building Bioc 3.10, new devel.

### Friday April 19

* TxDb, OrgDb packages built and posted to devel annotation data repository.

### Wednesday April 24

* No API changes to Bioc 3.9.

* Deadline to add new packages to the BiocC 3.9 manifest.

* Contributed annotation packages posted to devel annotation data repository.

### Friday April 26

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without errors or warnings. This includes software, data experiment
  and workflow packages. Some warnings will be accepted, clarification
  on the bioc-devel mailing list.

* Bioconductor release candidate.  Package maintainers should limit
  changes to "show-stopper" bugs and documentation improvements.

### Wednesday May 1

* Last day to commit changes to the Bioc 3.9 branch. NEWS files
  must be updated before the builds start at 4:45pm EST or they will
  not be included in the release announcement.

  The branch will be frozen prior to creating the release branch on Monday,
  April 29.  Committing last minute changes could break your package in both
  release and devel! Be sure to run 'R CMD build' and 'R CMD check' locally
  before committing any changes.

### Thursday May 2

* Creation of the Bioc 3.9 release branch. Development can resume on
  trunk, but changes there will not be part of the release.

### Friday May 3

* Bioconductor Release 3.9.


## Post-release

* Build AMIs for release and devel

* Build Dockers for release and devel

* Packages marked as deprecated in Bioc 3.9 are now removed from the
  Bioc 3.9/3.10 nightly builds.

* Identify packages to be deprecated in the new devel, Bioc 3.10.
  This includes packages with errors and unresponsive maintainers.
