# Bioconductor 3.10 Release Schedule

This release will use the latest R-3.6.1 ("Action of the Toes").
The official release date is schedule for Wednsday October 30th.
The following highlights important deadlines for the release.


### Friday October 4

* Deadline for new package submissions.

### Monday October 7

* Announce deprecated packages for Bioc 3.10. They'll be removed from Bioc 3.11.
  See [End of Life](/developers/package-end-of-life) for details.

### Monday October 14

* Stop building Bioc 3.9, current release. Commits to this branch will be
  disabled.

* Start building Bioc 3.11, new devel.

### Wednesday October 23

* No API changes to Bioc 3.10.

* Deadline to add new packages to the BiocC 3.10 manifest. Package submitted to
  tracker must have completed the review processes and been accepted to be added
  to the manifest

* Contributed annotation packages posted to devel annotation data repository.

### Friday October 25

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without errors or warnings. This includes software, data experiment
  and workflow packages. Some warnings will be accepted, clarification
  on the bioc-devel mailing list.

* Bioconductor release candidate.  Package maintainers should limit
  changes to "show-stopper" bugs and documentation improvements.

### Monday October 28

* Last day to commit changes to the Bioc 3.10 branch. NEWS files
  must be updated before the builds start at 4:45pm EST or they will
  not be included in the release announcement.

  The branch will be frozen prior to creating the release branch on Tuesday,
  October 29.  Committing last minute changes could break your package in both
  release and devel! Be sure to run 'R CMD build' and 'R CMD check' locally
  before committing any changes.

### Tuesday October 29

* Creation of the Bioc 3.10 release branch. Development can resume on
  trunk, but changes there will not be part of the release.

### Wednsday October 30

* Bioconductor Release 3.10.


## Post-release

* Build AMIs for release and devel

* Build Dockers for release and devel

* Packages marked as deprecated in Bioc 3.10 are now removed from the
  Bioc 3.10/3.11 nightly builds.

* Identify packages to be deprecated in the new devel, Bioc 3.11.
  This includes packages with errors and unresponsive maintainers.
