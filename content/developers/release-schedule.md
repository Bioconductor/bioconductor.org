# Bioconductor 3.11 Release Schedule

This release will use the latest R-4.0.0 (scheduled release Friday April 24th)
The official release date is schedule for Tuesday April 28th.
The following highlights important deadlines for the release.


### Friday April 3

* Deadline for new package submissions.

### Monday April 6

* Announce deprecated packages for Bioc 3.11. They'll be removed from Bioc 3.12.
  See [End of Life](/developers/package-end-of-life) for details.

### Tuesday April 14

* Stop building Bioc 3.10, current release. Commits to this branch will be
  disabled.

* Start building Bioc 3.12, new devel.

### Wednesday April 22

* No API changes to Bioc 3.11.

* Deadline to add new packages to the BiocC 3.11 manifest. Package submitted to
  tracker must have completed the review processes and been accepted to be added
  to the manifest

* Contributed annotation packages posted to devel annotation data repository.

### Friday April 24

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without errors or warnings. This includes software, data experiment
  and workflow packages. Some warnings will be accepted, clarification
  on the bioc-devel mailing list.

* Bioconductor release candidate.  Package maintainers should limit
  changes to "show-stopper" bugs and documentation improvements.

### Sunday April 26

* Last day to commit changes to the Bioc 3.11 branch. NEWS files
  must be updated before the builds start at 4:45pm EST or they will
  not be included in the release announcement.

  The branch will be frozen prior to creating the release branch on Monday,
  April 27.  Committing last minute changes could break your package in both
  release and devel! Be sure to run 'R CMD build' and 'R CMD check' locally
  before committing any changes.

### Monday April 27

* Creation of the Bioc 3.11 release branch. Development can resume on
  trunk, but changes there will not be part of the release.

### Tuesday April 28

* Bioconductor Release 3.11.


## Post-release

* Build AMIs for release and devel

* Build Dockers for release and devel

* Packages marked as deprecated in Bioc 3.11 are now removed from the
  Bioc 3.11/3.12 nightly builds.

* Identify packages to be deprecated in the new devel, Bioc 3.12.
  This includes packages with errors and unresponsive maintainers.
