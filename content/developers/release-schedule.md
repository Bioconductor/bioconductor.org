# Bioconductor 3.14 Release Schedule

This release will use  R-4.2.0
The official release date is tentatively schedule for Wednesday October 27th.
The following highlights important deadlines for the release.


### Friday October 1

* Deadline for new package submissions.

### Monday October 6

* Announce deprecated packages for Bioc 3.14. They'll be removed from Bioc 3.15.
  See [End of Life](/developers/package-end-of-life) for details.

### Monday October 11

* Stop building Bioc 3.13, current release. Commits to this branch will be
  disabled.

* Start building Bioc 3.15, new devel (using R-devel).

### Wednesday October 13

* No API changes to Bioc 3.14.

* Contributed annotation packages posted to devel annotation data repository.

### Wednesday October 20

* Deadline to add new packages to the Bioc 3.14 manifest. Package submitted to
  tracker must have completed the review processes and been accepted to be added
  to the manifest

### Friday October 22

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without errors or warnings. This includes software, data experiment
  and workflow packages. Some warnings will be accepted, clarification
  on the bioc-devel mailing list.

* Bioconductor release candidate.  Package maintainers should limit
  changes to "show-stopper" bugs and documentation improvements.

### Monday October 25

* Last day to commit changes to the Bioc 3.14 branch. NEWS files
  must be updated before the builds start at 4:45pm EST or they will
  not be included in the release announcement.

  The branch will be frozen prior to creating the release branch on Tuesday,
  October 26.  Committing last minute changes could break your package in both
  release and devel! Be sure to run 'R CMD build' and 'R CMD check' locally
  before committing any changes.

### Tuesday October 26

* Creation of the Bioc 3.14 release branch. Development can resume on
  trunk, but changes there will not be part of the release.

### Wednesday October 27

* Bioconductor Release 3.14.


## Post-release

* Build Dockers for release and devel

* Packages marked as deprecated in Bioc 3.14 are now removed from the
  Bioc 3.15 nightly builds.

* Identify packages to be deprecated in the new devel, Bioc 3.15.
  This includes packages with errors and unresponsive maintainers.
