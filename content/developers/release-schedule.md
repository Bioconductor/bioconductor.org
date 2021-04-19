# Bioconductor 3.13 Release Schedule

This release will use  R-4.1.0.
The official release date is schedule for Wednesday May 19th.
The following highlights important deadlines for the release.

<b>Note:</b> As R has not officially announced the release of R-4.1.0 there is a possibility
that the release will be pushed back. 


### Friday April 23

* Deadline for new package submissions.

### Monday April 26

* Announce deprecated packages for Bioc 3.13. They'll be removed from Bioc 3.14.
  See [End of Life](/developers/package-end-of-life) for details.

### Tuesday May 4

* Stop building Bioc 3.12, current release. Commits to this branch will be
  disabled.

* Start building Bioc 3.14, new devel.

### Wednesday May 12

* No API changes to Bioc 3.13.

* Deadline to add new packages to the Bioc 3.13 manifest. Package submitted to
  tracker must have completed the review processes and been accepted to be added
  to the manifest

* Contributed annotation packages posted to devel annotation data repository.

### Friday May 14

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without errors or warnings. This includes software, data experiment
  and workflow packages. Some warnings will be accepted, clarification
  on the bioc-devel mailing list.

* Bioconductor release candidate.  Package maintainers should limit
  changes to "show-stopper" bugs and documentation improvements.

### Tuesday May 18

* Last day to commit changes to the Bioc 3.13 branch. NEWS files
  must be updated before the builds start at 4:45pm EST or they will
  not be included in the release announcement.

  The branch will be frozen prior to creating the release branch on Tuesday,
  May 18.  Committing last minute changes could break your package in both
  release and devel! Be sure to run 'R CMD build' and 'R CMD check' locally
  before committing any changes.

### Wednesday May 19

* Creation of the Bioc 3.13 release branch. Development can resume on
  trunk, but changes there will not be part of the release.

### Thursday May 20

* Bioconductor Release 3.13.


## Post-release

* Build AMIs for release and devel

* Build Dockers for release and devel

* Packages marked as deprecated in Bioc 3.13 are now removed from the
  Bioc 3.14 nightly builds.

* Identify packages to be deprecated in the new devel, Bioc 3.14.
  This includes packages with errors and unresponsive maintainers.
