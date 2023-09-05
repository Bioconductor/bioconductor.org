# Bioconductor 3.18 Release Schedule

This release will use R-devel. The release date for Bioc 3.18 is tentatively
schedule for Wednesday October 25th. The following highlights important
deadlines for the release.


### Friday September 29

* Deadline for new package submissions.

### Monday October 2

* Announce deprecated packages for Bioc 3.18. They'll be removed from Bioc 3.19.
  See [End of Life](http://contributions.bioconductor.org/package-end-of-life-policy.html) for details.

### Wednesday October 4

*  No major API changes to Bioc 3.18.

### Monday October 9

* Stop building Bioc 3.17, current release. Commits to this branch will be
  disabled. Start configuring builders for 3.19.

### Wednesday October 11

* Contributed annotation packages posted to devel annotation data repository.

* Bioconductor 3.18 release candidate.  Package maintainers should limit
  changes to "show-stopper" bugs and documentation improvements.

### Monday October 16

* Start building Bioc 3.19.

### Wednesday October 18

* Deadline to add new packages to the Bioc 3.18 manifest. Package submitted to
  tracker must have completed the review processes and been accepted to be added
  to the manifest

### Friday October 20

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without errors or warnings. This includes software, data experiment
  and workflow packages. Some warnings will be accepted, clarification
  on the bioc-devel mailing list.

### Monday October 23

* Last day to commit changes to the Bioc 3.18 branch. NEWS files
  must be updated before the builds start at 1:30 pm EST or they will
  not be included in the release announcement.

  The branch will be frozen prior to creating the release branch on Tuesday,
  October 24.  Committing last minute changes could break your package in both
  release and devel! Be sure to run 'R CMD build' and 'R CMD check' locally
  before committing any changes.

### Tuesday October 24

* Creation of the Bioc 3.18 release branch. Development can resume on
  devel branch, but changes there will not be part of the release.

### Wednesday October 25

* Bioconductor Release 3.18.


## Post-release

* Build Dockers for release and devel

* Packages marked as deprecated in Bioc 3.18 are now removed from the
  Bioc 3.19 nightly builds.

* Identify packages to be deprecated in the new devel, Bioc 3.20.
  This includes packages with errors and unresponsive maintainers.

* Move 3.17 products to archive. Mirrors should adjust accordingly to not have
  3.17 mirror deleted.

* Update New Submission to use 3.19 devel.
