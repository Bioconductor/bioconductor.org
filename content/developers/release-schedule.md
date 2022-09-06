# Bioconductor 3.16 Release Schedule

This release will use R-4.2.1. The release date for Bioc 3.16 is tentatively
schedule for Wednesday November 2nd.
The following highlights important deadlines for the release.


### Friday October 7

* Deadline for new package submissions.

### Monday October 10

* Announce deprecated packages for Bioc 3.16. They'll be removed from Bioc 3.17.
  See [End of Life](http://contributions.bioconductor.org/package-end-of-life-policy.html) for details.

### Wednesday October 12

*  No API changes to Bioc 3.16.

### Monday October 17

* Stop building Bioc 3.15, current release. Commits to this branch will be
  disabled.

* Start building Bioc 3.17.

### Wednesday October 19

* Contributed annotation packages posted to devel annotation data repository.

* Bioconductor 3.16 release candidate.  Package maintainers should limit
  changes to "show-stopper" bugs and documentation improvements.

### Wednesday October 26

* Deadline to add new packages to the Bioc 3.16 manifest. Package submitted to
  tracker must have completed the review processes and been accepted to be added
  to the manifest

### Friday October 28

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without errors or warnings. This includes software, data experiment
  and workflow packages. Some warnings will be accepted, clarification
  on the bioc-devel mailing list.

### Monday October 31

* Last day to commit changes to the Bioc 3.16 branch. NEWS files
  must be updated before the builds start at 2:30 pm EST or they will
  not be included in the release announcement.

  The branch will be frozen prior to creating the release branch on Tuesday,
  November 1.  Committing last minute changes could break your package in both
  release and devel! Be sure to run 'R CMD build' and 'R CMD check' locally
  before committing any changes.

### Tuesday November 1

* Creation of the Bioc 3.16 release branch. Development can resume on
  trunk, but changes there will not be part of the release.

### Wednesday November 2

* Bioconductor Release 3.16.


## Post-release

* Build Dockers for release and devel

* Packages marked as deprecated in Bioc 3.16 are now removed from the
  Bioc 3.17 nightly builds.

* Identify packages to be deprecated in the new devel, Bioc 3.17.
  This includes packages with errors and unresponsive maintainers.
