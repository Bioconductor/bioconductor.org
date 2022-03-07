# Bioconductor 3.15 Release Schedule

This release will use R-4.2.0. This release schedule is based off the tentative
April 22nd release of R-4.2.0; If the R release is rescheduled, we will also
adjust dates accordingly. 
The release date for Bioc 3.15 is tentatively schedule for Wednesday April 27th.
The following highlights important deadlines for the release.


### Friday April 1

* Deadline for new package submissions.

### Monday April 4

* Announce deprecated packages for Bioc 3.14. They'll be removed from Bioc 3.15.
  See [End of Life](http://contributions.bioconductor.org/package-end-of-life-policy.html) for details.

### Monday April 11

* Stop building Bioc 3.14, current release. Commits to this branch will be
  disabled.

* Start building Bioc 3.16.

### Wednesday April 13

* No API changes to Bioc 3.15.

* Contributed annotation packages posted to devel annotation data repository.

### Wednesday April 20

* Deadline to add new packages to the Bioc 3.15 manifest. Package submitted to
  tracker must have completed the review processes and been accepted to be added
  to the manifest

### Friday April 22

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without errors or warnings. This includes software, data experiment
  and workflow packages. Some warnings will be accepted, clarification
  on the bioc-devel mailing list.

* Bioconductor release candidate.  Package maintainers should limit
  changes to "show-stopper" bugs and documentation improvements.

### Monday April 25

* Last day to commit changes to the Bioc 3.15 branch. NEWS files
  must be updated before the builds start at 2:30 pm EST or they will
  not be included in the release announcement.

  The branch will be frozen prior to creating the release branch on Tuesday,
  April 26.  Committing last minute changes could break your package in both
  release and devel! Be sure to run 'R CMD build' and 'R CMD check' locally
  before committing any changes.

### Tuesday April 26

* Creation of the Bioc 3.15 release branch. Development can resume on
  trunk, but changes there will not be part of the release.

### Wednesday April 27

* Bioconductor Release 3.15.


## Post-release

* Build Dockers for release and devel

* Packages marked as deprecated in Bioc 3.15 are now removed from the
  Bioc 3.16 nightly builds.

* Identify packages to be deprecated in the new devel, Bioc 3.16.
  This includes packages with errors and unresponsive maintainers.
