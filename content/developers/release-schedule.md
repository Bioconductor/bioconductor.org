# Bioconductor 3.17 Release Schedule

This release will use R-4.3. The release date for Bioc 3.17 is tentatively
schedule for Wednesday April 26th pending the scheduled release of R 4.3 on
Friday April 21st. The following highlights important deadlines for the release.


### Friday March 31 

* Deadline for new package submissions.

### Monday April 3

* Announce deprecated packages for Bioc 3.17. They'll be removed from Bioc 3.18.
  See [End of Life](http://contributions.bioconductor.org/package-end-of-life-policy.html) for details.

### Wednesday April 5

*  No major API changes to Bioc 3.17.

### Monday April 10

* Stop building Bioc 3.16, current release. Commits to this branch will be
  disabled. Start configuring builders for 3.18.

### Wednesday April 12

* Contributed annotation packages posted to devel annotation data repository.

* Bioconductor 3.17 release candidate.  Package maintainers should limit
  changes to "show-stopper" bugs and documentation improvements.

### Monday April 17

* Start building Bioc 3.18.

### Wednesday April 19

* Deadline to add new packages to the Bioc 3.17 manifest. Package submitted to
  tracker must have completed the review processes and been accepted to be added
  to the manifest

### Friday April 21

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without errors or warnings. This includes software, data experiment
  and workflow packages. Some warnings will be accepted, clarification
  on the bioc-devel mailing list.

### Monday April 24

* Last day to commit changes to the Bioc 3.17 branch. NEWS files
  must be updated before the builds start at 1:30 pm EST or they will
  not be included in the release announcement.

  The branch will be frozen prior to creating the release branch on Tuesday,
  April 25.  Committing last minute changes could break your package in both
  release and devel! Be sure to run 'R CMD build' and 'R CMD check' locally
  before committing any changes.

### Tuesday April 25

* Creation of the Bioc 3.17 release branch. Development can resume on
  trunk, but changes there will not be part of the release.

### Wednesday April 26

* Bioconductor Release 3.17.


## Post-release

* Build Dockers for release and devel

* Packages marked as deprecated in Bioc 3.17 are now removed from the
  Bioc 3.18 nightly builds.

* Identify packages to be deprecated in the new devel, Bioc 3.19.
  This includes packages with errors and unresponsive maintainers.
