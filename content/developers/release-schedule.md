# Bioconductor 3.19 Release Schedule

The release date for Bioc 3.19 is schedule for Wednesday May 1st. The 3.19
release will use R-4.4. R-4.4 is schedule to be released April 24th.

The following highlights important deadlines for the release:

### Friday March 22

* Deadline for new package submissions.

### Monday April 8

* Announce deprecated packages for Bioc 3.19. They'll be removed from Bioc 3.20.
  See [End of Life](http://contributions.bioconductor.org/package-end-of-life-policy.html) for details.

### Wednesday April 10

*  No major API changes to Bioc 3.19.

### Monday April 15

* Stop building Bioc 3.18, current release. Commits to this branch will be
  disabled. Start configuring builders for 3.20 and start 3.20 builds as soon as
  possible.  

### Wednesday April 17

* Contributed annotation packages posted to devel annotation data repository.

* Bioconductor 3.19 release candidate.  Package maintainers should limit
  changes to "show-stopper" bugs and documentation improvements.

### Wednesday April 24

* Deadline to add new packages to the Bioc 3.19 manifest. Packages submitted to
  Bioconductor new package submission process must have completed the review
  processes and been accepted to be added to the manifest.

### Friday April 26

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without errors or warnings. This includes software, data experiment
  and workflow packages. Some warnings will be accepted, clarification
  on the bioc-devel mailing list.

### Monday April 29

* Last day to commit changes to the Bioc 3.19 branch. NEWS files
  must be updated before the builds start at 1:30 pm EST or they will
  not be included in the release announcement.

  The branch will be frozen prior to creating the release branch on Tuesday,
  April 23.  Committing last minute changes could break your package in both
  release and devel! Be sure to run 'R CMD build' and 'R CMD check' locally
  before committing any changes.

### Tuesday April 30

* Creation of the Bioc 3.19 release branch. Development can resume on
  devel branch, but changes there will not be part of the release.

### Wednesday May 1

* Bioconductor Release 3.19.


## Post-release

* Build Dockers for release and devel

* Packages marked as deprecated in Bioc 3.19 are now removed from the
  Bioc 3.20 nightly builds.

* Identify packages to be deprecated in the new devel, Bioc 3.20.
  This includes packages with errors and unresponsive maintainers.

* Move 3.18 products to archive. Mirrors should adjust accordingly to not have
  3.18 mirror deleted.

* Update New Submission to use 3.20 devel.
