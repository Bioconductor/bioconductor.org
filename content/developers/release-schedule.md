# Bioconductor 3.23 Release Schedule

The release date for Bioc 3.23 is schedule for Wednesday April 29th. The 3.23
release will use R-4.6. With the R dependency for 3.23 being R-4.6, these
deadlines are subject to change based on CRAN official R-4.6 release date.

The following highlights important deadlines for the release:

### Friday March 20

* Deadline for new package submissions. Packages will still need to pass the
  formal review process to be included in the release. Package reviewers are
  volunteer and have a limited capacity for review. Packages submitted by this
  date we will try to have at least an initial review of the pacakge but is
  still subject to reviewer availability; packages submitted after this date are
  not guaranteed to be reviewed and will so only as reviewer time permits.
    
### Monday March 30

* Announce deprecated packages for Bioc 3.23. They'll be removed from Bioc 3.24.
  See [End of Life](http://contributions.bioconductor.org/package-end-of-life-policy.html) for details.

### Wednesday April 1

*  No major API changes to Bioc 3.23.

### Monday April 6

* Stop building Bioc 3.22, current release. Commits to this branch will be
  disabled. Start configuring builders for 3.24 and start 3.24 builds as soon as
  possible.  

### Wednesday April 8

* Contributed annotation packages posted to devel annotation data repository.

* Bioconductor 3.23 release candidate.  Package maintainers should limit
  changes to "show-stopper" bugs and documentation improvements.

### Wednesday April 22

* Deadline to add new packages to the Bioc 3.23 manifest. Packages submitted to
  Bioconductor new package submission process must have completed the review
  processes and been accepted to be added to the manifest.

### Friday April 24

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without errors or warnings. This includes software, data experiment
  and workflow packages. Some warnings will be accepted, clarification
  on the bioc-devel mailing list.

### Monday April 27

* Last day to commit changes to the Bioc 3.23 branch. NEWS files
  must be updated before the builds start at 1:30 pm EST or they will
  not be included in the release announcement.

  The branch will be frozen prior to creating the release branch on Tuesday,
  April 28.  Committing last minute changes could break your package in both
  release and devel! Be sure to run 'R CMD build' and 'R CMD check' locally
  before committing any changes.

### Tuesday April 28

* Creation of the Bioc 3.23 release branch. Development can resume on
  devel branch, but changes there will not be part of the release.

### Wednesday April 29

* Bioconductor Release 3.23.


## Post-release

* Build Dockers for release and devel

* Packages marked as deprecated in Bioc 3.23 are now removed from the
  Bioc 3.24 nightly builds.

* Identify packages to be deprecated in the new devel, Bioc 3.24.
  This includes packages with errors and unresponsive maintainers.

* Move 3.22 products to archive. Mirrors should adjust accordingly to not have
  3.22 mirror deleted.

* Update New Submission to use 3.24 devel.
