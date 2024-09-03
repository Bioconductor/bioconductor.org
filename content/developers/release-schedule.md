# Bioconductor 3.20 Release Schedule

The release date for Bioc 3.20 is schedule for Wednesday October 30th. The 3.20
release will use R-4.4.

The following highlights important deadlines for the release:

### Friday September 20

* Deadline for new package submissions. Packages will still need to pass the
  formal review process to be included in the release. Package reviewers are
  volunteer and have a limited capacity for review. Packages submitted by this
  date we will try to have at least an initial review of the pacakge; packages
  submitted after this date are not guaranteed to be reviewed.
  
### Monday October 7

* Announce deprecated packages for Bioc 3.20. They'll be removed from Bioc 3.21.
  See [End of Life](http://contributions.bioconductor.org/package-end-of-life-policy.html) for details.

### Wednesday October 9

*  No major API changes to Bioc 3.20.

### Tuesday October 15

* Stop building Bioc 3.19, current release. Commits to this branch will be
  disabled. Start configuring builders for 3.21 and start 3.21 builds as soon as
  possible.  

### Wednesday October 16

* Contributed annotation packages posted to devel annotation data repository.

* Bioconductor 3.20 release candidate.  Package maintainers should limit
  changes to "show-stopper" bugs and documentation improvements.

### Wednesday October 23

* Deadline to add new packages to the Bioc 3.20 manifest. Packages submitted to
  Bioconductor new package submission process must have completed the review
  processes and been accepted to be added to the manifest.

### Friday October 25

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without errors or warnings. This includes software, data experiment
  and workflow packages. Some warnings will be accepted, clarification
  on the bioc-devel mailing list.

### Monday October 28

* Last day to commit changes to the Bioc 3.20 branch. NEWS files
  must be updated before the builds start at 1:30 pm EST or they will
  not be included in the release announcement.

  The branch will be frozen prior to creating the release branch on Tuesday,
  October 29.  Committing last minute changes could break your package in both
  release and devel! Be sure to run 'R CMD build' and 'R CMD check' locally
  before committing any changes.

### Tuesday October 29

* Creation of the Bioc 3.20 release branch. Development can resume on
  devel branch, but changes there will not be part of the release.

### Wednesday October 30

* Bioconductor Release 3.20.


## Post-release

* Build Dockers for release and devel

* Packages marked as deprecated in Bioc 3.20 are now removed from the
  Bioc 3.21 nightly builds.

* Identify packages to be deprecated in the new devel, Bioc 3.21.
  This includes packages with errors and unresponsive maintainers.

* Move 3.19 products to archive. Mirrors should adjust accordingly to not have
  3.19 mirror deleted.

* Update New Submission to use 3.21 devel.
