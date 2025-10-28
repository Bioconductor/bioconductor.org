# Bioconductor 3.22 Release Schedule

The release date for Bioc 3.22 is schedule for <s>Wednesday October 29th</s>
Thursday October 30th. The 3.22 release will use R-4.5.

<!--
With the R dependency for 3.22 being R-4.5, these
deadlines are subject to change based on CRAN official R-4.5 release date.
-->

The following highlights important deadlines for the release:

### Friday September 26

* Deadline for new package submissions. Packages will still need to pass the
  formal review process to be included in the release. Package reviewers are
  volunteer and have a limited capacity for review. Packages submitted by this
  date we will try to have at least an initial review of the pacakge but is
  still subject to reviewer availability; packages submitted after this date are
  not guaranteed to be reviewed and will so only as reviewer time permits.
    
### Monday October 6

* Announce deprecated packages for Bioc 3.22. They'll be removed from Bioc 3.23.
  See [End of Life](http://contributions.bioconductor.org/package-end-of-life-policy.html) for details.

### Wednesday October 8

*  No major API changes to Bioc 3.22.

### Monday October 13

* Stop building Bioc 3.21, current release. Commits to this branch will be
  disabled. Start configuring builders for 3.23 and start 3.23 builds as soon as
  possible.  

### Wednesday October 15

* Contributed annotation packages posted to devel annotation data repository.

* Bioconductor 3.22 release candidate.  Package maintainers should limit
  changes to "show-stopper" bugs and documentation improvements.

### Wednesday October 22

* Deadline to add new packages to the Bioc 3.22 manifest. Packages submitted to
  Bioconductor new package submission process must have completed the review
  processes and been accepted to be added to the manifest.

### Friday October 24

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without errors or warnings. This includes software, data experiment
  and workflow packages. Some warnings will be accepted, clarification
  on the bioc-devel mailing list.

### Monday October 27

* Last day to commit changes to the Bioc 3.22 branch. NEWS files
  must be updated before the builds start at 1:30 pm EST or they will
  not be included in the release announcement.

  The branch will be frozen prior to creating the release branch on Tuesday,
  October 28.  Committing last minute changes could break your package in both
  release and devel! Be sure to run 'R CMD build' and 'R CMD check' locally
  before committing any changes.

### <s>Tuesday October 28</s> Wednesday October 29

* Creation of the Bioc 3.22 release branch. Development can resume on
  devel branch, but changes there will not be part of the release.

### <s>Wednesday October 29</s> Thursday October 30

* Bioconductor Release 3.22.


## Post-release

* Build Dockers for release and devel

* Packages marked as deprecated in Bioc 3.22 are now removed from the
  Bioc 3.23 nightly builds.

* Identify packages to be deprecated in the new devel, Bioc 3.23.
  This includes packages with errors and unresponsive maintainers.

* Move 3.21 products to archive. Mirrors should adjust accordingly to not have
  3.21 mirror deleted.

* Update New Submission to use 3.23 devel.
