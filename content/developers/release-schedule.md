# Bioconductor 3.24 Release Schedule

The release date for Bioc 3.24 is schedule for Wednesday November 25th. The 3.24
release will use R-4.6.

<!--
With the R dependency for 3.23 being R-4.6, these
deadlines are subject to change based on CRAN official R-4.6 release date.
-->

The following highlights important deadlines for the release:
  
### Monday November 2

* Announce deprecated packages for Bioc 3.24. They'll be removed from Bioc 3.25.
  See [End of Life](http://contributions.bioconductor.org/package-end-of-life-policy.html) for details.

### Wednesday November 4

*  No major API changes to Bioc 3.24.

### Monday November 9

* Stop building Bioc 3.23, current release. Commits to this branch will be
  disabled. Start configuring builders for 3.25 and start 3.25 builds as soon as
  possible.  

### Wednesday November 11

* Contributed annotation packages posted to devel annotation data repository.

* Bioconductor 3.24 release candidate.  Package maintainers should limit
  changes to "show-stopper" bugs and documentation improvements.

### Wednesday November 18

* Deadline to add new packages to the Bioc 3.24 manifest. Packages submitted to
  Bioconductor new package submission process must have completed the review
  processes and been accepted to be added to the manifest.

### Friday Novemeber 20

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without errors or warnings. This includes software, data experiment
  and workflow packages. Some warnings will be accepted, clarification
  on the bioc-devel mailing list.

### Monday November 23

* Last day to commit changes to the Bioc 3.24 branch. NEWS files
  must be updated before the builds start at 1:30 pm EST or they will
  not be included in the release announcement.

  The branch will be frozen prior to creating the release branch on Tuesday,
  November 24.  Committing last minute changes could break your package in both
  release and devel! Be sure to run 'R CMD build' and 'R CMD check' locally
  before committing any changes.

### Tuesday November 24

* Creation of the Bioc 3.24 release branch. Development can resume on
  devel branch, but changes there will not be part of the release.

### Wednesday November 25

* Bioconductor Release 3.24.


## Post-release

* Build Dockers for release and devel

* Packages marked as deprecated in Bioc 3.24 are now removed from the
  Bioc 3.25 nightly builds.

* Identify packages to be deprecated in the new devel, Bioc 3.25.
  This includes packages with errors and unresponsive maintainers.

* Move 3.23 products to archive. Mirrors should adjust accordingly to not have
  3.23 mirror deleted.

* Update New Submission to use 3.25 devel.
