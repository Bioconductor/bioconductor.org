# ![](/images/icons/magnifier.gif)Bioconductor 3.4 Release Schedule

This release will use R-3.3.1 (Bug in Your Hair) unless 3.3.2 comes
out in time (would need to be at least 1 week before release day
though).


## Tentative schedule


### Tuesday September 20, 2016

* Packages for deprecation identified.


### Monday September 26

* Deadline for new package submissions.

* All packages conform to Bioconductor [guidelines][guidelines].

* .db0 packages available for developers.

[guidelines]: /developers/package-guidelines


### Monday October 3

* Release candidate of annotation packages built and posted to devel
  annotation data repository.


### Friday October 7

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without error.

* Bioconductor Feature Freeze: no API changes to BioC.


### Wednesday October 12

* No new packages added to BioC release roster.

* Deadline for annotation package contributors to upload updated
  packages.

* Stop release (BioC 3.3) builds.  Commits to this branch will be
  disabled.


### Thursday October 13

* Deadline for packages passing ''R CMD build'' and ''R CMD check''
  without warning.  Some warnings will be accepted, clarification on
  the bioc-devel mailing list.

* Annotation data packages finalized.

* Start rebuilding/checking workflows.


### Friday October 14

* Bioconductor release candidate.  Package maintainers should limit
  changes to "show-stopper" bugs and documentation improvements.

* Package NEWS files updated.  Latest NEWS will be collated and included
  in release announcement.


### Monday October 17

* Creation of the BioC 3.4 release branch.  Development can resume on
  trunk, but changes there will not be part of the release.

* Build final release repositories.

* Test install scripts, GUI installation.


### Tuesday October 18

* Bioconductor release.

