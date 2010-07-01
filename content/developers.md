![](/images/icons/magnifier.gif)Developers
==========================================

Bioconductor is an open development project, meaning that all developers from
the scientific community are able to contribute software. In order to foster
coherence within the project, developers are encouraged to follow a common
set of [package guidelines](/developers/guidelines/) to make it easier for
other to use and extend the software. New software is added to the project
through a [package submission process](/developers/submission/).

Additionally, there is a
[Developer's Wiki](http://wiki.fhcrc.org/bioc/DeveloperPage/) that provides
information on the current direction of the core team and a
[developer mailing list](https://stat.ethz.ch/mailman/listinfo/bioc-devel)
that facilitates communication amongst the Bioconductor developers.  

* [Developer Wiki](http://wiki.fhcrc.org/bioc/DeveloperPage/)
* [Developer Mailing List](https://stat.ethz.ch/mailman/listinfo/bioc-devel)

Once a package has been submitted, the rapid development of R and
Bioconductor means that there is some danger it could become
disfunctional over time.  Therefore you will want to check on it
periodically to ensure that it is still working as expected.  To
facilitate this, Bioconductor maintains a [build
system](http://bioconductor.org/checkResults/) to check that every
package in the release and development branches can be built and
checked every night.  If you have incorporated unit tests in the right
way, these will also be run every night during check.  An output of
these efforts is available here:

* [Build/Check Results](http://bioconductor.org/checkResults/)
