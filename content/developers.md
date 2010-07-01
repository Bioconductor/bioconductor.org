![](/images/icons/magnifier.gif)Developers
==========================================

Bioconductor is an open development project, meaning that all developers from
the scientific community are able to contribute software. In order to foster
coherence within the project, developers are encouraged to follow a common
set of [package guidelines](/developers/package-guidelines/) to make it easier for
other to use and extend the software. New software is added to the project
through a [package submission process](/developers/package-submission/).

<h2 id="wiki">Developer's wiki and mailing list</h2>

Additionally, there is a
[Developer's Wiki](http://wiki.fhcrc.org/bioc/DeveloperPage/) that provides
information on the current direction of the core team and a
[developer mailing list](https://stat.ethz.ch/mailman/listinfo/bioc-devel)
that facilitates communication amongst the Bioconductor developers.

<h2 id="svn">Source control system</h2>

Packages that are contributed to Bioconductor are added to the project
[subversion
repository](https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks).
Contributing developers are then provided with password protected
write access to their portion of the codebase. Anonymous access to
this repository is also available and only requires the following:
User: `readonly`, password: `readonly`.

<h2 id="maintenance">Package maintenance</h2>

Changes in R and Bioconductor could result in the malfunction of
software packages. Therefore, package maintainers should periodically
check their packages to ensure that they are still working as
expected. To facilitate this, Bioconductor maintains a [daily build
system](http://bioconductor.org/checkResults/) to check that every
package in the release and development branches can be built via
`R CMD build` and checked via `R CMD check`. Packages that fail to
pass these checks will regretfully be dropped from the next release
of Bioconductor.

<h2 id="howto">How to Wiki page</h2>

There is also a [wiki how to](http://wiki.fhcrc.org/bioc/HowTo) page that
explains how to do many common things.
