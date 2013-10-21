Bioconductor is an open development project, meaning that all
developers from the scientific community are able to contribute
software. In order to foster coherence within the project, developers
are encouraged to follow [package
guidelines](/developers/package-guidelines/) to make it easier for
others to use and extend the software. New software is added to the
project through a [package submission
process](/developers/package-submission/).

<h2 id="mailing_list">Mailing List</h2>

The [developer mailing
list](https://stat.ethz.ch/mailman/listinfo/bioc-devel) facilitates
communication amongst Bioconductor developers.

<h2 id="svn">Source Control System</h2>

Packages contributed to Bioconductor are added to the project
[subversion
repository](https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks).
Contributing developers are then provided with password protected
write access to their portion of the codebase. Anonymous access to
this repository is also available; the [source
control](/developers/how-to/source-control) page provides more information.

<h2 id="maintenance">Package Maintenance</h2>

Changes in R and Bioconductor could result in the malfunction of
software packages. Therefore, package maintainers should periodically
check their packages to ensure that they are still working as
expected. To facilitate this, Bioconductor maintains a [daily build
system](/checkResults/) to check that every
package in the release and development branches can be built via
`R CMD build` and checked via `R CMD check`. Packages that fail to
pass these checks will regretfully be dropped from the next release
of Bioconductor.


