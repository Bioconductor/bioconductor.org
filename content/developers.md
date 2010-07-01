![](/images/icons/magnifier.gif)Developers
==========================================

Bioconductor is an open development project, meaning that all developers from
the scientific community are able to contribute software. In order to foster
coherence within the project, developers are encouraged to follow a common
set of [package guidelines](/developers/guidelines/) to make it easier for
other to use and extend the software. New software is added to the project
through a [package submission process](/developers/submission/).


<h2 id="wiki">The Wiki and Developers Mailing List</h2>

Additionally, there is a
[Developer's Wiki](http://wiki.fhcrc.org/bioc/DeveloperPage/) that provides
information on the current direction of the core team and a
[developer mailing list](https://stat.ethz.ch/mailman/listinfo/bioc-devel)
that facilitates communication amongst the Bioconductor developers.  

* [Developer Wiki](http://wiki.fhcrc.org/bioc/DeveloperPage/)
* [Developer Mailing List](https://stat.ethz.ch/mailman/listinfo/bioc-devel)


<h2 id="svn">The svn Repository</h2>

Packages that are contributed to Bioconductor are added to the project
[subversion
repository](https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks).
Contributing developers are then provided with password protected
write access to their portion of the codebase.  Anonymous access to
this repository is also available and only requires the following:
User = readonly, password = readonly.

* [svn respository](https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks)


<h2 id="maintenance">Package Maintenance</h2>

Once a package has been submitted, the rapid development of R and
Bioconductor means that there is some danger it could become
dysfunctional over time.  Therefore you will want to check it
periodically to ensure that it is still working as expected.  To
facilitate this, Bioconductor maintains a [build
system](http://bioconductor.org/checkResults/) to check that every
package in the release and development branches can be built and
checked every night.  If you have incorporated unit tests in the right
way, these will also be run every night during check.  An output of
these efforts is available here:

* [Build/Check Results](http://bioconductor.org/checkResults/)

Packages that are not maintained, and which fail to pass these checks,
will regretfully have to be dropped from release.

