Using the 'Devel' Version of _Bioconductor_
===========================================

Which version of R?
-------------------

Package authors should develop against the version of _R_ that will be
available to users when the _Bioconductor_ devel branch becomes the
_Bioconductor_ release branch.

_R_ has a '.y' release in x.y.z every year (typically mid-April), but
_Bioconductor_ has a .y release (where current devel becomes release)
every 6 months (mid-April and mid-October).

This means that, from mid-October through mid-April, _Bioconductor_
developers should be developing against R-devel. From mid-April to
mid-October, developers should use R-release (actually, the R snapshot
from the R-x-y-branch) for _Bioconductor_ development.

See the [BiocManager][] vignette for instructions on using multiple versions
of _Bioconductor_.

[BiocManager]: https://CRAN.R-project.org/package=BiocManager

Using 'bioc-devel' during mid-April to mid-October
--------------------------------------------------

In order to use the 'bioc-devel' version of _Bioconductor_ during the
mid-April to mid-October release cycle, use the release version of _R_
and invoke the function `install(version="devel")` (from the
[BiocManager][] package):

    if (!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager")
    BiocManager::install(version = "devel")
    BiocManager::valid()              # checks for out of date packages

After doing this, all packages will be installed from the 'bioc-devel'
repository.

Using 'bioc-devel' during mid-October to mid-April
--------------------------------------------------

In order to use the 'bioc-devel' version of _Bioconductor_ during the
mid-October to mid-April release cycle, you must install the devel
version of _R_.

* [Source](https://stat.ethz.ch/R/daily/)
* [macOS](https://mac.r-project.org/)
* [Windows](https://cran.r-project.org/bin/windows/base/rdevel.html)

Then, make sure that your version of [BiocManager][] is current and
your packages up-to-date.

    if (!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager")
    BiocManager::install(version="devel")
    BiocManager::valid()
