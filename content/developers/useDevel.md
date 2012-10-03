Using The Devel Version of Bioconductor
=======================================

Beginning on October 2, 2012, with the release of Bioconductor 2.11, the
way to use the development (devel) version of Bioconductor (2.12) is to install
R-devel (R-2.16). Packages can then be installed normally; for example, this will install the devel version of IRanges and its dependencies:

    source("http://bioconductor.org/biocLite.R")
    biocLite("IRanges")

R-2.15 users who want to upgrade from Bioconductor 2.10 to 2.11 should
read the [install page](/install).

Users who were already using the devel version of Bioconductor (2.11) prior
to the release of 2.11 simply need to update their packages to the latest
version, as follows:

    source("http://bioconductor.org/biocLite.R")
    biocLite(character())


