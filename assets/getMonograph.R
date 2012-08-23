
thisRVer <- paste(R.Version()[c("major", "minor")], collapse=".")
if (thisRVer > "2.14")
{
    getMonograph <- function(...) 
    {
        source("http://bioconductor.org/biocLite.R")
        biocLite(monograph_group())
    }
} else {
    source("http://bioconductor.org/getBioC.R")

    getMonograph <- function(pkgs, groupName="monograph", ...)
    {
        if (missing(pkgs))
            biocinstall(groupName=groupName, ...)
        else
            biocinstall(pkgs=pkgs, groupName=groupName, ...)
    }

}
