source("http://bioconductor.org/getBioC.R")

getMonograph <- function(pkgs, groupName="monograph", ...)
{
    if (missing(pkgs))
        biocinstall(groupName=groupName, ...)
    else
        biocinstall(pkgs=pkgs, groupName=groupName, ...)
}

