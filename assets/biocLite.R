source("http://bioconductor.org/getBioC.R")

biocLite <- function(pkgs, groupName="lite", ...)
{
    if (missing(pkgs))
        biocinstall(groupName=groupName, ...)
    else
        biocinstall(pkgs=pkgs, groupName=groupName, ...)
}

