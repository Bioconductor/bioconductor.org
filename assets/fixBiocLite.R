biocinstall <- function(pkgs, groupName="default", repos, ...)
{
    message <- paste("You have an outdated biocLite() function.",
    "Run 'rm(biocLite)' and try again.")
    stop(message)
}
