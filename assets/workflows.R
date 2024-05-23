source("http://bioconductor.org/biocLite.R")

workflowInstall <- function(pkg, ...)
{
    vers <- getRversion()
    if (vers >= "3.6"){
        stop(
            "With R version 3.6 or greater, install Bioconductor ",
            "packages using BiocManager; see https://bioconductor.org/install"
            )
    }else{
        warning("Sourcing this script is defunct!\nThis script is schedule to be removed from this location.\nPlease obtain from https://github.com/Bioconductor/LegacyInstall")
    }
       
    repos <- c(biocinstallRepos(), 
               sprintf("%s//bioconductor.org/packages/%s/workflows", 
                       BiocInstaller:::.protocol(), 
                       BiocInstaller::biocVersion()))
    
    install.packages(pkg, repos=repos, ...)
}
