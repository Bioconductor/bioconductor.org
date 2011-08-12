biocinstall <- function(pkgs, lib, develOK=FALSE, groupName="default", ...) {
    ## Install Bioconductor packages using CRAN-style repositories
    ##
    ## pkgs: character vector of Bioconductor packages to install.
    ##       The groupName argument will be ignored if pkgs is specified.
    ##
    ## lib: character vector giving the library directories where to
    ##      install the packages.  Recycled as needed.
    ##
    ## groupName: character matching one of "affy", "default",
    ##            "graph", "all"
    ##
    ## ...: extra arguments passed to install.packages.  Of particular interest
    ##      may be the dependencies argument.  This script defaults to TRUE for
    ##      R 2.0.x and to "Depends" for R 2.1.x.  For R 2.1.x, you can set
    ##      dependencies=c("Depends", "Suggests", "Imports").
    ##
    ## develOK: Deprecated.  To receive developmental software and packages,
    ##          running the development version of R should be practiced.


    ## !!! Always change version number when updating this file !!!
    VERSION <- "1.0"

    if (develOK) {
        stop("develOK has been deprecated.  If you would like to use this script for installing developmental packages, you must install the devel version of R.")
    }

    
    stopifnot(require("utils"))
    ## Verify we're running a recent enough version of R
    requiredRVer <- "2.1.0"
    develVer <- "2.2.0"
    thisRVer <- paste(R.Version()[c("major", "minor")], collapse=".")
    if (compareVersion(thisRVer, requiredRVer) < 0) {
        stop(paste("\nYou are currently running R version ", thisRVer,
                   ", however R version ", requiredRVer,
                   " is required.",sep=""))
    }
    if (compareVersion(thisRVer, develVer) >= 0)
      RDEVEL <- TRUE
    else
      RDEVEL <- FALSE
    
    cat(paste("\nRunning biocinstall version", VERSION," with R version ",
              thisRVer, "\n"))

    biocRepRoot = "http://www.bioconductor.org/packages/bioc"
    if (!missing(pkgs)) {
        if (is.null(pkgs) || length(pkgs) == 0)
          stop("pkgs was set to ", sQuote(pkgs),
               ".  I was expecting a non-empty character vector.") 
        pkgsToInstall <- pkgs
    } else { ## pkgs missing
        affyPkgs <- c("affy", "affycomp", "affydata", "affyPLM", "annaffy",
                      "gcrma", "makecdfenv", "matchprobes", "marray")
        
        graphPkgs <- c("graph", "Rgraphviz", "RBGL")

        defaultPkgs <- c("affy", "affydata", "affyPLM", "annaffy", "annotate",
                         "Biobase", "Biostrings", "DynDoc", "edd", "gcrma",
                         "genefilter", "geneplotter", "hgu95av2", "limma",
                         "makecdfenv", "marray", "matchprobes", "multtest",
                         "pamr", "reposTools", "ROC", "siggenes", "sma",
                         "statmod", "tkWidgets", "vsn", "widgetTools", "xtable")

        if (RDEVEL)
          biocRep = paste(biocRepRoot, "devel", sep="/")
        else
          biocRep = paste(biocRepRoot, "stable", sep="/") # XXX: FIXME!!!!!!!!!!!!!!!!!!!!!
        allPkgs = CRAN.packages(CRAN=biocRep)[, "Package"]
        pkgsToInstall <- switch(groupName,
                                "affy"=affyPkgs,
                                "graph"=graphPkgs,
                                "default"=defaultPkgs,
                                "all"=allPkgs,
                                stop("unknown groupName ", sQuote(groupName)))
    }
    if (missing(lib)) {
      lib <- .libPaths()[1]
    } 

    args <- list(...)
    nms <- names(args)
    if (! "dependencies" %in% nms)
      args[["dependencies"]] <- c("Depends", "Imports")
    ## CRAN-style Repositories where we'll look for packages
    reposList = c(
      "http://www.bioconductor.org/packages/bioc/devel",
      "http://www.bioconductor.org/packages/data/annotation/devel",
      "http://www.bioconductor.org/packages/data/experiment/devel",
      "http://www.bioconductor.org/packages/omegahat/devel",
      "http://www.bioconductor.org/packages/lindsey/devel",
      "http://cran.fhcrc.org")
    args <- c(list(pkgs=pkgsToInstall, repos=reposList, lib=lib), args)
    do.call("install.packages", args)
}
