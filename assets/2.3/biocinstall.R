biocinstall.defaultPkgs <- function()
{
        c("affy", "affydata", "affyPLM", "annaffy", "annotate",
          "Biobase", "Biostrings", "DynDoc", "edd", "gcrma",
          "genefilter", "geneplotter", "globaltest", "hgu95av2", "limma",
          "makecdfenv", "marray", "matchprobes", "multtest",
          "pamr", "ROC", "siggenes", "sma",
          "statmod", "tkWidgets", "vsn", "widgetTools", "xtable")
}

biocinstall.affyPkgs <- function()
{
        c("affy", "affycomp", "affydata", "affyPLM", "annaffy",
          "gcrma", "makecdfenv", "matchprobes", "marray")
}

biocinstall.graphPkgs <- function()
{
        c("graph", "Rgraphviz", "RBGL")
}

biocinstall.litePkgs <- function()
{
        c("affy", "affydata", "affyPLM", "annaffy", "annotate",
          "Biobase", "Biostrings", "DynDoc", "gcrma",
          "genefilter", "geneplotter", "hgu95av2", "limma",
          "marray", "matchprobes", "multtest",
          "ROC", "vsn", "xtable")
}

biocinstall.monographPkgs <- function()
{
        c("affycomp", "affydata", "affypdnn", "affyPLM",
          "ALL", "ALLMLL", "AmpAffyExample", "annaffy",
          "AnnBuilder", "annotate", "arrayMagic",
          "arrayQuality", "beta7", "Biobase", "bioDist",
          "Biostrings", "cMAP", "CoCiteStats", "convert",
          "e1071", "edd", "estrogen", "exactRankTests",
          "facsDorit", "factDesign", "gbm", "gcrma",
          "geneplotter", "golubEsets", "GOstats", "gpls",
          "graph", "hexbin", "hgu133a", "hgu133atagcdf",
          "hgu133acdf", "hgu133bcdf", "hgu2beta7", "hgu95av2",
          "hgu95av2cdf", "hgu95av2probe", "hopach",
          "hsahomology", "hu6800cdf", "hu6800probe",
          "humanLLMappings", "ipred", "KEGG", "KEGGSOAP",
          "kidpack", "limma", "locfit", "LogitBoost",
          "matchprobes", "mclust", "mlbench",
          "MLInterfaces", "multtest", "pamr", "prada",
          "PROcess", "ProData", "randomForest", "rat2302",
          "RbcBook1", "RBGL", "RColorBrewer", "RCurl",
          "reposTools", "Rgraphviz", "rrcov", "simpleaffy",
          "sma", "SpikeInSubset", "SSOAP", "statmod",
          "vsn", "XML", "xtable", "YEAST", "yeastExpData")
}

biocinstall.allPkgs <- function()
{
        contriburl = "http://bioconductor.org/packages/1.8/bioc/src/contrib"
        available.packages(contriburl)[, "Package"]
}

## biocinstall() version 1.8.5
## Called by biocLite() when R version 2.3.z is detected.
##
## Install Bioconductor packages using CRAN-style repositories.
## Arguments:
##
##   pkgs: character vector of Bioconductor packages to install.
##         The groupName argument will be ignored if pkgs is specified.
##
##   groupName: character matching one of "default", "affy", "graph", "lite",
##              "monograph", "all".
##
##   repos: the user should not try to pass this argument. It's in the argument list
##          of biocinstall() just because we wan't to catch the situation when
##          the user tries to pass it.
##
##   dependencies: passed to install.packages (see install.packages documentation).
##
##   type: passed to install.packages (see install.packages documentation).
##
##   ...: extra arguments passed to install.packages (see install.packages documentation).

biocinstall <- function(pkgs, groupName="default", repos,
                        dependencies=c("Depends", "Imports"), type=getOption("pkgType"), ...)
{
    # !!! Always change version number when updating this file !!!
    VERSION <- "1.8.5"

    # R version 2.3.z should have been detected by biocLite()
    thisRVer <- paste(R.Version()[c("major", "minor")], collapse=".")
    cat(paste("Running biocinstall version", VERSION, "with R version", thisRVer, "\n"))
    cat("Your version of R requires version 1.8 of Bioconductor.\n")

    if (missing(pkgs)) {
        pkgs <- switch(groupName,
                       "default"=biocinstall.defaultPkgs(),
                       "affy"=biocinstall.affyPkgs(),
                       "graph"=biocinstall.graphPkgs(),
                       "lite"=biocinstall.litePkgs(),
                       "monograph"=biocinstall.monographPkgs(),
                       "all"=biocinstall.allPkgs(),
                       stop("unknown groupName ", sQuote(groupName)))
        cat("Will install the following packages:\n")
        print(pkgs)
        cat("Please wait...\n\n")
    }

    if (!missing(repos))
        stop("You can't pass a 'repos' argument to the 'biocinstall' function\n")

    ## CRAN-style Repositories where we'll look for packages
    repos <- c(
        "bioc",
        "data/annotation",
        "data/experiment",
        "omegahat",
        "lindsey"
    )
    repos <- paste("http://bioconductor.org/packages/1.8", repos, sep="/")
    repos <- c(repos, "http://cran.fhcrc.org")

    install.packages(pkgs=pkgs, repos=repos, dependencies=dependencies, type=type, ...)
}
