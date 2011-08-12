biocinstallRepos <- function()
{
    ## CRAN-style Repositories where we'll look for packages.
    ## The software repo (bioc) _must_ be the first element.
    repos <- c(
        "bioc",
        "data/annotation",
        "data/experiment",
        "extra"
    )
    repos <- paste("http://bioconductor.org/packages/2.0", repos, sep="/")
    ## Temporary only (H. 11/10/2006)
    #repos <- c(repos, "http://cran.us.r-project.org")
    repos <- c(repos, "http://cran.fhcrc.org")
    repos
}

biocinstallPkgGroups <- function(groupName="default")
{
    groupName <- match.arg(groupName,
        c("default", "affy", "graph", "lite", "monograph", "all"))

    lite_group <- c(
        "affy", "affydata", "affyPLM", "annaffy", "annotate",
        "Biobase", "Biostrings", "DynDoc", "gcrma",
        "genefilter", "geneplotter", "hgu95av2", "limma",
        "marray", "matchprobes", "multtest",
        "ROC", "vsn", "xtable", "affyQCReport")
    if (groupName == "lite")
        return(lite_group)

    default_group <- c(
        lite_group,
        "edd", "globaltest", "makecdfenv", "pamr", "siggenes",
        "sma", "statmod", "tkWidgets", "widgetTools")
    if (groupName == "default")
        return(default_group)

    affy_group <- c(
        "affy", "affycomp", "affydata", "affyPLM", "annaffy",
        "gcrma", "makecdfenv", "matchprobes", "marray", "affyQCReport")
    if (groupName == "affy")
        return(affy_group)

    graph_group <- c("graph", "Rgraphviz", "RBGL")
    if (groupName == "graph")
        return(graph_group)

    monograph_group <- c(
        "affycomp", "affydata", "affypdnn", "affyPLM",
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
        "Rgraphviz", "rrcov", "simpleaffy",
        "sma", "SpikeInSubset", "SSOAP", "statmod",
        "vsn", "XML", "xtable", "YEAST", "yeastExpData")
    if (groupName == "monograph")
        return(monograph_group)

    bioc_url <- biocinstallRepos()[[1]]
    contriburl <- paste(bioc_url, "src/contrib", sep="/")
    all_group <- available.packages(contriburl)[, "Package"]
    names(all_group) <- NULL
    if (groupName == "all")
        return(all_group)

    ## Should never happen
    stop("unknown groupName ", sQuote(groupName))
}

### biocinstall() version 2.0.8
### Called by biocLite() when R version 2.5 is detected.
###
### Install Bioconductor packages using CRAN-style repositories.
### Arguments:
###
###   pkgs: character vector of Bioconductor packages to install.
###         The groupName argument will be ignored if pkgs is specified.
###
###   groupName: character matching one of "default", "affy", "graph", "lite",
###              "monograph", "all".
###
###   repos: the user should not try to pass this argument. It's in the argument list
###          of biocinstall() only because we want to raise an error when
###          the user tries to pass it.
###
###   dependencies: passed to install.packages (see install.packages documentation).
###
###   type: passed to install.packages (see install.packages documentation).
###
###   ...: extra arguments passed to install.packages (see install.packages documentation).

biocinstall <- function(pkgs, groupName="default", repos,
                        dependencies=c("Depends", "Imports"), type=getOption("pkgType"), ...)
{
    ## !!! Always change version number when updating this file !!!
    VERSION <- "2.0.8"
    ## R version 2.5.z should have been detected by biocLite()
    thisRVer <- paste(R.Version()[c("major", "minor")], collapse=".")
    cat(paste("Running biocinstall version", VERSION, "with R version", thisRVer, "\n"))
    cat("Your version of R requires version 2.0 of Bioconductor.\n")

    #if (type == "mac.binary") {
    #    line0 <- "Sorry, Mac OS X binaries are not available yet\n"
    #    line1 <- "Mac OS X binaries for Bioconductor 2.0 will be available soon."
    #    line2 <- "If your system is set up to build binaries from source packages"
    #    line3 <- "(which is probably the case if you've compiled R yourself)"
    #    line4 <- "then you can retry with 'biocLite(..., type=\"source\")'."
    #    line5 <- "Otherwise, please check back later and sorry for the inconvenience..."
    #    stop(paste(line0, line1, line2, line3, line4, line5, sep="\n  "))
    #}

    if (missing(pkgs)) {
        pkgs <- biocinstallPkgGroups(groupName)
        cat("Will install the following packages:\n")
        print(pkgs)
        cat("Please wait...\n\n")
    }

    if (!missing(repos))
        stop("You can't pass a 'repos' argument to the 'biocinstall' function\n")

    repos <- biocinstallRepos()

    install.packages(pkgs=pkgs, repos=repos, dependencies=dependencies, type=type, ...)
}

