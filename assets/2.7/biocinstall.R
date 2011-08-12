biocinstallRepos <- function()
{
    ## BioConductor CRAN-style repositories.
    ## The software repo (bioc) _must_ be the first element.
    repos <- c(
        "bioc",
        "data/annotation",
        "data/experiment",
        "extra"
    )
    repos <- paste("http://bioconductor.org/packages/2.2", repos, sep="/")
	repos <- c(repos, "http://cran.fhcrc.org")
    repos
}

biocinstallPkgGroups <- function(groupName="default")
{
    groupName <- match.arg(groupName,
        c("default", "affy", "graph", "lite", "monograph", "RBioinf", "biocases", "all"))

    lite_group <- c(
        "affy", "affydata", "affyPLM", "annaffy", "annotate",
        "Biobase", "Biostrings", "DynDoc", "gcrma",
        "genefilter", "geneplotter", "hgu95av2.db", "limma",
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
        "AnnBuilder", "annotate", "arrayQuality",
        "beta7", "Biobase", "bioDist", "Biostrings",
        "cMAP", "CoCiteStats", "convert",
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
    if (groupName == "monograph") {
        warning("'arrayMagic' package is no longer available")
        return(monograph_group)
    }

    RBioinf_group <- c(
        lite_group, graph_group,
        "RBioinf", "BiocCaseStudies",
        "XML", "RCurl",
        "biomaRt", "GEOquery",
        "KEGG", "KEGGSOAP",
        "hgu95av2", "hgu95av2probe", "hgu95av2cdf", "human.db0",
        "BSgenome.Hsapiens.UCSC.hg18")
    if (groupName == "RBioinf")
        return(RBioinf_group)

    biocases_group <- c(
        lite_group, graph_group,
        "ALL",
        "apComplex",
        "bioDist",
        "biocGraph",
        "biomaRt",
        "CCl4",
        "CLL",
        "Category",
        "class",
        "convert",
        "GO.db",
        "GOstats",
        "GSEABase",
        "hgu133a.db",
        "hgu95av2cdf",
        "hgu95av2probe",
        "hopach",
        "KEGG.db",
        "kohonen",
        "lattice",
        "latticeExtra",
        "MASS",
        "MLInterfaces",
        "org.Hs.eg.db",
        "ppiStats",
        "randomForest",
        "RColorBrewer",
        "Rintact",
        "sma",
        "weaver",
        "yeastExpData"
    )
    if (groupName == "biocases")
        return(biocases_group)

    bioc_url <- biocinstallRepos()[[1]]
    contriburl <- paste(bioc_url, "src/contrib", sep="/")
    all_group <- available.packages(contriburl)[, "Package"]
    names(all_group) <- NULL
    if (groupName == "all")
        return(all_group)

    ## Should never happen
    stop("unknown groupName ", sQuote(groupName))
}

### biocinstall() version 2.2.11
### Called by biocLite() when R version 2.7 is detected.
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
    VERSION <- "2.2.11"
    ## R version 2.7.0 should have been detected by biocLite()
    thisRVer <- paste(R.Version()[c("major", "minor")], collapse=".")
    cat(paste("Running biocinstall version", VERSION, "with R version", thisRVer, "\n"))

    stopifnot(require("utils"))
    cat("Your version of R requires version 2.2 of BioConductor.\n")

	if (type == "mac.binary" && substring(R.Version()[["os"]], 1, 8) != "darwin8.") {
		cat("WARNING:  The Mac binary versions of BioConductor 2.2 packages are built for Mac OS X Tiger.\n",
			"         Use 'biocLite(<<pkgs>>, type = \"source\")' to ensure packages have compatible binary libraries.\n")
	}

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
