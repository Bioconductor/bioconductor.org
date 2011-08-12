message('BioC_mirror = ',
        getOption("BioC_mirror", "http://bioconductor.org"),
        '\nChange using chooseBioCmirror().')

biocinstallRepos <- function()
{
    ## Bioconductor CRAN-style repositories.
    ## The software repo (bioc) _must_ be the first element.
    repos <- c(
        "2.7/bioc",
        "2.7/data/annotation",
        "2.7/data/experiment",
        "2.7/extra"
    )
    mirror <-
      as.character(getOption("BioC_mirror", "http://bioconductor.org"))
    repos <- paste(mirror, "packages", repos, sep="/")
    # repos <- c(repos, "http://brainarray.mbni.med.umich.edu/bioc")
    repos <- c(repos, ifelse(getOption("repos")=="@CRAN@", "http://cran.fhcrc.org", getOption("repos")))
    return(repos)
}

biocinstallPkgGroups <- function(groupName="default")
{
    groupName <- match.arg(groupName,
        c("default", "affy", "graph", "lite", "monograph", "RBioinf",
          "biocases", "all"))

    lite_group <- c(
        "affy", "affydata", "affyPLM", "affyQCReport", "annaffy",
        "annotate", "Biobase", "biomaRt", "Biostrings", "DynDoc",
        "gcrma", "genefilter", "geneplotter", "GenomicRanges",
        "hgu95av2.db", "limma", "marray", "multtest", "vsn", "xtable")
    if (groupName == "lite")
        return(lite_group)

    default_group <- c(
        lite_group,
        "globaltest", "makecdfenv", "pamr", "siggenes",
        "sma", "statmod", "tkWidgets", "widgetTools")
    if (groupName == "default")
        return(default_group)

    affy_group <- c(
        "affy", "affycomp", "affydata", "affyPLM", "affyQCReport",
        "annaffy", "gcrma", "makecdfenv", "marray")
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
        "BiocCaseStudies",
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
        "matchprobes",
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

### biocinstall() version 2.7.7
### Called by biocLite() when R 2.12 is detected.
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
###   ...: extra arguments passed to install.packages (see install.packages documentation).

biocinstall <- function(pkgs, groupName="default", repos, ...)
{
    ## !!! Always change VERSION when updating this file !!!
    VERSION <- "2.7.7"               # this script version
    BIOC_VERSION <- "2.7"            # this version of Bioconductor

    thisRVer <- paste(R.Version()[c("major", "minor")], collapse=".")
    cat("Using R version ", thisRVer, ", ",
        "biocinstall version ", VERSION, ".\n", sep="")

    stopifnot(require("utils"))

    if (missing(pkgs))
        pkgs <- biocinstallPkgGroups(groupName)

    cat("Installing Bioconductor version", BIOC_VERSION, "packages:\n")
    print(pkgs)
    cat("Please wait...\n\n")

    if (!missing(repos))
        stop("'repos' argument to the 'biocinstall' function not allowed\n")
    repos <- biocinstallRepos()

    type <- list(...)[["type"]]
    if (is.null(type))
        type <- getOption("pkgType")
    if (type == "win64.binary") {
        cat("Please note that support for 64-bit Windows binary packages in\n",
            "Bioconductor is experimental and that some packages are not\n",
            "available for this platform at the moment (but are for 32-bit\n",
            "Windows). Thank you for your comprehension.\n\n",
            sep="")
    } else if (type %in% c("mac.binary", "mac.binary.leopard")) {
        mbniPkgs <-
          intersect(pkgs,
                    row.names(available.packages(
                       "http://brainarray.mbni.med.umich.edu/bioc/src/contrib")))
        if (length(mbniPkgs) > 0) {
            mbniList <- paste(mbniPkgs, collapse = ", ")
            warning("MBNI Brain Array packages (", mbniList,
                    ") are not available as Mac binaries. ",
                    "Please use biocLite again with type='source'.")
        }
    }
    
    install.packages(pkgs=pkgs, repos=repos, ...)
}

