## Install Monograph packages using CRAN-style repositories
## ...: arguments passed to install.packages.

getBCBS <- function(...)
{
    # !!! Always change version number when updating this file !!!
    VERSION <- "1.0.0"
    thisRVer <- paste(R.Version()[c("major", "minor")], collapse=".")
    cat(paste("Running getBCBS version", VERSION, "with R version", thisRVer, "...\n"))

    choppedRVer <- gsub("(\\w+).(\\w+).(\\w+)", "\\1.\\2", thisRVer)
    biocVer <- switch(choppedRVer,
                      "2.2" = "1.7",
                      "2.3" = "1.8",
                      "2.4" = "1.9"
               )
    if (is.null(biocVer))
        stop("The Monograph packages are available for ",
             "R >= 2.2 and R < 2.5 only.")
    cat(paste("Your version of R requires version", biocVer,
              "of the Monograph. Please wait...\n"))

    # Fails gracefully for users trying to install Monograph 1.9 binaries
    if (biocVer == "1.9") {
        args <- list(...)
        if ("type" %in% names(args))
            type <- args$type
        else
            type <- getOption("pkgType")
        if (type == "win.binary") {
            line0 <- "Windows binaries are not available for Monograph 1.9 packages.\n"
            stop(line0)
            #line0 <- "Windows binaries are not available.\n"
            #line1 <- "Windows binaries for Monograph 1.9 will be available soon."
            #line2 <- "If your system is set up to build binaries from source packages"
            #line3 <- "(which is probably the case if you've compiled R yourself)"
            #line4 <- "then you can retry with 'getMonograph(type=\"source\")'."
            #line5 <- "Otherwise, please check back later and sorry for the inconvenience..."
            #stop(paste(line0, line1, line2, line3, line4, line5, sep="\n  "))
        }
        if (type == "mac.binary") {
            line0 <- "Mac OS X binaries are not available for Monograph 1.9 packages.\n"
            stop(line0)
            #line0 <- "Mac OS X binaries are not available.\n"
            #line1 <- "Mac OS X binaries for Monograph 1.9 will be available soon."
            #line2 <- "If your system is set up to build binaries from source packages"
            #line3 <- "(which is probably the case if you've compiled R yourself)"
            #line4 <- "then you can retry with 'getMonograph(type=\"source\")'."
            #line5 <- "Otherwise, please check back later and sorry for the inconvenience..."
            #stop(paste(line0, line1, line2, line3, line4, line5, sep="\n  "))
        }
    }

    # CRAN-style Repositories where we'll look for packages
    if (biocVer == "1.7") {
        repos <- c(
            "monograph/1.2",
            "bioc/1.7",
            "data/annotation/1.7",
            "data/experiment/1.7",
            "omegahat/1.7",
            "lindsey/1.7")
    } else {
        repos <- c(
            "monograph",
            "bioc",
            "data/annotation",
            "data/experiment",
            "omegahat",
            "lindsey")
        repos <- paste(biocVer, repos, sep="/")
    }
    repos <- c(paste("http://bioconductor.org/packages", repos, sep="/"),
                   "http://cran.fhcrc.org")

    pkgManifest <- c("affycomp", "affydata", "affypdnn", "affyPLM",
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

    havePkgs <- installed.packages()[, "Package"]
    instPkgs <- pkgManifest[!pkgManifest %in% havePkgs]
    args <- list(...)
    nms <- names(args)
    # This needs to be improved to allow the user to use an abbreviated
    # name for arg "dependencies" e.g. "dep"
    if (! "dependencies" %in% nms)
        args$dependencies <- TRUE
    args <- c(list(pkgs=instPkgs, repos=repos), args)
    if (length(instPkgs)) {
        cat("Installing packages:\n",
            paste(instPkgs, collapse=", "), "\n")
        do.call("install.packages", args)
        cat("DONE\n\n")
    }
    args <- list(repos=repos)
    cat("Running update.packages...\n")
    do.call("update.packages", args)
}
