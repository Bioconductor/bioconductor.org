## This function gets and installs the required Bioconductor libraries.
##
## libName: a single character string (no vectors) for the name of the library to be
##          installed. Valid names include (the default is "default"):
##
##              "affy" - packages "affy", "affydata",
##                       "affyPLM", "annaffy", "matchprobes",
##                       "gcrma", "makecdfenv", plus exprs.
##              "cdna" - packages marray, vsn, plus exprs
##              "default" - All packages from "affy", "cdna" and "exprs".
##              "exprs" - packages "Biobase" "annotate", "genefilter",
##                        "geneploter", "edd", "ROC", "limma",
##                        "globaltest", "RMAGEML", "siggenes",
##                        "pamr" "vsn", and "multtest".
##              "prog" - packages "graph", "hexbin", "Ruuid", and
##                         "externalVector".
##              "graph" - packages "graph", "Rgraphviz", and "RBGL"
##              "widgets" - packages "tkWidgets", "widgetTools", and "DynDoc"
##              "database" - packages AnnBuilder, SAGElyzer, Rdbi,
##                           RdbiPgSql
##              "design" - packages daMa and factDesign
##               annotation - packages annotate, AnnBuilder,
##                            humanLLMappings, KEGG, GO, SNPtools,
##                            makecdfenv, ontoTools
##              "analyses" - Biobase, ctc, daMA, edd, factDesign,
##                           genefilter, geneplotter, globaltest,
##                           gpls, limma, RMAGMEL, multtest, pamr,
##                           ROC, siggenes, splicegear, qvalue
##             "externalData" - externalVector
##             "proteomics" - gpls, PROcess, apComplex
##             "arrayCGH" - DNAcopy, aCGH

## develOK: Deprecated.  To receive developmental software and packages,
##          running the development version of R should be practiced.

## destdir: The installation directory.  If no directory is provided,
##          getBioC() will default to using your standard library path.

## verbose: a boolean indicating whether any error related to the
##          downloading process will be (TRUE) printed. Error messages will
##          still be returned but invisible if verbose is set to FALSE.

## versForce: by default, if a binary package is being downloaded (e.g. a
##        Win32 zip package) and was built using a different version
##        of R then the user is running, it will install.  Setting
##        this to FALSE will disable this behaviour and explicitly check.

## force: by default, if all of the dependencies for a package are not
##        available, the package itself will be downloaded.  If
##        force is FALSE, the system will not download a package for
##        which it can not obtain all of the dependencies

## getAllDeps: Will not prompt the user to download depenedencies, but
##        will get them automatically.

## method: Method to be used for downloading files.  Currently download
##         methods `"internal"', `"wget"' and `"lynx"' are available.
##         The default is to choose the first of these which will be
##         `"internal"'.

biocinstall <- function (libName = "default",  develOK=FALSE,
                     destdir, versForce=TRUE,
                     verbose = TRUE,
                     force=TRUE, getAllDeps=TRUE, method="auto") {
    if (develOK) {
        stop("develOK has been deprecated.  If you would like to use this script for installing developmental packages, you must install the devel version of R.")
    }

    ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ## !!! Always change version number when updating this file
    getBioCVersion <- "1.2.70"
    writeLines(paste("Running getBioC version ",getBioCVersion,"....\n",
                     "If you encounter problems, first make sure that\n",
                     "you are running the latest version of getBioC()\n",
                     "which can be found at: www.bioconductor.org/getBioC.R",
                     "\n\n",
                     "Please direct any concerns or questions to",
                     " bioconductor@stat.math.ethz.ch.\n",sep=""))
    ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    MINIMUMRel <- "2.0.0"
    MINIMUMDev <- "2.0.0"

    rInfo <- R.Version()

    rVers <- paste(rInfo$major,rInfo$minor,sep=".")
    ## !!! Again, using compareVersion() here as we
    ## !!! have not bootstrapped w/ reposTools for versionNumber
    ## !!! class
    if (develOK)
        MINIMUM <- MINIMUMDev
    else
        MINIMUM <- MINIMUMRel


    if (compareVersion(rVers,MINIMUM) < 0)
    {
        stop(paste("\nYou are currently running R version ",rVers,
                   ", however R version ",MINIMUM,
                   " is required.",sep=""))
    }

    ## Check the specified libName.  If it is 'all' we want to warn
    ## that they're about ready to get a metric ton of packages.
    if (length(libName) > 1)
        stop("The 'libName' parameter does not accept vectors ",
             "please select a single 'libName' option.")

    if (libName == "all") {
        ## Make sure they want to get all packages
        msg <- paste("\nYou are downloading all of the Bioconductor",
                     " packages and any dependencies.\n",
                     "Depending on your system this will be about ",
                     "160 packages and be roughly 1GB in size.\n",
                     "\nAre you sure that you want to do this?", sep="")
        out <- GBCuserQuery(msg,c("y","n"))
        if (out == "n") {
            cat("\nNot downloading. if you wish to see other download options,\n",
                "please go to the URL:",
                " http://www.bioconductor.org/faq.html#getBioC\n", sep="")
            return(invisible(NULL))
        }
    }
    curLibPaths <- .libPaths()
    on.exit(.libPaths(curLibPaths), add=TRUE)


     ## make sure to expand out the destdir param
    if (!missing(destdir))
        destdir <- path.expand(destdir)

    ## Stifle the "connected to www.... garbage output
    curNetOpt <- getOption("internet.info")
    on.exit(options(internet.info=curNetOpt), add=TRUE)
    options(internet.info=3)

    ## First check to make sure they have HTTP capability.  If they do
    ## not, there is no point to this exercise.
    http <- as.logical(capabilities(what="http/ftp"))
    if (http == FALSE) {
        stop(paste("Your R is not currently configured to allow HTTP",
                   "\nconnections, which is required for getBioC to",
                   "work properly."))
    }

    ## find out where we think that bioC is
    bioCoption <- getOption("BIOC")
    if (is.null(bioCoption))
        bioCoption <- "http://www.bioconductor.org"

    ## Now check to see if we can connect to the BioC website
    biocURL <- url(paste(bioCoption,"/main.html",sep=""))
    options(show.error.messages=FALSE)
    test <- try(readLines(biocURL)[1])
    options(show.error.messages=TRUE)
    if (inherits(test,"try-error"))
        stop(paste("Your R can not connect to the Bioconductor",
                   "website, which is required for getBioC to",
                   "work properly.  The most likely cause of this",
                   "is the internet configuration of R"))
    else
        close(biocURL)

    ## Get the destination directory
    if (missing(destdir)) {
        lP <- .libPaths()
        if (length(lP) == 1)
            destdir <- lP
        else {
            dDval <- menu(lP,
                          title="Please select an installation directory:")
            if (dDval == 0)
                stop("No installation directory selected")
            else
                destdir <- lP[dDval]
        }
    }
    else
        .libPaths(destdir)

    if (length(destdir) > 1)
        stop("Invalid destdir parameter, must be of length 1")

    PLATFORM <- .Platform$OS.type
    if (file.access(destdir,mode=0) < 0)
        stop(paste("Directory",destdir,"does not seem to exist.\n",
                   "Please check your 'destdir' parameter and try again."))

    if (file.access(destdir,mode=2) < 0)
        stop(paste("You do not have write access to",destdir,
                   "\nPlease check your permissions or provide",
                   "a different 'destdir' parameter"))

    if (file.exists(file.path(destdir, "00LOCK")))
        stop("You have a 00LOCK directory in\n", destdir,
             ", presumably from a failed installation prior to running this",
             ".  Please remove this directory before running getBioC.")

    messages <- paste("Your packages are up to date.",
                      "No downloading/installation will be performed.",
                      sep="\n")


    ## Download and install reposTools and Biobase first

    ## Get the package description file from Bioconductor
    getReposTools <- getReposTools(develOK, PLATFORM, destdir, method=method,
                                   bioCoption=bioCoption)
    require(reposTools) || stop("Needs reposTools to continue")

    ## Get Repository entries from Bioconductor
    urlPath <- switch(PLATFORM,
                      "unix"="/Source",
                      "/Win32")

    bioCRepURL <- getReposURL(develOK, urlPath, bioCoption)
    bioCEntries <- getReposEntry(bioCRepURL)

    curOps <- getOption("repositories2")
    on.exit(options(repositories2=curOps), add=TRUE)
    repNames <- names(curOps)
    curOps <- gsub("http://www.bioconductor.org",bioCoption,curOps)
    names(curOps) <- repNames

    ## FIXME:
    ## This can be removed for the 1.5 release
    if (develOK == TRUE)
        optReps <- curOps[c("BIOCDevel","BIOCDevData",
                            "BIOCCourses","BIOCcdf",
                            "BIOCExdata",
                            "BIOCprobes","CRAN",
                            "BIOCOmegahat", "BIOCLinds")]
    else
        optReps <- curOps[c("BIOCRel1.5","BIOCData",
                            "BIOCCourses","BIOCcdf",
                            "BIOCExdata", "BIOCprobes",
                            "CRAN", "BIOCOmegahat", "BIOCLinds")]

    optReps <- optReps[!is.na(optReps)]
    options(repositories2=optReps)

    ## Sync lib list
    syncLocalLibList(destdir)

    ## 'packs' might be NULL, implying everything in the
    ## main repository (release/devel)

    if (libName == "all")
        packs <- repPkgs(bioCEntries)
    else {
        ## FIXME:
        ## Despite using themes, still need to go through
        ## this headache right now because we're doing an
        ## update and then an install to pick up any new
        ## values.  Need to have functionality like the old
        ## install.packages2() which only updates old
        ## packages and installs new packages built into
        ## reposTools.  So we're only pseudo-using themes
        ## right now.

        repThemes <- reposThemes(bioCEntries)
        themeNames <- unlist(lapply(repThemes, repThemeName))
        whichTheme <- match(libName, themeNames)
        if (is.na(whichTheme))
            stop(packNameOutput())
        packs <- unlist(lapply(repThemePkgs(repThemes[[whichTheme]]),
                               pkgName))
    }

    ## Need to determine which 'packs' are alreaedy
    ## installed and which are not.  Call install on the latter
    ## and update on the former.

    if (load.locLib(destdir)) {
        locPkgs <- unlist(lapply(locLibList, Package))
        havePkgs <- packs %in% locPkgs
        installPkgs <- packs[! havePkgs]
        updatePkgs <- packs[havePkgs]
    }
    else {
        installPkgs <- packs
        updatePkgs <- NULL
    }


    out <- new("pkgStatusList", statusList=list())
    if (length(updatePkgs) > 0) {
        updateList <- update.packages2(updatePkgs, bioCEntries,
                                       libs=destdir,
                                       type = ifelse(PLATFORM == "unix", "Source",
                                       "Win32"), versForce=versForce, recurse=FALSE,
                                       getAllDeps=getAllDeps, method=method,
                                       force=force, prevRepos=FALSE,
                                       searchOptions=TRUE,
                                       develOK=develOK)

        statusList(out) <- updateList[[destdir]]
    }

    syncLocalLibList(destdir, quiet=TRUE)
    if (length(installPkgs) > 0)
        statusList(out) <- install.packages2(installPkgs, bioCEntries, lib=destdir,
                                             type = ifelse(PLATFORM == "unix", "Source",
                                             "Win32"), versForce=versForce,
                                             recurse=FALSE,
                                             getAllDeps=getAllDeps,
                                             method=method,
                                             force=force,
                                             searchOptions=TRUE,
                                             develOK=develOK)

    if (length(updated(out)) == 0) {
        print("All requested packages are up to date")
    } else {
        print(out)
    }

    ## Windows doesn't currently have Rgraphviz
    if (PLATFORM != "windows") {
        if (libName %in% c("all","prog","graph")) {
            otherPkgsOut <- paste("Package Rgraphviz requires",
                                  " special libraries to be installed.\n",
                                  "Please see the URL ",
                                  "http://www.bioconductor.org/faq.html#Other Notes",
                                  " for\n",
                                  "more details on installing these packages",
                                  " if they fail\nto install properly\n\n",
                                  sep=""
                                  )
            cat(otherPkgsOut)
        }
    }

    ## If they are using 'default', alert the user that they have not
    ## gotten all packages
    if (libName == "default") {
        out <- paste("You have downloaded a default set of packages.\n",
                     "If you wish to see other download options, please",
                     " go to the URL:\n",
                     "http://www.bioconductor.org/faq.html#getBioC\n",
                     sep="")
        cat(out)
    }
}

getReposTools <- function(develOK, platform, destdir=NULL,
                          method="auto", bioCoption) {
    ## This funciton will check to see if reposTools needs to be
    ## updated, and if so will download/install it

    PACKAGES <- getPACKAGES(develOK, bioCoption)

    ### check reposTools ala checkLibs
    if (checkReposTools(PACKAGES)) {
        sourceUrl <- getDLURL("reposTools", PACKAGES, platform)
        sourceUrl <- paste(bioCoption, sourceUrl, sep="/")
        ## Get the package file name for reposTools
        fileName <- getFileName(sourceUrl, destdir)
        ## Try the connection first before downloading
        options(show.error.messages = FALSE)
        tryMe <- try(url(sourceUrl, "r"))
        options(show.error.messages = TRUE)
        if(inherits(tryMe, "try-error"))
            stop("Could not get the required package reposTools")
        else {
            ## Close the connection for checking
            close(tryMe)
            ## Download and install
            print("Installing reposTools ...")
            download.file(sourceUrl, fileName,
                          mode = getMode(platform), quiet = TRUE, method=method)
            installPack(platform, fileName, destdir)
            if (!("reposTools" %in% installed.packages(lib.loc=destdir)[,"Package"]))
                stop("Failed to install package reposTools")
            unlink(fileName)
        }
        return(invisible(NULL))
    }
}

packNameOutput <- function() {
    out <- paste("\ndefault:\ttargets affy, cdna and exprs.\n",
                   "exprs:\t\tpackages Biobase, annotate, genefilter, ",
                   "geneploter, edd, \n\t\tROC, multtest, pamr vsn, and limma.\n",
                 "affy:\t\tpackages affy, affydata, ",
                 "annaffy, affyPLM, makecdfenv,\n\t\t",
                 "and matchprobes plus 'exprs'.\n",
                 "cdna:\t\tpackages marray, vsn,",
                 " plus 'exprs'.\nprog:\t\tpackages graph, hexbin, ",
                 "externalVector.\n",
                 "graph:\t\tpackages graph, Rgraphviz, RBGL",
                 "\nwidgets:\tpackages tkWidgets, widgetTools,",
                 " DynDoc.\ndesign:\t\tpackages daMA and factDesign\n",
                 "externalData:\tpackage externalVector.\n",
                 "database:\tAnnBuilder, SAGElyzer, Rdbi and ",
                 "RdbiPgSQL.\n",
                 "analyses:\tpackages Biobase, ctc, daMA, edd, ",
                 "factDesign,\n\t\tgenefilter, geneplotter, globaltest, ",
                 "gpls, limma,\n\t\tRMAGEML, multtest, pamr, wvalue, ROC, ",
                 "siggenes and splicegear.\n",
                 "annotation:\tpackages annotate, AnnBuilder, ",
                 "humanLLMappings\n\t\tKEGG, GO, SNPtools, ",
                 "makecdfenv and ontoTools.",
                 "\nproteomics:\tpackages gpls, PROcess and apComplex.",
                 "\narrayCGH:\tpackages aCGH, DNAcopy, repeated, ",
                 "and rmutil.",
                 "\nall:\t\tAll of the Bioconductor packages.\n",
                 sep="")
    out
}

## This function put together a vector containing Bioconductor's
## packages based on a defined libName

## Returns the mode that is going to be used to call download.file
## depending on the platform
getMode <- function(platform){
    switch(platform,
           "unix" = return("w"),
           "windows" = return("wb"),
           stop(paste(platform,"is not currently supported")))
}

## Installs a given package
installPack <- function(platform, fileName, destdir=NULL){
    if(platform == "unix"){
        cmd <- paste(file.path(R.home(), "bin", "R"),
                     "CMD INSTALL")
        if (!is.null(destdir))
            cmd <- paste(cmd, "-l", destdir)
        cmd <- paste(cmd, fileName)
        system(cmd)
    }else{
        if(platform == "windows"){
            zip.unpack(fileName, .libPaths()[1])
        }else{
            stop(paste(platform,"is not currently supported"))
        }
    }
}

## Returns the surce url for a given package
getDLURL <- function(pakName, rep, platform){
    temp <- rep[rep[, "Package"] == pakName]
    names(temp) <- colnames(rep)
    path <- switch(platform,
                   "unix" = temp[names(temp) == "SourceURL"],
                   "windows" = temp[names(temp) == "WIN32URL"],
                   stop(paste(platform,"is not currently supported")))
}

## Returns the description file (PACKAGE) that contains the name,
## version number, url, ... of Bioconductor packages.
getPACKAGES <- function (develOK, bioCoption){
    URL <- getReposURL(develOK,"/PACKAGES", bioCoption)
    con <- url(URL)
    options(show.error.messages = FALSE)
    tryMe <- try(read.dcf(con))
    options(show.error.messages = TRUE)

    if(inherits(tryMe, "try-error"))
       stop(paste("The url:",URL,
                  "does not seem to have a valid PACKAGES file."))

    close(con)
    return(tryMe)
}

## Returns the url for some files that are needed to perform the
## functions. name is added to teh end of the URL
getReposURL <- function(develOK, name="", bioCoption){
    if (develOK)
        URL <- paste(bioCoption, "repository/devel/package",
                     name, sep ="/")
    else
        URL <- paste(bioCoption, "repository/release1.5",
                     "/package",name,sep="/")
    URL
}

## Returns the file name with the destination path (if any) attached
getFileName <- function(url, destdir){
    temp <- unlist(strsplit(url, "/"))
    if(is.null(destdir))
        return(temp[length(temp)])
    else
        return(file.path(destdir, temp[length(temp)]))
}

## getBioC has to check to see if "reposTools" has
## already been loaded and generates a message if any has.
checkReposTools <- function(PACKAGES){
    pkgVers <- PACKAGES[,"Version"]

    ## First get package version
    ## !!! Not yet using VersionNumber classes here
    ## !!! bootstrapping issue as this comes from reposTools
    ## !!! use compareVersion for now
    if ("reposTools" %in% installed.packages()[,"Package"]) {
        curVers <- packageDescription("reposTools",fields="Version")
        if (compareVersion(curVers,pkgVers) < 0) {
            if ("package:reposTools" %in% search()) {
                error <- paste("reposTools is out of date but",
                               " currently loaded in your R session.",
                               "\nIf you would like to continue,",
                               " please either detach this package",
                               " or restart\nyour R seesion before",
                               " running getBioC.",sep="")
                stop(error)
            }
        }
        else
            return(FALSE)
    }

    return(TRUE)
}

## From reposTools
GBCuserQuery <- function(msg, allowed=c("yes","y","no","n")) {
    ## Prompts the user with a string and for an answer
    ## repeats until it gets allowable input
    repeat {
        allowMsg <- paste("[",paste(allowed,collapse="/"),
                          "] ", sep="")
        outMsg <- paste(msg,allowMsg)
        cat(outMsg)
        ans <- readLines(n=1)
        if (ans %in% allowed)
            break
        else
            cat(paste(ans,"is not a valid response, try again.\n"))
    }
    ans
}
