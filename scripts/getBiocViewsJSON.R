
# todo  - add code that makes sure biocViews and rjson are installed

# run me like this:
# R CMD BATCH -q --vanilla --no-save --no-restore \
# '--args versions=c("2.5","2.6","2.7"); outdir="directory/to/write/json"; devel_repos=(c("bioc", "data/experiment", "data/anotation")); devel_version=2.7;' \
# scripts/getBiocViewsJSON.R log.txt
# where versions are all bioC versions for which we want to retrieve JSON.


args <- (commandArgs(TRUE))
if (length(args) == 0) {
    print("No arguments supplied.") #todo - more useful usage message
    q("no")
}

cmds = args[1:4] # modify if # of args changes

for (i in 1:length(cmds)) {
    eval(parse(text=cmds[i]))
}


if (!exists("versions") && !exists("outdir") && !exists("devel_repos")
  && !exists("devel_version")) {
    print("'outdir', 'versions', 'devel_repos', and 'devel_version' arguments required")
    q("no")
}

if (!require(biocViews)) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("biocViews")
} 

if (!require(rjson)) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("rjson")
}


library(biocViews)
library(rjson)

data(biocViewsVocab)

getPackageAsList = function(packageDetail) {
    
    slotNames <- slotNames(packageDetail)
    
    packageAsList <- lapply(slotNames, function(slotName){
        l <- list()
        slot <- slot(packageDetail, slotName);
        if (length(slot)) {
            l[slotName] <- slot
        } else {
            l[slotName] = list()#"<unknown>"
        }
        l
    })
    
    
    
    packageAsList
}

getPackageAsEnv <- function(packageDetail) {
    e <- new.env()
    slotNames <- slotNames(packageDetail)
    for (slotName in slotNames) {
        slot <- slot(packageDetail, slotName)
        assign(slotName, slot, envir=e)
    }
    return(e);
}


getPackages <- function(y) {
    if (length(y@packageList)) {
        packagesForThisNode <- lapply(y@packageList, function(package){
            p <- get("packages", envir=.GlobalEnv)
            #p[package@Package] = package
            assign(package@Package, getPackageAsEnv(package), envir=p)
            assign("packages", p, envir=.GlobalEnv) # i know this is evil
        })
        invisible(NULL)
    }
    invisible(NULL)
}



getPackageName <- function(packageDetail) {
    packageDetail@Package
}


getItem <- function(y) {
	packageList <- list()
	if (length(y@packageList)) {
		packageList <- lapply(y@packageList, getPackageName)
	} else {
		packageList <- list()
	}
	item <- list(name=y@name, subViews=y@subViews, parentViews=y@parentViews, packageList=packageList)
	all[[count]] <- item
	count <- count + 1
	all
	#cat(paste("json = ", toJSON(item), "\n"))
}



for (version in versions) {
    
    repos = c("bioc", "data/annotation", "data/experiment")
    if (version == devel_version)
        repos <- devel_repos
        
    for (repo in repos) {
        if (repo == "bioc") defaultView = "Software"
        if (repo == "data/annotation") defaultView = "AnnotationData"
        if (repo == "data/experiment") defaultView = "ExperimentData"
      
        reposUrl = paste("http://bioconductor.org/packages/", version, "/", repo, sep="")

            biocViews <- getBiocViews(reposUrl, biocViewsVocab, defaultView)
            count <- 1
            all <- list()
            result <- lapply(biocViews, getItem)


            outputDirectory <- paste(outdir, version, repo, sep="/")
            dir.create(outputDirectory, recursive=T)

            biocViewsFile <- paste(outputDirectory, "biocViews.json", sep="/")


            conn <- file(biocViewsFile, "w")
            json <- toJSON(result)
            cat(json, file=conn)
            close(conn)

            packages <- new.env()
        #    packages <- list()

            lapply(biocViews, getPackages)
            packagesFile <- paste(outputDirectory, "packages.json", sep="/")
            conn <- file(packagesFile, "w")
            json <- toJSON(packages)
            cat(json, file=conn)
            close(conn)
            warnings()

        
    }
    
    
    
warnings()    
    
}


print("-------")
warnings()
