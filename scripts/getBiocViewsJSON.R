
# todo  - add code that makes sure biocViews and rjson are installed

# run me like this:
# R CMD BATCH -q --vanilla --no-save --no-restore '--args versions=c("2.6","2.7") outdir="directory/to/write/json"' scripts/getBiocViewsJSON.R log.txt
# where versions are the RELEASE and DEVEL versions respectively


args <- (commandArgs(T))
if (length(args) == 0) {
    print("No arguments supplied.") #todo - more useful usage message
    q("no")
}

for (i in 1:length(args)) {
    eval(parse(text=args[[i]]))
}

if (!exists("versions") && !exists("outdir")) {
    print("'outdir' and 'versions' arguments required")
    q("no")
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
    reposUrl = paste("http://bioconductor.org/packages/", version, "/bioc", sep="")
    
    biocViews <- getBiocViews(reposUrl, biocViewsVocab, "No View Provided")
    outputDirectory <- paste(outdir, "/", version, sep="")
    dir.create(outputDirectory, recursive=T)
    
    biocViewsFile <- paste(outputDirectory, "biocViews.json", sep="/")
    count <- 1
    all <- list()
    result <- lapply(biocViews, getItem)

    
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
    
    
    
}




