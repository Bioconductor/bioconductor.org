# Writing Robust Web-Query Functions

Web resources (files) have been known to move location and many urls can be
temporarily unavailable. Functions that query web resources, for reading or
download, should handle such situations and avoid the infinite loop, using up
all available _R_ connections and unclear error messages.

## Guiding Principles

1. Download a file of reasonable size. Use `system.time()` to estimate the
   download time. Remember the package should require less than 5 minutes to
   run `R CMD check --no-build-vignettes`.

2. Set a limit on the number of times the function tries a url. Avoid
   `while()` statements that have no guaranteed termination. These
   become infinite loops and eventually `timeout` on the build report.

3. Supply an informative error message.

## Sample Functions

These functions try a url 3 times and then fail. The last error message is 
reported and warnings are allowed through.

    readURL <- function(URL, n=-1L) {
        tries <- 0L
        msg <- character()
        while (tries < 3L) {
            URLdata <- tryCatch(readLines(URL, n), error=identity)
            if (!inherits(URLdata, "error"))
                break
            tries <- tries + 1L
        }
        if (tries == 3L)
            stop("failed to get URL after 3 tries:",
                 "\n  url: ", URL,
                 "\n  error: ", conditionMessage(URLdata))
        URLdata
    }

This function uses RCurl::getURL() but easily could have been written with
download.file(), httr::GET() or curl::curl() etc.

    getURL <- function(URL) {
        tries <- 0L
        msg <- character()
        while (tries < 3L) {
            URLdata <- tryCatch(getURL(URL, dirlistonly = TRUE), error=identity)
            if (!inherits(URLdata, "error"))
                break
            tries <- tries + 1L
        }
        if (tries == 3L)
            stop("failed to get URL after 3 tries:",
                 "\n  url: ", URL,
                 "\n  error: ", conditionMessage(URLdata))
        URLdata
    }
