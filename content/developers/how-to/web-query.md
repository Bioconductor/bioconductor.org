# Querying Web Resources

Packages that rely on access to web resources need to be written
carefully. Web resources can change location, can be temporarily
unavailable, or can be very slow to access and retrieve. Functions
that query web resources, should anticipate and handle such situations
gracefully -- failing quickly and clearly when the resource is not
available in a reasonable time frame. Some avoidable problems seen in
_Bioconductor_ package code include infinite loops, use of all
available _R_ connections, and unclear error messages.

## Guiding Principles

Remember the _Bioconductor_ packages are built nightly across multiple
operating systems, and that users benefit from easy-to-run vignettes
and examples.

1. Download files of reasonable size. Use `system.time()` to estimate the
   download time. Remember the package should require less than 5 minutes to
   run `R CMD check --no-build-vignettes`.

2. Set a limit on the number of times the function tries a URL. Avoid
   `while()` statements that have no guaranteed termination. These
   become infinite loops and eventually result in build-system `TIMEOUT`s.

3. Supply an informative error message.

## Template for Resource Queries

This function can serve as a template for appropriate resource
retrieval. It tries to retrieve the resource one or several times before
failing, and takes as arguments:

- `URL`, the resource to be queried, typically `character(1)` or
  `url()`.
- `FUN`, the function to be used to query the resource. Examples might
  include `readLines()`, `download.file()`, `httr::GET()`,
  `RCurl::getURL()`.
- `...`: additional arguments used by `FUN`.
- `N.TRIES`: the number of times the URL will be attempted; only under
  exceptional circumstances might this differ from its default value.

The return value is the retrieved resource. If resource retrieval
fails, the function indicates the failure, including the condition
(error) message on the last attempt. Warnings propagate to the user in
the normal way.

    getURL <- function(URL, FUN, ..., N.TRIES=1L) {
        N.TRIES <- as.integer(N.TRIES)
        stopifnot(length(N.TRIES) == 1L, !is.na(N.TRIES))

        while (N.TRIES > 0L) {
            result <- tryCatch(FUN(URL, ...), error=identity)
            if (!inherits(result, "error"))
                break
            N.TRIES <- N.TRIES - 1L
        }

        if (N.TRIES == 0L) {
            stop("'getURL()' failed:",
                 "\n  URL: ", URL,
                 "\n  error: ", conditionMessage(result))
        }

        result
    }


Base _R_ functions using `url()` connections respect
`getOption("timeout")`; see `?url` for details.

`FUN` might be implemented to retrieve the resource and test for
status, e.g.,

    FUN <- function(URL, ...) {
        response <- httr::GET(URL, timeout(getOption("timeout")), ...)
        stop_for_status(response)
        response
    }
