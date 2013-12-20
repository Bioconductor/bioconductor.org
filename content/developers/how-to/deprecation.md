![](/images/icons/magnifier.gif)Function Deprecation Guidelines
==================================================

** About Deprecation **

In the normal course of software development, a function, method, or 
other object may be removed. Here are some guidelines for ensuring
that this process is minimally disruptive to your users.

** What to Deprecate **

Any function, method, or data that is no longer used.

** When to Follow These Guidelines **

If you introduce a function into the devel branch of your package, 
then soon after decide not to use it, you may simply remove the function
without following these guidelines. It is expected that the devel
branch is unstable and subject to API changes without notice (though
you may decide to communicate these changes to your users via the
[bioconductor](/help/mailing-list) mailing list).

However, if a function has existed in at least one released version
of Bioconductor, these guidelines must be followed.

** How To Deprecate **

The full process for removing a function takes three release cycles
(18 months).

*** Step 1: Deprecate the function ***

When you first decide to eliminate a function, you should mark it as
deprecated in the devel branch. Do this by calling <code>.Deprecated()</code>
inside the function. To do this, you must provide a replacement function 
which should be used in place of the old function. Example:

    myOldFunc <- function()
    {
        .Deprecated("myNewFunc")
        ## use new function, or remainder of myOldFunc
    }

This causes a warning to be emitted whenever a user calls 
<code>myOldFunc()</code>. See <code>?.Deprecated</code> for more information.

Indicate in the man page of the old function that it has been
deprecated, and suggest a replacement function. Be sure the old
function is not called in man page examples or vignette code chunks; R
CMD check should report this.

    \name{MyPkg-deprecated}
    \alias{MyPkg-deprecated}
    \title{Deprecated functions in package \sQuote{MyPkg}}

    \description{
      These functions are provided for compatibility with older versions
      of \sQuote{MyPkg} only, and will be defunct at the next release.
    }

    \details{
    
      The following functions are deprecated and will be made defunct; use
      the replacement indicated below:
      \itemize{
    
        \item{myOldFunc: \code{\link{newFunc}}}
    
      }
    }
    
*** Step 2: Mark the function as defunct ***

In the next release cycle, after your function has been deprecated, it
must be made defunct in the devel branch.  This means a call to the
old function will return an informative error but not run any
additional code. Example:

    myOldFunc <- function()
    {
        .Defunct("myNewFunc")
    }

See <code>?Defunct</code> for more information.

Remove the documentation of the defunct function, and add to a man
page such as the following:

    \name{MyPkg-defunct}
    \alias{myOldFunc}
    \title{Defunct functions in package \sQuote{pkg}}
    \description{These functions are defunct and no longer available.}

    \details{
      Defunct functions are: \code{myOldFunc}
    }

*** Step 3: Remove the function ***

In the next release cycle, after your function has been marked as defunct,
remove it entirely from your package R code and NAMESPACE in the devel
branch. Also remove any man page content that documents the function.

Leave the man page from the previous step in place so that 

    help("MyPkg-defunct")
    
still shows the list of defunct functions and their appropriate replacements.

