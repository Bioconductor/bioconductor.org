![](/images/icons/magnifier.gif)Function Deprecation Guidelines
==================================================

** About Deprecation **

In the normal course of software development, a function, method, or 
other object may be removed. Here are some guidelines for ensuring
that this process is minimally disruptive to your users.

** What to Deprecate **

Any function or method that is no longer used.

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

Indicate in the man page of the old function that it has been deprecated, and
suggest a replacement function. Be sure the old function is not called in 
man page examples or vignette code chunks; <code>R CMD check</code> should
report this.


Ensure that if the user enters:

    help("myOldFunc-deprecated")
    
that they are shown the manual page for the original old function,
which should mention that the function is deprecated and point users towards
the replacement function.

Additionally, a man page reachable by

    help("myPackageName-deprecated")
    
should list the functions in your package that are deprecated and
discuss the appropriate replacement functions. See 

    help("base-deprecated")
    
for an example.

*** Step 2: Mark the function as defunct ***

In the next release cycle, after your function has been deprecated,
it must be made defunct in the devel branch.
This means a call to the old function will
return an informative error but not run any additional code. Example:

    myOldFunc <- function()
    {
        .Defunct("myNewFunc")
    }

See <code>?Defunct</code> for more information.

Additionally, ensure that a call to:

    help("myPackageName-defunct")
    
shows a manual page that lists the defunct functions in your package, 
and discusses the appropriate replacements. See

    help("base-defunct")
    
for an example. Make sure that 

    help("myPackageName-deprecated")
    
<b>no longer</b> lists your function.

*** Step 3: Remove the function ***

In the next release cycle, after your function has been marked as defunct,
remove it entirely from your package R code and NAMESPACE in the devel
branch. Also remove any man page content that documents the function.

Leave the man page from the previous step in place so that 

    help("myPackageName-defunct")
    
still shows the list of defunct functions and their appropriate replacements.

