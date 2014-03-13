Using the `devel` Version of Bioconductor
=========================================

In this release cycle, use R-devel with Bioconductor 2.14
(`devel`). From within R-alpha, simply invoke

    source("http://bioconductor.org/biocLite.R")
    biocLite()
    
Separate libraries
------------------

Keep Bioconductor `release` and `devel` packages in
separate libraries. This is the default when using R-3.0.2 for
Bioc-2.13 (`release`), and R-alpha for Bioc-2.14 (`devel`).

<!--
An easy way to do this is to have two separate installation of R-3.0.

A more complicated way is to use the `R_LIBS_USER` environment
variable.  First, create two separate directories. Suggested directory
names are Linux:
    
    ~/R/x86_64-unknown-linux-gnu-library/3.0-bioc-release
    ~/R/x86_64-unknown-linux-gnu-library/3.0-bioc-devel

Mac OS:
    
    ~/Library/R/3.0-bioc-release/library
    ~/Library/R/3.0-bioc-devel/library

and Windows:
    
    C:\Users\YOUR_NAME\Documents\R\win-library\3.0-bioc-release
    C:\Users\YOUR_NAME\Documents\R\win-library\3.0-bioc-devel
    
(change `YOUR_NAME` to your user name)
    
Invoke "R for bioc-devel" or "R for bioc-release" from the command
line on Linux:
    
    R_LIBS_USER=~/R/x86_64-unknown-linux-gnu-library/3.0-bioc-release R
    R_LIBS_USER=~/R/x86_64-unknown-linux-gnu-library/3.0-bioc-devel R
    
Mac OS X:
    
    R_LIBS_USER=~~/Library/R/3.0-bioc-release/library R
    R_LIBS_USER=~~/Library/R/3.0-bioc-devel/library R

and Windows (assuming that R.exe is in PATH):
    
    cmd /C "set R_LIBS_USER=C:\Users\YOUR_NAME\Documents\R\win-library\3.0-bioc-release &&  R"
    cmd /C "set R_LIBS_USER=C:\Users\YOUR_NAME\Documents\R\win-library\3.0-bioc-devel &&  R"
    
When correctly configured, R's `.libPaths()` function will return the
`release` or `devel` directory as its first entry. Packages are
installed to that directory, and that is the first place that
`library()` looks for them.  <code>biocLite()</code> and
<code>install.packages()</code> respect this setting;
<code>update.packages()</code> attempts to update packages in the
directory where the current package is installed.

Aliases
-------

On Linux and Mac OS X, you can create a bash alias to save typing. Add the
following to your ~/bash_profile on Linux:
    
    alias Rdevel='R_LIBS_USER=~/R/x86_64-unknown-linux-gnu-library/3.0-bioc-devel R'
    alias Rrelease='R_LIBS_USER=~/R/x86_64-unknown-linux-gnu-library/3.0-bioc-release R'
    
or Mac OS X
    
    alias Rdevel='R_LIBS_USER=~/Library/R/3.0-bioc-devel/library R'
    alias Rrelease='R_LIBS_USER=~/Library/R/3.0-bioc-release/library R'
    
Invoke R from the command line as `Rdevel` or `Rrelease`.

On Windows, create two shortcuts, one for `release` and one for
`devel`. For `devel` (do similar steps for `release`) go to My
Computer and navigate to a directory that is in your PATH. Then
right-click and choose New->Shortcut.  In the "type the location of
the item" box, put:

    cmd /C "set R_LIBS_USER=C:\Users\YOUR_NAME\Documents\R\win-library\3.0-bioc-devel &&  R"

(again, it's assumed R.exe is in your PATH) Click "Next", and in the
"Type a name for this shortcut" box, type

    Rdevel
    
Invoke these from the command line as `Rdevel.lnk`.
    
Because `R_LIBS_USER` is an environment variable, its value should be
inherited by any subprocesses started by R, so they should do the
right thing as well.
-->
