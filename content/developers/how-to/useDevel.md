Using the `devel` Version of Bioconductor
=========================================

In order to use the `devel` (3.1) version of Bioconductor, you
must install `R-devel`, available from:

* [Source](ftp://ftp.stat.math.ethz.ch/Software/R/)
* [Mac](http://r.research.att.com/)
* [Windows](http://cran.r-project.org/bin/windows/base/rdevel.html)

Then, in R-devel, do the following:

    source("http://bioconductor.org/biocLite.R")

Now, any package you install will be the `devel` version.


<!--
In order to use the `devel` version of Bioconductor, simply call
the function `useDevel()` (from the `BiocInstaller`) package:

    ## In R-3.1.0
    library(BiocInstaller) 
    useDevel()

After doing this, all packages will be installed from the `devel`
(BioC 3.0) repository.

If you also want to work with the `release` version of Bioconductor
(2.14), we recommend maintaining two separate installations of R
3.1.0, one to be used with Bioconductor 2.14 (BioC-release) and the
other to be used with Bioconductor 3.0 (BioC-devel). Run `useDevel()`
as described above in this latter installation.
-->

<!--
An easy way to do this is to have two separate installations of R-3.1.

A more complicated way is to use the `R_LIBS_USER` environment
variable.  First, create two separate directories. Suggested directory
names are Linux:
    
    ~/R/x86_64-unknown-linux-gnu-library/3.1-bioc-release
    ~/R/x86_64-unknown-linux-gnu-library/3.1-bioc-devel

Mac OS:
    
    ~/Library/R/3.1-bioc-release/library
    ~/Library/R/3.1-bioc-devel/library

and Windows:
    
    C:\Users\YOUR_NAME\Documents\R\win-library\3.1-bioc-release
    C:\Users\YOUR_NAME\Documents\R\win-library\3.1-bioc-devel
    
(change `YOUR_NAME` to your user name)
    
Invoke "R for bioc-devel" or "R for bioc-release" from the command
line on Linux:
    
    R_LIBS_USER=~/R/x86_64-unknown-linux-gnu-library/3.1-bioc-release R
    R_LIBS_USER=~/R/x86_64-unknown-linux-gnu-library/3.1-bioc-devel R
    
Mac OS X:
    
    R_LIBS_USER=~~/Library/R/3.1-bioc-release/library R
    R_LIBS_USER=~~/Library/R/3.1-bioc-devel/library R

and Windows (assuming that R.exe is in PATH):
    
    cmd /C "set R_LIBS_USER=C:\Users\YOUR_NAME\Documents\R\win-library\3.1-bioc-release &&  R"
    cmd /C "set R_LIBS_USER=C:\Users\YOUR_NAME\Documents\R\win-library\3.1-bioc-devel &&  R"
    
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
    
    alias Rdevel='R_LIBS_USER=~/R/x86_64-unknown-linux-gnu-library/3.1-bioc-devel R'
    alias Rrelease='R_LIBS_USER=~/R/x86_64-unknown-linux-gnu-library/3.1-bioc-release R'
    
or Mac OS X
    
    alias Rdevel='R_LIBS_USER=~/Library/R/3.1-bioc-devel/library R'
    alias Rrelease='R_LIBS_USER=~/Library/R/3.1-bioc-release/library R'
    
Invoke R from the command line as `Rdevel` or `Rrelease`.

On Windows, create two shortcuts, one for `release` and one for
`devel`. For `devel` (do similar steps for `release`) go to My
Computer and navigate to a directory that is in your PATH. Then
right-click and choose New->Shortcut.  In the "type the location of
the item" box, put:

    cmd /C "set R_LIBS_USER=C:\Users\YOUR_NAME\Documents\R\win-library\3.1-bioc-devel &&  R"

(again, it's assumed R.exe is in your PATH) Click "Next", and in the
"Type a name for this shortcut" box, type

    Rdevel
    
Invoke these from the command line as `Rdevel.lnk`.
    
Because `R_LIBS_USER` is an environment variable, its value should be
inherited by any subprocesses started by R, so they should do the
right thing as well.
-->
