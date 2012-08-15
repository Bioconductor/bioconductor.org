Using The Devel Version of Bioconductor
=======================================

With R going to a yearly release schedule and Bioconductor keeping its
twice-yearly release schedule, the same version of R (2.15) can be used
with two different versions of Bioconductor (2.10, release, and 2.11, devel).
The version number of the BiocInstaller package is what is used to determine
whether to download packages from the BioC 2.10 or 2.11 repositories.
In keeping with
[Bioconductor versioning conventions](/developers/version-numbering),
if the middle number
(y in x.y.z) is even, the package is part of a release version; if odd, it's
part of a devel version.
    
By default, when BiocInstaller is first installed on R-2.15, it will be
set up to download BioC 2.10 packages. 
    
If you want to change this, you can run the <code>useDevel()</code>
function. With argument <code>TRUE</code> (the default), it will download
the devel version of BiocInstaller and subsequently all packages
downloaded with <code>biocLite()</code> will be from the BioC 2.11
(devel) repository. You should run <code>useDevel()</code> only once.

Note: The information on this page can also be accessed from within R
as follows:

    library(BiocInstaller)
    ?useDevel
    
It is possible to keep BioC 2.10 and 2.11 libraries separate, within
the same installation of R.
    
The trick is to use the <code>R_LIBS_USER</code> environment variable.
First, create two separate directories for your BioC release
and devel packages. Suggested directory names are as follows:
    
Linux:
    
    ~/R/x86_64-unknown-linux-gnu-library/2.15-bioc-release
    
    ~/R/x86_64-unknown-linux-gnu-library/2.15-bioc-devel

Mac OS:
    
    ~/Library/R/2.15-bioc-release/library
    
    ~/Library/R/2.15-bioc-devel/library

Windows:
    
    C:\Users\YOUR_USER_NAME\Documents\R\win-library\2.15-bioc-release
    
    C:\Users\YOUR_USER_NAME\Documents\R\win-library\2.15-bioc-devel
    
(change YOUR_USER_NAME to your user name)
    

You can then invoke "R for bioc-devel" or "R for bioc-release" from
the command line as follows:

Linux:
    
    R_LIBS_USER=~/R/x86_64-unknown-linux-gnu-library/2.15-bioc-release R
    
    R_LIBS_USER=~/R/x86_64-unknown-linux-gnu-library/2.15-bioc-devel R
    
    
Mac OS X:
    
    R_LIBS_USER=~~/Library/R/2.15-bioc-release/library R
    R_LIBS_USER=~~/Library/R/2.15-bioc-devel/library R

Windows:
    
    cmd /C "set R_LIBS_USER=C:\Users\YOUR_USER_NAME\Documents\R\win-library\2.15-bioc-release &&  R"
    
    cmd /C "set R_LIBS_USER=C:\Users\YOUR_USER_NAME\Documents\R\win-library\2.15-bioc-devel &&  R"
    
(Note: this assumes that R.exe is in your PATH.)


If you launch R in this way and then invoke <code>.libPaths()</code>,
you'll see that the first item is your special release or
devel directory. Packages will be installed to that directory and that
is the first place that <code>library()</code> will look for them. 
<code>biocLite()</code> and <code>install.packages()</code>
respect this setting; <code>update.packages()</code> attempts to update
packages in the directory where the current package is installed.


On Linux and Mac OS X, you can create a bash alias to save typing. Add the
following to your ~/bash_profile:
    

Linux
    
    alias Rdevel='R_LIBS_USER=~/R/x86_64-unknown-linux-gnu-library/2.15-bioc-devel R'

    alias Rrelease='R_LIBS_USER=~/R/x86_64-unknown-linux-gnu-library/2.15-bioc-release R'
    
Mac OS X
    
    alias Rdevel='R_LIBS_USER=~/Library/R/2.15-bioc-devel/library R'
    alias Rrelease='R_LIBS_USER=~/Library/R/2.15-bioc-release/library R'

    
You can then invoke these from the command line as
    
    Rdevel
    
...and...
    
    Rrelease


On Windows, you can create two shortcuts, one for devel and one for
release. Go to My Computer and navigate to a directory that is in your
PATH. Then right-click and choose New->Shortcut.
    
in the "type the location of the item" box, put:
    
    cmd /C "set R_LIBS_USER=C:\Users\YOUR_USER_NAME\Documents\R\win-library\2.15-bioc-release &&  R"
    
...for release and
    
    cmd /C "set R_LIBS_USER=C:\Users\YOUR_USER_NAME\Documents\R\win-library\2.15-bioc-devel &&  R"

...for devel.

(again, it's assumed R.exe is in your PATH)

Click "Next".

In the "Type a name for this shortcut" box, type

    Rdevel

or

    Rrelease
    
You can invoke these from the command line as
    
    Rdevel.lnk
    
...and...
    
    Rrelease.lnk
    
(You must type in the .lnk extension.)
    
Because <code>R_LIBS_USER</code> is an environment variable, its value should be
inherited by any subprocesses started by R, so they should do the
right thing as well.
    
    

