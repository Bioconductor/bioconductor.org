Using the 'Devel' Version of _Bioconductor_
===========================================

Which version of R?
-------------------

Package developers want to develop against the version of _R_ that
will be available to users when the _Bioconductor_ devel branch
becomes the _Bioconductor_ release branch.

_R_ has a '.y' release in x.y.z every year, but Bioconductor has a .y
release (where current devel becomes release) every 6 months.

When the next (typically mid-April) .y releases of _R_ and
_Bioconductor_ coincide, Bioc-devel is based on R-devel.

When the next (typically mid-October) .y release of _Bioconductor_
corresponds to no change in _R_'s y, Bioc-devel is based on release R.

This means that, from mid-October through mid-April, _Bioconductor_
developers should be developing against R-devel. From mid-April to
mid-October, developers should use R-release (actually, the R snapshot
from the R-x-y-branch) for _Bioconductor_ development.

<% rvers = config[:r_version_associated_with_release].split(".")[0..1].join(".") %>

[BiocInstaller]: /packages/BiocInstaller

Using Bioc-devel
----------------

<!--

In order to use the `devel` version of _Bioconductor_ during the current
release cycle, you must install `R-devel`:

* [Source](https://stat.ethz.ch/R/daily/)
* [Mac OS X](http://r.research.att.com/)
* [Windows](https://cran.r-project.org/bin/windows/base/rdevel.html)

Then, make sure that your version of [BiocInstaller][] is current and
your packages up-to-date. Do this by removing all versions of
[BiocInstaller][]

    remove.packages("BiocInstaller")  # repeat until R says there is no
                                      # package 'BiocInstaller' to remove
    source("https://bioconductor.org/biocLite.R")  # install correct version
    BiocInstaller::biocValid()

-->

In order to use the `devel` version of _Bioconductor_ during the current
release cycle, simply call the function `useDevel()` (from the
`BiocInstaller`) package:

    ## In R-<%= rvers %>
    library(BiocInstaller)
    useDevel()
    biocValid()              # checks for out of date packages
    biocLite()               # (optional) updates out of date packages

After doing this, all packages will be installed from the `devel`
(BioC <%=config[:devel_version]%>) repository.

If you also want to work with the `release` version of _Bioconductor_
(<%=config[:release_version]%>), we recommend maintaining two separate installations of R
<%= rvers %>, one to be used with _Bioconductor_ <%=config[:release_version]%> (BioC-release) and the
other to be used with _Bioconductor_ <%=config[:devel_version]%> (BioC-devel). Run `useDevel()`
as described above in this latter installation.

An easy way to do this is to have two separate installations of R-<%=config[:release_version]%>.

A more complicated way is to use the `R_LIBS_USER` environment
variable.  First, create two separate directories. Suggested directory
names are Linux:

    ~/R/x86_64-unknown-linux-gnu-library/<%=config[:release_version]%>-bioc-release
    ~/R/x86_64-unknown-linux-gnu-library/<%=config[:release_version]%>-bioc-devel

Mac OS:

    ~/Library/R/<%=config[:release_version]%>-bioc-release/library
    ~/Library/R/<%=config[:release_version]%>-bioc-devel/library

and Windows:

    C:\Users\YOUR_NAME\Documents\R\win-library\<%=config[:release_version]%>-bioc-release
    C:\Users\YOUR_NAME\Documents\R\win-library\<%=config[:release_version]%>-bioc-devel

(change `YOUR_NAME` to your user name)

Invoke "R for bioc-devel" or "R for bioc-release" from the command
line on Linux:

    R_LIBS_USER=~/R/x86_64-unknown-linux-gnu-library/<%=config[:release_version]%>-bioc-release R
    R_LIBS_USER=~/R/x86_64-unknown-linux-gnu-library/<%=config[:release_version]%>-bioc-devel R

Mac OS X:

    R_LIBS_USER=~~/Library/R/<%=config[:release_version]%>-bioc-release/library R
    R_LIBS_USER=~~/Library/R/<%=config[:release_version]%>-bioc-devel/library R

and Windows (assuming that R.exe is in PATH):

    cmd /C "set R_LIBS_USER=C:\Users\YOUR_NAME\Documents\R\win-library\<%=config[:release_version]%>-bioc-release &&  R"
    cmd /C "set R_LIBS_USER=C:\Users\YOUR_NAME\Documents\R\win-library\<%=config[:release_version]%>-bioc-devel &&  R"

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

    alias Rdevel='R_LIBS_USER=~/R/x86_64-unknown-linux-gnu-library/<%=config[:release_version]%>-bioc-devel R'
    alias Rrelease='R_LIBS_USER=~/R/x86_64-unknown-linux-gnu-library/<%=config[:release_version]%>-bioc-release R'

or Mac OS X

    alias Rdevel='R_LIBS_USER=~/Library/R/<%=config[:release_version]%>-bioc-devel/library R'
    alias Rrelease='R_LIBS_USER=~/Library/R/<%=config[:release_version]%>-bioc-release/library R'

Invoke R from the command line as `Rdevel` or `Rrelease`.

On Windows, create two shortcuts, one for `release` and one for
`devel`. For `devel` (do similar steps for `release`) go to My
Computer and navigate to a directory that is in your PATH. Then
right-click and choose New->Shortcut.  In the "type the location of
the item" box, put:

    cmd /C "set R_LIBS_USER=C:\Users\YOUR_NAME\Documents\R\win-library\<%=config[:release_version]%>-bioc-devel &&  R"

(again, it's assumed R.exe is in your PATH) Click "Next", and in the
"Type a name for this shortcut" box, type

    Rdevel

Invoke these from the command line as `Rdevel.lnk`.

Because `R_LIBS_USER` is an environment variable, its value should be
inherited by any subprocesses started by R, so they should do the
right thing as well.
