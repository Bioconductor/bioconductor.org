Using the 'Devel' Version of _Bioconductor_
===========================================

Which version of R?
-------------------

Package authors should develop against the version of _R_ that will be
available to users when the _Bioconductor_ devel branch becomes the
_Bioconductor_ release branch.

_R_ has a '.y' release in x.y.z every year (typically mid-April), but
_Bioconductor_ has a .y release (where current devel becomes release)
every 6 months (mid-April and mid-October).

This means that, from mid-October through mid-April, _Bioconductor_
developers should be developing against R-devel. From mid-April to
mid-October, developers should use R-release (actually, the R snapshot
from the R-x-y-branch) for _Bioconductor_ development.

[BiocInstaller]: /packages/BiocInstaller

Using 'bioc-devel' during mid-April to mid-October
--------------------------------------------------

In order to use the 'bioc-devel' version of _Bioconductor_ during the
mid-April to mid-October release cycle, use the release version of _R_
and invoke the function `useDevel()` (from the `BiocInstaller`)
package:

    library(BiocInstaller)
    useDevel()
    biocValid()              # checks for out of date packages
    biocLite()               # (optional) updates out of date packages

After doing this, all packages will be installed from the 'bioc-devel'
repository.

One way to work with 'bioc-release' and 'bioc-devel' versions of
_Bioconductor_, is to have two separate installations of the release
version of _R_, one to be used with 'bioc-release' and the other to be
used with 'bioc-devel'. Run `useDevel()` as described above in the
'bioc-devel' installation.

Another way of working with release and devel versions of
_Bioconductor_ is with a single _R_ installation, using the
`R_LIBS_USER` environment variable to create separate 'bioc-release'
and 'bioc-devel' libraries.  Before installing any packages, create
two separate directories. Suggested directory names are Linux:

    ~/R/x86_64-unknown-linux-gnu-library/<R version>-bioc-release
    ~/R/x86_64-unknown-linux-gnu-library/<R version>-bioc-devel

where `<R version>` is replaced by the x.y release version of R. On
macOS:

    ~/Library/R/<R version>-bioc-release/library
    ~/Library/R/<R Version>-bioc-devel/library

and Windows:

    C:\Users\<User>\Documents\R\win-library\<R version>-bioc-release
    C:\Users\<User>\Documents\R\win-library\<R version>-bioc-devel

(change `<User>` to your user name)

Invoke 'R for bioc-devel' or 'R for bioc-release' from the command
line on Linux:

    R_LIBS_USER=~/R/x86_64-unknown-linux-gnu-library/<R version>-bioc-release R
    R_LIBS_USER=~/R/x86_64-unknown-linux-gnu-library/<R version>-bioc-devel R

macOS:

    R_LIBS_USER=~~/Library/R/<R version>-bioc-release/library R
    R_LIBS_USER=~~/Library/R/<R version>-bioc-devel/library R

and Windows (assuming that R.exe is in PATH):

    cmd /C "set R_LIBS_USER=C:\Users\<User>\Documents\R\win-library\<R version>-bioc-release &&  R"
    cmd /C "set R_LIBS_USER=C:\Users\<User>\Documents\R\win-library\<R version>-bioc-devel &&  R"

When correctly configured, _R_'s `.libPaths()` function will return
the 'bioc-release' or bioc-devel' directory as its first
entry. Packages are installed to that directory, and that is the first
place that `library()` looks for them.  `biocLite()` and
`install.packages()` respect this setting; `update.packages()`
attempts to update packages in the directory where the current package
is installed.

Simplify invokation using an alias (Linux, macOS) or shortcut
(Windows). Add the following to your `~/.bash_profile` on Linux:

    alias Rdevel='R_LIBS_USER=~/R/x86_64-unknown-linux-gnu-library/<R version>-bioc-devel R'
    alias Rrelease='R_LIBS_USER=~/R/x86_64-unknown-linux-gnu-library/<R version>-bioc-release R'

or macOS

    alias Rdevel='R_LIBS_USER=~/Library/R/<R version>-bioc-devel/library R'
    alias Rrelease='R_LIBS_USER=~/Library/R/<R version>-bioc-release/library R'

Invoke _R_ from the command line as `Rdevel` or `Rrelease`.

On Windows, create two shortcuts, one for 'bioc-release' and one for
'bioc-devel'. For 'bioc-devel' (do similar steps for `release`) go to My
Computer and navigate to a directory that is in your PATH. Then
right-click and choose New->Shortcut.  In the "type the location of
the item" box, put:

    cmd /C "set R_LIBS_USER=C:\Users\<User>\Documents\R\win-library\<R version>-bioc-devel &&  R"

(this assumes that `R.exe` is in your PATH). Click "Next", and in the
"Type a name for this shortcut" box, type `Rdevel`.  Invoke these from
the command line as `Rdevel.lnk`.

Using 'bioc-devel' during mid-October to mid-April
--------------------------------------------------

In order to use the 'bioc-devel' version of _Bioconductor_ during the
mid-October to mid-April release cycle, you must install the devel
version of _R_.

* [Source](https://stat.ethz.ch/R/daily/)
* [macOS](http://r.research.att.com/)
* [Windows](https://cran.r-project.org/bin/windows/base/rdevel.html)

Then, make sure that your version of [BiocInstaller][] is current and
your packages up-to-date. Do this by removing all versions of
[BiocInstaller][]

    remove.packages("BiocInstaller")  # repeat until R says there is no
                                      # package 'BiocInstaller' to remove
    source("https://bioconductor.org/biocLite.R")  # install correct version
    BiocInstaller::biocValid()
