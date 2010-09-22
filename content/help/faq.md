# ![](/images/icons/magnifier.gif)FAQ

This is the FAQ, currently under construction.

* [Package Installation](#install-packages)
* [Developing Packages](#developer-faq)

<h2 id="install-packages">Package Installation</h2>

* Package XXX fails to install

A common reason for a package to fail to install is that `R` or
`Bioconductor` software dependencies are not satisfied, as shown here
from a user trying to install the `affyPLM` package:

    ...
    ** inst
    ** preparing package for lazy loading
    Error: package 'affy' required by 'affyPLM' could not be found
    Execution halted
    ERROR: lazy loading failed for package 'affyPLM'

Be sure to use `biocLite` to [install packages][2] that are
appropriate for your system and version of `R`.  Be sure that your
installed packages are up-to-date by following [update packages][1].

Less commonly, packages may install but then fail to load, as here
with the `Rsamtools` package:

    Error in dyn.load(file, DLLpath = DLLpath, ...) :
    unable to load shared library 
    '/usr/local/lib64/R/library/Rsamtools/libs/Rsamtools.so':
      /usr/local/lib64/R/library/Rsamtools/libs/Rsamtools.so: undefined symbol: ecrc32

This is likely a system configuration issue, e.g., a Windows `PATH`
environment variable is not set correctly, or the Linux `ldconfig`
program or `LD_LIBRARY_PATH` environment variable is incorrect.

Packages may also fail to install because third party software is not
available. This often happens during the `configure` part of the
package installation process, as illustrated here with the `XML`
package:

    * Installing *source* package 'XML' ...
    ...
    checking for xml2-config... no
    Cannot find xml2-config
    ERROR: configuration failed for package 'XML'

These types of errors can sometimes be easily solved (installing
necessary libraries or other software, perhaps referenced on the
[package home page][4]). It will often be necessary to understand your
system more thoroughly than you'd like, perhaps with the assistance of
the Bioconductor [mailing list][3].

<h2 id="developer-faq">Developing Packages</h2>

* What packages belong in the Depends:, Imports:, or Suggests: fields?

Two relevant mailing list posts
([a](https://stat.ethz.ch/pipermail/r-devel/2008-December/051602.html),
[b](https://stat.ethz.ch/pipermail/bioc-devel/2010-September/002310.html))
address this. Generally, packages whose code you use in your own
package should where ever possible be Import:'ed. Packages required
for vignettes are often Suggest:'ed. Depends: is appropriate for
packages that cannot be Import:'ed (e.g., because they do not have a
NAMESPACE) or for packages that provide essential functionality needed
by the user of your package, e.g., your functions always return
`GRanges` objects, so the user will necessarily need `GenomicRanges`
on their search path.

[1]: /install/index.html#update-bioconductor-packages
[2]: /install/index.html#install-bioconductor-packages
[3]: /help/mailing-list/
[4]: /packages/release/bioc/
