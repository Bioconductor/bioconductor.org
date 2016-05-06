# Package End of Life Policy

Creation Date: 3 March, 2015.<br />
Last Edit Date: 3 May, 2016

The Bioconductor project strives to provide a software repository that is stable
and relevant to users across the community. Each year 100-150 new software
packages are added to the repository; as of Spring 2015 over 1000 packages are
hosted. The Bioconductor build system builds and checks each software package
every 24 hours. Regular execution of example, vignette and unit test code
ensures the package is operating as expected and all required dependencies are
available. The
[build system](http://www.bioconductor.org/checkResults/devel/bioc-LATEST/)
provides a detailed report for each package across three platforms: Linux, 
Windows and Mac.

In an effort to maintain a high quality repository we have adopted a
one-year end of life (EOL) process for packages that no longer pass
build or check and do not have an active maintainer. Packages are
assessed for EOL deprecation prior to each Bioconductor release; the
EOL policies apply to software, annotation, and experiment data
packages.

## Criteria for package deprecation

1. R CMD build or check errors on one or more platforms

   The package must build and check without error on all platforms
   (exceptions to cross-platform builds are available under limited
   circumstances) at each Bioconductor release. All efforts will be
   made to keep a package in the repository if the maintainer is
   actively attempting a fix.

2. Inactive maintainer

   The maintainer listed in the DESCRIPTION file must be responsive to
   questions on the support site, package-related email from users and
   Bioconductor team members, package-related errors in the build
   system, and requests for bug fixes.

Alternatively, a package maintainer may decide that they no longer
wish to maintain their package.

## End of Life process

**Step I**: Deprecation

Packages to be deprecated will be marked with a deprecation warning and the
package name will have a strikethrough on the build report. The warning is
emitted when the package is loaded, and is reported on the package 'landing
page'. The message alerts users that the package currently fails the minimal
build and check criteria, and that the package will likely be removed from
Bioconductor in the next release.

If at any time in this 6 month period the required criteria are met (e.g., the
package returns to active maintenance, perhaps after 'adoption' by a third
party) the warning is removed.

**Step II**: Defunct

When a package has gone through one development cycle as 'deprecated' without
remedial action, the package is marked as 'Defunct'. The package is removed
from the nightly build system, is no longer available through 'biocLite()', and
does not have a current 'landing page'.

The package remains available in the subversion archive, and in previous
versions of Bioconductor.

Defunct packages cannot re-enter the Bioconductor repository except through
review as a 'new package'.

**Example**:

A package identified for deprecation during Bioconductor release 3.2 / 
devel 3.3 would make the following transitions:

    - deprecate in release 3.2
    - deprecate in devel 3.3
    ** Bioc Release **
    - deprecate in release 3.3 (i.e., no change in release status)
    - defunct in devel 3.4
    ** Bioc Release **
    - defunct in release 3.4
    - gone in devel 3.5
    ** Bioc Release **
    - gone in release 3.5


## Implementation detail

1. Notify the bioc-devel mailing list and maintainers of packages
   Depending, Importing, or Suggesting the package that the package
   will be deprecated. If appropriate, indicate that a new maintainer
   is welcome to take over.

2. Add the following code chunk to the 'devel' version of the package
   in a file `R/zzz.R`, adjusting the Bioconductor version to the
   version _after_ the current devel version.

       .onAttach <- function(libname, pkgname) {
           msg <- sprintf(
               "Package '%s' is deprecated and will be removed from Bioconductor
                version %s", pkgname, "3.4")
           .Deprecated(msg=paste(strwrap(msg, exdent=2), collapse="\n"))
       }

3. Add the following annotation to the package DESCRIPTION file.

      PackageStatus: Deprecated

4. The package remains deprecated in the 'devel' branch for up to 6
   months, after which time Bioconductor core
   team members remove the package from the 'devel' package manifest.
