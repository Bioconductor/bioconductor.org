# Troubleshooting Build Report


## How and When does the builder pull? When will my changes propagate?

Please remember the daily builder pulls, installs, builds, and checks packages
only once per day.  This process starts around 1:45 PM (13:45) EST everyday
(i.e., [UTC−05:00](https://everytimezone.com/s/b47302b3)).
Changes pushed to Bioconductor before 1:45 PM will be reflected in the
following day's build report that is posted around 12:45 PM EST.  The
build report has a timestamp at the top of its main page when it was generated.
Changes after 1:45 PM EST will not be reflected until the day after tomorrow,
therefore possibly taking up to 36-48 hours.

The build reports for
[devel](https://bioconductor.org/checkResults/devel/bioc-LATEST/) and
[release](https://bioconductor.org/checkResults/release/bioc-LATEST/) show
the package version and commit id that is being reflected for that build.
The package landing pages (example:
[Biobase](https://bioconductor.org/packages/Biobase) in release and
[Biobase](https://bioconductor.org/packages/devel/Biobase) in devel)
will not display the new version of a package until the package
installs/builds/checks without ERROR; We do not propagate broken packages.
Please also remember a package ALWAYS needs a valid version bump to
propagate to users.


## How do I reproduce the build system ERROR?

In order to reproduce the ERROR's accuately locally, remember to be using the
correct version of R and Bioconductor.  The version of R used for the build
report can be found at the top of the
[release](https://bioconductor.org/checkResults/release/bioc-LATEST/) and
[devel](https://bioconductor.org/checkResults/devel/bioc-LATEST/) build
reports. Once you are using the correct version of R make sure all your packages
are up-to-date with `BiocManager::valid()` and `BiocManager::install()`. There
are some additional environment variables the daily builder uses during `R CMD
check`. Those are the following:
```
export _R_CHECK_EXECUTABLES_=false
export _R_CHECK_EXECUTABLES_EXCLUSIONS_=false
export _R_CHECK_LENGTH_1_CONDITION_=package:_R_CHECK_PACKAGE_NAME_,abort,verbose
export _R_CHECK_LENGTH_1_LOGIC2_=package:_R_CHECK_PACKAGE_NAME_,abort,verbose
export _R_CHECK_S3_METHODS_NOT_REGISTERED_=true
```

The Single Package build has some extra
[documentation](https://github.com/Bioconductor/Contributions#r-cmd-check-environment)
about how to set up your local system to use optional environment
variables. Please note that if you look at the file listed on this page, it has
many additional variables; _Bioconductor_ does a much more stringent check on
incoming packages than on the daily builder (for now). You are welcome to use
this file if you wish as it is a more comprehensive check but the above listed
environment variables should be included minimally.

Another option to debug and test is to utilize the _Bioconductor_ docker
image. The documentation for using docker images can be found
[here](https://bioconductor.org/help/docker/). The docker image does contain the
environment variable setting found on the daily builder.


## Common Build Report Errors

Often common Error's will arise as R develops and matures or as Bioconductor
packages are modified and advance. This document provides some guidance on Error's
and potential solutions.


<a name="top"></a>

- [Bioconductor 3.11 with R 4.0](#Bioc3.11R4.0)


<a name="Bioc3.11R4.0"></a>

## Bioconductor 3.11 with R 4.0

R switched from 3.x to 4.0 which generally means some significant changes.

- [S3 method registration](#s3method)
- [Removed Settings in R CMD config](#rcmdconfig)
- [Conditional length > 1](#condLen)
- [Scalar / Vector Logic](#scalarvec)
- [Class ==  vs  is/inherits](#classEq)
  - [Matrix is now Array](#matarr)
- [data.frame stringAsFactors](#stringsAsFactors)
- [stats::smoorthEnds](#statsSmoothEnds)
- [grid package changes](#grid)
- [plot generic moved](#plot)
- [Partial Argument Matching](#partMatch)
- [Invalid UTF-8](#invalidUTF)
- [Dependency Issues](#dep311)
- [Deprecated Functions](#depFun311)

<p class="back_to_top">[ <a href="#top">Back to Bioconductor / R List </a> ]</p>

<a name="s3method"></a>

### S3 method registration

Many packages are currently failing because of undeclared S3 methods in the
NAMESPACE. There is some background information found on the R developers blog
post: [S3 Method
Lookup](https://developer.r-project.org/Blog/public/2019/08/19/s3-method-lookup/index.html)

This ERROR takes many different forms on the build report. Some of the more
common forms include

* Cannot coerce class \<structure\> to a data.frame
* Cannot coerce type ‘S4’ to vector of type ‘double’
* No applicable method for \<foo\> applied to an object of
  class \<bar\>
* ‘X’ is a list, but does not have components ‘x’ and ‘y’
* Error in colMeans(x, na.rm = TRUE) : 'x' must be numeric
* Error in RG\[1:2, \] : incorrect number of dimensions
* formal argument \<foo\> matched by multiple actual arguments
* object \<foo\> of mode 'function' was not found
* 'x' and 'y' lengths differ

<b>Solution:</b> Register the S3 method in the NAMESPACE

	S3method(<function>, <dispatch>)

A simple example which effects many packages is a S3 plotting method.
The following line would be added to the package NAMESPACE.

    S3method(plot, TCC)  # example from TCC package


<p class="back_to_top">[ <a href="#Bioc3.11R4.0">Back to Bioc 3.11 R 4.0</a> ]</p>


<a name="rcmdconfig"></a>

### Removed Setting in R CMD config

The source of the ERROR is utilizing settings in package configure script that have been removed or replaced.
There is a section of R NEWS "R CMD config no longer knows about the unused
settings F77 and FCPIFCPLAGS, nor CXX98 and similar."
Executing the configuration script when installing the package fails, and the
output contains lot of messages along the lines of the following:

* configure: WARNING: The flags FFLAGS="" do not work
* checking whether the ERROR: no information for variable 'F77'
* configure: WARNING: This value for FFLAGS does not work.

<b>Solution:</b> Replace instances of "`${R} CMD config F77`" with "`${R} CMD config FC`"


<p class="back_to_top">[ <a href="#Bioc3.11R4.0">Back to Bioc 3.11 R 4.0</a> ]</p>


<a name="condLen"></a>

### Conditional Length > 1

In R 4.0 a conditional with a length greater than 1 will produce a WARNING. On
the Bioconductor Daily Builder and Single Package Builder this is increased to
an ERROR.

Traditionally `if / while` statements could accept vectors using the first
element as the conditional value and ignoring the remaining values.  This now
produces a WARNING as seem in this dummy example and documented at
[Conditions of Length Greater Than One](https://developer.r-project.org/Blog/public/2018/10/12/conditions-of-length-greater-than-one/index.html)

    > if (c(TRUE, FALSE)) {}
    NULL
    Warning message:
    In if (c(TRUE, FALSE)) { :
        the condition has length > 1 and only the first element will be used

<b>Solution:</b>
Bioconductor increased the severity as in most cases this is a misjudgment in
the length of the argument rather than intentional.  The code should be reviewed
to see if argument is being assigned correctly.  In most cases it might be
appropriate to use an `any( )` or `all( )` surrounding the vector.

See also [mailing list post](https://stat.ethz.ch/pipermail/bioc-devel/2020-January/016081.html)

<p class="back_to_top">[ <a href="#Bioc3.11R4.0">Back to Bioc 3.11 R 4.0</a> ]</p>

<a name="scalarvec"></a>

### Scalar / Vector Logic

This is not a change in R yet but we have been notified that it is forth coming
and have escalated to an ERROR on our daily builders in preparation. This type
of ERROR occurs with the misuse of `&&` and `||`.  The double `&&` and `||`
imply a scalar comparison rather than a vector comparison that the singular `&`
and `|` expect. See the dummy example below:

    > c(TRUE, TRUE) && TRUE
    Error in c(TRUE, TRUE) && TRUE :
       'length(x) = 2 > 1' in coercion to 'logical(1)'

<b>Solution:</b>
Most cases are misjudgment and misunderstanding of the use of a scalar
comparison from a vector comparison. Changing the double `&&` / `||` to a singular
`&` / `|` will generally be sufficient if a vector comparison is intended or having the vector argument use an
appropriate `any( )` or `all( )` surrounding the vector will result in the
appropriate scalar comparison. **Note:** If this comparison is in a conditional
please see the section above; `any( )` or `all( )` will most likely be a better alternative.

See also [mailing list post](https://stat.ethz.ch/pipermail/bioc-devel/2020-January/016081.html)

<p class="back_to_top">[ <a href="#Bioc3.11R4.0">Back to Bioc 3.11 R 4.0</a> ]</p>

<a name="classEq"></a>

### Class ==  vs  is/inherits

While this isn't a change in R / Bioconductor as of yet, there is strong discussion
about the affects and consequences of this code structure.  A better discussion
and explanation can be found [When you think `class(.) == *`, think again!](https://developer.r-project.org/Blog/public/2019/11/09/when-you-think-class.-think-again/index.html)

The sum up is `class( x ) == "foo"` should be avoided. It can be misleading if classes extend other classes. The
better option is to use `is( x , "foo")` or `inherits(x, "foo")`.

This is also advised in [Bioconductor best practices](https://bioconductor.org/developers/package-guidelines/#rcode)

<a name="matarr"></a>

Starting in R 4.0, a matrix is considered an extension of array.

```
> m = matrix()
> class(m)
[1] "matrix" "array"

> class(m) == "matrix"
[1]  TRUE FALSE
> if ( class(m) == "matrix"){}
Error in if (class(m) == "matrix") { : the condition has length > 1
```
This change along with the previous section regarding conditional length results
in many errors where users were doing something along the lines of  `if
(class(m) == "matrix")`; This is an excellent example where the following is
the appropriate change `if(is(m, "matrix"))` or `if(inherits(m, "matrix"))` or
`if(is.matrix(m))`.

Another common ERROR now occurring because of this change is something similar
to the following:

```
Error in vapply(experiments(object), class, character(1)) :
  values must be length 1,
 but FUN(X[[4]]) result is length 2
```


<p class="back_to_top">[ <a href="#Bioc3.11R4.0">Back to Bioc 3.11 R 4.0</a> ]</p>

<a name="stringsAsFactors"></a>

### data.frame stringsAsFactors

In R 4.0, the default for data.frame argument `stringsAsFactors` changed from
TRUE to FALSE.  This changes is causing the most breakage in tests where there
are checks for particular factor levels or constructing factor levels.
The ERROR’s take many different forms. The simple solution is to change or add
the `stringAsFactors=TRUE` to the data.frame call,  however maintainers may want
to re-evaluate code for potential restructuring or ease of use.


<p class="back_to_top">[ <a href="#Bioc3.11R4.0">Back to Bioc 3.11 R 4.0</a> ]</p>

<a name="statsSmoothEnds"></a>

### stats::smoothEnds

A recent change to `stats::smoothEnds()`, now returns an integer vector with the
input is an integer vector. Previously it could return a number vector.

Example R 3.6.3
```
> class(smoothEnds(c(401:403)))
[1] "integer"
> class(smoothEnds(c(401:403, 555L)))
[1] "numeric"
```

Example from 4.0.0

```
> class(smoothEnds(c(401:403)))
[1] "integer"
> class(smoothEnds(c(401:403, 555L)))
[1] "integer"
```
This has the potential to cause ERROR's if a class type was checked.

<p class="back_to_top">[ <a href="#Bioc3.11R4.0">Back to Bioc 3.11 R 4.0</a> ]</p>

<a name="grid"></a>

### Grid package changes

We do not have a lot of specifics on what has changed but were notified by
email.  Important sections of the email as follows:

```
  	I am about to commit some internal changes to 'grid' units
        (for, in some cases, 100x speed-up of unit operations).
    	A number of packages have already been fixed to work with these
        changes, but, according to my testing, the following CRAN
        packages will still fail R CMD check.
	Some of those are cascades ('armada', 'countToFPKM', and 'wilson'
        from 'ComplexHeatmap' - see below - and 'fingertipscharts' from 'lemon'
        and 'xpose' is actually a 'ggforce' problem), but all of the other
        package authors have been notified and several are already working on
        fixes.
    	The most serious of those is 'ComplexHeatmap' because it causes multiple
        follow-on failures, the CRAN ones above and others on BioConductor:
	Again, the main package authors have been notified and the
  	'ComplexHeatmap' author is working on an update.
```

<p class="back_to_top">[ <a href="#Bioc3.11R4.0">Back to Bioc 3.11 R 4.0</a> ]</p>

<a name="plot"></a>

### plot generic moved

The plot generic has moved from graphics to base.  The ERROR's seen from this
change are non specific and can take many forms.  Some of the ERROR's we see
are

```
Error in getGeneric(f, TRUE, envir, package) :
 no generic function found for 'plot'
```
or
```
Error in as.double(y) :
      cannot coerce type 'S4' to vector of type 'double'

```

The explanation given:

      “The namespace controls the search strategy for variables used by
      functions in the package. If not found locally, R searches the package
      namespace first, then the imports, then the base namespace and then  the
      normal search path." as per
      https://cran.r-project.org/doc/manuals/r-devel/R-exts.html#Package-namespaces):

      CRAN and Bioconductor  had a few packages that "worked" because the right
      plot() was found in the normal search path, but now fail because it's
      calling the one in base instead.



<p class="back_to_top">[ <a href="#Bioc3.11R4.0">Back to Bioc 3.11 R 4.0</a> ]</p>

<a name="partMatch"></a>

### Partial Argument Matching

There is now more strict checking of argument matching with regards to partial
argument matching. Best described with the following example

    setGeneric(“mycoolfunction”,   function(object,  breaks)
	    standardGeneric(“mycoolfunction”)
    setMethod(“mycoolfunction”,
        signature=c(object=”GRanges”, break=”GRanges”),
        <code>)

Notice the generic uses `breaks` while the setMethod uses `break`; This is an
example of a partial argument match that will no longer be valid.

Partial argument matching when envoking functions should also be avoided as
part of best practices. For example

    mycoolfunction <- function( x, myargum, secondarg ) { code }

    mycoolfunction(x=2, myar=1:2, second=3)          # BAD Coding!

    mycoolfunction(x=2, myargum=1:2, secondarg=3)    # Good Practice!


<p class="back_to_top">[ <a href="#Bioc3.11R4.0">Back to Bioc 3.11 R 4.0</a> ]</p>

<a name="invalidUTF"></a>

### Package inputenc Error: Invalid UTF-8

This ERROR started to appear on tokay2 (windows) in Spring 2020. We are not sure
the exact source of the ERROR (change in MiKTek, Change in R, other?) but the
solution is simple:

Please place `\usepackage[utf8]{inputenc}` in the beginning of your Sweave
vignette right after the `\documentclass` line.

<p class="back_to_top">[ <a href="#Bioc3.11R4.0">Back to Bioc 3.11 R 4.0</a> ]</p>

<a name="dep311"></a>

### Dependency Issues

Dependency Issues can fall into a few sub-categories:

#### CRAN binaries not available for R 4.0

The fall cycle of Bioconductor uses R-devel in preparation for the new release of
R in the spring. This is always a slightly more disruptive cycle with regards to
package dependencies from CRAN. CRAN over the next 6 months leading up to the new release of R
will make binaries for Windows and MacOS available. As they become available the
Bioconductor builders will automatically add these binaries. If the binaries have
not been created yet, they will be unavailable and result in a `package not available
error`.  Bioconductor will not go to extra efforts to find work around to install
these packages; when they are available, they will be added.
<b>Solution:</b> Please be patient!

#### Package have been removed from CRAN

CRAN packages are occasionally removed. Unfortunately, Bioconductor will only
allow package dependencies to be actively maintained packages on CRAN or Bioconductor.
A package will have to alter their package to not utilize code and not rely on this
dependency. You may of course try to pentition CRAN for reinstatement or reach
out to the package maintainer to fix and submit to CRAN. Good Luck!

#### Package has been removed from Bioconductor

We try to be more aware of orphaned packages and packages that remain broken for
extended periods of time. Package deprecation and removal occurs and packages
will have to alter to not utilize code from these packages or could potential
offer to take over maintenance of broken packages but that would require original
maintainers permission. Bioconductor Package deprecation is announced throughout
the release cycle on the mailing list and support site to try and allow dependent
packages time to adjust code before removal. This release the most
notable maintainer requested deprecation from 3.10 (therefore removed in 3.11) are
SNPchip and GenomeGraphs. A full list of deprecated packages can be found
[List of Deprecated Packages 3.10](https://support.bioconductor.org/p/125352/).
We also documented removed packages on our [Removed Package Page](/about/removed-packages/)


<p class="back_to_top">[ <a href="#Bioc3.11R4.0">Back to Bioc 3.11 R 4.0</a> ]</p>


<a name="depFun311"></a>

### Deprecated Functions

Functions can be deprecated, defunct, and eventually removed.  Bioconductor tries
to enforce this progression to allow maintainers to adjust code. Most deprecated
or defunct functions will (should) suggest the alternative. The following are
noted in Bioconductor 3.11

#### RangedData

    Error : RangedData objects are defunct. Please migrate your code to use GRanges
      or GRangesList objects instead. See IMPORTANT NOTE in ?RangedData

#### Normalize

    Error:
    'normalize' is defunct.
    Use ''normalize,SingleCellExperiment-method' is defunct.
    Use 'logNormCounts' instead' instead.


#### calculateQCMetrics

    Error:
    'calculateQCMetrics' is defunct.
    Use 'perCellQCMetrics' instead.


<p class="back_to_top">[ <a href="#Bioc3.11R4.0">Back to Bioc 3.11 R 4.0</a> ]</p>
