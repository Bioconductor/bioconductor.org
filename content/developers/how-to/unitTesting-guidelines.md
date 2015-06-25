![](/images/icons/magnifier.gif)Unit Testing Guidelines
==================================================


* [Introduction](#introduction)
* [Motivation](#motivation)
* [Deciding Which Test Framework To Use](#choosingTestFramework)
* [RUnit Usage](#RUnitUsage)
  * [Adding Tests to Your Code](#RUnitAddingTests)
  * [Conventions for the Build Process](#RUnitConventions)
  * [Using Tests During Development](#RUnitDuringDevelopment)
  * [File Summary](#RUnitFileSummary)
* [Testthat Usage](#testthatUsage)
  * [Conversion from RUnit to testthat](#testthatConversion)
* [Test Coverage](#coverage)
* [Additional Resources](#resources)

<h2 id="introduction">Introduction</h2>

Unit tests are simple to write, easily invoked, and confer large
benefits throughout the software development process, from early stage
exploratory code, to late stage maintenance of a long-established
project.  Unit testing often becomes indispensable to those who give
it a try. Here we explain how to write unit tests, how to run them,
and how they are woven into the standard Bioconductor build process.
We hope that unit tests will become a standard part of your software
development, and an integral part of your Bioconductor package.

We recommend either the [RUnit] or [testthat] packages from CRAN to write unit
tests. RUnit is an _R_ implementation of the [agile] software development 'XUnit'
tools (see also [JUnit], [PyUnit]) each of which tries to encourage, in their
respective language, the rapid development of robust useful
software.  Testthat also draws inspiration from the xUnit family of testing
packages, as well as from many of the innovative ruby testing libraries, like
[rspec](https://rspec.info/), [testy](https://github.com/ahoward/testy),
[bacon](https://github.com/chneukirchen/bacon) and
[cucumber](https://cucumber.io).

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="motivation">Motivation</h2>

Why bother with unit testing?

Imagine that you need a function `divideBy` taking two arguments,
which you might define like this:

    divideBy <- function(dividend, divisor) {
        if (divisor == 0)
           return(NA)
        dividend / divisor
    }

As you develop this function you would very likely test it out in a
variety of ways, using different arguments, checking the results,
until eventually you are satisfied that it performs properly.  Unless
you adopt some sort of software testing protocol, however, your tests
are unlikely to become an integral part of your code.  They may be
scattered across different files, or they may not exist as re-runnable
code in a file at all, just as ad hoc command-line function calls you
sometimes remember to make.

A far better approach, we propose, is to use **lightweight,
formalized** unit testing.  This requires only a very few conventions
and practices:

* Store the test functions in a standard directory.
* Use simple functions from the *RUnit* or *testthat* packages to check your results.
* Run the tests as a routine part of your development process.

Here is a RUnit test for `divideBy`:

    test_divideBy <- function() {
        checkEquals(divideBy(4, 2), 2)
        checkTrue(is.na(divideBy(4, 0)))
        checkEqualsNumeric(divideBy(4, 1.2345), 3.24, tolerance=1.0e-4)
    }

And the equivalent test suing testthat:

    test_that("divideBy works properly", {
      expect_equal(divideBy(4, 2), 2)
      expect_true(is.na(divideBy(4, 0)))
      expect_equal(divideBy(4, 1.2345), 3.24, tolerance = 1.0e-4)
    })

Adopting these practices will cost you very little.  Most developers
find that these practices simplify and shorten development time.  In
addition, they create an **executable contract** &mdash; a concise and
verifiable description of what your code is supposed to do. The
experienced unit-testing programmer will create such a test function
to accompany every function, method and class they write.  (But don't
let this scare you off.  Even adding a single test to your package is
worthwhile, for reasons explained below.)

Developers often rebel when unit tests are recommended to them,
calculating that creating unit tests for existing code would be a
lengthy and tedious job, and that their productivity will suffer.

Unit tests, however, are best written **as you develop** code, rather
than after your package is written.  Replace your informal testing
with a few lightweight formal practices, and you will see both your
immediate and long-term productivity increase.

Consider that every unit of software (every function, method, or
class) is designed to do a job, to return specific outputs for
specific inputs, or to cause some specific side effects.  A unit test
specifies these behaviors, and provides a single mechanism &mdash; one
or more test functions residing in one or more files, within a
standard directory structure &mdash; to ensure that the target
function, method or class does its job.  With that assurance, the
programmer (and their collaborators) can then, with confidence, proceed
to use it in a larger program.  When a bug appears, or new features
are needed and added, one adds new tests to the existing collection.
Your code becomes progressively more powerful, more robust, and yet
remains easily and automatically validated.

Some proponents suggest that the benefits of unit testing extend
further: that code design itself improves.  They argue that the
operational definition of a function through its tests encourages
clean design, the 'separation of concerns', and sensible handling of
edge cases.

Finally, unit testing can be **adopted piecemeal**.  Add a single test
to your package, even if only a test for a minor feature, and both you
and your users will benefit.  Add more tests as you go, as bugs arise,
as new features are added, when you find yourself puzzling over code
your wrote some months before.  Soon, unit testing will be part of
your standard practice, and your package will have an increasingly
complete set of tests.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="choosingTestFramework">Deciding Which Test Framework To Use</h2>

RUnit and testthat are both robust testing solutions that are great tools for
package development, which you choose to use for your package largely comes
down to personal preference.  However here is a brief list of strengths and
weaknesses of each.

### RUnit Strengths ###
- Longer history (first release 2005)
- Direct analog to other xUnit projects in other languages.
- Only need to learn a small set of check functions.
- Used extensively in Bioconductor (210 Bioconductor packages, overall 339 circa May 2015), particularly in
  the core packages.

### RUnit Weaknesses ###
- No RUnit development activity since 2010, and has no active maintainer.
- Need to manually source package and test code to run interactively.
- More difficult to setup and run natively (although see
  `BiocGenerics:::testPackage()` below which handles some of this).

### Testthat Strengths ###
- Active development with over 39 contributors.
- Greater variety of test functions available, including partial matching and
  catching errors, warnings and messages.
- Easy to setup with `devtools::use_testthat()`.
- Integrates with `devtools::test()` to automatically reload package source and
  run tests during development.
- Test failures and errors are more informative than RUnit.
- A number of different reporting functions available, including visual
  real-time test results.
- Used extensively in CRAN (546 CRAN packages, overall 598 circa May 2015).

### Testthat Weaknesses ###
- Test code is slightly more verbose than the equivalent RUnit tests.
- Has been available for less time (only since 2009).

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="RUnitUsage">RUnit Usage</h2>

<h3 id="RUnitAddingTests">Adding Tests For Your Code</h3>

Three things are required:

1. Create a file containing functions in the style of `test_dividesBy`
   for each function you want to test, using *RUnit*-provided check
   functions.
2. Add a few small (and idiosyncratic) files in other directories.
3. Make sure the *RUnit* and [BiocGenerics] packages are available.

Steps two and three are explained in [conventions for the build
process](#conventions).

These are the *RUnit* check methods:

    checkEquals(expression-A, expression-B)
    checkTrue(condition)
    checkEqualsNumeric(a, b, tolerance)

In a typical test function, as you can see in `test_divideBy`, you
invoke one of your program's functions or methods, then call an
appropriate *RUnit* check function to make sure that the result is
correct.  *RUnit* reports failures, if there are any, with enough
context to track down the error.

*RUnit* can test that an exception (error) occurs with

    checkException(expr, msg)

but it is often convenient to test specific exceptions, e.g., that a
warning "unusual condition" is generated in the function `f <- function()
{ warning("unusual condition"); 1 }` with

    obs <- tryCatch(f(), warning=conditionMessage)
    checkIdentical("unusual condition", obs)

use `error=...` to test for specific errors.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h3 id="RUnitConventions">Conventions for the Build Process</h3>

Writing unit tests is easy, though your Bioconductor package must be
set up properly so that `R CMD check MyPackage` finds and run your
tests.  We take some pains to describe exactly how things should be
set up, and what is going on behind the scenes.  (See the [next
section](#RUnitDuringDeveloment) for the simple technique to use when you
want to test only a small part of your code).

The standard command `R CMD check MyPackage` sources and runs all R
files found in your `MyPackage/tests/` directory.  Historically, and
sometimes still, *R* package developers place test code of their own
invention and style into one or more files in this `tests` directory.

*RUnit* was added to this already-existing structure and practice
about 2005, and the additions can be confusing, beginning with the
indirect way in which your test functions are found and executed. (But
follow these steps and all should be well.  Post to [bioc-devel] if
you run into any difficulty.)

There are two steps:

1. Create the file `MyPackage/tests/runTests.R` with these contents:

       BiocGenerics:::testPackage("MyPackage")

2. Create any number of files in `MyPackage/inst/unitTests/` for your
   unit test functions.  You can put your tests all in one file in
   that directory, or distributed among multiple files.  All files
   must follow the naming convention specified in this regular
   expression:

       pattern="^test_.*\\.R$"

   For our example, therefore, a good choice would be
   `MyPackage/inst/unitTests/test_divideBy.R` or if the `dividesBy`
   function was one of several home-brewed arithmetic functions you
   wrote, and for which you provide tests, a more descriptive filename
   (a practice we always recommend) might be
   `MyPackage/inst/unitTests/test_homeBrewArithmetic.R`

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h3 id="RUnitDuringDevelopment">Using Tests During Development</h3>

    R CMD check MyPackage

will run all of your tests.  But when developing a class, or debugging
a method or function, you will probably want to run just one test at a
time, and to do so when an earlier version of the package is
installed, against which you are making local exploratory
changes. Assuming you have followed the directory structure and naming
conventions recommended above, that your current working directory is
inst, here is what you would do:

    library(RUnit)
    library(MyPackage)

    source('../R/divideBy.R')
    source('unitTests/test_divideBy.R')
    test_divideBy()
    [1] TRUE

A failed test is reported like this:

    Error in checkEquals(divideBy(4, 2), 3) : Mean relative difference: 0.5

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h3 id="RUnitFileSummary">Summary: the minimal setup</h3>

A minimal Bioconductor **unitTest** setup requires only this one-line addition to
the `MyPackage/DESCRIPTION` file

    Suggests: RUnit, BiocGenerics

and two files, `MyPackage/tests/runTests.R`:

    BiocGenerics:::testPackage("MyPackage")

and `MyPackage/inst/unitTests/test_divideBy.R`:

    test_divideBy <- function() {
        checkEquals(divideBy(4, 2), 2)
        checkTrue(is.na(divideBy(4, 0)))
        checkEqualsNumeric(divideBy(4, 1.2345), 3.24, tolerance=1.0e-4)
    }

Remember that your `unitTests/test_XXXX.R` file, or files, can have any
name(s), as long as they start with `test_`.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="testthatUsage">Testthat Usage</h2>
Hadley Wickham, the primary author of testthat has a comprehensive chapter on
[Testing with testthat] in his R packages book.  There is also an article
[testthat: Get Started with Testing] in the R-Journal.

The easiest way to setup the testthat infrastructure for a package is using
`devtools::use_testthat()`.

You can then automatically reload your code and tests and re-run them using
`devtools::test()`.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h3 id="testthatConversion">Conversion from RUnit to testthat</h3>

If you have an existing RUnit project you would like to convert to using
testthat you will need to change the following things in your package
structure.

1. `devtools::use_testthat()` can be used to setup the testthat testing structure.
2. Test files are stored in `tests/testthat` rather than `inst/unitTests` and
   should start with `test`.  Richard Cotton's
   [runittotesthat](https://github.com/richierocks/runittotestthat]) package
   can be used to programmatically convert RUnit tests to testthat format.
3. You need to add `Suggests: testthat` to your `DESCRIPTION` file rather than
   `Suggests: RUnit, BiocGenerics`.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="coverage">Test Coverage</h2>

[Test coverage](https://en.wikipedia.org/wiki/Code_coverage)
refers to the percentage of your package code
that is tested by your unit tests. Packages with higher coverage
have a lower chance of containing bugs. 

Our build system calculates test coverage for every software
package (using the [covr](https://github.com/jimhester/covr)
package) and reports the results in a "test coverage" shield
on the package landing page. The value shown in this shield is
either a percentage (a number from 0 to 100) or 'unknown', which
could mean:

* The package has no unit tests, or the unit tests are not
  properly configured. Read this page from
  [the beginning](#top) to learn how to add them.
* The unit tests in the package failed.
* There was a problem calculating the test coverage.
  Inquire on the [bioc-devel list](/help/support/#bioc-devel)
  if you have questions about this.


<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>


<h2 id="resources">Additional Resources</h2>

Some web resources worth reading:

* [An Overview](http://en.wikipedia.org/wiki/Unit_testing)
* [An informal account](http://www.daedtech.com/addicted-to-unit-testing)
* [Test-driven development][tdd]
* [Agile software development][agile]
* [Testing with testthat]
* [testthat: Get Started with Testing]

[BiocGenerics]: /packages/release/bioc/html/BiocGenerics.html
[RUnit]: http://cran.r-project.org/web/packages/RUnit/index.html
[agile]: http://en.wikipedia.org/wiki/Agile_software_development
[JUnit]: http://www.junit.org
[PyUnit]: http://pyunit.sourceforge.net
[tdd]: http://en.wikipedia.org/wiki/Test-driven_development
[bioc-devel]: /help/mailing-list/#bioc-devel
[testthat]: http://cran.r-project.org/web/packages/testthat/index.html
[Testing with testthat]: http://r-pkgs.had.co.nz/tests.html
[testthat: Get Started with Testing]: http://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf
<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>
