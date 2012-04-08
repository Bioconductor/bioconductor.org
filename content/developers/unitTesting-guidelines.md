![](/images/icons/magnifier.gif)Unit Testing Guidelines
==================================================


* [Introduction](#introduction)
* [Motivation](#motivation)
* [Adding Tests to Your Code](#addingTests)
* [Conventions for the Build Process](#conventions)
* [Using Tests During Development](#duringDevelopment)
* [File Summary](#fileSummary)
* [Additional Resources](#resources)

<h2 id="introduction">Introduction</h2>

Unit tests are simple to write, easily invoked, and confer large benefits throughout the software development process,
from early stage exploratory code, to late stage maintenance of a long-established project.  Unit testing often  becomes
indispensable to those who give it a try. Here we explain how to write unit tests, 
how to run them, and how they are woven into the standard Biocondcutor build process.  We hope that unit tests will become
a standard part of your software development, and an integral part of your Bioconductor package.
<p><br>

We use and recommend the <a href="http://cran.r-project.org/web/packages/RUnit/index.html">RUnit</a> package from CRAN to write unit tests &mdash;
an <b><i>R</i></b> implementation of the <a href="http://en.wikipedia.org/wiki/Agile_software_development">agile</a>
software development 'XUnit' </i></b> tools  (see also <a href="http://www.junit.org">JUnit</a>, <a href="http://pyunit.sourceforge.net">PyUnit</a>)
each of which tries to encourage, in their respective language,  the rapid development of robust useful software.

(<p>An alternative approach, used by some Bioconductor packages, is Hadley Wickham's 
<a href="http://cran.r-project.org/web/packages/testthat/index.html">'testthat'</a> package.)

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="motivation">Motivation</h2>


Why bother with unit testing?  
<br><br>
Imagine that you need a function <code>divideBy</code> taking two arguments, which you might define like this:

<code><pre>
divideBy <- function(dividend, divisor) { 
    if(divisor==0) 
       return(NA)
    return(dividend/divisor) 
}
</pre></code>

As you develop this function you would very likely test it out in a variety of ways, using different arguments, checking
the results, until eventually you are satisfied that it performs properly.  Unless you adopt some sort of software
testing protocol, however, your tests are unlikely to become an integral part of your code.  They may be scattered
across different files, or they may not exist as re-runnable code in a file at all, just as ad hoc command-line function
calls you sometimes remember to make.

<br><br>

A far better approach, we propose, is to use <b>lightweight, formalized</b> unit testing.  This requires only a very few conventions and practices:
<br><br>
<ul>
  <li> Store the test functions in a standard directory.
  <li> Use simple functions from the RUnit package to check your results.
  <li> Run the tests as a routine part of your development process.
</ul>


<p>
<br>
Here is such a unit test for <code>divideBy</code>:

<code><pre>
test_divideBy <- function() {
    checkEquals(divideBy(4, 2), 2)
    checkTrue(is.na(divideBy(4, 0)))
    checkEqualsNumeric(divideBy(4, 1.2345), 3.24, tolerance=1.0e-4)
}
</pre></code>

Adopting these practices will cost you very little.  Most developers find that these practices simplify and shorten
development time.  In addition, they create an <i><b> executable contract</b></i> &mdash; a concise and verifiable
description of what your code is supposed to do. The experienced unit-testing programmer will create such a test
function to accompany every function, method and class they write.  (But don't let this scare you off.  Even adding a single test 
to your package is worthwhile, for reasons explained below.)
<p>
<br>
Developers often rebel when unit tests are recommended to them, calculating that creating unit tests
for existing code would be a lengthy and tedious job, and that their productivity will suffer.

<p>
Unit tests, however, are best written <b><i>AS</i></b> you develop code, rather than after your package is written.  Replace your informal testing
with a few lightweight formal practices, and you will see both your immediate and long-term productivity increase.
<p>

Consider that every unit of software (every function,  method, or class) is designed to do a job, to return specific
outputs for specific inputs, or to cause some specific side effects.  A unit test specifies these behaviors, and
provides a single mechanism &mdash; one or more test functions residing in one or more files, within a standard directory structure &mdash; to
ensure that the target function, method or class does its job.  With that assurance, the programmer (and her
collaborators) can then, with confidence, proceed to use it in a larger program.  When a bug appears, or new features are
needed and added, one adds new tests to the existing collection.  Your code becomes progressively more powerful, more
robust, and yet remains easily and automatically validated.

<p>
Some proponents suggest that the benefits of unit testing extend further:  that code design itself improves.   They argue that the
operational definition of a function through its tests encourages clean design, the 'separation of concerns', and sensible
handling of edge cases.
<p>

<code>test_dividesBy</code> illustrates all the crucial features of a good test function.  It uses the
simple RUnit <b><i>check</i></b> functions.  It makes sure that reasonable values are returned across a range
of normal and pathological conditions.  Its name begins with <code>test_</code> so that it is recognized and run by the
Bioc build process.  It would reside (more about this below) in the <code>inst/unitTests</code> directory.

<p>

Finally, unit testing can be <b><i>adopted piecemeal</i></b>.   Add a single test to your package, even if only a test for a minor feature, and 
both you and your users will benefit.  Add more tests as you go, as bugs arise, as new features are added, when you find yourself puzzling
over code your wrote some months before.  Soon, unit testing will be part of your standard practice, and your package will have
an increasingly complete set of tests.


<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="addingTests">Adding Tests For Your Code</h2>

Three things are required:
<br>
<ol>
   <li> Create a file containing functions in the style of  <code>test_dividesBy</code> for each function you want to test, using RUnit-provided check functions.
   <li> Add a few small (and idiosyncratic) files in other directories.
   <li> Make sure the RUnit and BiocGenerics packages are available.  
</ol>
Steps two and three are explained in <a href=#conventions>conventions for the build process</a>, below. 

These are the  <b>RUnit</b> check methods:

<pre><code>checkEquals(expression-A, expression-B)
checkTrue(condition)
checkEqualsNumeric(a, b, tolerance)
</code></pre>

Less commonly used, but recommended if your code throws exceptions:

<pre><code>checkException(expr, msg)
</code></pre>

In a typical test function, as you can see in <code>test_divideBy</code>, you invoke one of your program's functions or methods, then call an 
appropriate RUNit check function to make sure that the result is correct.   RUnit reports failures, if there are any, with enough context to track down the error.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="conventions">Conventions for the Build Process</h2>

Though writing unit tests is easy, and though they make software development much easier than doing without, it can be
quite confusing to set up your Bioconductor package properly so that 
<pre><code>R CMD check <b>YOURPACKAGE</b></code></pre> finds and run your tests.  This can be confusing, so we
take some pains here to describe exactly how things should be set up, and what is going on behind the scenes.  (See the
<a href='#duringDeveloment'>next section</a> for the much simpler technique to use when you want to test only a small part of your code.)
<br><br>
<ul>

 <li> The standard command <code>R CMD check MyPackage</code> sources and runs all R files found in your MyPackage/tests/  directory.

 <li> Historically, and sometimes still, R package developers placed test code of their own invention and style into one or more files in this 'tests' directory.

 <li> RUnit support was added to this already-existing structure and practice about 2005, and the additions can be confusing,
    beginning with the indirect way in which your test functions are found and executed. (But follow these steps and all
    should be well.  Post to bioc-devel if you run into any difficulty.) 

    There are three steps:<br><br>

   </ul>  <ol>

     <li> Create the file  'MyPackage/tests/runTests.R' with these contents:

   require("MyPackage") || stop("unable to load MyPackage")
   BiocGenerics:::testPackage("MyPackage")

     <li> Create any number of files in MyPackage/inst/unitTests/ for your unit test functions.  You can put your tests
        all in one file in that directory, or distributed  among multiple files.  In either case, all files must follow the nameing convention specified in this regular expression:
  <code><pre>pattern="^test_.*\\.R$"</pre></code>

Thus, for our example, a good choice would be:

  <code><pre>MyPackage/inst/unitTests/test_divideBy.R</pre></code>

   </ol>


<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="duringDeveloment">Using Tests During Development</h2>

<pre><code>R CMD check <b>MyPackage</b></code></pre> will run all of your tests.  But when developing a class, or
debugging a method or function, you will probably want to run just one test at a time, and to do so before the package
is complete.  Assuming you have followed the directory structure and naming conventions recommended above, here is what
you would do:
<pre><code>library (RUnit)
library (MyPackage)

source ('MyPackage/inst/unitTests/myTests.R')
test_divideBy ()
 [1] TRUE
</code></pre>
A failed test is reported like this:

<pre><code> Error in checkEquals(divideBy(4, 2), 3) : Mean relative difference: 0.5
</code></pre>

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="fileSummary">Summary: the minimal setup</h2>
A minimal Bioconductor <b>unitTest</b> setup requires two additions to the <code>MyPackage/DESCRIPTION</code> file, then three specific files, each residing in a 
different directory.  (You can use other names for the three <b><i>R</i></b> source files shown below if you wish.)

<h4>MyPackage/DESCRIPTION</h4>
<pre><code>Imports: BiocGenerics 
Suggests: RUnit</code></pre>



<h4> MyPackage/tests/runTests.R </h4>
<pre><code>require("MyPackage") || stop("unable to load MyPackage")
MyPackage:::.test()</code></pre>

<h4> MyPackage/R/testPackage.R </h4>
<pre><code>.test <- function() BiocGenerics:::testPackage("MyPackage")</code></pre>

<h4> MyPackage/inst/unitTests/myTests.R</h4>
<pre><code>test_divideBy <- function () {
    checkEquals(divideBy (4, 2), 2)
    checkTrue(is.na (divideBy (4, 0)))
    checkEqualsNumeric(divideBy (4, 1.2345), 3.24, tolerance=1.0e-4)
}
</code></pre>


<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="resources">Additional Resources</h2>

Some web resources worth reading:
<br><br>
<ul>
  <li> <a href="http://en.wikipedia.org/wiki/Unit_testing">An Overview</a>
  <li> <a href="http://www.daedtech.com/addicted-to-unit-testing">An informal account</a>
  <li> <a href="http://en.wikipedia.org/wiki/Test-driven_development">Test-driven development</a>
  <li> <a href="http://en.wikipedia.org/wiki/Agile_software_development">Agile software development</a>
</ul>

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>
