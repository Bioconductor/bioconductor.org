<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{A primer for building packages}
-->

<!-- TODO: -->


![](/images/icons/help.gif)Building packages for Bioconductor
=======================================================================

The aim of this guide is to be a primer for new package authors who aspire to write packages for Bioconductor.  It will take you through the steps necessary to do this using tools available in the widely available IDE RStudio.  And it will also attempt to explain the most fundamental aspects of package anatomy with the aim of helping your package to eventually meet or exceed the [Bioconductor Package Guidelines](/developers/package-guidelines/).  We also have other resources that can help developers with various aspects of package development.  For those not interested in this very introductory guide, there are other developer resources on the Bioconductor website [here](/developers/).  And there are also other [how to documents](/developers/how-to/) that can be explored to answer various questions about developing for Bioconductor.

Bioconductor packages are really just R packages.  So the final word on what
can be done with them can be found in this document
[here](http://cran.fhcrc.org/doc/manuals/R-exts.html).

But because that document is exhaustively detailed, it can be difficult
for a beginner to know what they need to get started with writing R
packages for Bioconductor.  This document will attempt to
simplify the discussion by focusing on the bare minimum of what is
needed to make R code run as a Bioconductor package.

* [Creating a package using RStudio](#creating-packages-rstudio)  
* [Basic package Anatomy](#basic-package-anatomy)  
* [The DESCRIPTION file](#description)  
* [The 'R' directory](#the-r-dir)
* [The NAMESPACE](#namespace)  
* [Documenting manual pages](#manual-pages)  
* [Writing vignettes](#vignettes)  
* [Storing data in your package](#external-data-dirs)  
* [BiocCheck](#bioccheck)  
* [Tips for writing good packages](#good-practices)  

<!--  -->



<h2 id="package-type">Package Types</h2> 

The Bioconductor project recognizes that there are many different
general uses for R packages. Speaking very broadly, some packages are
used primarily for holding annotations or examples of data sets while
other packages are primarily for distributing software.  For the
purposes of this document, we are focusing on software packages since
they are by far the most complex to create and maintain.  In other
words, these are packages that wrap some functionality up so that
other R users can also easily reproduce and build on your work.




<h2 id="creating-packages-rstudio">Creating a package using RStudio</h2> 

Today a lot of people use RStudio, and this is a great place to get
started because RStudio has already added tools to make package
creation easier for end users.  But before we get started lets configure a couple of things so that your code will be formatted the way that we prefer it for Bioconductor.  

* Set the column width marker to 80 columns.  You can find this in the 'Tools' menu if you select the 'Global Options...' and then look at the 'Code Editing' panel.  Then make sure that you click the 'show margin' option and that it is set to 80 columns.  This will help you see if your lines are too long.

* Second set up the tab to be 4 spaces.  You can do this right above where you set the column width marker for 80.

To make your first package with RStudio, go to the
'File' menu and choose 'New Project'.  From here you probably want to
choose 'New Directory' and then you definitely want o choose 'R
Package' to set up a new R package.  This will take you to the final
step to get started which will just ask you what you want to call your
new package and where to put it.  It is worth taking a moment to
choose a name that you think you will like and to put the package in a
location where you are comfortable working.  These things can be
easily changed later but it's still best to get them correct right
now.

For this example lets call our package 'MyPkg'.  Now go ahead an fill
this information in and push the 'Create Project' button.

If you look at the screen now you will see several files listed in the
files tab.  These are files that are needed to make a package (and
which RStudio conveniently creates for you).  A couple of these files
are 'extras' that are put there by RStudio to help it manage your
package building ('.Rbuildignore', 'Read-and-delete-me' and
'MyPkg.Rproj') You can just ignore these files, but I wanted you to
understand that they are not part of 'required' R package anatomy.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>


<h2 id="basic-package-anatomy">Basic package Anatomy</h2>

Now lets look at the core files and directories that RStudio has
placed in this directory and describe what they are each for:

*  The DESCRIPTION file both describes and helps to configure the R package.
*  The NAMESPACE file determines which objects that the package defines will be visible to end users and which symbols from other packages need to be imported for the code in the package to function.
*  The 'R' directory is where your R code will go.  You can put
multiple files in here and by default they will load in alphanumeric
order.
* The 'man' directory is where you put .Rd files that will be
rendered by R into manual pages when they are needed by users.

These are just the parts that RStudio will put in for you. Some other important things you will need to know about are:

* The 'vignettes' directory which is in the base directory and is used to store '.Rnw' or '.Rmd' files.  These files will contain vignettes that explain how to use your package to new users.  
* The 'data' directory.  This is where you can put R objects stored in .Rda files that are intended to be used for examples.
* The 'inst/extdata' directory.  The 'inst' directory is just a generic place where you can put stuff that will be stored in your installed package.  And the inst/extdata sub directory is a very popular place to store other kinds of data files that you might need for examples or unit tests.


In the next several sections we will go over most of these components in a little bit more detail.  But remember, this is meant to be a high level guide for just getting you started.  So I am not going to put everything in here.  For that kind of guide I really think you should see the official guide at CRAN called [Writing R extensions](http://cran.fhcrc.org/doc/manuals/R-exts.html).

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>


<h2 id="description">The DESCRIPTION file</h2>

The DESCRIPTION file has a LOT of roles in an R package.  And technically it is a '.dcf' file.  It does its job by detailing a series of fields to tell R about what is in a package and how to use that information.  The one that was auto generated for you should already look like this:

```
    Package: MyPkg
    Type: Package
    Title: What the package does (short line)
    Version: 1.0
    Date: 2014-06-04
    Author: Who wrote it
    Maintainer: Who to complain to <yourfault@somewhere.net>
    Description: More about what it does (maybe more than one line)
    License: What license is it under?
```

Which is a good start really for almost any DESCRIPTION file.  You will see that for some fields they are giving you hints about what to put there and for other fields they have just filled them in already.  You will need to fill in at least that many fields (and a few more) in order to make a package for Bioconductor. Lets start by describing what goes into the fields that are present already in our example, and then I will talk a bit about adding some of the fields that are not.

* Package - This string is just the official package name.
* Type - This will basically always be 'Package'.  If you are wondering why this field even exists the short answer is: history.  The long answer is already explained in the [Writing R extensions manual](http://cran.fhcrc.org/doc/manuals/R-exts.html)
* Title - This is a slightly longer and more descriptive name for your package.
* Version - This is not filled in correctly for a Bioconductor package in the default example above.  Bioconductor packages need a three part version number X.Y.Z. The 'X' value is the major version number (normally 0 until there has been a release) The 'Y' value is normally bumped automatically for you with each release such that the development branch is an odd number while the release is an even number, so these usually skip a number each release.  And the 'Z' number is the minor release number which can be anything as long as it's bigger than the last number you used...
* Date - This is fine as is.  A lot of people leave this field out.
* Author - List of people (preferably with email addresses who wrote the original package)
* Maintainer - Who should people contact if they have a problem with this package.?  You must have someones name and email address here.  You must make sure this is up to date and accurate.
* Description - This is where you can write a nice verbose description of what your package is about.
* License - Here you need to fill in the official license.  Normally you will want to just use one of the 'official' licenses.  Which basically means choosing from one of the following strings for this field: 'GPL-2', 'GPL-3', 'LGPL-2', 'LGPL-2.1', 'LGPL-3', 'AGPL-3', 'Artistic-2.0', 'BSD_2_clause', 'BSD_3_clause', 'MIT'.

Those are the fields that RStudio will have created for you.  Here are some others that you should also think about:

* Depends - This field is for listing hard dependencies.  If you have another package and your package is really meant to be used along with that package then you might want to list it here.  Doing so will basically always load that other package before loading yours thus guaranteeing that all the contents exported by both packages are available to end users.
* Imports - This is what you should use when you want to access the contents of another package from within the code of your own package.  If that is your intention, then you should list those packages here.  This field is normally used in conjunction with import statements inside the NAMESPACE file that will spell out which functions and objects need to be available for your package.
* Suggests - This is for dependencies that are not needed in order to load and use your package.  If you have a function somewhere that is rarely used and that needs some specific package resource then Suggests is a good option.  Suggests is normally used along with require() to ensure that the package in question is loaded when and where it is needed. 
* Collate - Remember when I mentioned that the contents of the 'R' directory will be read in at load time in alphanumeric order.  Some of you probably cringed at that point.  Now is your chance to correct this.  Simple list the '.R' source files here in the order that you want them to load and R will do so. Failing to use a Collate field can cause peculiar behaviors later on so if you have more than on '.R' file then this is recommended.
* BiocViews - These are unique to Bioconductor and they describe each package using a controlled vocabulary. You can see an example of how they are used and also browse existing ones by clicking [here](/packages/devel/BiocViews.html#___Software). While you are browsing, be sure to check the box that says 'Developers: check this box to toggle the visibility of childless biocViews.'  That way you will be able to see the full ontology.  You can also make use of the recommendBiocViews(pkgdir) function from the biocViews package to suggest potential views for your package.


Now that we have discussed the DESCRIPTION file lets fill ours out so that it can actually work for our test package:

```
    Package: MyPkg
    Type: Package
    Title: A hello world package
    Version: 0.99.1
    Date: 2014-06-04
    Author: Some person <somePerson@someplace.com>
    Maintainer: Some person <somePerson@someplace.com>
    Description: We learn best how to program by practicing it.  
        So in order to practice creating packages, we are going to 
        copy this stuff into an actual DESCRITION file and test it 
        out for ourselves.
    License: Artistic-2.0
```

Once you have done this.  Build the package for yourself by hitting 'Ctrl+Shift+B' (or you can find the command for 'build and reload' on the 'build' menu).

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>


<h2 id="the-r-dir">The 'R' directory</h2>

The R directory is used for holding your '.R' source files.  The example we have has created a little file for you called 'MyPkg.R'.

Lets create a little test function in that package like this:

```
    myFun <- function(arg1){
        arg1 + 1
    }
```
now that you have code for a simple R function in your 'R' directory, hit 'Ctrl+Shift+B' to build and reload.  You will see now that not only will your package load, but that the function you just wrote is available if you list the contents of your package like this:

```
    library("MyPkg")
    ls("package:MyPkg")
```

And you can even test call your new function just like this:

```
    myFun(5)
```

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>


<h2 id="namespace">The NAMESPACE file</h2>

The NAMESPACE file is used to control which objects are exposed by a package (to end users) as well as which ones are available to it (as imported from other packages).  Unlike the DESCRIPTION file the NAMESPACE file is not a '.dcf' file.  It is actually meant to contain R code.  You can see a simple example of this by looking at the NAMESPACE file that was generated for you.  It contains the following line of code:

```
    exportPattern("^[[:alpha:]]+")
```

This basically tells R to export every single thing in your package. And this is usually a bad idea. In fact, this is such a consistently bad idea that if you send us a package that does this we will probably ask you to change it.  

Doing this will expose all the contents of your package to end users.  Including any helper functions or utilities that you may have meant to only be used internally.  This means that your users may have to sift through many of your functions that they don't need to know about just in order to learn about the ones that they do.  It also means that you will have to spend a lot of time documenting all these utilities (even ones that you don't intend for anyone else to ever use). This results in man pages that say things like 'for internal use only'.  There are times when documenting something that is used internally can be a good idea.  But you have to weight the benefits of having everything documented against the consequences of everyone having to read documentation that they don't really need.  In general, it is best to only document the stuff that end users intend to make use of.

There are a lot of commands that are normally put into a NAMESPACE file.  The 1st big class are commands to selectively export things that you want to document for end users.  The most common of these is export() which is used like this:

```
    export('myFun')
```

Now try this for yourself.  Replace the exportPattern() command with the call to export() that only exposes your function.  Then hit 'Ctrl+Shift+B' and make sure that it works as expected.

The other big use for NAMESPACE files is to import functions from other packages.  When you import a function you don't need to attach that package to the search path.  This is good as having a lot of things on the search path can slow R down. So by importing you allow your package to have access to this code without paying this performance penalty.  Lets suppose that you wanted to make use of the genomeStyles function from the GenomeInfoDb package.  In this case you need to do two things.  First you need to add GenomeInfoDb to the Imports field in the DESCRIPTION file.  This means you need to add a line like this to your DESCRIPTION file:

```
    Imports: GenomeInfoDb
```
    
And then add something like this to your NAMESPACE file.

```
    import('GenomeInfoDb')
```

Alternatively you could also add this to your NAMESPACE file:

```
    importFrom('GenomeInfoDb', 'genomeStyles')
```

Either of these approaches will get you access to the genomeStyles function from the GenomeInfoDb package and without fully loading the GenomeInfoDb package onto the search path.  The 1st case is more general and will make all the variables in the GenomeInfoDb package available.  But the 2nd case will only get the 'genomeStyles' function from the 'GenomeInfoDb' package.

For more details on how to use the NAMESPACE file, you should see [this](http://cran.r-project.org/doc/manuals/R-exts.html#Package-namespaces) section of the 'Writing R Extensions' manual.

<h3>Exercise 1: </h3>
Now create another new function.  This time, have your function return just the names of your imported genomeStyles function.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>




<h2 id="manual-pages">Documenting manual pages</h2> 

The R manual pages are written in 'Rd format'.  The style is mostly self explanatory and the example file 'MyPkg-package.Rd' that is provided for you is worth taking a moment to look at.  The topic of how to write these kinds of files, what the different fields mean and how to fill them out has already been very well documented [here](http://cran.r-project.org/doc/manuals/R-exts.html#Writing-R-documentation-files), in its own specific section of the 'Writing R Extensions' manual.  There are a few specific things worth mentioning though.

* Use aliases.  The alias tags are how the manual pages will know when to bring a user to your manual page.  So use a lot of appropriate aliases.  Also it is OK to document a few related functions on the same page.  Not only is this efficient, but it will make it more likely that you will come back and do a better job of maintaining your documentation if it is less repetitive.

* When possible fill in the see also field.  If you have other related functions you should reference them here.

* Always take credit for your work.  Nobody is going to think you are vain.  But they might wonder about your code if you are shy about owning up to having created it...  For the same reasons you should also be generous with giving people credit if they have contributed code.

* Always be thinking about the examples section.  The most important thing on your entire manual page is the example section.  If you don't have good examples, it will be hard for many users to sort out how they should use your code.  So plan ahead and write examples that demonstrate clearly how the functions should be called and that run quickly.  It is not important that the example data produces a meaningful result, but it is important that the code clearly demonstrate how it is intended to be used.

Now look again at that manual page that RStudio generated for your new package ('Rd format').  You will notice that (among other things) the example section is not filled in.  If you were to try and run check on your package right now you would get an error from this.  Go ahead and try that to see for yourself what will happen.  From the 'Build' menu, choose 'Check Package' (CTRL-SHIFT-E).  You should get an error.

To correct this put a simple example into your manual page.  This is just the page that describes the package so there is not much for it to demonstrate.  For now just add the following simple example section as a replacement for what is currently there: 

```
    myFun(10)
```

Then run check again.  You will still have warnings, but it should at least complete the process.

Now lets look at another one of RStudios nicer features.  Lets add a manual page for our new function.

To do this go to the 'File' menu, choose 'New File' and 'Rd File'.  Be sure to choose a man page that is pre-configured for a 'Function' and be sure to name it 'myFun' after our new function. This should drop a new 'pre-fabricated' manual page into your packages existing 'man' directory.  

Now obviously you need to fill in your manual pages with actual descriptions and working examples etc.  But for now fill in a value for the new manual page title so that R will be able to build and check the package again.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>




<h2 id="vignettes">Writing vignettes</h2> 

You might be wondering why you would need a vignette if you already have a series of manual pages. The answer is that end users need a high level guide to show them how to connect the dots.  It's not enough to just show how all the pieces work. You also need to know how they are supposed to work together.  A good vignette should first explain what the package is for, and then it should provide a common workflow for how that package could be used to accomplish that.  In doing so, it puts all the functions and classes into context so that new users can find their way.

Structurally, vignettes go in the 'vignettes' directory. These consist of either an '.Rnw' or a .Rmd' file. That R will render into a .pdf or html document later.  It used to be that most vignettes were written in Sweave as described [here](http://cran.r-project.org/doc/manuals/R-exts.html#Writing-package-vignettes), but now more and more users are choosing to write them using the knitr package as described [here](http://cran.r-project.org/doc/manuals/r-devel/R-exts.html#Non_002dSweave-vignettes).  Either form of vignette is acceptable, but writing a '.Rmd' document is a lot easier to do.

To get RStudio to create a vignette for you simply go to the 'File' menu, and look under 'New File'. From here you could choose to create either a new 'R Sweave' or 'R Markdown' file.  But for our purposes in this demo, choose an 'R Markdown' file.  Now at this point RStudio will create a sample R Markdown file for you.  But unlike the case with the man pages, it will not put it in the correct place for you (at least not as of the time of this writing).  So go again to the 'File' menu and choose 'Save As'.  This will open up a little widget for you to specify how you want this markdown file to be saved.  From here click the button at the bottom left called 'New Folder' and create a sub directory in your MyPkg directory called 'vignettes'.  Then tell the widget to save the markdown file there under the name 'MyPkg.Rmd'.

Finally you must add the following lines to your DESCRIPTION file:

```
    Suggests: knitr
    VignetteBuilder: knitr
```

These are needed because knitr is needed to support non-Sweave based markdown files as vignettes and knitr is not part of base R.


You now have the beginnings of a vignette.  Now there are really two things that a vignette must always do: 1) it must explain sufficient background information so that new users can determine if a package is going to help them solve a problem or not (even if they are new to the topic) and 2) it must demonstrate how to use the various functions from the package in the expected way.  In order to do the latter you will sometimes need to repeat examples that you may have already shown in the example sections of your manual pages.  Running these examples twice can be inefficient.  How can you avoid this problem?  Well you can always mark your code chunks as eval=FALSE.  Here is an example of a markdown code chunk that is marked to not be evaluated:


<pre>
    &#96;&#96;&#96;&#123;r chunkName, eval&#61;FALSE&#125;
        myFun(3)
    &#96;&#96;&#96;
</pre>


If you were using Sweave, it would have looked like this instead:

<pre>
    &#60;&#60;chunkName, eval&#61;FALSE&#62;&#62;&#61;
        myFun(3)
    &#64;
</pre>


<h3>Exercise 2: </h3>
Now Document your genomeStyleNames function with a manual page and also update your vignette.  Be sure to include examples in each case that demonstrate how to use it.  When you put the example into your vignette mark the code chunk as FALSE.


<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>






<h2 id="external-data-dirs">Storing data in your package</h2> 

As mentioned above, there are basically two places where R packages can store helpful user data.  The 1st case is for when you have your data already in an R object.  In that case you probably want to just save the R object in the 'data' directory.  This directory is where you can put '.Rda' objects that you want to document and use in examples by calling load().  

The second case it the 'inst/extdata' directory.  The 'inst' directory is just a place where you can put stuff that will be stored in your installed package for retrieval later.  You can think of it a little bit like a storage place for miscellaneous items.  For example, some users will put an 'inst/scripts' directory in their package as a place to store scripts.  But the 'inst/extdata' sub-directory combination gets special mention here because so many packages use it as a place to store extra data like a database or '.csv' files etc.  Basically you should think of using it for anything data oriented that is not a '.Rda' file.  This comes up a lot because not everything that you may want to load into your R session as an example of your packages abilities will necessarily start as a '.Rda' file.

Finally, we put a maximum value on how big the built tarball for a software package can be (five megabytes).  If you need more data than this, it is time to consider making an companion package. A companion package can be either an experiment data package or an annotation package. In either case its primary use will be to hold data for an end user.  The biggest pragmatic difference between an experiment data package and an annotation package is primarily how often they are updated.  Annotations are data that people need to rely upon to be current, so these data need to be updated (or checked to make sure they are still current) for every release.  In contrast most experimental data never really changes once the experiment has been run and so it is used primarily so that people can have an authentic data set for running their algorithms on.

<h3>Exercise 3: </h3>
Part 1: GenomeStyles() returns a list object.  Save this list object out to the file system and add it to a 'data' directory so that it can be loaded by end users.
  
Part 2: Now call genomeStyles('Homo_sapiens') and trap the result.
Then use write.table() to create a tab delimited file with the contents of that data.frame on disc.   Finally, store that information in inst/extdata and write an extractor function using system.file that pulls that data back into your R session.



<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>



<h2 id="bioccheck">BiocCheck</h2> 

If you look in our [package guidelines](/developers/package-guidelines/) you will see that not only must a package pass R CMD build and R CMD check without ERRORs or WARNINGs, but it must also pass a bunch of other requirements intended to make sure that your code is well documented and that it will integrate well into the rest of the project.  Fortunately we have a tool to help you with this.  From RStudio you can run build and check from the 'Build' menu.  But you can also test your own package for BiocCheck compliance like this:

```
    library(BiocCheck)
    pathToPkg <- file.path("..","MyPkg")
    BiocCheck(pathToPkg)
```

If you have followed the instructions up to this point, then calling the function above will produce an list object as output to tell you which things are needed for your package to be compliant with the [Bioconductor package guidelines](/developers/package-guidelines/).  The BiocCheck output is broken into three distinct categories:

* Requirements are things that we are confident will always be necessary.  if you have failed some of these, you should aim to fix it right away.
* Recommendations are things that are usually really important to change, but their interpretation can sometimes depend on context.  If you are failing some of these you should try to fix it if at all possible.
* Considerations are usually a good idea.  But they might or might not apply in your particular case.
    

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>



<h2 id="good-practices">Tips for writing good packages</h2> 

Entire books have been written about good coding practices.  But here are just a few tips for writing good Bioconductor packages.

* Write small functions.  It's best if each function has only one job that needs to do.  And it's also best if it does that job in as few lines of code as possible.  If you find yourself writing great big functions that wrap on for more than a screen then you should probably take a moment to split it up into smaller helper functions.  Smaller functions are easier to read, debug and to reuse.
* Don't trust users to provide correct argument values. For functions that you export, check to make sure that the values provided are the correct kind and are in the expected range etc..  It's not uncommon in a well written function for half a function to be dedicated to this kind of checking.
* Write live examples for your manual pages and vignettes.  This may sometimes means taking a minute to make sure that your code will run fast enough to get through check in less than five minutes.  But it is worth it to have your code checked when the package is checked.
* Use toy data sets.  Unless your data is abnormally small, genomic data tends to be large.  This means that using 'real' data in your examples will likely take too long.  Instead we recommend that you use 'toy' data sets.  It doesn't matter if the data is real as long as it allows you to demonstrate that the function still works and to show someone else how to use it.
* Write unit tests that quickly test internal steps or to verify important functionality in your code.  We have an entire [workflow](/developers/how-to/unitTesting-guidelines/) written about how to implement unit tests.  It is well worth your time to read it.  Using unit tests is a great way to make sure that your code is robust and continues to do what it initially advertised.  Later on, when some upstream change breaks your code, unit tests will allow you to quickly and easily identify and fix the problem.
* Choose a good package name.  Make sure that your package name is not already in use and try to make sure that it gives the end user some idea about what the package will actually do for them.
* Make sure that your package idea is actually new.  A lot of packages get written that provide functionality that either already exists or is only marginally different from an existing package.  Before you spend a lot of time on a package look at the existing [website](/packages/release/BiocViews.html#___Software) to make sure that someone has not already done all this work for you.  Your contribution will be appreciated by the community much more if it builds on and extends the existing work rather than reinventing the wheel.


<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>


<h2 id="exercise-answers">Answers for exercises:</h2> 

<h3>Exercise 1: </h3>
You should have created a function that looks like this:
```
    genomeStyleNames <- function(){
        names(genomeStyles())
    }
```

And you should have exported it so that your NAMESPACE should now look something like this:

```
    export("myFun","genomeStyleNames")  
    import('GenomeInfoDb')
```

<h3>Exercise 2: </h3>
You can either choose to document this function in the manual pages by either adding a new alias to your older manual page or by creating a new page entirely.  Either way there is not really a 'correct' answer for how to document things (or at least not one that R will not already tell you about).



<h3>Exercise 3: </h3>
Part 1: 
You should have used simple code like this to save the R object to file
``` 
    gs <- genomeStyles()
    save(gs, file='gs.rda')
```

Then once the file 'gs.Rda' was stored in 'data' you could do this once your package was loaded. 

```
    data(package='MyPkg')
    data(gs)
```


  
Part 2: Now call genomeStyles('Homo_sapiens') and trap the result.
Then use write.table() to create a tab delimited file with the contents of that data.frame on disc.   Finally, store that information in inst/extdata and write an extractor function using system.file that pulls that data back into your R session.

```
    res <- genomeStyles('Homo_sapiens')
    write.table(res, sep="\t",row.names=FALSE, col.names=FALSE, 
                file='human.txt')
```

Then if you move this to inst/extdata you can write a simple function to extract it like this:


```
  getData <- function(){
    read.delim(file=
               system.file('extdata','human.txt', package="MyPkg"), 
               header=FALSE)
  }
  getData()
```





