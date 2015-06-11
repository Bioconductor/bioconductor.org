<script type="text/javascript"
  src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
</script>

# Creating Workflow Vignettes

## What is a workflow vignette?

Workflow vignettes are documents which describe a bioinformatics workflow that involves 
multiple Bioconductor packages. These workflows are usually more extensive than 
the vignettes that accompany individual Bioconductor packages.

[Existing Workflow Vignettes](/help/workflows/)

Workflow vignettes may deal with larger data sets and/or be more computationally intensive
than typical Bioconductor package vignettes. For this reason, the automated builder that
produces these vignettes does not have a time limit (in contrast to the Bioconductor package 
building system which will time out if package building takes too long).

## Who should write a workflow vignette?

Anyone who is a bioinformatics domain expert.

## How do I write a workflow vignette?

* Request that a new directory be created in our SVN repository
  under "/trunk/madman/workflows" and that you be given read/write access to this
  directory. You can request this access by emailing 
  "maintainer at bioconductor dot org". (You can view existing workflow sources
  [here](https://hedgehog.fhcrc.org/bioconductor/trunk/madman/workflows/), username
  and password is **readonly**.)

* Write a vignette in LaTeX or Markdown, using the 
 [knitr](http://yihui.name/knitr/) package. Commit it to the 
 svn location above. Alternatively, you can write a full
 R package in this location (in this case, it's not required to use
 knitr for your vignette).

 * Go to the [DocBuilder Web App](https://docbuilder.bioconductor.org/app/).
   Log in with your SVN username and password.

* Click on "Create New Jenkins Project".
  Fill in the directory name from the first step above, and your email address.

* The workflow builder will now try and build your vignette on Mac, Windows, and Linux.
  You'll receive an email if there were any errors. You can monitor the progress of
  builds [here](http://docbuilder.bioconductor.org:8080/).

* Every time you commit a change to your workflow directory, another 
  build will be triggered and you will receive email if it fails.

* When you are ready for your workflow to appear on the 
  Bioconductor web site, contact "maintainer at bioconductor dot org"
  and we will allow the workflow to propagate to our web site where it
  will be listed alongside [the other workflows](/help/workflows/).
  It will be updated every time there is an SVN commit and 
  a successful build.

## Using Math Symbols in a Markdown workflow vignette

If you want to include math symbols in a workflow vignette, put the following 
snippet at the beginning of your .Rmd file:

    <script type="text/javascript"
      src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
    </script>

Then you can use math the same way you would in `LaTeX`, except the symbols for escaping it 
are different. For inline formulae, use <span>\\</span>(N<span>\\)</span>, and for displayed 
equations, use <span>$</span>$N<span>$</span>$. 

The first will render as \\(N\\) and the second as $$N$$ .

See the [Mathjax](https://www.mathjax.org/) documentation for
more information.

## Adding metadata to a Markdown vignette

With a markdown vignette (.Rmd file) you can add a YAML metadata preamble
at the beginning of the file, for example:

    ---
    Author: My Name
    Date: March 16, 2015
    ---

This metadata will be rendered in the "About This Document" box that appears at the
top of the web page.

## Tidying package loading output

Most workflows load a number of packages and you do not want
the output of loading those packages to clutter your workflow 
document. Here's how you would solve this in markdown; you can
do something similar in Latex.

First, set up a code chunk that is evaluated but not echoed, and whose
results are hidden. We also set warning=FALSE to be sure that 
no output from this chunk ends up in the document:

    ```{r, echo=FALSE, results="hide", warning=FALSE}
    suppressPackageStartupMessages({
        library(GenomicRanges)
        library(GenomicAlignments) 
        library(Biostrings)
        library(Rsamtools)
        library(ShortRead)
        library(BiocParallel)
        library(rtracklayer)
        library(VariantAnnotation)
        library(AnnotationHub)
        library(BSgenome.Hsapiens.UCSC.hg19)
        library(RNAseqData.HNRNPC.bam.chr14)
    })
    ```

Then you can set up another code chunk that *is* echoed,
which has almost the same contents. The second invocation of `library()`
will not produce any output since the package has already been loaded:


    ```{r}
    library(GenomicRanges)
    library(GenomicAlignments) 
    library(Biostrings)
    library(Rsamtools)
    library(ShortRead)
    library(BiocParallel)
    library(rtracklayer)
    library(VariantAnnotation)
    library(AnnotationHub)
    library(BSgenome.Hsapiens.UCSC.hg19)
    library(RNAseqData.HNRNPC.bam.chr14)
    ```

## Questions

If you have any questions, please ask on the bioc-devel
[mailing list](/help/mailing-list).



