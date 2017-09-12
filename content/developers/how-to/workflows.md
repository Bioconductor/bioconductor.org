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

## How do I write and submit a workflow vignette?

* Write a package with the same name as the workflow. The workflow vignette
 written in Markdown, using the [rmarkdown](http://rmarkdown.rstudio.com/)
 package should be included in the vignette directory. The package does not need
 man/ or R/ directories nor a data/ directory as ideally workflows make use of
 existing data in a Bioconductor repository or on the web; the workflow package
 itself should not contain large data files. Please include a detailed
 Description in the DESCRIPTION file as this will be used as the abstract on the
 Bioconductor web site alongside [the other workflows](/help/workflows/).

* The package should pass R CMD build in <= 20 minutes and should require
  <= 4GB RAM.

* In the DESCRIPTION file, include the line "Workflow: True"

* Submit the package to the [GitHub submission
  tracker](https://github.com/Bioconductor/Contributions) for a formal
  review. Please also indicate in the tracker issue that this package is a
  workflow. 
 
* Once the package is approved a new directory will be created in our SVN repository
  under "/trunk/madman/workflows/\<YourPackage\>" that you be given read/write
  access to. Your SVN credentials will be sent to you via email after
  approval. The SVN repository is where you will make any future commits/updates
  to your workflow package. (You can view existing workflow sources
  [here](https://hedgehog.fhcrc.org/bioconductor/trunk/madman/workflows/),
  username and password is **readonly**.)

* A Bioconductor Team member will "Create New Jenkins Project" for your
  package. Jenkins is an Amazon Web Serice (AWS) instance that manages the
  Biocondcutor workflows. You can access Jenkins at [DocBuilder Web
  App](https://docbuilder.bioconductor.org/app/). Log in with your SVN username
  and password. Every time you commit a change to your workflow directory, the
  workflow builder will try and build your vignette on Mac, Windows, and
  Linux. You will receive an email if there were any errors. You can monitor the
  progress of builds [here](http://docbuilder.bioconductor.org:8080/). After a
  successful build, the Bioconductor website will be updated accordingly. 

## Using Math Symbols in a Markdown workflow vignette

If you want to include math symbols in a workflow vignette, put the following 
snippet at the beginning of your .Rmd file:

    <script type="text/javascript"
      src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
    </script>

Then you can use math the same way you would in `LaTeX`, except the symbols for escaping it 
are different. For inline formulas, use <span>\\</span>(N<span>\\)</span>, and for displayed 
equations, use <span>$</span>$N<span>$</span>$. 

The first will render as \\(N\\) and the second as $$N$$ .

See the [Mathjax](https://www.mathjax.org/) documentation for
more information.

## Tidying package loading output

Most workflows load a number of packages and you do not want
the output of loading those packages to clutter your workflow 
document. Here's how you would solve this in markdown; you can
do something similar in Latex.

First, set up a code chunk that is evaluated but not echoed, and whose
results are hidden. We also set `warning=FALSE` to be sure that 
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

## Citations

To manage citations in your workflow document,
specify the bibliography file in the document metadata header.

    bibliography: references.bib
    
You can then use citation keys in the form of [@label] to cite an entry with an identifier "label".

Normally, you will want to end your document with a section header "References" or similar, after which the bibliography will be appended.

For more details see the [rmarkdown documentation](http://rmarkdown.rstudio.com/authoring_pandoc_markdown.html#citations).

## Questions

If you have any questions, please ask on the bioc-devel
[mailing list](/help/mailing-list).



