<script type="text/javascript"
  src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
</script>

# Creating Workflow Package

The main focus of a workflow package is the vignette!

## What is a workflow vignette?

Workflow vignettes are documents which describe a bioinformatics workflow that involves
multiple Bioconductor packages. These workflows are usually more extensive than
the vignettes that accompany individual Bioconductor packages.

[Existing Workflows](http://bioconductor.org/packages/devel/BiocViews.html#___Workflow)

Workflow vignettes may deal with larger data sets and/or be more computationally intensive
than typical Bioconductor package vignettes. For this reason, the automated builder that
produces these vignettes does not have a time limit (in contrast to the Bioconductor package
building system which will time out if package building takes too long). It is
expected the majority of vignette code chunks are evaluated.

## Who should write a workflow vignette?

Anyone who is a bioinformatics domain expert.

## How do I write and submit a workflow vignette?

* Write a package with the same name as the workflow. The workflow vignette
 written in Markdown, using the [rmarkdown](http://rmarkdown.rstudio.com/)
 package should be included in the vignette directory. You may include more than
 one vignette but please use useful identifying names.

* The package does not need man/ or R/ directories nor a data/ directory as
 ideally workflows make use of existing data in a Bioconductor repository or on
 the web; the workflow package itself should not contain large data files.

* In the DESCRIPTION file, include the line "BiocType: Workflow". Please also
  include a detailed Description field in the DESCRIPTION file. The DESCRIPTION
  file should contain biocViews which should be from the [Workflow
  branch](http://bioconductor.org/packages/devel/BiocViews.html#___Workflow). If
  you think a new term is relevant please reach out to
  <lori.shepherd@roswellpark.org>.

* Submit the package to the [GitHub submission
  tracker](https://github.com/Bioconductor/Contributions) for a formal
  review. Please also indicate in the tracker issue that this package is a
  workflow.

* Workflows are git version controlled. Once the package is accepted it will be
  added to our git repository at git@git.bioconductor.org and instructions will
  be sent for gaining access for maintainence.


## Consistent formatting

* In an effort to standardize the workflow vignette format, it is strongly
  encouraged to use either BiocStyle for formatting or utilize
  [BiocWorkflowTools](http://bioconductor.org/packages/BiocWorkflowTools/). 
  The following header shows how to use BiocStyle in the vignette:

	```
	output:
		BiocStyle::html_document
	```

* The following should also be include

      - author affiliations
      - a date representing when the workflow vignette has been modified

* The first section should have some versioning information. The  R version,
  Bioconductor version, and package version should be visible. 
  The following is an example of how this could be achieved:

	```
	<p>
	**R version**: `r R.version.string`
	<br />
	**Bioconductor version**: `r BiocManager::version()`
	<br />
	**Package**: `r packageVersion("annotation")`
	<\p>
	```
* An example start to a workflow vignette:

The following is taken as an example header from the variants workflow package:

<pre>
    <code>
&ndash; &ndash; &ndash;
title&#58; Annotating Genomic Variants
author&#58; 
&ndash;name&#58; Valerie Obenchain
  affiliation&#58; Fred Hutchinson Cancer Research Center, 1100 Fairview Ave. N., P.O. Box 19024, Seattle, WA, USA 98109&ndash;1024
date&#58; 11 April 2018
vignette&#58; &#62;
  &#37;&#92;VignetteIndexEntry&#123;Annotating Genomic Variants&#125;
  &#37;&#92;VignetteEngine&#123;knitr&#58;&#58;rmarkdown&#125;
output&#58; 
    BiocStyle&#58;&#58;html&#95;document
&ndash; &ndash; &ndash;


&#35; Version Info
&#96;&#96;&#96;&#123;r, echo=FALSE, results=&quot;hide&quot;, warning=FALSE&#125;
suppressPackageStartupMessages&#40;&#123;library&#40;&#39;variants&#39;&#41;&#125;&#41;
&#96;&#96;&#96;
&#60;p>
&#42;&#42;R version&#42;&#42;&#58; &#96;r R.version.string&#96;
&#60;br &#47;&#62;
&#42;&#42;Bioconductor version&#42;&#42;&#58; &#96;r BiocInstaller::biocVersion&#40;&#41;&#96;
&#60;br &#47;&#62;
&#42;&#42;Package version&#42;&#42;&#58; &#96;r packageVersion&#40;&#34;variants&#34;&#41;&#96;
&#60;&#47;p&#62;


    </code>
</pre>

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

You can then use citation keys in the form of &#91;@label&#93; to cite an entry with an identifier "label".

Normally, you will want to end your document with a section header "References" or similar, after which the bibliography will be appended.

For more details see the [rmarkdown documentation](http://rmarkdown.rstudio.com/authoring_pandoc_markdown.html#citations).

## Questions

If you have any questions, please ask on the bioc-devel
[mailing list](/help/mailing-list).
