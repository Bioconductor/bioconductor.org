<!-- %\VignetteEngine{knitr::knitr} %\VignetteIndexEntry{} -->

<!-- To deploy this to the website, do the following: library(knitr)
render_markdown(strict=TRUE) knit("2015_April.Rmd")

This will produce an .md file that can be used  on the website. The {:toc} and
{:no_toc} tags will be rendered by nanoc, the tool that produces the website.
-->

{::options parse_block_html="true" /}

# _Bioconductor_ Newsletter
{:.no_toc}

posted by [Valerie Obenchain](mailto:valerie.obenchain@roswellpark.org),
October 2015

The _Bioconductor_ newsletter is a quarterly review of core infrastructure
developments, community projects and future directions. We cover topics of
general interest as well as those with the greatest impact on the software.
This quarter the _Bioconductor_ core team relocated from Seattle, Washington to
Buffalo, New York and welcomed two new team members.  Two new forums for
workflows and analysis pipelines were introduced: the recently launched
_Bioconductor_ F1000Research Channel and the new section for
_R / Bioconductor_ Workflows in the BMC Source Code for Biology and Medicine
journal. New features and functions are summarized for the infrastructure
packages as well as a handful of contributed packages that have been
especially active over the past devel cycle.

## Contents
{:.no_toc}

* Table of contents will replace this text.
{:toc}


## _Bioconductor_ relocates to Roswell Park Cancer Institute

In September the _Bioconductor_ core team relocated from Seattle, Washington to
Buffalo, New York. The new home institution is
[Roswell Park Cancer Institute](https://www.roswellpark.org/about-us)
(RPCI) which is part of the Buffalo Niagara Medical Complex. The Complex covers
about 120 acres in downtown Buffalo and is a consortium of hospitals, health care facilities and
educational institutions. The campus is growing rapidly with current employee
numbers at 12000 estimated to reach 17000 by 2017. RPCI has close ties with
SUNY University at Buffalo.

The change in physical location also brought some staffing changes. Marc
Carlson, Sonail Arora and Paul Shannon left the project at the end of August.
All were instrumental in the design and implementation of [AnnotationHub][AH]
(and annotations in general), tools for biological network analysis,
educational material and many other areas. A big thanks to them for their many
contributions.

In Buffalo we are pleased to welcome new team members Jim Java and Brian Long.
Jim formerly worked as a biostatistician for the Gynecologic Oncology
Group at RPCI where he analyzed clinical trial data. He also worked as a
software engineer for several companies on projects ranging from point-of-sale
software to embedded systems. He has a Ph.D. in computer science (focused on
natural language processing), and MA degrees in biostatistics and English
literature. His research interests include scientific and statistical
computing, natural language processing (especially authorship identification),
linguistics, climate science, gun violence, ethics, and literature.

Brian's work in software development has focused on relational database & NoSQL
storage, UI/UX design, security and scale. Most recently, he's been working
with a team to develop a Platform-as-a-Service (similar to Google App Engine).
The platform offering allows application developers to deploy components to
on-demand virtual infrastructure and supports a large number of programming
languages via Apache Thrift. Most of the work was not done as open source but
the [code](https://github.com/ezbake) has been released on GitHub. Feel free
to check it out!


## Reproducible Research

### F1000Research channel launched

The _Bioconductor_ [home page](http://www.bioconductor.org/) has a new link
for the
[F1000 Research Channel](http://f1000research.com/channels/bioconductor).

[F1000Research](http://f1000research.com/) is an open science publishing
platform with the following goals:

- Fast publication:
  All scientific research is published with a few days making new results
  immediately visible. Review is conducted after publication.

- Open peer review:
  Open peer review removes the potential bias that can be present with anonymous
  pre-publication review.

- Publication of all findings:
  All results are published including null/negative results, small findings,
  case reports, data notes and observation articles.

- Data provided:
  All articles are accompanied by the data used to generate the results
  enabling reproducible research.

The motivation for a _Bioconductor_ F1000Research channel was a forum for
task-oriented workflows that cover a current, genome-scale analysis problem
from start to finish. In contrast to package vignettes, these workflows combine
resources from several different packages, demonstrating integrative analysis
and modeling techniques. The end goal is to make it easier for _Bioconductor_
users to navigate the software offerings by adapting these workflows to
quickly arrive at a solution for their specific problem.

This effort has been led by Wolfgang Huber, Vince Carey, Sean Davis, Kasper
Daniel Hansen and Martin Morgan.

### Source Code for Biology and Medicine hosts new _R / Bioconductor_ Workflows section

Levi Waldron is an Assistant Professor of Biostatistics at the CUNY School of
Public Health in New York City. He serves on the _Bioconductor_ technical
advisory board and is currently heads up the group developing new
infrastructure for the analysis of
[multi-assay genomic experiments](https://github.com/vjcitn/biocMultiAssay).
Levi recently became a section editor for the
[BMC Source Code for Biology and Medicine Journal](http://www.scfbm.org/)
and earlier this month posted an
[announcement](https://support.bioconductor.org/p/71811/) calling for
_R / Bioconducor_ workflows submissions to the new section. The journal is
open access and aims to publish source code for distribution and use
to advance biological and medical research. Manuscripts are considered on all
aspects of workflow for information systems, decision support systems, client
user networks, database management, and data mining.

The workflows should address a problem of general interest and be easily
adaptable by other researchers. Implementation can be as a workflow on the
[Bioconductor web site](http://www.bioconductor.org/help/workflows/), an
[Amazon](http://www.bioconductor.org/help/bioconductor-cloud-ami/) or
[docker](http://www.bioconductor.org/help/docker/) image or other
cross-platform supported approach. More details are available in the
[author instructions](http://www.scfbm.org/authors/instructions/workflow).


## Infrastructure Updates

### Build machines to cloud

Over the past weeks Dan has been busy transferring package builds, the
single package builder and the web site to the cloud. The move provides more
flexible login access as well as machine configuration. Everything should look
the same from the outside so user interaction with these resources will
not change.

This has been a huge effort and is almost complete. Thanks Dan!

### HTS core package stack

_Bioconductor_ encourages software reuse and aims for a flexible, integrated
set of infrastructure packages. One consequence of interrelationship is
that a change in a low-level package can affect packages downstream.

Throughout a devel cycle new classes are added and code is reorganized which may
cause the dependency hierarchy to change. We do our best to identify and fix
packages affected by these changes. Modified packages and their dependencies
are committed to svn/git and should propagate together through the
nightly builds and become available via `biocLite()` the following day.

Herv&eacute; recently added a graphic to S4Vectors/inst/doc/ that depicts the
High Throughput Sequencing (HTS) core package stack. Knowledge of this
hierarchy is useful when developing new S4 classes and methods. It also
benefits leading-edge developers working directly from svn (vs `biocLite()`). A package
installed from svn may have updated dependencies that have not yet propagated
through the build system and must be installed (in order) by hand.

The package stack file will be updated as the core packages change. This is the
current snapshot.

  as of August 2015

                   VariantAnnotation
                        |     |
                        v     v
           GenomicFeatures   BSgenome
                        |     |
                        v     v
                      rtracklayer
                           |
                           v
                   GenomicAlignments
                      |         |
                      v         v
     SummarizedExperiment   Rsamtools
                    |       |      |
                    v       v      v
                GenomicRanges   Biostrings
                        |          |
                        v          v
               GenomeInfoDb   XVector
                        |     |
                        v     v
                        IRanges
                           |
                           v
                       S4Vectors


### New functions

All are available in the devel branch, _Bioconductor_ 3.2:

*   *GenomicFeatures::coverageByTranscripts()*

    Computes transcripts (of CDS) coverage of a set of ranges.

    (contributed by Herv&eacute; Pages)

*   *improvements to rtracklayer::import() for GFF files*

    Reads data from a GFF file into a data.frame or DataFrame object. The
    function auto-detects the GFF version but has a 'version' argument if
    needed. All columns are loaded by default or individual columns can be
    specified in the 'columns' argument. Additional flexibility is provided
    by `rtracklayer::readGFF()`.

    (contributed by Herv&eacute; Pages)

*   *coercion between GRanges object and character vector*

    A `GRanges` object can be coerced to a character vector and back. See
    ?GRanges in the [GenomicRanges][GR] package for details.

    (contributed by Herv&eacute; Pages)

        ## From GRanges to character:
        > gr <- GRanges("chr1", IRanges(1, 5), "-")
        > x <- as.character(gr)
        > x
        [1] "chr1:1-5:-"

        ## From character to GRanges:
        > GRanges(x)
        GRanges object with 1 range and 0 metadata columns:
              seqnames    ranges strand
                 <Rle> <IRanges>  <Rle>
          [1]     chr1    [1, 5]      -
          -------
          seqinfo: 1 sequence from an unspecified genome; no seqlengths

*   *coercion to and from ExpressionSet and SummarizedExperiment*

    The venerable `ExpressionSet` can be coerced to and from the more
    modern `SummarizedExperiment`. Coercion using `as()` supports
    mapping identifiers from common probe- or gene-based annotation
    labels to genomic ranges;
    `makeSummarizedExperimentFromExpressionSet()` allows
    user-specified identifier conversions.

    (contributed by Jim Hester)

## Activity in contributed packages

As complement to the section on new features and functions in the
infrastructure packages we want to highlight significant changes made in
contributed packages. Not all packages keep current NEWS files so it can be
tricky to determine which packages have added new features. One way of gauging
active development is the number of svn/git commits over a period of time.

As of September 22, these packages all had 50+ commits since the April 2015
release: [CopywriteR][CopywriteR], [systemPipeR][systemPipeR], 
[ComplexHeatmap][ComplexHeatmap], [derfinderHelper][derfinderHelper], 
[ggtree][ggtree], [RnBeads][RnBeads] and [cogena][cognea].

A few authors said the commits were due to maintenance and the package had
not changed much since the last release. Other packages did change
significantly and the authors have summarized the changes below. Comments have
been lightly edited for length.

  [ComplexHeatmap][ComplexHeatmap]
  Author: Zuguang Gu

  This package provides a framework to combine and visualize multiple heat
  maps. Combined maps can be flexibly annotated and decorated (add graphics
  post-generation). The `oncoPrint()` function offers a compact means of
  visualizing genomic alteration events.

  - vignettes are separated into separate topics and are extensively improved.

  - support text rotation for heatmap titles.

  - rows can be split if `cluster_rows` is a clustering object.

  - legends for continuous values can be set as continuous color bars.

  - row title and column title as well as legend title support math expressions.

  - add `densityHeatmap()` which visualizes density distribution in a
  matrix/list through a heatmap.

  - add `decorate*` functions which make it easy and straightforward to add
  more graphics on the plot.

  - add `oncoPrint()` which makes it easy to make oncoprints

  - add `select()` function to interactively select sub-region in the heatmap
  and retrieve row/column index in the selected sub region.

  - add `row*` and `column*` helper functions (e.g. `rowAnnotation()`,
  `row_anno_barplot()`)


  [ggtree][ggtree]
  Authors: Guangchuang Yu and Tommy Tsan-Yuk Lam

  A phylogenetic tree viewer with extensive capabilities for adding node-level
  (taxa) annotation layers. Extends the graph grammar and infrastructure in
  `ggplot` and `ggplot2`.

  - `read.raxml()` and `read.r8s()` to support RAxML and r8s input;
  data are stored in `raxml` and `r8s` objects

  - `merge_tree()` which combines statistical evidences inferred from
  different software making it possible to compare results

  - `gheatmap()` which annotates a tree with associated numerical matrix
  (e.g. genotype table)

  - `msaplot()` which annotates a tree with multiple sequence alignment

  Both `gheatmap()` and `msaplot()` add a new layer to tree view and can be
  transformed to circular form by adding `+coord_polar(theta="y")` to the
  grammar.


  [cogena][cogena]
  Authors: Zhilong Jia and Michael Barnes

  cogena is a workflow for co-expressed gene-set enrichment analysis.

  - Add pipeline for drug discovery and drug repositioning based on the cogena
  workflow. Candidate drugs can be predicted based on the gene expression of
  disease-related data, or other similar drugs can be identified based on the
  gene expression of drug-related data. Moreover, the drug mode of action can
  be disclosed by the associated pathway analysis.

  - add functions `coExp()` and `clEnrich()` used in the pipeline

  - add gene sets CmapDn100.gmt and CmapUp100.gmt, based on
  [Connectivity
  Map](https://www.broadinstitute.org/genome_bio/connectivitymap.html)
  to enable drug repositioning analysis

  - add new gene set MsigDB 5.0


  [systemPipeR][systemPipeR]
  Author: Thomas Girke

  [systemPipeR][systemPipeR] provides infrastructure for building and running
  automated analysis workflows for a wide range of next generation sequence
  (NGS) applications.

  NGS WORKFLOWS

  - added new end-to-end workflows for Ribo-Seq and polyRibo-Seq, ChIP-Seq and
    VAR-Seq

  - added the data package 'systemPipeRdata' to generate systemPipeR workflow
    environments with a single command (genWorkenvir) containing all parameter
    files and sample data required to quickly test and run workflows

  - about 20 new functions have been added to the package. Some examples are:
    - Read pre-processor function with support for SE and PE reads
    - Parallelization option of detailed FASTQ quality reports
    - Read distribution plots across all features available in a
      genome annotation (see ?featuretypeCounts)
    - Visualization of coverage trends along transcripts summarized
      for any number of transcripts (see ?featureCoverage)
    - Differential expression/binding analysis includes now DESeq2 as
      well as edgeR

  - adaption of R Markdown for main vignette

  WORKFLOW FRAMEWORK

  - simplified design of complex analysis workflows

  - improvements to workflow automation and parallelization on single
    machines and computer clusters


## Developers

### Why vectorize?

When looking for ways to optimize _R_ code one of the most common suggestions
is to 'vectorize'. A 'vectorized' function in _R_ is one that has a
loop-like construct written in a compiled language with a light _R_ wrapper.
There are a couple of qualities that give these functions their performance
advantage.

Data are passed to a 'vectorized' function as a vector instead of individual
elements. Vectors in _R_ are 'typed' meaning all elements must be of the same
data type. Passing a whole vector to the complied code reduces the amount of
work _R_ has to do to interpret data type. Second, the loop computation is done
in a complied language such as C, C++ or FORTRAN.

Vectorized functions call .C, .Call, .Primitive or .Internal in the source
code. One example is `base::which()`

    > which
    function (x, arr.ind = FALSE, useNames = TRUE)
    {
        wh <- .Internal(which(x))
        if (arr.ind && !is.null(d <- dim(x)))
            arrayInd(wh, d, dimnames(x), useNames = useNames)
        else wh
    }
    <bytecode: 0x3699c90>
    <environment: namespace:base>

An concrete example of how these functions can improve performance is
this
[post on the support site](https://support.bioconductor.org/p/70432/#70545)
where Herv&eacute; helped the author of the
[ChIPseeker](https://bioconductor.org/packages/ChIPseeker) package
identify a bottleneck in `ChIPseeker:::getFirstHitIndex()`. The
solution was to replace a call to `sapply()` (i.e., _R_-level
iteration) with two vectorized functions, `duplicated()` and
`which()`. This elegant one-liner reduced the algorithm from quadratic
in time to linear in time. Nice!

Herv&eacute; took a classic approach to improving the code:
simplifying the original support site example (thanks!)  to one that
illustrated the problem in a timely fashion; stepping through the code
a line at a time and noticing that one function call took inordinately
long; understanding the small section of code that caused problems;
and identifying a vectorized alternative.

We may not all see such drastic improvements but take home is that
these functions are worth knowing about and implementing when
possible.

### The ellipsis

In _R_, the ellipsis (...) is used to pass a variable number of arguments to
a function. A common question from developers writing their own S4
generics is when (and why) to include the ellipsis in the function signature.
Herv&eacute; answered a recent
[post on
Bioc-devel](https://stat.ethz.ch/pipermail/bioc-devel/2015-September/008014.html)
about this topic. Main ideas are summarized here but for the full details
you'll want to read the post.

Defining generics and methods instead of ordinary functions for the getters and
setters of an S4 object is generally considered good practice. When introducing
a new generic, the recommendation is to keep the signature as "generic" as
possible so methods can add arguments that are specific to the objects they
deal with.

     setGeneric("parameters",
         function(object, ...) standardGeneric("parameters")
     )

Another case is that of having extra arguments such as a modifier or toggle
precede the ellipsis in the generic.

     setGeneric("enrichment",
         function(object, method="auto", ...)
             standardGeneric("enrichment"),
         signature="object"
     )

The purpose of `signature="object"` is to limit dispatch to `object`
only. Defining sensible default values (e.g., method="auto" or
verbose=FALSE) is also important to encourage consistent,
user-friendly behavior across methods.

Going the direction of more restrictive, an example where the ellipsis
was intentionally omitted is the `organism()` generic recently added
to
[BiocGenerics](https://bioconductor.org/packages/BiocGenerics). This
is a straightforward getter and we want to enforce this and encourage
methods to stick to that very simple contract. If someone comes up
with a use-case where extra arguments are needed for their method then
the ellipsis can be added to the generic.

While `...` is often appropriate in the signature of generics, methods
are often best written to include only the arguments actually used. In
this way, invalid arguments are not silently ignored:

    .A = setClass("A", representation(x="integer"))
    setMethod("enrichment", "A", function(object, method="auto", ...) {
        cat("enrichment:", class(object), "; method:", method, "\n")
    })
    enrichment(.A(), mehods="special")    # typo silently ignored
    ## enrichment: A ; method: auto

    setMethod("enrichment", "A", function(object, method="auto") {
        cat("enrichment:", class(object), "\nmethod:", method, "\n")
    })
    enrichment(.A(), mehods="special")
    ## Error in .local(object, method, ...) :
    ## unused argument (mehods = "special")

The key when designing generics is to identify or anticipate the
greatest common factor across methods and formalize that at the level
of the generic. When in doubt, better to underestimate than
overestimate. For methods, use the most restrictive signature
possible.


## Project Statistics

### Website traffic

The following compares the number of sessions and new users from the third
quarter of 2015 (July 1 - September 25) with the third quarter of 2015. Sessions
are broken down by new and returning visitors. New visitors correspond to the
total new users.

<table border="0" cellpadding="5" cellspacing="0">
 <caption><b>Third Quarter Website Traffic 2015 vs 2014</b></caption>
  <tbody valign="top">
    <tr>
        <td><b>Sessions: Total</b></td>
        <td>28.05% increase</td>
        <td>(339,991 vs 265,507)</td>
    </tr>
    <tr>
        <td><b>Sessions: Returning Visitor</b></td>
        <td>34.03% increase</td>
        <td>(227,171 vs 169,498)</td>
    </tr>
    <tr>
        <td><b>Sessions: New Visitor</b></td>
        <td>17.50% increase</td>
        <td>(112,820 vs 96,009)</td>
    </tr>
  </tbody>
</table>

<br/>
Statistics generated with [Google Analytics](http://www.google.com/analytics/).

### Package downloads and new submissions

The number of unique IP downloads of software packages for July, August and
September of 2015 were 35349, 33600 and 29953 respectively.  For the same time
period in 2014, numbers were 36362, 36118 and 47918. Numbers must be
compared by month (vs sum) because some IPs are the same between months.
See the web site for a full summary of [download
stats](http://bioconductor.org/packages/stats/).

As of September 25, a total of 54 software packages have been added in the
third quarter of 2015 bringing counts to 1078 in devel (_Bioconductor_ 3.2)
and 1024 in release (_Bioconductor_ 3.1).


## News, Events and Courses

See the [events page](http://www.bioconductor.org/help/events/) for a listing
of all courses and conferences.

* [European Bioconductor Developers Conference](https://sites.google.com/site/eurobioc2015/)
07 - 08 December 2015 — Cambridge, UK

* ["Bioconductor for Genomic Data Science" Coursera course](http://kasperdanielhansen.github.io/genbioconductor/)
Launched September 7 2015. This series of month-long courses is part of the
JHU Genomic Data Science Specialization. All 6 classes will run every month.

* [_Bioconductor_ 3.2 Release](http://www.bioconductor.org/developers/release-schedule/)
14 October 2015 - Worldwide!

* [A Short Introduction to Bioconductor](http://blog.revolutionanalytics.com/2015/08/a-short-introduction-to-bioconductor.html)
A brief summary of the project by Pete Hickey written for the Revolutions blog.


Send comments or questions to 
[Valerie](mailto:valerie.obenchain@roswellpark.org).


[AH]: //bioconductor.org/packages/AnnotationHub
[GR]: //bioconductor.org/packages/GenomicRanges
[CopywriteR]: //bioconductor.org/packages/CopywriteR
[systemPipeR]: //bioconductor.org/packages/systemPipeR
[ComplexHeatmap]: //bioconductor.org/packages/ComplexHeatmap
[derfinderHelper]: //bioconductor.org/packages/derfinderHelper
[ggtree]: //bioconductor.org/packages/ggtree
[RnBeads]: //bioconductor.org/packages/RnBeads
[cogena]: //bioconductor.org/packages/cogena
