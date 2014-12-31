<!-- %\VignetteEngine{knitr::knitr} %\VignetteIndexEntry{} -->

<!-- To deploy this to the website, do the following: library(knitr)
render_markdown(strict=TRUE) knit("2015_January.Rmd")

This will produce an .md file that can be used  on the website. The {:toc} and
{:no_toc} tags will be rendered by nanoc, the tool that produces the website.
-->

{::options parse_block_html="true" /}

# *Bioconductor* Newsletter *DRAFT*
{:.no_toc}

posted by [Valerie Obenchain](mailto:vobencha@fhcrc.org), January 2015

## Contents 
{:.no_toc}

* Table of contents will replace this text. 
{:toc}


## Software Infrastructure

### NCList

One of Herv&eacute;'s recent projects was to find a replacement for the interval
tree-based algorithm used in finding and counting overlaps.  The change was
motivated by performance problems especially when the comparison involved many
different chromosomes or many overlapping ranges. He decided on an approach 
based on the Nested Containment List algorithm by Alexander V. Alekseyenko and
Christopher J. Lee. The change has been implemented in BioC 3.1 (current devel). 

In BioC 3.1 overlaps operations on GRanges and/or GRangesList objects are
approximately 3x to 10x faster than in BioC 3.0 for a medium size data set 
(e.g. 25 million reads). Memory usage is also reduced by ~ 25% or more.
Numbers will vary depending on the size of the data; the larger the data the
greater the improvement.

The user-visible change is the 'algorithm' argument added to most overlap-based 
operations. This allows a choice between the new (algorithm="nclist") and the 
old (algorithm="intervaltree") algorithm. Except for the 3 situations below, 
choosing one algorithm or the other does not affect the output, only the 
performance.

Situations where 'nclist' and 'intervaltree' produce different output:

* With 'nclist', zero-width ranges are interpreted as insertion points and are
  considered to overlap with ranges that contain them.

* When using 'select="arbitrary"', the 'nclist' algorithm will generally
  not select the same hits as the old algorithm.

* When using 'intervaltree', the 'maxgap' argument has a special
  meaning if 'type' is "start", "end", or "within". This is not yet
  the case with the new 'nclist'.

For a complete description of changes and future activity please see this
[post](https://stat.ethz.ch/pipermail/bioc-devel/2014-December/006749.html)
on the bioc-devel mailing list.

### Docker

Anyone who has attended the annual BioC conference or participated in a
class at the Hutchinson Center is familiar with the `Bioconductor` Amazon 
Machine Images (AMI). Dan configures these images with the current
version of `R` and all necessary package dependencies and non-R software.
Providing a pre-configured environment to participants has eliminated 
many "day of" problems such as internet overload due to concurrent downloads, 
installation of the wrong package version and the inability to access common 
sample data used in the course.

Recently Dan has been looking into the Docker software as another approach to
providing pre-configured environments. The Docker containers operate at the
single process level and are therefore more lightweight than a virtual machine.
Additionally the containers have access to the full physical machine.

These containers are isolated, reproducible environments useful for development
or production. In contrast to the AMI, a Docker container could be used on a
desktop or laptop instead of the cloud or other remote hardware. We envision
them being useful for 

* reproducibility: providing identical, reproducible environments

* convenience: deploying on any (Linux, Mac, or Windows) machine, or in 
  the cloud via Amazon's Elastic Container Service (ECS)

* time saver; instead of waiting for packages to compile on Linux, you can 
  download a container in which packages are already installed

Examples of potential applications are available at the 
[bioc_docker](https://github.com/Bioconductor/bioc_docker) GitHub repository. 


## Coordinate Mapping

Translation (mapping) of coordinates between genome, transcript and protein
space is a common task in bioinformatics. Over the past months we've been
working to expand and harmonize the mapping capabilities in our infrastructure.
High-level functions for mapping between genome and transcript space, 
`mapToGenome` and `mapToTranscript`, have been added to `GenomicRanges` and
`GenomicAlignments` in the devel branch.  Element-wise (aka "parallel") versions
of the functions, `pmapToGenome` and `pmapToTranscript`, map the i-th element
of 'x' to the i-th element of 'alignment'. Still on the TODO is a function for
protein space.

Others involved in the project are  Michael, Herve, Laurent, Robert and
Martin. Discussions and progress can be followed at the
[biocCoordinateMapping](https://groups.google.com/forum/#!forum/bioccoordinatemapping) Google Group.

A related mapping task is migrating data from one assembly to another
either via alignment tool or by converting assembly coordinates. The 
implementation of the UCSC LifOver tool in `rtracklayer` 
(thanks Michael) and the availability of the UCSC chain files in 
`AnnotationHub` (thanks Sonali) make this a straightforward operation. 

    library(AnnotationHub)
    hub <- AnnotationHub()

The chain file format describes a pairwise alignment that allows gaps in both 
sequences simultaneously. `AnnotationHub` currently hosts 1113 chain files.

    allChains <- query(hub, 'chain')

    ## > length(allChains)
    ## [1] 1113

Search for the conversion from hg38 to hg19:

    query(hub, 'hg38ToHg19')

    ## > query(hub, 'hg38ToHg19')
    ## class: AnnotationHub 
    ## hub: https://annotationhub.bioconductor.org 
    ## cache: /home/vobencha/.AnnotationHub 
    ## display()ing 1 of 1 records on 6 mcols()
    ##                            title            dataprovider      species
    ## AH14108 hg38ToHg19.over.chain.gz hgdownload.cse.ucsc.edu Homo sapiens
    ##         taxonomyid genome                                description
    ## AH14108       9606   hg38 UCSC liftOver chain file from hg38 to hg19
    ##                                            tags
    ## AH14108 liftOver, chain, UCSC, genome, homology

Methods in `AnnotationHub` use `import.chain()` from `rtracklayer` to read 
data into `R`. The return object is a `Chain` class where data for each 
chromosome are parsed into a `ChainBlock` class.

    chain <- query(hub, 'hg38ToHg19')[[1]] 
    chain

    ## > chain
    ## Chain of length 25
    ## names(25): chr1 chr2 chr3 chr4 chr5 chr6 ... chr20 chr21 chr22 chrX 
    ## chrY chrM

    class(chain$chr6)

    ## > class(chain$chr6)
    ## [1] "ChainBlock"
    ## attr(,"package")
    ## [1] "rtracklayer"

`liftOver()` translates coordinates and outputs a `GRangesList` .

    gr <- GRanges(c("chr7", "chr2"), IRanges(c(75625897, 68010781), width=1)) 
    res <- liftOver(gr, chain)

    ## > res
    ## GRangesList object of length 2:
    ## $1 
    ## GRanges object with 1 range and 0 metadata columns:
    ##       seqnames               ranges strand
    ##          <Rle>            <IRanges>  <Rle>
    ##   [1]     chr7 [75255215, 75255215]      *
    ## 
    ## $2 
    ## GRanges object with 1 range and 0 metadata columns:
    ##       seqnames               ranges strand
    ##   [1]     chr2 [68237913, 68237913]      *


Another example of using `liftOver()` is in the
[Changing genomic coordinate systems](http://www.bioconductor.org/help/workflows/liftOver/) workflow. SNPS from the NHGRI GWAS catalogue are mapped from hg38 
to hg19 resulting in a few lost loci.


## Overview of the `csaw` package

Gordon Smyth is a professor at Walter and Eliza Hall Institute of Medical 
Research in Victoria and has been active member in the `Bioconductor` project 
since inception. Many know Gordon from his detailed responses on the
[support site](https://support.bioconductor.org/) where he provides
statistical guidance and thoughtful discussion to the novice and advanced
user alike. The Smyth group has authored a number of
`Bioconductor` packages including cornerstone contributions such as `limma` and
`edgeR`. `limma` is consistently in the top 10 and `edgeR` in the
top 25 of the 
[package download stats](http://bioconductor.org/packages/stats/index.html).

A trademark of a Smyth group package is a well-written vignette with detailed
scientific and statistical background. These documents are an excellent starting
point for anyone new to microarray or RNA-seq analysis. The most recent 
contribution from the group is the `csaw` package (ChIP-seq analysis with 
windows) by Aaron Lun. I asked Aaron and Gordon if they would answer a few 
questions about the package and they graciously agreed.


Aaron Lun and Gordon Smyth on the csaw package:

**Q: What is differential binding (DB) in ChIP-seq and how does this
differ from differential expression (DE) in RNA-seq?**

As most readers will know, ChIP-seq sequences genomic DNA that is bound to a
target protein whereas RNA-seq sequences RNA transcripts. From a scientific
point of view, RNA-seq measures gene expression whereas ChIP-seq is used to
examine the regulation of gene expression. ChIP-seq is often used for example
to find the binding sites of a transcription factor (TF) or to examine the
positioning of an epigenetic histone mark across the genome.  Many people
analyze ChIP-seq data by calling peaks in each individual ChIP-seq library.
These peaks are then taken to identify where the TF binds or where epigenetic
marking is active. In the csaw package, however, we envisage ChIP-seq
experiments with multiple experimental conditions and with multiple biological
replicates within each condition, in other words with a structure much like a
typical gene expression experiment. We focus on testing for DB in a quantitative
way between the experimental conditions rather than on making present/absent
calls.  ChIP-seq and RNA-seq experiments both produce short sequence reads that
can be aligned to a reference genome. In principle, the two types of experiments
can be analyzed in broadly similar ways. In each case, we can choose a set of
genomic locations of interest, count the number of reads mapping to those
regions, and then test for differential coverage between the experimental
conditions relative to biological variability. The key difference between
ChIP-seq and RNA-seq is how the genomic locations are chosen and combined.

**Q: What scientific question does a DB ChIP-seq analysis help answer?**

A conventional ChIP-seq analysis produces a list of peaks in each library, with
the implication that the protein is bound to DNA within the peak regions and not
bound elsewhere. We take an alternative view that binding coverage is
quantitative, not just present or absent.  DB analyses identify a list of
genomic regions where the strength of binding changes between biological
conditions. We believe that the regions where the binding is subject to change
are the most likely to be biologically important, even if they are not
necessarily the regions with highest coverage. DB regions are necessarily linked
with the phenotype or treatment of interest being studied in the experiment. We
can therefore start investigating mechanisms through which DB can lead to
observed differences in biology. One of our favorite analyses is to relate DB to
differential expression for the same set of genes. These insights are harder to
obtain with the conventional analyses, as the identified regions are only
associated with each condition in isolation.

**Q: What was the motivation for creating the package? Has your group recently
started doing DB ChIP-seq or have you been doing it for some time?**

Our group has been doing DB analyses for a while. However, these analyses have
mostly been gene-orientated. We would count reads into pre-defined intervals
like the promoters or gene bodies and then test for whether these intervals were
DB using the edgeR package. This type of analysis can be done very similarly for
ChIP-seq as for RNA-seq. The csaw package is our first attempt at performing de
novo DB analyses where the regions of interest are not known in advance. We
wanted to be able to identify novel DNA elements like new enhancers or promoters
and we wanted to avoid potential headaches regarding the regions of interest,
for example misspecifying the promoter region of a gene.  We were motivated to
take a new approach to de novo DB analysis by our observation that some commonly
used approaches do not control the error rates correctly. A simple method is to
call peaks in each condition separately, then simply compare the regions to
identify unique peaks in each condition. This ad hoc approach does not provide
any statistical error rate control and tends to overstate the differences
between the conditions. A more sophisticated version of this approach is to
conduct statistical tests between conditions for each of the regions called as
peaks in either condition. Again, this over-estimates the differences between
the conditions.  Our approach gives similar flexibility to peak calling but with
rigorous error rate control. The idea is to slide windows across the genome,
count reads into those windows, and then use those counts to test for
significant differences between conditions. Adjacent regions are then merged and
the p-values combined in a statistically rigorous way. This provides the same
level of statistical rigor as for our previous gene-orientated analysis but
without the need to specify regions of interest beforehand. Using small windows
also provides excellent spatial resolution.  There are a number of other Bioc
packages that can do DB analyses, diffBind and DBChIP for example, but these
require a set of peaks identified with external software like MACS. The
motivation for the csaw package is to avoid having to specific the genomic
regions externally . We wanted to generate and discover the DB completely de
novo as an integral part of the analysis. csaw is the first Bioc package to take
a windowing approach.

**Q: How are the statistical methods from edgeR leveraged or extended in
csaw? Please describe from a statistical perspective, why methods used in
RNA-seq are appropriate for DB ChIP-seq?**

Once a set of genomic regions has been chosen, a DB analysis is very similar to
an RNA-seq experiment. So it was logical to leverage the statistical methods
provided by edgeR. edgeR accounts for biological variability between replicates,
which is critical for DB analyses where the counts for each region are typically
over-dispersed. Failure to account for such variability will lead to spurious DB
calls. The generalized linear model functionality of edgeR provides a rich
framework with which to analyze complex experiments with multiple experimental
factors or covariate and batch effects.  We chose edgeR over limma and voom
because the ChIP-seq counts for each window can be quite small. The reads can be
sparsely distributed across the whole genome. limma and voom are very effective
for moderate to large counts, but edgeR is able to more accurately model the
distribution of very small integer counts. csaw uses the quasi-likelihood
functionality of edgeR because of its rigorous error rate control and because it
gives access to the adaptive empirical Bayes functionality of limma within the
edgeR negative binomial framework.  We made many improvements to edgeR as part
of the csaw project. Our sliding window approach to ChIP-seq generates a much
larger number of regions than is typical for an RNA-seq DE analysis, so it was
important  to be as computationally efficient as possible. One of us (AL)
converted many of the edgeR functions into C++ for speed and memory efficiency.
There are a number of statistical extensions specific to csaw. These include an
implementation of Simes' method for combining the p-values of adjacent windows
within a region; a non-linear normalization procedure adapted to low counts; and
a method to calculate the average abundance of a region scaled to its width. 


## Project Statistics 

### Publications and Citations

This past year a number of `Bioconductor` core and community members contributed
to the manuscript, Orchestrating high-throughput genomic analysis with 
`Bioconductor`, scheduled to appear in Nature Methods early 2015. The article
provides an overview of the project for potential users and developers
highlighting reproducibility and flexibility in current applications.

Other recent project-level (vs package-level) publications are Scalable Genomics
with `R` and `Bioconductor` by Michael Lawrence and Martin, and a review of `R`
and `Bioconductor` as applied to genomics, oceanography and ecology in by Sylvia
Tippman. Links to all articles are available on the 
[publications page](http://www.bioconductor.org/help/publications/).

A PubMed query for title or PubMed ID returned 22838 citation references
for software packages with CITATION files. For packages with no CITATION, a
PubMed title search returned 1281 references.

Full text citations of `Bioconductor` are available on the web site from 
[PubMed](http://www.ncbi.nlm.nih.gov/pubmed/?term=bioconductor&sort=ePubDate),
[PubMedCentral](http://www.ncbi.nlm.nih.gov/pmc/?term=bioconductor&sort=ePubDate)and 
[Google Scholar](http://scholar.google.com/scholar?q=bioconductor&btnG=Search).

### Website traffic

The following tables compare traffic from the fourth quarter of 2014 with 
the fourth quarter of 2013 (October 1 - December 28).

In the fourth quarter we saw an increase in total sessions, new sessions and
the number of total users.

<table>
 <caption>Website traffic Q4 2014 vs Q4 2013</caption>
  <tr>
    <th></th><th></th>
  </tr>
  <tr>
   <td width="55%">Sessions</td><td><b>23.28%</b> (311,731 vs 252,873)</td>
  </tr>
  <tr>
   <td width="55%">% New Sessions</td><td><b>0.62%</b> (36.01% vs 35.79%)</td>
  </tr>
  <tr>
   <td width="55%">Users</td><td><b>24.32%</b> (133,839 vs 107,655)</td>
  </tr>
</table>
<p></p>

The greatest increase (percent change) in total sessions was seen in China
followed by Spain then the United States and Italy.

<table>
 <caption>Total Sessions by Location Q4 2014 vs Q4 2013</caption>
  <tr>
    <th></th><th></th>
  </tr>
  <tr>
    <td width="55%">United States</td><td><b>24.95%</b> (101011 vs 80840)</td>
  </tr>
  <tr>
    <td width="55%">China</td><td><b>28.43%</b> (27593 vs 21485)</td>
  </tr>
  <tr>
    <td width="55%">United Kingdom</td><td><b>11.19%</b> (22627 vs 20350)</td>
  </tr>
  <tr>
    <td width="55%">Germany</td><td><b>20.75%</b> (21138 vs 17506)</td>
  </tr>
  <tr>
    <td width="55%">France</td><td><b>20.47%</b> (10446 vs 8671)</td>
  </tr>
  <tr>
   <td width="55%">Canada</td><td><b>21.96%</b> (9808 vs 8042)</td>
  </tr>
  <tr>
   <td width="55%">Japan</td><td><b>19.21%</b> (9406 vs 7890)</td>
  </tr>
  <tr>
   <td width="55%">Spain</td><td><b>27.30%</b> (8444 vs 6633)</td>
  </tr>
  <tr>
   <td width="55%">India</td><td><b>4.15%</b> (8048 vs 7727)</td>
  </tr>
  <tr>
   <td width="55%">Italy</td><td><b>24.68%</b> (7494 vs 6010)</td>
  </tr>
</table>
<p></p>

Statistics were generated with Google Analytics.

### Package downloads and new submissions 

There was a 9% increase in the number of (distinct IP) downloads of software 
packages in the fourth quarter of 2015 (106212 total) as compared to the fourth
quarter of 2014 (97462 total). See the website for a full summary of package
[download stats](http://www.bioconductor.org/packages/stats/).

During the period of October 1 to December 31 a total of 63 software packages
were added bringing the current counts to 954 packages in devel (`Bioconductor`
3.1) and 934 in release (`Bioconductor` 3.0).

## Resources and Upcoming Events

### Videos and Webinars

On average, 20-30 software packages are submitted to `Bioconductor` each 
quarter resulting in 100+ new packages each year. New package guidelines are
posted on the
[website](http://www.bioconductor.org/developers/package-guidelines/) but 
developers often have additional questions or special case situations. We
thought a more interactive exchange might help avoid common errors and
encourage questions before a package was submitted for review.

In mid December Marc and Sonail hosted a Google Hangout to share some key 
development tips and solicit questions from a wider audience. Topics covered
were package organization, documentation and code reuse. Questions were posted 
via a YouTube chat window and answered after Marc finished the slide 
presentation. If you are considering submitting a package you may want to
check out the Webinar posted on the `Bioconductor` 
[video page](https://www.youtube.com/watch?v=QfqaK_BHebU).

### Conferences and Courses 

The [events page](http://www.bioconductor.org/help/events/) is updated regularly 
with new courses and conference announcements. In January 2015, EMBL is hosting 
both the European Developer's conference and an Advanced `R` Programming course.

[European Bioconductor Developer's
Meeting](http://www-huber.embl.de/BiocEurope/)
European Molecular Biology Laboratory
Heidelberg, Germany
January 12-13, 2015

[Advanced R Programming and
Development](http://www.dataprogrammers.net/embl_jan2015/)
European Molecular Biology Laboratory
Heidelberg, Germany
January 15-16, 2015

Send comments or questions to Valerie at 
[vobencha@fhcrc.org](vobencha@fhcrc.org).
