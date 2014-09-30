<!-- %\VignetteEngine{knitr::knitr} %\VignetteIndexEntry{} -->

<!-- To deploy this to the website, do the following: library(knitr)
render_markdown(strict=TRUE) knit("2014_October.Rmd")

This will produce an .md file that can be used  on the website. The {:toc} and
{:no_toc} tags will be rendered by nanoc, the tool that produces the website.
-->

{::options parse_block_html="true" /}

# *Bioconductor* Newsletter  
{:.no_toc}

posted by [Valerie Obenchain](mailto:vobencha@fhcrc.org), October 2014

## Contents 
{:.no_toc}

* Table of contents will replace this text. 
{:toc}


## Software Infrastructure

### GRCh38 assembly

The GRCh38 human genome assembly is available in `Bioconductor` as
[BSgenome](http://www.bioconductor.org/packages/devel/data/annotation/html/BSgenome.Hsapiens.NCBI.GRCh38.html)
[TranscriptDb](http://www.bioconductor.org/packages/devel/data/annotation/html/TxDb.Hsapiens.UCSC.hg38.knownGene.html)
and
[SNPloc](http://www.bioconductor.org/packages/devel/data/annotation/html/SNPlocs.Hsapiens.dbSNP141.GRCh38.html)
packages.

The GRCh38 assembly includes both a primary assembly (non-redundant haploid
assembly) and alternate sequences (alt loci). Alt loci are provided for regions
of the genome where variation prevents representation by a single sequence.
These regions are not new but have become more prominent as tools for variant
detection have matured.

The previous GRCh37 assembly included patch releases tagged as &lsquo;fix&rsquo; 
or &lsquo;novel&rsquo;. The &lsquo;fix&rsquo; patches were incorporated in the 
primary assembly of GRCh38 while the &lsquo;novel&rsquo; patches were moved 
into the alt loci units. The &lquo;multi-sequence&rsquo; nature of GRCh38 
raises questions about how to best work with these alternate sequences
with respect to alignment and downstream analysis.

### htslib

The samtools library and associated sub-tools play and integral role in the
analysis of HTS data. The htslib is the successor of libbam which is currently
provided by samtools. Specifically, htslib is a C library for handling 
high-throughput sequencing data, providing APIs for manipulating SAM, BAM, and 
CRAM sequence files (similar to but more flexible than the old Samtools API) 
and for manipulating VCF or BCF variant files.

An implementation of htslib is in the works for `Bioconductor` and will likely 
be implemented as stand-alone package. Follow the development on Martin's 
[GitHub site](https://github.com/mtmorgan/Rhtslib).

### `S4Vectors` and `IRanges` split completed

In September Herv&eacute; completed the move of non-range based code from
`IRanges` to `S4Vectors`. The virtual `Vector` and `List` classes moved as
well as `DataFrame`, `Rle` and `Hits`. Developers using or building on these
classes should now import from `S4Vectors`.


## `Bioconductor` support site 

In September the `Bioconductor` mailing list was replaced with a fork of
Biostars and renamed the 
[Bioconductor Support Site](https://support.bioconductor.org/). This affects 
the bioconductor list only; bioc-devel remains unchanged (bioc-devel@r-project.org).

The move was motivated by the volume of list traffic which highlighted
the need for advanced searching, tagging and real-time editing of posts. 
Ideally the new interface will encourage participation from first-time users 
and simplify topic management.

Marc has imported the last 11+ years of posts to create continuity in the
new environment. A [FAQ](https://support.bioconductor.org/info/faq/) is
available to help with site navigation and common tasks such as posting,
merging, or tracking topics. 

Thanks to Marc and Dan for their work on this.


## Developer's corner 

### Style markdown documents with `BiocStyle`

The `BiocStyle` package provides a fast and easy approach to styling markdown
documents in `Bioconductor` fashion. It includes all standard formatting 
styles for creating PDF and HTML documents of vignettes, workflows 
or other project documents. 
The [html](http://www.bioconductor.org/packages/3.0/bioc/vignettes/BiocStyle/inst/doc/HtmlStyle.html)
version of the package vignette is a demo of the styling and color theme.

The package offers formatting advantages over standard markdown such as
automatic centering of figures, improved table display and Latex-compatible math
symbols. Custom style sheets can be included by wrapping them in `
BiocStyle:::markdown`:

  BiocStyle::markdown(css.files = c('my.css'))

### Reuse and recycle: The power of `import`

The `Bioconductor` infrastructure contains a wealth of tools for HTS analysis.
Because these methods and containers exist across a number of packages they can
be difficult for new users to discover and for developers to remember when
adding new functionality.

The `import` generic in rtracklayer is one such tool. `import` reads and parses
large file formats such as BED, BAM, BigWig, GFF, Fasta, and Chain files. The
methods operate on *File objects (e.g., `BamFile`) and param (e.g.,
`ScanBamParam`) objects, which allow flexible control over the parsing and
subsetting of data. Data returned from `import` are parsed into useful
downstream containers such as `GRanges`  or `Rle`.

`import` should be the tool of choice when interacting with large files, such as
those available in `AnnotationHub`, or when developing new reading/parsing
functions.

### biocMultiAssay

During Developer Day at BioC 2014, Levi Waldron's discussion of his
biocMultiAssay project generated a good deal of interest in the community. This
effort, lead by Levi, Vince, Kasper and Martin, aims to create `Bioconductor`
tools for the efficient manipulation and analysis of multi-assay omics
experiments.

The primary motivation is to combine data across multiple experiments for a
common group of samples or patients. Goals are to develop classes and methods
for the extraction of data subsets defined by indices such as genomic position
or gene ID, and to streamline analyses that span multiple genomic data types.
The data are high-dimensional assays such as gene and protein expression, copy
number, methylation,  somatic mutation, or microRNA.

The project has a both a
[GitHub](https://github.com/vjcitn/biocMultiAssay) site with code prototypes
and a
[Google Group](https://groups.google.com/forum/#!forum/biocmultiassay) for
conference call announcements and tracking progress.


## Interview with Dr. Janet Young, FHCRC

This section of the newsletter highlights the work of an individual or group in
the `Bioconductor` community. This month we spoke to Janet Young from the  Fred
Hutchinson Cancer Research Center. Janet is originally from the UK with an
undergraduate degree in Natural Sciences from the University of Cambridge, and a
PhD in Genetics from University College London. She is currently a Staff
Scientist in the Malik lab in the Basic Sciences Division.

**Q: To begin would you tell us a bit about yourself?**

I joined Fred Hutch in 2000 and worked in the Trask lab first as a post-doc then
as a staff scientist. My own research focused on the evolution and
transcriptional regulation of mammalian olfactory receptor gene families, but I
also helped others with projects to measure genomic copy-number gain and loss in
prostate cancer and measurement of methylation levels in healthy human tissues.
When Barb (Trask) retired I spent time in the Tapscott lab where we studied how
transposable elements might be involved in a form of muscular dystrophy.
Currently I provide bioinformatics support to a variety of projects in the Malik
lab. The group studies evolutionary biology and genetic conflict, primarily in
drosophila, primates and yeast.

**Q: How did you get started with `Bioconductor`?**

I started working with `Bioconductor` when helping others in the Trask and
Tapscott labs with various microarray projects. Initially I used `R` /
`Bioconductor` simply for creating diagnostic plots of microarray data, but soon
started using limma and lumi for the analysis steps.

**Q: How does `Bioconductor` fit into your current workflows?**

I'm largely using it for analysis of deep sequencing data these days. We use a
variety of upstream software such as TopHat, BWA, and GATK. I use `Bioconductor`
for things like differential expression analysis, comparing coverage to look for
genomic copy number changes, filtering SNPs, or retrieving and analyzing gene/
annotations. Often I use rtracklayer to export the data for viewing in IGV or
the UCSC genome browser. As well as being a great analysis tool itself,
`Bioconductor` acts as the glue to help me integrate results from other tools.

**Q: Are there any `Bioconductor` resources you find particularly useful?**

The local classes offered at the Hutch were very helpful. I also like the
responsive Q and A on the mailing list. All software has bugs; knowing that the
bugs get fixed in a timely manner makes you keep using it. The package vignettes 
are a valuable 'stand alone' resource that help get you going with a specific
package or task right away.

**Thanks for talking with us and sharing your insights.**


## Quarterly Project Statistics 

The `Bioconductor` project continues to expand globally. Over the next quarter 
there are [course offerings](http://www.bioconductor.org/help/events/) in 
Japan, Germany the UK and US. In August 2014, the Latin American Bioconductor 
[LAB](http://lab.foundation/) foundation held its official inauguration
in Ribeirao Preto, Brazil. LAB is a non-profit scientific initiative created 
to represent and expand `Bioconductor` to the research community in Latin America
and is headed up by Benilton Carvalho and Houtan Noushmehr.

### Website traffic

Google analytics reports the following new visitors to the website
for the period of July 1 to September 28, 2014:

<table>
 <caption> New vs Returning Users </caption>
  <tr>
   <th></th><th>Sessions</th>
  </tr>
  <tr>
   <td width="55%">Returning Visitor</td><td><b>179,242</b> (63.81%)</td>
  </tr>
  <tr>
   <td width="55%">New Visitor</td><td><b>101,668</b> (36.19%)</td>
  </tr>
</table>
<p></p>

Overall website traffic by country:

<table>
 <caption>Website Visits by Location</caption>
  <tr>
    <th></th><th>Sessions</th>
  </tr>
  <tr>
    <td width="55%">Americas</td><td><b>118,398</b> (42.15%)</td>
  </tr>
  <tr>
    <td width="55%">Europe</td><td><b>93,973</b> (33.45%)</td>
  </tr>
  <tr>
    <td width="55%">Asia</td><td><b>57,744</b> (20.56%)</td>
  </tr>
  <tr>
    <td width="55%">Oceania</td><td><b>8,015</b> (2.85%)</td>
  </tr>
  <tr>
    <td width="55%">Africa</td><td><b>2,403</b> (0.86%)</td>
  </tr>
  <tr>
   <td width="55%">Other</td><td><b>377</b> (0.13%)</td>
  </tr>
</table>


### Package downloads

The number of distinct IP downloads of `Bioconductor` software packages  for
July, August and September were 36900, 36749, and 36618 respectively for an
average of 36756. A full summary of package download stats is available
[here](http://www.bioconductor.org/packages/stats/).


## Resources, Courses and Conferences

### Search `Bioconductor` materials by topic

Materials from past courses and conferences have long been available on the
`Bioconductor` web site categorized by conference name and date. At BioC 2014
this year we had several requests for a more refined search of these
materials by topic area or key word.

In response, Sonali and Dan have categorized all 2014 materials and implemented
a new key word(s) [search table](http://www.bioconductor.org/help/course-materials/)
interface. The plan is to index all future materials while years prior to 2014
will be available in the old format (see &lsquo;Courses by year&rsquo; below 
the search table).

### Publications 

If you are looking for resources to enhance your knowledge of working with
genomic ranges and sequences in `Bioconductor` the following publications maybe
of  interest:

[Software for Computing and Annotating GenomicRanges](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3738458/)  
This manuscript describes data structures available in the `Bioconductor`
infrastructure for representing and annotating ranges on the genome. Focus is
on the `IRanges`, `GenomicRanges` and `GeomicFeatures` packages which provide
support for transcript structures, read alignments and coverage vectors.

[Scalable Genomics with R and Bioconductor](http://arxiv.org/abs/1409.2864)
Strategies for analyzing large genomic data are described and implemented in `R`
and `Bioconductor`. Topics include scalable processing, summarization and
visualization.

### Bioconductor 3.0 release

The release of `Bioconductor` 3.0 is scheduled for October 14. This version
will continue to use the current version of `R` (3.1.1). Visit the website for 
help [updating packages](http://www.bioconductor.org/install/) and for a 
look at the 
[release schedule](http://www.bioconductor.org/developers/release-schedule/).

### Upcoming events

[Practical Course on Analysis of High-Throughput Sequencing Data](https://www.ebi.ac.uk/training/course/HTS2014)
EBI, Hinxton, UK   
October 20-25, 2014

[Learning R / Bioconductor for Sequence Analysis](https://register.bioconductor.org/Seattle-Oct-2014/)
FHCRC, Seattle WA, USA
October 27-29, 2014

[BioC Europe 2015](http://www.bioconductor.org/help/course-materials/2015/BiocEurope2015/)
EMBL, Heidelberg, Germany
January 12-15, 2015


Please send comments or questions to Valerie at 
[vobencha@fhcrc.org](vobencha@fhcrc.org).
