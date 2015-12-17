<!-- %\VignetteEngine{knitr::knitr} %\VignetteIndexEntry{} -->

{::options parse_block_html="true" /}

# *Bioconductor* Newsletter
{:.no_toc}

posted by [Valerie Obenchain](mailto:vobencha@roswellpark.org), January 2016


## Contents 
{:.no_toc}

* Table of contents will replace this text. 
{:toc}

TODO: intro


## October release

`Bioconductor` 3.2 was release on October 14, consisting of 1104 software
packages, 257 experiment data packages, and 917 annotation packages. There are
80 new software packages.

This is the last version of `Bioconductor` to be supported on Snow Leopard. Snow
Leopard users should plan to migrate to Mavericks or newer before the next
release in Spring 2016.

There are 80 new software packages included in this release. Package summaries
and the official release schedule can be found on the 
[web site](http://www.bioconductor.org/news/bioc_3_2_release/).


## Design matrices for differential gene expression

Mike Love is an author of the
[DESeq2](http://www.bioconductor.org/packages/3.3/bioc/html/DESeq2.html)
package for differential gene expression of RNASeq data.  A visit to the
[support site](https://support.bioconductor.org/) shows the number of questions
he answers on a daily basis related not only to using the [DESeq2
package](http://www.bioconductor.org/packages/3.3/bioc/html/DESeq2.html) but
about analysis of gene expression data in general. In addition to supporting
[DESeq2](), he is a postdoc in [Rafael Irizarry's
lab](http://rafalab.dfci.harvard.edu/) in the Department of Biostatistics and
Computational Biology at the Dana Farber Cancer Institute and Harvard School of
Public Health where he develops quantitative methods for genomics and
epigenetics, teaches
[edX
courses](https://www.edx.org/course/data-analysis-life-sciences-1-statistics-harvardx-ph525-1x)
and and occasionally [blogs](https://mikelove.wordpress.com/) about statistics and
`R`.

Of the many DESeq2-related posts on the support site,
creating an appropriate design matrix is a regular, and one that causes a
fair bit of confusion. I asked Mike if he had observed any patterns in these
questions or could share any insights as to what users were struggling with.
Below he explains 

### What is 'study design'?

TODO:
Explain what study 'design' is eg, 'design of a study
explains the inter-relationship between the sample, essentially describes how
samples are distributed between groups' ...

and why it's important, eg, the design matrix affects (dictates?) how the
scientific question is asked/answered ...? I'm guessing the design matrix
carries with it assumptions that if not true or specified correctly will
invalidation results. Right?


### Case vs control

Simple designs don't seem to pose much issue. For example, control and treated
samples, or control, treatment 1 and treatment 2. These are easily modeled
using R's built-in `formula` and `model.matrix` functions, and then input to
[limma]()http://www.bioconductor.org/packages/3.3/bioc/html/limma.html), 
[edgeR](http://www.bioconductor.org/packages/3.3/bioc/html/edgeR.html), 
[DESeq2](http://www.bioconductor.org/packages/3.3/bioc/html/DESeq2.html) or other 
`Bioconductor` packages.

### Confounding and batch effects

Sometimes, quantitative/computational problems arise in the form of
error messages which indicate inherent problems in the experimental
design. One of these is when comparisons of interest are *confounded*
with technical factors, such as the sample preparation batch. There
is the canonical case of confounding when control and treatment
samples are prepared in their own batches, but I also see cases of bad
experiment design such as:

| condition | batch |
|:---------:|:-----:|
| control   | 1     |
| control   | 1     |
| treat. A  | 1     |
| treat. A  | 1     |
| treat. B  | 2     |
| treat. B  | 2     |
| treat. C  | 2     |
| treat. C  | 2     |

While treatment A can be compared against control, and treatment C can
be compared against treatment B, no comparisons can be made across the
batches. The reason some comparisons cannot be made is that the
difference in gene expression due to, for example the effect of
treatment B compared to control, cannot separated from the differences
which could arise between different sample preparation batches. The
most effective solution here is to use a blocked design, where the
batches each include all of the possible conditions. At the least,
control samples should be included in each batch, so that the batch
effect can be measured using these samples.

For more on why batch effects pose a big problem for high-throughput
experiments (or any experiments) here are two links:

* http://simplystatistics.org/2015/05/20/is-it-species-or-is-it-batch-they-are-confounded-so-we-cant-know/
* http://www.nature.com/nrg/journal/v11/n10/abs/nrg2825.html

### Blocking, interactions and nested designs

Blocked experimental designs, and others, such as those where the
significance of interactions between conditions is tested, or nested
interactions, can be read about in the excellent limma User's Guide,
in the section on  
[Single-Channel Experimental
Designs](https://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf)

The Guide describes in detail how the design matrix can be formulated in
different ways to answer the same question and explains how the different
parametrizations affect interpretation of the results.  The approaches
recommended by the limma authors can typically be applied as well to other
Bioconductor packages.

### Advanced designs 

Then there are some very complicated designs with myriad technical and
biological factors, where the investigator has many comparisons to
make and not a solid sense how to make them. In these cases I highly
recommend, for people who find themselves not knowing what
design to use or how to interpret the coefficients, that they consider
partnering with a local statistician or someone with a
background in linear modeling or quantitative analysis.

Interpreting quantitative analyses is hard stuff, and while
Bioconductor simplifies the analysis of high-throughput assays to a
large degree, and it's not necessarily reasonable to expect that
complicated results can be compiled or interpreted by someone without
a quantitative background. I like to compare this expectation
to a person not trained in laboratory procedures expecting to walk into a lab and
perform an experiment that would be complicated even for an experience
technician. I think it's safer and more reasonable to find a
collaborator who can assist, although it's preferable to include such
collaborators on projects from the beginning.

## Annotation Tour

### The primary packages

`Bioconductor` hosts many different types of annotation packages. This
section highlights the most heavily used and the most common applications.

* `OrgDb`:

  The OrgDb packages encapsulate all the information we know about a
  given organism's genes as of a given date. GO terms, ontology, Entrez IDs,
  RefSeq ID, Ensembl IDs and many others. Because these data have nothing
  to do with where a gene is found they are not related to a genome build.

* `ChipDb`:

  The `ChipDb` packages contain a single mapping: probe Id to Entrez gene ID.
  The Entrez gene ID is also found in the `OrgDb` package but for those
  interested in just this single mapping, the `ChipDb` is a lighter weight
  option.

* `TxDb`:

  `TxDb` packages contain location information of transcripts, genes,
  exons and other gene-related features for a specific organism
  based on a given build of the genome.

* `BSgenome`:

  These packages contain full genome sequences for a specific organism
  based on a given build of the genome.

* `SNPlocs`:

  SNP locations and alleles for a specific organism extracted from a
  particular dbSNP build which is based on genome build. 

* `AnnotationHub`:

  This package provides an interface to browse and download a wide collection
  of annotation packages and individual resources. Much of the data are
  pre-parsed into `R` / `Bioconductor` objects.


@Jim: I know this next chunk is long but it such excellent background I feel 
we must include it (or parts of it).

The NCBI databases are a hierarchy of sorts, where people submit sequences they
think were expressed in a particular species.  These sequences come from lots
of different groups, and start out as sort of provisional transcripts at NCBI
or Ensembl. As more people find the same things and evidence accrues for a
particular sequence being real, the provisional transcripts get collapsed into
a single sequence, and are given a RefSeq or GenBank ID (or maybe both).  Then
if a set of RefSeq or GenBank transcripts appear to be variants from a single
locus, they might be collapsed into a single UniGene ID.  But the general idea
is that a jumble of submitted (and predicted) transcripts are slowly collapsed
from a bunch of hypothetical gene like sequences, into a smaller set of
transcripts or DNA sequences that we think are 'for real'. 

This collapsing process goes on all the time. In addition, duplicates are
constantly being found in RefSeq or Gene or whatever, and one ID is deprecated
in favor of the other. RefSeq and GenBank have weekly releases so this is a
fast-moving target.

Genome builds, on the other hand, have to do with trying to piece
together the actual structure of the genome. Since the build has to do
with the entire genome, it isn't possible or reasonable to update
weekly, so they do releases every so often. So in mm9, GeneX may have
been thought to have been found on Chr1 at 1234567-1237654, but when
they did mm10, they may have shuffled things around so that gene may
now be in a different place on the genome. But that doesn't change the
sequence of the gene, nor what it does, nor what transcripts it is
thought to make, nor the gene ontology terms appended to it, or
anything else.


### Common tasks

One common task is to annotate a microarray experiment by mapping the
manufacturer's IDs to something more general, such as a HUGO gene
symbol, or an NCBI (Gene, GenBank, RefSeq, UniGene) or Ensembl (Ensembl
gene, Ensembl transcript) ID. As an example, we can map an Affymetrix
ID from the Human Gene 1.0 ST array to the corresponding HUGO symbol.

	library(hugene10sttranscriptcluster.db)
	hugene <- hugene10sttranscriptcluster.db ## minimize typing
	select(hugene, "8012257", "SYMBOL")

    'select()' returned 1:1 mapping between keys and columns
     PROBEID SYMBOL
	1 8012257   TP53

This is a very simple example, and probably not that useful, except
as a quick interactive query. Note that we did not specify the
`keytype` argument. The default `keytype` argument for `select` is the
central ID for the package being used. In this case, that is the
PROBEID, and since we are using probeset IDs as `keys`, it is not
necessary to specify the `keytype`.

A more common use case is to annotate a vector of `keys`, returning
one or more output IDs (or `columns`). As an example, we will use just
five `keys` from the hugene10sttranscriptcluster.db package, and query
for the HUGO symbol and Entrez Gene IDs. 

	ids <- keys(hugene)[15000:15005]
	ids
	[1] "8005171" "8005191" "8005200" "8005202" "8005204"
	
	annot <- c("SYMBOL","ENTREZID")
	select(hugene, ids, annot)
	'select()' returned 1:many mapping between keys and columns
	   PROBEID       SYMBOL  ENTREZID
	1  8005171        TRPV2     51393
	2  8005191  LRRC75A-AS1    125144
	3  8005191     SNORD49A     26800
	4  8005191     SNORD49B    692087
	5  8005191      SNORD65    692106
	6  8005200  LRRC75A-AS1    125144
	7  8005200     SNORD49A     26800
	8  8005200     SNORD49B    692087
	9  8005200      SNORD65    692106
	10 8005202     SNORD49A     26800
	11 8005202  LRRC75A-AS1    125144
	12 8005202     SNORD49B    692087
	13 8005202      SNORD65    692106
	14 8005204     CCDC144A      9720
	15 8005204    CCDC144CP    348254
	16 8005204     CCDC144B    284047
	17 8005204    CCDC144NL    339184
	18 8005204 LOC101929141 101929141
	19 8005221         <NA>      <NA>

Please note two things about the above results. First, the PROBEID
column in the returned `data.frame` has the same order as the input
ids. Second, some of the Affymetrix IDs map to more than one
gene. All of the mappings are returned, with a message that there
was a 1:many mapping for some of the `keys`. Because of the 1:many
mappings, the dimensions of the returned `data.frame` do not match the
dimensions of the data we would like to annotate (e.g., we wanted
information for five IDs, and got 19 rows of data returned).

If we want to guarantee that the returned data are in the same order
*and* are the same length as the input `keys` vector, we can use
`mapIds` instead. However, `mapIds` can only do one `keytype` at a
time, and returns a `vector` rather than a `data.frame`. Unlike
`select`, which has a default value for the keytype, `mapIds` requires
a fourth argument, specifying the `keytype` of the `keys` we are
using.

	> mapIds(hugene, ids, "SYMBOL", "PROBEID")
      8005171       8005191       8005200       8005202       8005204 
      "TRPV2" "LRRC75A-AS1" "LRRC75A-AS1"    "SNORD49A"    "CCDC144A" 
      8005221 
      NA

We can easily wrap this in a small script to return a `data.frame`
with just one row per `key`.

	> d.f <- as.data.frame(lapply(annot, function(x) mapIds(hugene, ids, x, "PROBEID")))
	> names(d.f) <- annot
	> d.f
		         SYMBOL ENTREZID
	8005171       TRPV2    51393
	8005191 LRRC75A-AS1   125144
	8005200 LRRC75A-AS1   125144
	8005202    SNORD49A    26800
	8005204    CCDC144A     9720
	8005221        <NA>     <NA>

The default for `mapIds` is to take the first instance for any 1:many
mappings. We can use the `multiVals` argument to control what is
returned. Please note that this argument comes after an ellipsis
(`...`) argument, so you cannot use positional arguments, and must
instead specify the `multiVals` argument directly.

	> mapIds(hugene, ids, "SYMBOL", "PROBEID", multiVals = "list")
	$`8005171`
	[1] "TRPV2"

	$`8005191`
	[1] "LRRC75A-AS1" "SNORD49A"    "SNORD49B"    "SNORD65"    

	$`8005200`
	[1] "LRRC75A-AS1" "SNORD49A"    "SNORD49B"    "SNORD65"    

	$`8005202`
	[1] "SNORD49A"    "LRRC75A-AS1" "SNORD49B"    "SNORD65"    

	$`8005204`
	[1] "CCDC144A"     "CCDC144CP"    "CCDC144B"     "CCDC144NL"    "LOC101929141"

	$`8005221`
	[1] NA

If we want to have a rectangular format for our annotation, where we
keep all the 1:many mappings while ensuring that each row maps
directly to our array of expression values, we can use a `DataFrame`
instead, telling `mapIds` to return a `CharacterList`.

	> lst <- lapply(annot, function(x)
	mapIds(hugene, ids, x, "PROBEID", multiVals = "CharacterList")
	> d.f <- as(lst, "DataFrame")
	> names(d.f) <- annot
	> d.f
	DataFrame with 6 rows and 2 columns
                                   SYMBOL                ENTREZID
                          <CharacterList>         <CharacterList>
	8005171                             TRPV2                   51393
	8005191 LRRC75A-AS1,SNORD49A,SNORD49B,... 125144,26800,692087,...
	8005200 LRRC75A-AS1,SNORD49A,SNORD49B,... 125144,26800,692087,...
	8005202 SNORD49A,LRRC75A-AS1,SNORD49B,... 26800,125144,692087,...
	8005204   CCDC144A,CCDC144CP,CCDC144B,...  9720,348254,284047,...
	8005221                                NA                      NA


Given the above data, perhaps we are interested in TRPV2,
and want to know its chromosomal location. We can use the Homo.sapiens
package to get that information. While it is possible to use the HUGO
symbol for this gene to get the location, it is a better idea to use
the Entrez Gene ID, which is more likely to be unique.

	> select(Homo.sapiens, "51393", c("TXCHROM","TXSTART","TXEND"), "SYMBOL")
	'select()' returned 1:1 mapping between keys and columns
	   ENTREZID TXCHROM  TXSTART    TXEND
	1  51393     chr17 16318856 16340317

This just tells us the start and stop positions for the transcript. If
we want exonic locations, we can get those as well.

	> select(Homo.sapiens, "TRPV2", c("EXONCHROM","EXONSTART","EXONEND"), "SYMBOL")
	'select()' returned 1:many mapping between keys and columns
	    SYMBOL EXONCHROM EXONSTART  EXONEND
	1   TRPV2     chr17  16318856 16319147
	2   TRPV2     chr17  16320876 16321182
	3   TRPV2     chr17  16323429 16323562
	4   TRPV2     chr17  16325913 16326203
	5   TRPV2     chr17  16326783 16327081
	6   TRPV2     chr17  16329413 16329583
	7   TRPV2     chr17  16330036 16330191
	8   TRPV2     chr17  16330763 16330861
	9   TRPV2     chr17  16331631 16331701
	10  TRPV2     chr17  16332131 16332296
	11  TRPV2     chr17  16335098 16335164
	12  TRPV2     chr17  16335280 16335614
	13  TRPV2     chr17  16336888 16337012
	14  TRPV2     chr17  16338204 16338283
	15  TRPV2     chr17  16340103 16340317
	16  TRPV2     chr17  16325913 16326207
	17  TRPV2     chr17  16330177 16330191
	18  TRPV2     chr17  16330766 16330861

While this is useful for a single gene, it can get unweildy for large
numbers of genes. We can instead use `transcriptsBy` or `exonsBy` with
the `TxDb.Hsapiens.UCSC.hg19.knownGene` package, to get information
about all genes at once, and subset to those we care about.

	> trscpts <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene")
	> trscpts[["51393"]]
	GRanges object with 2 ranges and 2 metadata columns:
        seqnames               ranges strand |     tx_id     tx_name
           <Rle>            <IRanges>  <Rle> | <integer> <character>
	[1]    chr17 [16318856, 16340317]      + |     60527  uc002gpy.3
	[2]    chr17 [16318856, 16340317]      + |     60528  uc002gpz.4

	> exns <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene")
	> exns[["51393"]]
	GRanges object with 18 ranges and 2 metadata columns:
         seqnames               ranges strand   |   exon_id   exon_name
            <Rle>            <IRanges>  <Rle>   | <integer> <character>
	 [1]    chr17 [16318856, 16319147]      +   |    217024        <NA>
	 [2]    chr17 [16320876, 16321182]      +   |    217025        <NA>
	 [3]    chr17 [16323429, 16323562]      +   |    217026        <NA>
	 [4]    chr17 [16325913, 16326203]      +   |    217027        <NA>
	 [5]    chr17 [16325913, 16326207]      +   |    217028        <NA>
	...      ...                  ...    ... ...       ...         ...
	[14]    chr17 [16335098, 16335164]      +   |    217037        <NA>
	[15]    chr17 [16335280, 16335614]      +   |    217038        <NA>
	[16]    chr17 [16336888, 16337012]      +   |    217039        <NA>
	[17]    chr17 [16338204, 16338283]      +   |    217040        <NA>
	[18]    chr17 [16340103, 16340317]      +   |    217041        <NA>
	-------
	seqinfo: 93 sequences (1 circular) from hg19 genome

Using `*Ranges` objects is beyond the scope of this newsletter. Please
see the
[IRanges](http://bioconductor.org/packages/release/bioc/vignettes/IRanges/inst/doc/IRangesOverview.pdf)
vignette, as well as the
[GRanges](http://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
vignettes for more information.

@Val - is this getting too long, or should I add to the advanced
tasks?

### Advanced tasks 

Then if we aren't getting too long, we could move from simple
task-based examples to use cases (I have an RNA-Seq experiment and I
want to generate a list of DE genes with annotation for my PI, etc).



## Reproducible Research

### Managing package versions with biocLite() 

`Bioconductor` follows a biannual schedule with one release in Spring and one
in Fall. `R` has a single release per year, usually in the Fall. Because each
`Bioconductor` release is tied to a version of `R` this asymmetrical schedule
creates some confusion. 

When releases coincide in the Fall, the development branches become release
branches. For the next 6 months, packages in the `Bioconductor`
'devel' branch are built against the 'devel' version of `R` and packages in the
'release' branch are built against the 'release' version of `R`.

In Spring, `Bioconductor` has a release but `R` does not. The `Bioconductor`
'devel' branch becomes the current 'release' and both branches are developed
against the 'release' version of `R`. The purpose of building `Bioconductor`
'devel' against `R` release is to allow for a smooth transition in Fall,
specifically, it allows the `Bioconductor` 'release' branch to always be in
sync with the `R` 'release' branch.

The [BiocInstaller]() package has several functions to help manage clean
'release' and 'devel' package repositories. Below are a few troubleshooting
tips for common install and version mis-match problems. 

* Confirm a single, correct version of `BiocInstaller`:

Make sure you have only a single installation directory defined by 

    .libPaths()

If multiple paths are reported, remove one.

Check the version of `BiocInstaller`:

    packageVersion("BiocInstaller")

The 'correct' version will depend on whether you are using the 'devel' or
'release' branch of `Bioconductor`. You can check the current version of
[BiocInstaller]() on the 
[release]()
and
[devel]() 
build pages.

If you have the wrong package version (or multiple versions) installed, remove
them with repeated calls to 

    remove.packages("BiocInstaller")

until `R` says there is no package to remove. Restart `R`, verify there is no
[BiocInstaller]() package and install the correct version with

    source("http://bioconductor.org/biocLite.R")

Invoking biocLite(), with no arguments, will update all packages.
When asked whether to update old packages, choose 'a' for 'all.

    biocLite()

* Identify mis-matched package versions with `biocValid()`:

Use biocValid() to identify version mis-matches between packages:

    BiocInstaller::biocValid()

Resolve by calling `remove.packages()` on the offending package, confirm the
correct version of [BiocInstaller]() and reinstall with `biocLite()`.

* TODO: BiocInstaller::biocLite("BiocUpgrade")
* TODO: BiocInstaller::useDevel()
* TODO: BiocInstaller::chooseBioCmirror

More information on keeping your versions in sync can be found at the 
[Why use biocLite()?](http://www.bioconductor.org/install/#why-biocLite) 
section of the web site.

### Package version tags in svn / git


## Infrastructure

### 'generics' package

### ongoing SummarizedExperiment development

### ExperimentHub

### InteractionSet class
ChiA-PET and Hi-C and InteractionSet class (https://github.com/LTLA/InteractionSet)


## New functions in R / Bioconductor

New functions added to `R` (3.3) and `Bioconductor` (3.3) this quarter:

*   *SummarizedExperiment::readKalisto()*

    TODO

*   *GenomicFeatures::mapToIds()* and *GenomicFeatures::mapToRanges()*

    TODO


## Project Statistics 

### Website traffic

The following compares the number of sessions and new users from the first
quarter of 2015 (January 1 - March 30) with the first quarter of 2014. Sessions
are broken down by new and returning visitors. New visitors correspond to the
total new users.

<table border="0" cellpadding="5" cellspacing="0">
 <caption><b>First Quarter Website Traffic 2015 vs 2014</b></caption>
  <tbody valign="top">
    <tr>
        <td><b>Sessions: Total</b></td>
        <td>24.03%</td>
        <td>(339,283 vs 273,559)</td>
    </tr>
    <tr>
        <td><b>Sessions: Returning Visitor</b></td>
        <td>21.42%</td>
        <td>(213,848 vs 176,128)</td>
    </tr>
    <tr>
        <td><b>Sessions: New Visitor</b></td>
        <td>28.74%</td>
        <td>(125,435 vs 97,431)</td>
    </tr>
    <tr>
        <td><b>New Users</b></td>
        <td>28.74%</td>
        <td>(125,435 vs 97,431)</td>
    </tr>
  </tbody>
</table>

<br/>
Statistics generated with [Google Analytics](http://www.google.com/analytics/).

### Package downloads and new submissions 

The number of unique IP downloads of software packages for January, February and
March of 2015 were 31720, 31956, and 38379, respectively.  For the same time
period in 2014, numbers were 29690, 28993 and 34634. Numbers must be
compared by month (vs sum) because some IPs are the same between months.
See the web site for a full summary of [download
stats](http://bioconductor.org/packages/stats/).

A total of 55 software packages were added in the first quarter of 2015 bringing counts to 991 in devel (`Bioconductor` 3.2) and 936 in release 
(`Bioconductor` 3.1).


## Upcoming Events

See the [events page](http://www.bioconductor.org/help/events/) for a listing
of all courses and conferences.

* [Use R / Bioconductor for Sequence Analysis](https://register.bioconductor.org/Seattle-Apr-2015/):
This intermediate level course is offered 06 - 07 April in Seattle, WA
USA.

* [Advanced RNA-Seq and ChiP-Seq Data Analysis](http://www.ebi.ac.uk/training/course/RNA2015):
Held in Hinxton, UK at EMBL-EBI, 11 - 14 May. 

* [CSAMA 2015 - Statistics and Computing in Genome Data Science](http://www-huber.embl.de/csama/):
Held the 14 - 19 of June in Bressanone-Brixen, Italy.

* [BioC2015](http://www.bioconductor.org/help/course-materials/2015/BioC2015/):
This year the 20 - 22 of July in Seattle, WA, USA.


Send comments or questions to Valerie at 
[vobencha@roswellpark.org](vobencha@roswellpark.org).
