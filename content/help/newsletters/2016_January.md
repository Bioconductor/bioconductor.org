<!-- %\VignetteEngine{knitr::knitr} %\VignetteIndexEntry{} -->

{::options parse_block_html="true" /}

# *Bioconductor* Newsletter
{:.no_toc}

posted by [Valerie Obenchain](mailto:vobencha@roswellpark.org), January 2016

The _Bioconductor_ newsletter is a quarterly review of core infrastructure
developments, community projects and future directions. Topics are of general
interest as well as those with the greatest impact on the software.  This
quarter has seen substantial development on the `ExperimentHub` resource and the
`InteractionSet` class. We review some tips for managing package repositories
with `biocLite()` and introduce the new version tagging in the svn / git repos
which makes it possible to retrieve a specific version of a _Bioconductor_
package. Mike Love talks about constructing design matrices for gene expression
experiments and Jim MacDonald takes us on tour of _Bioconductor_ annotation
packages.

## <a name="Contents"></a> Contents 
{:.no_toc}

* Table of contents will replace this text. 
{:toc}

## F1000 Research Support Prize

At the 
[European _Bioconductor_ Developers
meeting](https://sites.google.com/site/eurobioc2015/) last December, a prize
was awarded to the individual(s) with the greatest contribution to the
_Bioconductor_ [support site](https://support.bioconductor.org/) forum. 

The prize was sponsored by F1000 research which recently recently launched a
dedicated
[_Bioconductor_ channel](http://f1000research.com/channels/bioconductor).
The [terms](https://support.bioconductor.org/p/75500/) of the award were 
'greatest contribution to the support site' and 'those attending the European 
developer conference'.

Congratulations to winners Aaron Lun and Michael Love! Each were awarded the
prize of waived publication costs for an article appearing in the F1000
_Bioconductor_ channel. Other contributors with substantial posts to their credit
are Jim MacDonald, Gordon Smyth, Ryan Thompson and Steve Lianoglou. Thanks to
everyone who takes the time to answer questions and share their experience on
the support site. 

Thanks to Mark Dunning and Laurent Gatto for suggesting the
prize (and organizing the conference!) and to Thomas Ingraham and F1000 Research
for sponsoring it.

[back to top](#Contents)

## October release

_Bioconductor_ 3.2 was release on October 14, consisting of 1104 software
packages, 257 experiment data packages, and 917 annotation packages. There are
80 new software packages.

This is the last version of _Bioconductor_ to be supported on Snow Leopard. Snow
Leopard users should plan to migrate to Mavericks or newer before the next
release in Spring 2016.

There are 80 new software packages included in this release. Package summaries
and the official release schedule can be found on the 
[web site](http://www.bioconductor.org/news/bioc_3_2_release/).

[back to top](#Contents)


## Design matrices for differential gene expression

Mike Love is a postdoc in [Rafael Irizarry's
lab](http://rafalab.dfci.harvard.edu/) in the Department of Biostatistics and
Computational Biology at the Dana Farber Cancer Institute and Harvard T.H. Chan
School of Public Health. He develops quantitative methods for genomics and
epigenetics, teaches [edX
courses](https://www.edx.org/course/data-analysis-life-sciences-1-statistics-harvardx-ph525-1x)
and occasionally [tweets](http://twitter.com/mikelove) about biostatistics and
_R_. Many know him as the author and primary supporter of the very popular
[DESeq2](http://www.bioconductor.org/packages/3.3/bioc/html/DESeq2.html)
package for differential gene expression of RNASeq data.

A visit to the
[support site](https://support.bioconductor.org/) shows the number of questions
he answers on a daily basis related not only to the [DESeq2
package](http://www.bioconductor.org/packages/3.3/bioc/html/DESeq2.html) but
about gene expression analysis in general. 

Of the many DESeq2-related posts on the support site, creating an appropriate
design matrix is a regular one and appears to cause a fair bit of confusion.
Below Mike shares some of his observations and thoughts about what key concepts
cause the most problems.

### A little background about 'experimental design'

**Experimental design** refers to the inter-relationships between samples,
including the biological and experiment information (clone 1, treatment B,
etc.) and technical information (the batch in which the sample was prepared and
processed). Getting the experimental design correct -- this happens before the
experiment takes place -- is very important as we will see below, because the
wrong experimental design can lead to uninterpretable data.

It is best practice to keep track of the experimental design (including
preparation batches) in a table, either a CSV or TSV file, or an Excel
spreadsheet, which can be exported to CSV when it comes time to do
bioinformatic analysis.  This is the easiest way to explain the experimental
design to someone for either a planned experiment or an experiment that has
already taken place.

**Design matrix** or **model matrix** is a matrix, typically represented in
statistics by an $X$, that will be used in statistical modeling.  Every row in
the design matrix describes a sample, and every column provides pieces of
information about that sample, such as, whether the sample was treated, whether
the sample was in batch 1, 2, or 3, etc.  For every column of the design
matrix, the model has a matching coefficient, usually denoted by $\beta$'s, to
describe differences in expression across samples.  These coefficients are
additive differences on the log scale, so multiplicative differences (fold
changes) in RNA-seq counts or microarray expression values, hence they are *log
fold changes*. These coefficients are then estimated using the experimental
data.

A single experimental design does not imply a single design matrix, but there
are often many choices involved. For example, one could include a coefficient
for technical batches (typically a good idea) or not, which would give a design
matrix with an addition column. The design matrix encodes assumptions the
investigator wants to make regarding the samples, and is formed with respect to
the biological question of interest. During significance testing, the
biological question is phrased as a *null hypothesis*, that one or more of the
coefficients are equal to zero, resulting in a *p value*.  The p value is a
meaningful estimate of a probability only if the assumptions are reasonable and
the model is well specified for the data.

I should say, I learned about these topics both from textbooks (John Rice's
*Mathematical Statistics and Data Analysis*, Sanford Weisberg's *Applied Linear
Regression*, and [Bioconductor
books](http://www.bioconductor.org/help/publications/#books)), as well as from
reading lots of posts on the _Bioconductor_ support forum from people like
Wolfgang Huber, Gordon Smyth, Simon Anders, James MacDonald, Aaron Lun and
others.

### Case vs control

Simple designs don't seem to pose much issue. For example, control and treated
samples, or control, treatment 1 and treatment 2. These are easily modeled
using _R_'s built-in `formula` and `model.matrix` functions, and then input to
[limma](http://www.bioconductor.org/packages/3.3/bioc/html/limma.html),
[edgeR](http://www.bioconductor.org/packages/3.3/bioc/html/edgeR.html),
or other _Bioconductor_ packages.
[DESeq2](http://www.bioconductor.org/packages/3.3/bioc/html/DESeq2.html)
directly takes `formula` expressions and converts to design matrices internally.

### Confounding and batch effects

Sometimes, quantitative/computational problems arise in the form of error
messages which indicate inherent problems in the experimental design. One of
these is when comparisons of interest are *confounded* with technical factors,
such as the sample preparation batch. There is the canonical case of
confounding when control and treatment samples are prepared in their own
batches, but also common are cases of bad experimental design such as:

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

While treatment A can be compared against control, and treatment C can be
compared against treatment B, no comparisons can be made across the batches.
The reason some comparisons cannot be made is that the difference in gene
expression due to, for example the effect of treatment B compared to control,
cannot separated from the differences which could arise between different
sample preparation batches. The most effective solution here is to use a block
design, where the batches each include all of the possible conditions. At the
least, control samples should be included in each batch, so that the batch
effect can be estimated using these samples.

These two links explain why batch effects pose a big problem for
high-throughput experiments (or any experiments):

* http://simplystatistics.org/2015/05/20/is-it-species-or-is-it-batch-they-are-confounded-so-we-cant-know/
* http://www.nature.com/nrg/journal/v11/n10/abs/nrg2825.html

### Blocking, interactions and nested designs

Block experimental designs, and others, such as those where the significance of
interactions between conditions is tested, or nested interactions, can be read
about in the excellent limma User's Guide, in the section on [Single-Channel
Experimental
Designs](https://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf)

The Guide describes in detail how the design matrix can be formulated in
different ways to answer the same question and explains how the different
parametrizations affect interpretation of the results.  The approaches
recommended by the limma authors can be applied to other _Bioconductor_
packages as well.

### Advanced designs 

Then there are some very complicated designs with many technical and biological
factors, where the investigator has many comparisons to make and not a solid
sense how to make them. In these cases I highly recommend, for people who find
themselves not knowing what design to use or how to interpret the coefficients,
that they consider partnering with a local statistician or someone with a
background in linear modeling or quantitative analysis.

Interpreting quantitative analyses is hard stuff, and while _Bioconductor_
simplifies the analysis of high-throughput assays to a large degree, it's not
necessarily reasonable to expect that complicated results can be compiled or
interpreted by someone without a quantitative background.  I think it's safer
and more reasonable to find a collaborator who can assist, and such
collaborators can help identify issues with experimental design if they are
included on projects from the outset.

[back to top](#Contents)


## Getting started with _Bioconductor_ annotation packages

Jim MacDonald is a biostatistician at the University of Washington Department
of Environmental and Occupational Health Sciences. He has analyzed the gamut of
HTS data from expression (microarray, RNA-Seq), to genomics (SNP arrays, DNA-Seq,
ChIP-Seq, methylation arrays, BS-Seq) and other 'omics' data. He has been heavily
involved in the direction of the _Bioconductor_ project since inception and
contributes and maintains a large number of annotation packages.

During the October 2015 release we were short-handed after loosing staff to the
Buffalo move. Jim stepped in and took responsibility for building all 
internal _Bioconductor_ annotation packages. Jim's comprehensive understanding 
of the annotation world is evident in his numerous posts on the 
[suport site](https://support.bioconductor.org/). In this section we've teamed
up (90% Jim, 10% Val) to give an overview of key packages and how they can be
used to answer some common analysis questions.

### The primary packages

This section highlights the most heavily used _Bioconductor_ annotation
packages.

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
  pre-parsed into _R_ / _Bioconductor_ objects.

It's worth noting that some annotation packages are tied to specific genome
builds and others are not. The `TxDb` family contain the location of
genes/transcripts/exons/etc. based on a given build. `BSgenome` and `SNPlocs`
are other examples of build-specific packages. Because genome assembly requires
piecing together the structure of the whole genome it follows that new builds
are only released every few years.  Data in the build-specific annotation
packages can be quite stable and stay current for years.

Other packages have nothing to do with where a gene is found and are therefore
not related to a genome build, e.g., the `OrgDb` family. These packages can be
thought of as encapsulating all information we have about the genes of a given
organism on a given date, knowing that it can become obsolete, at least
in part, the very next week. They contain such information as RefSeq, GenBank,
or UniGene IDs which represent provisional transcripts. These are a work in
progress and constantly being updated and modified based on public submissions.
_Bioconductor_ updates these packages every 6 months at realse time. The 
`AnnotationForge` package offers functions to build your own `OrgDb` 
(or other package) if you want something more current.

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

[back to top](#Contents)

## Reproducible Research

### Managing package versions with biocLite() 

_Bioconductor_ follows a biannual schedule with one release in Spring and one
in Fall. _R_ has a single release per year, usually in the Fall. Because each
_Bioconductor_ release is tied to a version of _R_ this asymmetrical schedule
creates some confusion. 

When releases coincide in the Fall, the development branches become release
branches. For the next 6 months, packages in the _Bioconductor_
'devel' branch are built against the 'devel' version of _R_ and packages in the
'release' branch are built against the 'release' version of _R_.

In Spring, _Bioconductor_ has a release but _R_ does not. The _Bioconductor_
'devel' branch becomes the current 'release' and both branches are developed
against the 'release' version of `R`. The purpose of building _Bioconductor_
'devel' against _R_ release is to allow for a smooth transition in Fall,
specifically, it allows the _Bioconductor_ 'release' branch to always be in
sync with the _R_ 'release' branch.

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
'release' branch of _Bioconductor_. You can check the current version of
`BiocInstaller` on the 
[devel]http://www.bioconductor.org/checkResults/devel/bioc-LATEST/()
and
[release](http://www.bioconductor.org/checkResults/release/bioc-LATEST/) 
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

* Upgrade to the most recent _Bioconductor_ for a version of _R_:

When _Bioconductor_ has a release but _R_ does not, the current _R_ release
supports both the release and devel versions of _Bioconductor_. You can upgrade
to the most current _Bioconductor_ (devel) with

    BiocInstaller::biocLite("BiocUpgrade")

This installs the most recent _Bioconductor_ packages without having to
reinstall _R_.

More information on keeping your versions in sync can be found at the 
[Why use biocLite()?](http://www.bioconductor.org/install/#why-biocLite) 
section of the web site.

[back to top](#Contents)


## `InteractionSet`

Aaron Lun, Liz Ing-Simmons and Malcolm Perry have been working on an
[InteractionSet](https://github.com/LTLA/InteractionSet) package to store
and manipulate data from ChIA-PET and Hi-C experiments.

ChIA-PET stands for Chromatin Interaction Analysis with Paired-End Tags. These
experiments probe for genome-wide interactions brought about or associated with
some protein of interest. An essential step in this technology that
differentiates it from Hi-C is the antibody-driven immunoprecipitation step
to enrich for chromatin bound by a specific protein. Chromatin interactions can
only be determined for parts of the genome that have a binding site for the
protein of interest. Interaction networks can be elucidated for transcription
factors, insulator proteins or transcription machinery. A ChIA-PET experiment
gives information about the potential role of proteins in structuring 3D genome
organization. 

In contrast, the Hi-C method provides information about 3D genome structure by
identifying long range chromatin interactions on a genome-wide scale. These
data are often used to study aspects of genome architecture such as chromosome
territories, topological domains, open/closed compartments and chromatin
structure.

Data from both technologies enable the study of physical interactions
between pairs of genomic regions. The InteractionSet package provides classes
to represent these interactions and store associated experimental data. The
aim is to provide package developers with stable class definitions that can be
manipulated through a large set of methods.

The package defines the following classes:

* `GInteractions` : represents pairwise interactions between genomic regions
* `InteractionSet`: contains experimental data relevant to each interaction
* `ContactMatrix` : stores a data matrix where each row and column represent
                    a genomic region.

The classes have methods for sorting and duplicate detection; for performing
one- or two-dimensional overlaps; for calculating distances between interacting
loci; and for calculating the minimum bounding box for groups of interactions.
Methods are also available to convert between classes, or to standard
_Bioconductor_ objects like a `RangedSummarizedExperiment` or `GRangesList`.


[back to top](#Contents)

## New functions_

New functions added to _R_ (3.3) and _Bioconductor_ (3.3) this quarter:

*   *SummarizedExperiment::readKalisto()*

    Reads kalisto data into a SummarizedExperiment.

    (contributed by Martin Morgan)

*   *GenomicFeatures::mapToIds()* and *GenomicFeatures::mapToRanges()*

    Maps between genomic identifiers (gene names, symbols etc.) and
    the genomic ranges they represent.

    (contributed by Jim Hester)

[back to top](#Contents)

## Project Statistics 

### Website traffic

The following compares the number of sessions and new users from the fourth
quarter of 2015 (November 1 - December 28) with the fourth quarter of 2014.

The fourth quarter of 2015 saw an increase over 2014 in the
number of Sessions (18.05%), Users (11.27%) and Pageviews (12.11%). Decreases
were seen in the Average number of pages viewed (-5.39%), Average
session duration (-0.54%), Bounce rate (aka single page visits; -3.30%) and
Percent of new sessions (-8.86%).

The majority of users are still on desktops (and laptops) but the number of
mobile users is steadily increasing. In the fourth quarter 2015, the number of
new users increased by 39.8% for mobile devices vs 6.7% for desktops.

<table border="0" cellpadding="5" cellspacing="0">
 <caption><b>Website Traffic: New Users by Device</b></caption>
  <tbody valign="top">
    <tr>
    <td><b>Desktop (includes laptops)</b></td>
    <tr>
        <td>Nov 1, 2015 - Dec 27, 2015</td>
        <td>71,111 (92.86%)</td>
    <tr>
        <td>Nov 1, 2014 - Dec 27, 2014</td>
        <td>66,622 (93.96%)</td>
    <tr>
        <td>% change</td>
        <td>6.74%</td>
    </tr>
    <tr>
        <td><b>Mobile</b></td>
    <tr>
        <td>Nov 1, 2015 - Dec 27, 2015</td>
        <td>4,230 (5.52%)</td>
    <tr>
        <td>Nov 1, 2014 - Dec 27, 2014</td>
        <td>3,026 (4.27%)</td>
    <tr>
        <td>% change</td>
        <td>39.79%</td>
    </tr>
    <tr>
        <td><b>Tablet</b></td>
    <tr>
        <td>Nov 1, 2015 - Dec 27, 2015</td>
        <td>1,241 (1.62%)</td>
    <tr>
        <td>Nov 1, 2014 - Dec 27, 2014</td>
        <td>1,260 (1.78%)</td>
    <tr>
        <td>% change</td>
        <td>-1.51%</td>
    </tr>
  </tbody>
</table>

<br/>


Statistics generated with [Google Analytics](http://www.google.com/analytics/).

[back to top](#Contents)

### Package downloads and new submissions 

The number of unique IP downloads of software packages for October, November and
December of 2015 were 40085, 41499 and 31946, respectively.  For the same time
period in 2014, numbers were **TODO**. Numbers must be
compared by month (vs sum) because some IPs are the same between months.
See the web site for a full summary of [download
stats](http://bioconductor.org/packages/stats/).

A total of 23 software packages were added in the fourth quarter of 2015
bringing counts to 1101 in devel (_Bioconductor_ 3.3) and 1104 in release
(_Bioconductor_ 3.2).

[back to top](#Contents)

## Upcoming Events

See the [events page](http://www.bioconductor.org/help/events/) for a listing
of all courses and conferences.

* [CSAMA 2016 (14th edition) - Statistics and Computing in Genome Data Science](http://www-huber.embl.de/csama/):
10-15 of July in Bressanone-Brixen, Italy.

## Acknowledgements 

Thanks to Jim MacDonald and Mike Love for contributing sections, Aaron Lun for
proofing the `InteractionSet` section and the _Bioconductor_ core team for
editorial review.

[back to top](#Contents)

Send comments or questions to Valerie at 
[vobencha@roswellpark.org](vobencha@roswellpark.org).
