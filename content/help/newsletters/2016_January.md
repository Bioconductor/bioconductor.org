<!-- %\VignetteEngine{knitr::knitr} %\VignetteIndexEntry{} -->

{::options parse_block_html="true" /}

# *Bioconductor* Newsletter
{:.no_toc}

posted by [Valerie Obenchain](mailto:vobencha@roswellpark.org), January 2016


## Contents 
{:.no_toc}

* Table of contents will replace this text. 
{:toc}


## Annotations (need better title)

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
constantly being found in RefSeq or Gene or whatever, and one ID is deprecate
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


Given the above data, perhaps we are interested in LRRC75A-AS1,
and want to know its chromosomal location. We can use the Homo.sapiens
package to get that information, starting with the 



* convert an Affy ID to Entrez Gene ID, 
* find the genomic region for a gene, etc. 
Cover some of the inherent difficulties of annotating things, like 1:many 
mappings, and the tradeoffs you have to make.

### Advanced tasks 

Then if we aren't getting too long, we could move from simple
task-based examples to use cases (I have an RNA-Seq experiment and I
want to generate a list of DE genes with annotation for my PI, etc).



## ExperimentHub

## Reproducible Research
- manage your install with
  biocLite(), biocValid(), chooseBioCmirror(), biocinstallRepos() ...
- locate specific version with (if done) Jim's tag svn / git repos 


## Infrastructure
- 'generics' package
- ongoing SummarizedExperiment development
- ExperimentHub
- ChiA-PET and Hi-C and InteractionSet class (https://github.com/LTLA/InteractionSet)


## New and Noteworthy

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
