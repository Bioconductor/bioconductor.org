![](/images/icons/help.gif)Using Bioconductor for Sequence Data
===============================================================

Bioconductor can import diverse sequence-related file types, including
fasta, fastq, ELAND, MAQ, BWA, Bowtie, SSOAP, BAM, gff, bed, and wig
files. Packages support common and advanced sequence manipulation
operations such as trimming, transformation, and alignment.
Domain-specific analyses include quality assessment, ChIP-seq,
differential expression, RNA-seq, and other approaches. Bioconductor
includes an interface to the Sequence Read Archive (via the
[SRAdb](/help/bioc-views/release/bioc/html/SRAdb.html) package).

* [Sample Workflow](#sample-workflow)  
* [Installation and Use](#install-and-use)
* [Exploring Package Content](#exploring-package-content)
* [Sequencing Resources](#sequencing-resources)

<h2 id="sample-workflow">Sample Workflow</h2>

The following psuedo-code illustrates a typical R / Bioconductor
session. It shows initial exploration of 454 resequencing of a 16S RNA
microbial community samples.

The workflow loads the `ShortRead` package and its dependencies. It
inputs about 250,000 reads of 200-250 bp each from a fastq
file. Flexible pattern matching (note the ambiguity letter `V')
removes a PCR primer artifact. The final lines plot the cumulative
number of trimmed reads as a function of their (log) abundance.

The result shows that most of the reads are from relatively few
sequences that each occur many times.

    ## Load packages; also loads Biostrings, IRanges, ...
    > library(ShortRead)
    > library(lattice) # for advanced plotting
    
    ## Input
    > seq <- readFastq("/path/to/file.fastq")
    
    ## Remove a PCR primer
    > pcrPrimer <- "GGACTACCVGGGTATCTAAT"
    > trimmed <- trimLRPatterns(Lpattern=pcrPrimer,
                                subject=sread(seq),
                                Lfixed="subject")
    
    ## Calculate and plot cumulative reads vs. occurrences
    > tbl <- tables(trimmed)[[2]]
    > xyplot(cumsum(nReads * nOccurrences) ~ nOccurrences, tbl, 
             scales=list(x=list(log=TRUE)), type="b", pch=20,
             xlab="Number of Occurrences", 
             ylab="Cumulative Number of Reads")

![cumulative reads](cumulative-reads.png)

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="install-and-use">Installation and Use</h2>

Follow [installation instructions](/install/) to start using these
packages.  To install the `ShortRead` package and all of its
dependencies, evaluate the commands

    > source("http://bioconductor.org/biocLite.R")
    > biocLite("ShortRead")

Package installation is required only once per R installation. View a
full list of
[available packages](/help/bioc-views/release/bioc/).

To use the `ShortRead` package, evaluate the command

    > library("ShortRead")

This instruction is required once in each R session.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="exploring-package-content">Exploring Package Content</h2>

Packages have extensive help pages, and include vignettes highlighting
common use cases. The help pages and vignettes are available from
within R. After loading a package, use syntax like

    > help(package="ShortRead")
    > ?readFastq

to obtain an overview of help on the `ShortRead` package, and the
`readFastq` function, and

    > browseVignettes(package="ShortRead")

to view vignettes (providing a more comprehensive introduction to
package functionality) in the `ShortRead` package. Use

    > help.start()

to open a web page containing comprehensive help resources.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="sequencing-resources">Sequencing Resources</h2>

The following packages illustrate the diversity of functionality
available; all are in the release version of Bioconductor.

* [IRanges](/help/bioc-views/release/bioc/html/IRanges.html),
  [GenomicRanges](/help/bioc-views/release/bioc/html/GenomicRanges.html)
  and
  [genomeIntervals](/help/bioc-views/release/bioc/html/genomeIntervals.html)
  for range-based (e.g., chromosomal regions) calculation, data
  manipulation, and general-purpose data
  representation. [Biostrings](/help/bioc-views/release/bioc/html/Biostrings.html)
  for alignment, pattern matching (e.g., primer removal), and data
  manipulation of large biological sequences or sets of
  sequences.

* [ShortRead](/help/bioc-views/release/bioc/html/ShortRead.html)
  and
  [Rsamtools](/help/bioc-views/release/bioc/html/Rsamtools.html)
  for file I/O, quality assessment, and high-level, general purpose
  data summary.
  [rtracklayer](/help/bioc-views/release/bioc/html/rtracklayer.html)
  for import and export of tracks on the UCSC genome browser.

* [BSgenome](/help/bioc-views/release/bioc/html/BSgenome.html)
  for accessing and manipulating curated whole-genome representations.
  [GenomicFeatures](/help/bioc-views/release/bioc/html/GenomicFeatures.html)
  for annotation of sequence features across common genomes,
  [biomaRt](/help/bioc-views/release/bioc/html/biomaRt.html)
  for access to Biomart databases.

* [SRAdb](/help/bioc-views/release/bioc/html/SRAdb.html)
  for querying and retrieving data from the Sequence Read Archive.

* ChIP-seq and related (e.g., motif discovery, identification of
  high-coverage segments) activities are facilitated by packages such
  as
  [CSAR](/help/bioc-views/release/bioc/html/CSAR.html),
  [chipseq](/help/bioc-views/release/bioc/html/chipseq.html),
  [ChIPseqR](/help/bioc-views/release/bioc/html/ChIPseqR.html),
  [ChIPsim](/help/bioc-views/release/bioc/html/ChIPsim.html),
  [ChIPpeakAnno](/help/bioc-views/release/bioc/html/ChIPpeakAnno.html),
  [rGADEM](/help/bioc-views/release/bioc/html/rGADEM.html),
  [segmentSeq](/help/bioc-views/release/bioc/html/segmentSeq.html),
  [segmentSeq](/help/bioc-views/release/bioc/html/segmentSeq.html).

* Differential expression and RNA-seq style analysis can be
  accomplished with
  [Genominator](/help/bioc-views/release/bioc/html/Genominator.html),
  [edgeR](/help/bioc-views/release/bioc/html/edgeR.html),
  [baySeq](/help/bioc-views/release/bioc/html/baySeq.html),
  [DESeq](/help/bioc-views/release/bioc/html/DESeq.html),
  and
  [DEGseq](/help/bioc-views/release/bioc/html/DEGseq.html).

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>
