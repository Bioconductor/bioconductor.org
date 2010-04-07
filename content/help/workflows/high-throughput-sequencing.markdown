Using Bioconductor for Sequence Data
====================================

Bioconductor can input diverse sequence-related file types, including
fasta, fastq, MAQ, BWA, Bowtie, BAM, gff, bed, and wig files. Packages
support common and advanced sequence manipulation operations such as
trimming, transformation, and alignment.  Domain-specific analyses
include quality assessment, ChIP-seq, differential expression,
RNA-seq, and other approaches. Bioconductor includes an interface to
the Sequence Read Archive.


## Sample Work Flow ##

The following psuedo-code illustrates a typical R / Bioconductor
session. It shows initial exploration of 454 resequencing of a 16S
RNA microbial community samples.

The work flow loads the ShortRead package and its dependencies. It
inputs about 250,000 reads of 200-250 bp each from a fastq
file. Flexible pattern matching (note the ambiguity letter `V')
removes a PCR primer artefact. The final lines plot the cumulative
number of trimmed reads as a function of their (log) abundance.

The result shows that most of the reads are from relatively few
sequences that each occur many times.

    ## Load packages; also loads Biostrings, IRanges, ...
    library(ShortRead)
    library(lattice) # for advanced plotting
    
    ## Input
    seq <- readFastq("/path/to/file.fastq")
    
    ## Remove a PCR primer
    pcrPrimer <- "GGACTACCVGGGTATCTAAT"
    trimmed <- trimLRPatterns(pcrPrimer, subject=sread(seq))
    
    ## Calculate and plot cumulative reads vs. occurrences
    tbl <- tables(trimmed)[[2]]
    xyplot(cumsum(nReads * nOccurrences) ~ nOccurrences, tbl, 
          scales=list(x=list(log=TRUE)), type="b", pch=20,
          xlab="Number of Occurrences", 
          ylab="Cumulative Number of Reads")


![cumulative reads](/images/help/workflows/cumulative-reads.png)


## Installation ##

Follow [installation instructions]("/install/"") to start using these
packages.  To install the `ShortRead` package and all of
its dependencies, evaluate the commands

    > source("http://bioconductor.org/biocLite.R")
    > biocLite("ShortRead")

Package installation is required only once per R installation. To use
the `ShortRead` package, evaluate the command

    > library("ShortRead")

This command is required once in each R session.


## Exploring Package Content ##

Packages have extensive help pages, and include vignettes highlighting
common use cases; the help pages and vignettes are available from
within R. After loading a package, use syntax like

    > help(package="ShortRead")
    > ?readFastq

to obtain an overview of help on the `ShortRead` package,
and the `readFastq` function, and

    > browseVignettes(package="ShortRead")

to view vignettes (providing a more comprehensive introduction to
package functionality) in the `ShortRead` package. Use

    > help.start()

to open a web page containing comprehensive help resources.


## Sequencing Resources ##

The following packages illustrate the diversity of functionality
  available; all are in the release version of Bioconductor.

<ul>

  <li><a href="http://bioconductor.org/packages/release/bioc/html/IRanges.html">IRanges</a>,
    <a href="http://bioconductor.org/packages/release/bioc/html/GenomicRanges.html">GenomicRanges</a>
    and <a href="http://bioconductor.org/packages/release/bioc/html/genomeIntervals.html">genomeIntervals</a>
    for range-based (e.g., chromosomal regions) calculation, data
    manipulation, and general-purpose data representation.
    <a href="http://bioconductor.org/packages/release/bioc/html/Biostrings.html">Biostrings</a>
    for alignment, pattern matching (e.g., primer removal), and data
    manipulation of large biological sequences or sets of
    sequences.</li>

  <li><a href="http://bioconductor.org/packages/release/bioc/html/ShortRead.html">ShortRead</a> and
    <a href="http://bioconductor.org/packages/release/bioc/html/Rsamtools.html">Rsamtools</a>
    for file I/O, quality assessment, and high-level, general purpose
    data summary.
    <a href="http://bioconductor.org/packages/release/bioc/html/rtracklayer.html">rtracklayer</a>
    for import and export of tracks on the UCSC genome browser.</li>

  <li><a href="http://bioconductor.org/packages/release/bioc/html/BSgenome.html">BSgenome</a>
    for accessing and manipulating curated whole-genome
    representations.
    <a href="http://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html">GenomicFeatures</a>
    for annotation of sequence features across common genomes,
    <a href="http://bioconductor.org/packages/release/bioc/html/biomaRt.html">biomaRt</a>
    for access to Biomart data bases.</li>
    
  <li>
    <a href="http://bioconductor.org/packages/release/bioc/html/SRAdb.html">SRAdb</a>
    for querying and retrieving data from the Sequence Read Archive.</li>
    
    
  <li>ChIP-seq and related (e.g., motif discovery, identification of
    high-coverage segments) activities are facilitated by packages such as
    <a href="http://bioconductor.org/packages/release/bioc/html/CSAR.html">CSAR</a>,
    <a href="http://bioconductor.org/packages/release/bioc/html/chipseq.html">chipseq</a>,
    <a href="http://bioconductor.org/packages/release/bioc/html/ChIPseqR.html">ChIPseqR</a>,
    <a href="http://bioconductor.org/packages/release/bioc/html/ChIPsim.html">ChIPsim</a>,
    <a href="http://bioconductor.org/packages/release/bioc/html/ChIPpeakAnno.html">ChIPpeakAnno</a>,
    <a href="http://bioconductor.org/packages/release/bioc/html/rGADEM.html">rGADEM</a>,
    <a href="http://bioconductor.org/packages/release/bioc/html/segmentSeq.html">segmentSeq</a>,
    <a href="http://bioconductor.org/packages/release/bioc/html/segmentSeq.html">segmentSeq</a>.</li>

  <li>Differential expression and RNA-seq style analysis can be accomplished with
    <a href="http://bioconductor.org/packages/release/bioc/html/Genominator.html">Genominator</a>,
    <a href="http://bioconductor.org/packages/release/bioc/html/edgeR.html">edgeR</a>,
    <a href="http://bioconductor.org/packages/release/bioc/html/baySeq.html">baySeq</a>,
    <a href="http://bioconductor.org/packages/release/bioc/html/DESeg.html">DESeg</a>, and
    <a href="http://bioconductor.org/packages/release/bioc/html/DEGseq.html">DEGseq</a>.</li>

</ul>
