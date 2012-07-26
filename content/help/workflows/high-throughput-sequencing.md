![](/images/icons/help.gif)Using Bioconductor for Sequence Data
===============================================================

Bioconductor can import diverse sequence-related file types, including
fasta, fastq, BAM, gff, bed, and wig files, among others. Packages
support common and advanced sequence manipulation operations such as
trimming, transformation, and alignment.  Domain-specific analyses
include quality assessment, ChIP-seq, differential expression,
RNA-seq, and other approaches. Bioconductor includes an interface to
the Sequence Read Archive (via the
[SRAdb](/packages/release/bioc/html/SRAdb.html) package).

* [Sample Workflow](#sample-workflow)  
* [Installation and Use](#install-and-use)
* [Exploring Package Content](#exploring-package-content)
* [Sequencing Resources](#sequencing-resources)

<h2 id="sample-workflow">Sample Workflow</h2>

This is a simple work flow for single-end RNA-seq looking at
differential representation of known genes. The data is in a directory
`dataDir`, as implied by the commands below (<...> indicates
information to be provide by the user).

    dataDir <- <...>
    fastqDir <- file.path(dataDir, "fastq")
    bamDir <- file.path(dataDir, "bam")
    outputDir <- file.path(dataDir, "output")

Sequencing is done outside R; at the start of the work flow we have
access to fastq files, in the `fastqDir` directory.  A first step
after sequencing might use [ShortRead][1] to produce a quality
assessment report. Down-sample fastq files if they are big.

    library(ShortRead)
    fls <- list.files(fastqDir, "fastq$", full=TRUE)
    names(fls) <- sub(".fastq", "", basename(fls))
    ## use FastqSampler if fastq files are large
    qas <- lapply(seq_along(fls),
                  function(i, fls) qa(readFastq(fls[i]), names(fls)[i]),
                  fls)
    qa <- do.call(rbind, qas)
    save(qa, file=file.path(outputDir, "qa.rda")
    browseURL(report(qa))

At this stage, one might use tools from [ShortRead][1] or
[Biostrings][2] to remediate low quality reads (e.g., trim low quality
tails or remove artifacts of sample preparation).

Alignment is usually done outside R. Output is a BAM file, one per
sample. `Biostrings::matchPDict` (in some circumsances) and the
[Rsubread][6] package might also be used.

Differential representation typically proceeds from a collection of
known gene locations, including information about gene
structure. Known gene information can come from a variety of sources,
conveniently from UCSC or biomaRt using the [GenomicFeatures][3]
package. These can be saved to disk as sqlite data bases for future
use or to be shared with lab mates. With R >= 2.14.0, there are
semi-annual snapshots available, as used here.

    library(TxDb.Scerevisiae.UCSC.sacCer2.sgdGene)
    txdb <- Scerevisiae_UCSC_sacCer2_ensGene_TxDb
    gnModel <- exonsBy(txdb, "gene")

There are different ways to count, some of which are implemented in
`GenomicRanges::summarizeOverlaps` (available with R >= 2.15). Here is a
simple function that counts any kind of overlap between known genes as
a 'hit', and that discards a read hitting more than one gene. This is
not always satisfactory, but illustrates the flexibility available.

    bamFls <- list.files(bamDir, "bam$", full=TRUE)
    names(bamFls) <- sub("\\..*", "", basename(bamFls))
    counter <- function(fl, gnModel)
    {
        aln <- readGappedAlignments(fl)
        strand(aln) <- "*" # for strand-blind sample prep protocol
        hits <- countOverlaps(aln, gnModel)
        counts <- countOverlaps(gnModel, aln[hits==1])
        names(counts) <- names(gnModel)
        counts
    }
    counts <- sapply(bamFls, counter, gnModel)
    save(counts, file=file.path(outputDir, "counts.rda"))

[edgeR][7] and [DESeq][8] are mature packages that take a
sophisticated statistical approach to the analysis of differential
representation; `DEXSeq` is a recent package (available in R >= 2.14)
to identify differential exon use. An [edgeR][7] work flow is

    library(edgeR)

    ## identify treatment groups
    grp <- factor(<...>)

    ## create data structures
    dge <- DGEList(counts, group=grp)
    dge <- calcNormFactors(dge)

    ## filter uniformative genes
    m <- 1e6 * t(t(dge$counts) / dge$samples$lib.size)
    ridx <- rowSums(m > 1) >= 2
    dge <- dge[ridx,]

    ## comparison between groups
    design <- model.matrix( ~ grp )
    dge <- estimateCommonDisp(dge, design)
    fit <- glmFit(dge, design, dispersion=dge$common.dispersion)
    lrTest <- glmLRT(dge, fit, coef=2)
    tt <- topTags(lrTest, Inf)
    save(tt, file=file.path(dataDir, "tt.rda"))

Two sanity checks along the way are that the treatment groups are
relatively distinct in an MDS (multi-dimensional scaling) plot, and
that the differentially representation tags really are.

    plotMDS.DGEList(dge)
    sapply(rownames(tt$table)[1:4],
           function(x, grp) tapply(counts[x,], grp, mean), 
           grp)

Some annotation information is available from the `org*` model
organism packages; the [biomaRt][13] package is an alternative. Here
we get the sgd ids from the top table, and use these to add gene name
annotations

    library(org.Sc.sgd.db)
    ttids <- rownames(tt$table)
    sgd2gnname <- mget(ttids, org.Sc.sgdGENENAME, ifnotfound=NA)
    ttAnno <- cbind(tt$table, genename=unlist(sgd2gnname))

A next step might involve gene set analyses, as in the [goseq][9]
package.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="install-and-use">Installation and Use</h2>

Follow [installation instructions](/install/) to start using these
packages.  To install the `ShortRead` package and all of its
dependencies, evaluate the commands

    > source("http://bioconductor.org/biocLite.R")
    > biocLite("ShortRead")

Package installation is required only once per R installation. View a
full list of
[available packages](/packages/release/bioc/).

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
package functionality; the [edgeR][7] vignette is an excellent example
of this) in the `ShortRead` package. Use

    > help.start()

to open a web page containing comprehensive help resources.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="sequencing-resources">Sequencing Resources</h2>

The following packages illustrate the diversity of functionality
available; all are in the release version of Bioconductor.

* [IRanges][5], [GenomicRanges][4] and
  [genomeIntervals](/packages/release/bioc/html/genomeIntervals.html)
  for range-based (e.g., chromosomal regions) calculation, data
  manipulation, and general-purpose data
  representation. [Biostrings][2] for alignment, pattern matching
  (e.g., primer removal), and data manipulation of large biological
  sequences or sets of sequences.

* [ShortRead][1] and [Rsamtools][10] for file I/O, quality assessment,
  and high-level, general purpose data summary.  [rtracklayer][11] for
  import and export of tracks on the UCSC genome browser.

* [BSgenome][12] for accessing and manipulating curated whole-genome
  representations.  [GenomicFeatures][3] for annotation of sequence
  features across common genomes, [biomaRt][13] for access to Biomart
  databases.

* [SRAdb](/packages/release/bioc/html/SRAdb.html)
  for querying and retrieving data from the Sequence Read Archive.

* ChIP-seq and related (e.g., motif discovery, identification of
  high-coverage segments) activities are facilitated by packages such
  as
  [CSAR](/packages/release/bioc/html/CSAR.html),
  [chipseq](/packages/release/bioc/html/chipseq.html),
  [ChIPseqR](/packages/release/bioc/html/ChIPseqR.html),
  [ChIPsim](/packages/release/bioc/html/ChIPsim.html),
  [ChIPpeakAnno](/packages/release/bioc/html/ChIPpeakAnno.html),
  [DiffBind](/packages/release/bioc/html/DiffBind.html),
  [iSeq](/packages/release/bioc/html/iSeq.html),
  [rGADEM](/packages/release/bioc/html/rGADEM.html),
  [segmentSeq](/packages/release/bioc/html/segmentSeq.html),
  [BayesPeak](/packages/release/bioc/html/BayesPeak.html),
  [PICS](/packages/release/bioc/html/PICS.html).

* Differential expression and RNA-seq style analysis can be
  accomplished with [edgeR][7], [DESeq][8],
  [baySeq](/packages/release/bioc/html/baySeq.html),
  [DEGseq](/packages/release/bioc/html/DEGseq.html),
  [Genominator](/packages/release/bioc/html/Genominator.html). The
  [DEXSeq][14] package, available with R >= 2.14, allows exon-level
  differential representation.

Additional packages are classified in the [biocViews][15] hierarchy,
for instance under Software : AssayTechnologies :
HighThroughputSequencing.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

[1]: /packages/release/bioc/html/ShortRead.html
[2]: /packages/release/bioc/html/Biostrings.html
[3]: /packages/release/bioc/html/GenomicFeatures.html
[4]: /packages/release/bioc/html/GenomicRanges.html
[5]: /packages/release/bioc/html/IRanges.html
[6]: /packages/release/bioc/html/Rsubread.html
[7]: /packages/release/bioc/html/edgeR.html
[8]: /packages/release/bioc/html/DESeq.html
[9]: /packages/release/bioc/html/goseq.html
[10]: /packages/release/bioc/html/Rsamtools.html
[11]: /packages/release/bioc/html/rtracklayer.html
[12]: /packages/release/bioc/html/BSgenome.html
[13]: /packages/release/bioc/html/biomaRt.html
[14]: /packages/devel/bioc/html/DEXSeq.html
[15]: /packages/release/BiocViews.html
