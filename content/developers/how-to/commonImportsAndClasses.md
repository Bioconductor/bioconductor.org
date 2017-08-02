# Re-using Bioconductor Import Methods and Classes

It is highly recommended to reuse existing methods for importing data as well as
reuse already established classes and functionality rather than creating 'new'
or repeating code. Here are some suggestions for importing different file types
and commonly used _Bioconductor_ classes. For more classes and functionality
also try searching in
[BiocViews](http://bioconductor.org/packages/release/BiocViews.html#___Software)
for your data type.

## Importing

+ GTF, GFF, BED, BigWig, ...<br><code>rtracklayer::import()</code><br><br>
+ FASTA<br><code>Biostrings::readDNAStringSet()</code><br><br>
+ SAM / BAM<br><code>Rsamtools::scanBam()</code><br><code>GenomicAlignments::readGAlignment*()</code><br><br>
+ VCF<br><code>VariantAnnotation::readVcf()</code><br><br>
+ FASTQ<br><code>ShortRead::readFastq()</code><br><br>
+ MS data (XML-based and mfg formats)<br><code>MSnbase::readMSData()</code><br><code>MSnbase::readMgfData()</code>


## Common Classes

+ DNA / RNA / AA sequences<br><code>Biostrings::*Stringset()</code><br><br>
+ GeneSets<br><code>GSEABase::GeneSet()</code><br><code>GSEABase::GeneSetCollection()</code><br><br>
+ 1-based, closed-interval genomic coordinates<br><code>GenomicRanges::GRanges()</code><br><br>
+ rectangular freature x sample data (RNAseq count matrix, microarray, ...)<br><code>SummarizedExperiment::SummarizedExperiment()</code><br><br>
+ multi omics data<br><code>MultiAssayExperiment::MultiAssayExperiment()</code><br><br>
+ single cell data<br><code>SingleCellExperiment::SingleCellExperiment()</code><br><br>
+ mass spec data<br><code>MSnbase::MSnExp()</code>
