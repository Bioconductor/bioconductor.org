# Common Bioconductor Methods and Classes

We strongly recommend reusing existing methods for importing data, and
reusing established classes for representing data. Here are some
suggestions for importing different file types and commonly used
_Bioconductor_ classes. For more classes and functionality also try
searching in [BiocViews](/packages) for your data type.

## Importing

+ GTF, GFF, BED, BigWig, etc., -- [rtracklayer][]`::import()`
+ VCF -- [VariantAnnotation][]`::readVcf()`
+ SAM / BAM -- [Rsamtools][]`::scanBam()`,
  [GenomicAlignments][]`::readGAlignment*()`
+ FASTA -- [Biostrings][]`::readDNAStringSet()`
+ FASTQ -- [ShortRead][]`::readFastq()`
+ MS data (XML-based and mgf formats) -- [Spectra][]`::Spectra()`, 
  [MSnbase][]`::readMSData()`, [Spectra][]`::Spectra(source = MsBackendMgf::MsBackendMgf())`, 
  [MSnbase][]`::readMgfData()`

## Common Classes

+ Rectangular feature x sample data --
  [SummarizedExperiment][]`::SummarizedExperiment()` (RNAseq count
  matrix, microarray, ...)
+ Genomic coordinates -- [GenomicRanges][]`::GRanges()` (1-based,
  closed interval)
+ Genomic coordinates from multiple samples --
  [GenomicRanges][]`::GRangesList()`
+ Ragged genomic coordinates -- [RaggedExperiment][]`::RaggedExperiment()`
+ DNA / RNA / AA sequences -- [Biostrings][]`::*StringSet()`
+ Gene sets -- [BiocSet][]`::BiocSet()`,
  [GSEABase][]`::GeneSet()`,
  [GSEABase][]`::GeneSetCollection()`
+ Multi-omics data --
  [MultiAssayExperiment][]`::MultiAssayExperiment()`
+ Single cell data --
  [SingleCellExperiment][]`::SingleCellExperiment()`
+ Mass spec data -- [Spectra][]`::Spectra()`, [MSnbase][]`::MSnExp()`

[rtracklayer]: https://bioconductor.org/packages/rtracklayer
[Biostrings]: https://bioconductor.org/packages/Biostrings
[Rsamtools]: https://bioconductor.org/packages/Rsamtools
[GenomicAlignments]: https://bioconductor.org/packages/GenomicAlignments
[VariantAnnotation]: https://bioconductor.org/packages/VariantAnnotation
[ShortRead]: https://bioconductor.org/packages/ShortRead
[MSnbase]: https://bioconductor.org/packages/MSnbase
[SummarizedExperiment]: https://bioconductor.org/packages/SummarizedExperiment
[GenomicRanges]: https://bioconductor.org/packages/GenomicRanges
[RaggedExperiment]: https://bioconductor.org/packages/RaggedExperiment
[BiocSet]: https://bioconductor.org/packages/BiocSet
[GSEABase]: https://bioconductor.org/packages/GSEABase
[MultiAssayExperiment]: https://bioconductor.org/packages/MultiAssayExperiment
[SingleCellExperiment]: https://bioconductor.org/packages/SingleCellExperiment
[Spectra]: https://bioconductor.org/packages/Spectra
[MsBackendMgf]: https://bioconductor.org/packages/MsBackendMgf

