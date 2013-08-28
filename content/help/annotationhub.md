# AnnotationHub

_AnnotationHub_ is a new approach to providing annotation resources to
the Bioconductor community. The initial plan is to make available
fasta, GFF / GTF, BED, VCF, and similar files from entities such as
Ensembl, UCSC, ENCODE, and 1000 Genomes projects as objects ready for
work in Bioconductor, e.g., _GRanges_ representations of BED
files.

As an example, the following lines discover and retrieve an ENCODE
'narrowPeaks' file as a GRanges representation.

    library(AnnotationHub)
	hub <- AnnotationHub()
	
	## data exploration
	length(names(hub))                  # resources available
	md <- metadata(hub)                 # DataFrame
	hub$goldenpath.hg19.encodeDCC<tab>  # tab completion

	## retrieval
	res <- hub$goldenpath.hg19.encodeDCC.wgEncodeUwTfbs.wgEncodeUwTfbsNhlfCtcfStdPkRep1.narrowPeak_0.0.1.RData
    class(res)         # GRanges representation of ENCODE bed file  

_AnnotationHub_ is currently (27 January, 2013) available from the
Bioconductor subversion
[repository](/developers/source-control/) for
use with the development version of [R](http://r-project.org); our
intention is to include _AnnotationHub_ in Bioconductor 2.12.
