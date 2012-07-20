# ![](/images/icons/help.gif)Using Bioconductor to Annotate Genetic Variants 

The VariantAnnotation package has facilities for reading in all or portions 
of Variant Call Format (VCF) files. Structural location information can be 
determined as well as amino acid coding changes for non-synonymous variants. 
Consequences of the coding changes can be investigated with the SIFT and 
PolyPhen database packages.


* [Sample Workflow](#sample-workflow) 
* [Installation and Use](#install-and-use)
* [Exploring Package Content](#exploring-package-content)
* [Resources](#resources)

<h2 id="sample-workflow">Sample Workflow</h2>

Performed with Bioconductor 2.11 and R >= 2.15; 
VariantAnnotation 1.3.9 in the devel branch.

This workflow annotates variants found in the Transient Receptor Potential 
Vanilloid (TRPV) gene family on chromosome 17. The VCF file is available in 
the cgdv17 data package and contains Complete Genomics data for population 
type CEU.

    > library(VariantAnnotation)
    > library(cgdv17)
    > file <- system.file("vcf", "NA06985_17.vcf.gz", package = "cgdv17")
     
    ## Explore the file header with scanVcfHeader
    > hdr <- scanVcfHeader(file)
     
    > info(hdr) 
    DataFrame with 3 rows and 3 columns
            Number        Type                 Description
       <character> <character>                 <character>
    NS           1     Integer Number of Samples With Data
    DP           1     Integer                 Total Depth
    DB           0     Flag    dbSNP membership, build 131
     
    > geno(hdr) 
    DataFrame with 12 rows and 3 columns
                Number        Type                         Description
           <character> <character>                         <character>
    GT               1      String                            Genotype
    GQ               1     Integer                    Genotype Quality
    DP               1     Integer                          Read Depth
    HDP              2     Integer                Haplotype Read Depth
    HQ               2     Integer                   Haplotype Quality
    PS               2     Integer                           Phase Set
    GENE             .      String                    Overlaping Genes
    mRNA             .      String                     Overlaping mRNA
    rmsk             .      String                  Overlaping Repeats
    segDup           .      String Overlaping segmentation duplication
    rCov             1       Float                   relative Coverage
    cPd              1      String                called Ploidy(level)
     
Convert the gene symbols to gene ids compatible with the 
TxDb.Hsapiens.UCSC.hg19.knownGene annotations. The annotaions 
are used to define the TRPV ranges that will be extracted from
the VCF file.

    > ## get entrez ids from gene symbols
    > library(org.Hs.eg.db)
    > genesym <- c("TRPV1", "TRPV2", "TRPV3")
    > geneid <- select(org.Hs.eg.db, keys=genesym, keytype="SYMBOL",
    +                  cols="ENTREZID")
    > geneid
      SYMBOL ENTREZID
    1  TRPV1     7442
    2  TRPV2    51393
    3  TRPV3   162514
     
Load the annotation package and create a list of transcripts
by gene.

    > library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    > txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    > txbygene = transcriptsBy(txdb, "gene")
     
Subset the annotations on chromosome 17 and adjust the seqlevels 
to match those in the VCF file.

    > tx_chr17 <- keepSeqlevels(txbygene, "chr17")
    > tx_17 <- renameSeqlevels(tx_chr17, c(chr17="17"))
    
    ## Create the gene ranges for the TRPV genes
    > rngs <- lapply(geneid$ENTREZID, 
                     function(id) 
                         range(tx_17[names(tx_17) %in% id]))
    > gnrng <- unlist(do.call(c, rngs), use.names=FALSE)
    > names(gnrng) <- geneid$SYMBOL 
     
To retrieve a subset of data from a VCF file, create a ScanVcfParam
object. This object can specify genomic coordinates (ranges) or 
individual VCF elements to be extracted. When ranges are extracted, 
a tabix index file must exist for the VCF. See ?indexTabix for details.

    > param <- ScanVcfParam(which = gnrng, info = "DP", geno = c("GT", "cPd"))
    > param
    class: ScanVcfParam 
    vcfWhich: 1 elements
    vcfFixed: character() [All] 
    vcfInfo: DP 
    vcfGeno: GT, cPd 
     
    ## Extract the TRPV ranges from the VCF file 
    > vcf <- readVcf(file, "hg19", param)
    ## Inspect the VCF object with the 'fixed', 'info' and 'geno' accessors
    > vcf
    class: VCF 
    dim: 387 1 
    genome: hg19 
    exptData(1): header
    fixed(4): REF ALT QUAL FILTER
    info(1): DP
    geno(2): GT cPd
    rownames(387): rs402369 17:3469483 ... rs322959 rs322958
    rowData values names(1): paramRangeID
    colnames(1): GS06985-1100-37-ASM
    colData names(1): Samples
     
    > head(fixed(vcf))
    GRanges with 6 ranges and 5 elementMetadata cols:
                 seqnames             ranges strand | paramRangeID            REF
                    <Rle>          <IRanges>  <Rle> |     <factor> <DNAStringSet>
        rs402369       17 [3469401, 3469401]      * |        TRPV1              A
      17:3469483       17 [3469483, 3469483]      * |        TRPV1              A
      17:3470672       17 [3470672, 3470676]      * |        TRPV1          AAAAA
      17:3471626       17 [3471626, 3471627]      * |        TRPV1             AA
        rs161363       17 [3472181, 3472181]      * |        TRPV1              C
        rs224546       17 [3472871, 3472871]      * |        TRPV1              T
                                ALT      QUAL      FILTER
                 <DNAStringSetList> <numeric> <character>
        rs402369           ########       120        PASS
      17:3469483           ########         0        PASS
      17:3470672           ########         0        PASS
      17:3471626           ########         0        PASS
        rs161363           ########        59        PASS
        rs224546           ########       157        PASS
      ---
      seqlengths:
       17
       NA
    > geno(vcf)
    SimpleList of length 2
    names(2): cPd GT
     
To find the structural location of the variants, use the locateVariants 
function with the TxDb.Hsapiens.UCSC.hg19.knownGene package that was
loaded eariler. The variants in the VCF object have chromosome name "17" 
while the annotation has "chr17". Adjust the seqlevels (chromosome names) of 
the VCF object to match that of the annotation.

    > seqlevels(vcf)
    [1] "17"
    > head(seqlevels(txdb))
    [1] "chr1" "chr2" "chr3" "chr4" "chr5" "chr6"
    > ## seqlevels do not match
    > intersect(seqlevels(vcf), seqlevels(txdb))
    character(0)
    > vcf_mod <- renameSeqlevels(vcf, c("17"="chr17"))
    > ## seqlevels now match
    > intersect(seqlevels(vcf_mod), seqlevels(txdb))
    [1] "chr17"
     
    ## Use the 'region' argument to define the region
    ## of interest. See ?locateVariants for details.
    > cds <- locateVariants(vcf_mod, txdb, CodingVariants())
    > five <- locateVariants(vcf_mod, txdb, FiveUTRVariants())
    > splice <- locateVariants(vcf_mod, txdb, SpliceSiteVariants())
    > intron <- locateVariants(vcf_mod, txdb, IntronVariants())
    > all <- locateVariants(vcf_mod, txdb, AllVariants())
     
Each row in cds represents a variant-transcript match so multiple rows
per variant are possible. If we are interested in gene-centric questions
the data can be summarized by gene reguardless of transcript.
     
    > ## Did any variants match more than one gene
    > table(sapply(split(values(all)[["GENEID"]], values(all)[["QUERYID"]]), 
              function(x)
                  length(unique(x)) > 1))
    FALSE  TRUE 
      379     8 
     
    > ## Summarize the number of variants by gene
    > idx <- sapply(split(values(all)[["QUERYID"]], values(all)[["GENEID"]]), 
               unique)
    > sapply(idx, length)
    162514  23729  51393   7442  84690 
       178      2     63    146      6 
     
    > ## Summarize variant location by gene
    > sapply(names(idx), 
          function(nm) {
              d <- all[values(all)[["GENEID"]] %in% nm, c("QUERYID", "LOCATION")]
              table(values(d)[["LOCATION"]][duplicated(d) == FALSE])
          })
               162514 23729 51393 7442 84690
    spliceSite      2     0     0    1     0
    intron        162     0    58  131     1
    fiveUTR         2     0     1    4     5
    threeUTR        6     2     1    2     0
    coding          6     0     3    8     0
    intergenic      0     0     0    0     0
     
Amino acid coding for non-synonymous variants can be computed
with the function predictCoding. The BSgenome.Hsapiens.UCSC.hg19 
package is used as the source of the reference alleles. Variant 
alleles are provided by the user.

    > library(BSgenome.Hsapiens.UCSC.hg19)
    > aa <- predictCoding(vcf_mod, txdb, Hsapiens)
    Warning messages:
    1: In .predictCoding(query, subject, seqSource, varAllele, ...) :
      records with missing 'varAllele' values will be ignored
    2: In .predictCoding(query, subject, seqSource, varAllele, ...) :
      varAllele values containing 'N' will not be translated
     
predictCoding returns results for coding variants only. As with 
locateVariants, the output has one row per variant-transcript match
so multiple rows per variant are possible.
     
    ## Did any variants match more than one gene
    > table(sapply(split(values(aa)[["GENEID"]], values(aa)[["QUERYID"]]), 
              function(x)
                  length(unique(x)) > 1))
    FALSE 
       17 
    
    > ## Summarize the number of variants by gene
    > idx <- sapply(split(values(aa)[["QUERYID"]], values(aa)[["GENEID"]], 
                    drop=TRUE), unique)
    > sapply(idx, length)
    162514  51393   7442 
         6      3      8 
    
    > ## Summarize variant consequence by gene
    > sapply(names(idx), 
             function(nm) {
                 d <- aa[values(aa)[["GENEID"]] %in% nm, c("QUERYID","CONSEQUENCE")]
                 table(values(d)[["CONSEQUENCE"]][duplicated(d) == FALSE])
             })
                   162514 51393 7442
    nonsynonymous       2     1    1
    not translated      1     0    5
    synonymous          3     2    2
    
The variants 'not translated' are explained by the warnings thrown when
predictCoding was called. Variants that have a missing varAllele or have an
'N' in the varAllele are not translated. If the varAllele substitution had
resulted in a frameshift the consequence would be 'frameshift'. See
?predictCoding for details.
    
The SIFT.Hsapiens.dbSNP132 and PolyPhen.Hsapiens.dbSNP131 packages 
provide predictions of how damaging amino acid coding changes may
be to protein structure and function. Both packages search on rsid.
    
The pre-computed predictions in the SIFT and PolyPhen packages are based 
on specific gene models. SIFT is based on Ensembl and PolyPhen on UCSC 
Known Gene. The TranscriptDb we used to identify coding variants was 
from UCSC Known Gene so we will use PolyPhen for predictions.
    
    ## Load the PolyPhen package and explore the available keys and columns
    > library(PolyPhen.Hsapiens.dbSNP131)
    > keys <- keys(PolyPhen.Hsapiens.dbSNP131)
    > cols <- cols(PolyPhen.Hsapiens.dbSNP131)
    ## column descriptions are found at ?PolyPhenDbColumns
    > cols(PolyPhen.Hsapiens.dbSNP131)
     [1] "RSID"        "TRAININGSET" "OSNPID"      "OACC"        "OPOS"       
     [6] "OAA1"        "OAA2"        "SNPID"       "ACC"         "POS"        
    [11] "AA1"         "AA2"         "NT1"         "NT2"         "PREDICTION" 
    [16] "BASEDON"     "EFFECT"      "PPH2CLASS"   "PPH2PROB"    "PPH2FPR"    
    [21] "PPH2TPR"     "PPH2FDR"     "SITE"        "REGION"      "PHAT"       
    [26] "DSCORE"      "SCORE1"      "SCORE2"      "NOBS"        "NSTRUCT"    
    [31] "NFILT"       "PDBID"       "PDBPOS"      "PDBCH"       "IDENT"      
    [36] "LENGTH"      "NORMACC"     "SECSTR"      "MAPREG"      "DVOL"       
    [41] "DPROP"       "BFACT"       "HBONDS"      "AVENHET"     "MINDHET"    
    [46] "AVENINT"     "MINDINT"     "AVENSIT"     "MINDSIT"     "TRANSV"     
    [51] "CODPOS"      "CPG"         "MINDJNC"     "PFAMHIT"     "IDPMAX"     
    [56] "IDPSNP"      "IDQMIN"      "COMMENTS"
    
    ## Get the rsids for the non-synonymous variants from the
    ## predictCoding results
    > rsid <- unique(names(aa)[values(aa)[["CONSEQUENCE"]] == "nonsynonymous"]) 
    
    ## Retrieve predictions for non-synonymous variants. Two of the six variants 
    ## are found in the PolyPhen database. 
    > select(PolyPhen.Hsapiens.dbSNP131, keys=rsid, 
             cols=c("AA1", "AA2", "PREDICTION"))
            RSID  AA1  AA2 PREDICTION
    1   rs222747    M    I     benign
    2   rs222748 <NA> <NA>       <NA>
    3  rs1129235 <NA> <NA>       <NA>
    4 rs11078458 <NA> <NA>       <NA>
    5  rs1039519 <NA> <NA>       <NA>
    6   rs322965    I    V     benign
    Warning message:
    In .formatPPDbSelect(raw, keys) : 
    key not found in database : rs222748
    key not found in database : rs1129235
    key not found in database : rs11078458
    key not found in database : rs1039519

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>


<h2 id="install-and-use">Installation and Use</h2>

Follow [installation instructions](/install/) to start using these
packages.  To install VariantAnnotation use 

    > library(BiocInstaller) 
    > biocLite("VariantAnnotation")

Package installation is required only once per R installation. View a
full list of available
[software](/packages/release/bioc/)
and 
[annotation](/packages/release/data/annotation/)
packages.

To use the `VariantAnnotation`, evaluate the commands

    > library(VariantAnnotation)

These commands are required once in each R session.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="exploring-package-content">Exploring Package Content</h2>

Packages have extensive help pages, and include vignettes highlighting
common use cases. The help pages and vignettes are available from
within R. After loading a package, use syntax like

    > help(package="VariantAnnotation")
    > ?predictCoding

to obtain an overview of help on the `VariantAnnotation` package, and 
the `predictCoding` function. View the package vignette with

    > browseVignettes(package="VariantAnnotation")

To view vignettes providing a more comprehensive introduction to
package functionality use

    > help.start()


<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

