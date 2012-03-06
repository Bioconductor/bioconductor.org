# ![](/images/icons/help.gif)Using Bioconductor to Annotate Genetic Variants 

The VariantAnnotation package has facilities for reading
in all or portions of Variant Call Format (VCF) files. Structural 
locatation information can be determined as well as amino acid 
coding changes for the non-synonymous variants. Consequences of 
the coding changes can be investigated with the SIFT and PolyPhen 
database packages.


* [Sample Workflow](#sample-workflow)  
* [Installation and Use](#install-and-use)
* [Exploring Package Content](#exploring-package-content)
* [Resources](#resources)

<h2 id="sample-workflow">Sample Workflow</h2>

This workflow walks through the annotation of variants located 
in the Transient Receptor Potential Vanilloid (TRPV) gene family. 
The VCF file is available in the cgdv17 data package
and contains Complete Genomics data for chromosome 17 from population
type CEU.

    > library(VariantAnnotation)
    > library(cgdv17)
    > file <- system.file("vcf", "NA06985_17.vcf.gz", package = "cgdv17")
     
    ## Explore the file header with scanVcfHeader
    > hdr <- scanVcfHeader(file)
     
    > hdr[[1]]$Header$INFO
    DataFrame with 3 rows and 3 columns
            Number        Type                 Description
       <character> <character>                 <character>
    NS           1     Integer Number of Samples With Data
    DP           1     Integer                 Total Depth
    DB           0        Flag dbSNP membership, build 131
     
    > hdr[[1]]$Header$FORMAT
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

    > library(org.Hs.eg.db)
    > genesym <- c("TRPV1", "TRPV2", "TRPV3")
    > ## get entrez ids from gene symbols
    > geneid <- select(org.Hs.eg.db, keys=genesym, keytype="SYMBOL",
    +                  cols="ENTREZID")
    > geneid
          SYMBOL ENTREZID
    6180   TRPV1     7442
    12054  TRPV2    51393
    19699  TRPV3   162514
     
Load the annotation package and create a list of transcripts
by gene.

    > library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    > txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    > txbygene = transcriptsBy(txdb, "gene")
     
Subset the annotations on chromosome 17 and adjust the seqlevels 
to match those in the VCF file.

    > tx_chr17 <- keepSeqlevels(txbygene, "chr17")
    > tx_17 <- renameSeqlevels(tx_chr17, c(chr17="17"))
     
    ## Create the gene ranges for the TRPV genes and collapse to
    ## a single range 
    > rngs <- unlist(tx_17[names(tx_17) %in% geneid$ENTREZID], use.names=FALSE)
    > gnrng <- reduce(rngs)
     
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
    vcfAsGRanges: FALSE 
     
    ## Extract the TRPV ranges from the VCF file 
    > vcf <- readVcf(file, "hg19", param)
    ## Inspect the VCF object with the 'fixed', 'info' and 'geno' accessors
    > vcf
    class: VCF 
    dim: 387 1 
    genome: hg19 
    exptData(1): HEADER
    fixed(4): REF ALT QUAL FILTER
    info(1): DP
    geno(2): cPd GT
    rownames(387): rs3813769 rs72838762 ... rs460716 rs3840876
    rowData values names(1): rangeID
    colnames(1): GS06985-1100-37-ASM
    colData names(1): Samples
     
    > head(fixed(vcf))
    GRanges with 6 ranges and 5 elementMetadata cols:
                  seqnames               ranges strand |              rangeID
                     <Rle>            <IRanges>  <Rle> |             <factor>
        rs3813769       17 [16318932, 16318932]      * | 17:16318856-16340317
       rs72838762       17 [16319395, 16319395]      * | 17:16318856-16340317
       rs72838763       17 [16319512, 16319512]      * | 17:16318856-16340317
      17:16319536       17 [16319536, 16319545]      * | 17:16318856-16340317
      17:16320586       17 [16320586, 16320589]      * | 17:16318856-16340317
       rs12951377       17 [16321959, 16321959]      * | 17:16318856-16340317
                             REF                ALT      QUAL      FILTER
                  <DNAStringSet> <DNAStringSetList> <numeric> <character>
        rs3813769              T           ########       118        PASS
       rs72838762              A           ########        90        PASS
       rs72838763              A           ########        86        PASS
      17:16319536     CCACCTCCCA           ########         0        PASS
      17:16320586           TTTT           ########         0        PASS
       rs12951377              A           ########       276        PASS
      ---
      seqlengths:
       17
       NA
    > geno(vcf)
    SimpleList of length 2
    names(2): cPd GT
     
To find the structural location of the variants we use the locateVariants 
function with the TxDb.Hsapiens.UCSC.hg19.knownGene package we loaded 
eariler. The variants in our VCF object have chromosome name "17" while 
the annotation has "chr17". Adjust the seqlevels (chromosome names) of 
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
     
    ## call locateVariants
    > loc <- locateVariants(vcf_mod, txdb)
     
    > ## For the 387 variants in vcf_mod 
    > dim(vcf_mod)
    [1] 387   1
     
    > ## we have 2514 rows of output
    > dim(loc)
    [1] 2514    7
     
Each row in loc represents a variant-transcript match. To 
answer gene-centric questions we summarize the data by gene
reguardless of transcript.
     
    > ## Did any variants match more than one gene
    > table(sapply(split(loc$geneID, loc$queryID), function(x)
    +     length(unique(x)) > 1
    + ))
    FALSE  TRUE 
      379     8 
     
    > ## Summarize the number of variants by gene
    > idx <- sapply(split(loc$queryID, loc$geneID), unique)
    > sapply(idx, length)
    162514  23729  51393   7442  84690 
       178      2     63    146      6 
     
    > ## Summarize variant location by gene
    > sapply(names(idx), function(nm) {
    +        d <- loc[loc$geneID %in% nm, c("queryID", "location")]
    +        table(d$location[duplicated(d) == FALSE])
    + })
                      162514 23729 51393 7442 84690
    transcript_region      0     0     0    0     0
    intron               163     0    58  132     0
    5'UTR                  2     0     1    2     0
    3'UTR                  7     2     1    4     6
    coding                 6     0     3    8     0
    intergenic             0     0     0    0     0
     
Amino acid coding for non-synonymous variants :
Load the BSgenome.Hsapiens.UCSC.hg19 package and call predictCoding
on the VCF object with modified seqlevels.

    > library(BSgenome.Hsapiens.UCSC.hg19)
    > aa <- predictCoding(vcf_mod, txdb, Hsapiens)
    Warning messages:
    1: In predictCoding(query = rdf, subject = subject, seqSource = seqSource,  :
      records with missing 'varAllele' values will be ignored
    2: In predictCoding(query = rdf, subject = subject, seqSource = seqSource,  :
      varAllele values containing 'N' will not be translated
     
As with locateVariants, the output from predictCoding has one row per 
variant-transcript match. Recall that predictCoding returns results for 
coding variants only so the variant per gene numbers are reduced from what 
was seen locateVariants.
     
    ## Did any variants match more than one gene
    > table(sapply(split(aa$geneID, aa$queryID), function(x)
    +     length(unique(x)) > 1
    + ))
    FALSE 
       12 
    
    > ## Summarize the number of variants by gene
    > idx <- sapply(split(aa$queryID, aa$geneID, drop=TRUE), unique)
    > sapply(idx, length)
    162514  51393   7442 
         5      3      4 
    
    > ## Summarize variant consequence by gene
    > sapply(names(idx), function(nm) {
    +        d <- aa[aa$geneID %in% nm, c("queryID", "consequence")]
    +        table(d$consequence[duplicated(d) == FALSE])
    + })
                  162514 51393 7442
    frameshift         0     0    1
    nonsynonymous      3     0    2
    synonymous         2     3    1
    
The SIFT.Hsapiens.dbSNP132 and PolyPhen.Hsapiens.dbSNP131 packages 
provide predictions of how damaging amino acid coding changes may
be on protein structure and function. Both packages search on rsid.
    
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
    > nonsyn <- aa$queryID[aa$consequence == "nonsynonymous"]
    > rsid <- unique(names(rowData(vcf_mod))[nonsyn])
    
    ## Choose a subset of columns 
    > subst <- c("AA1", "AA2", "PREDICTION")
    
    ## Retrieve predictions for non-synonymous variants
    > select(PolyPhen.Hsapiens.dbSNP131, keys=rsid, cols=subst)
           RSID AA1 AA2 PREDICTION
    3  rs322965   I   V     benign
    31 rs224534   T   I     benign
    Warning message:
    In .formatPPDbSelect(raw, keys = keys) :
      key not found in database : rs11078458
    key not found in database : rs1039519
    key not found in database : rs222748

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

