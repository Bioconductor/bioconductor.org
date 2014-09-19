<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{Annotate Genetic Variants}
-->

# ![](/images/icons/help.gif)Annotate Genetic Variants 

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

Performed with Bioconductor 3.0 and R >= 3.1; 
VariantAnnotation 1.11.35.

This workflow annotates variants found in the Transient Receptor Potential 
Vanilloid (TRPV) gene family on chromosome 17. The VCF file is available in 
the cgdv17 data package and contains Complete Genomics data for population 
type CEU.


    library(VariantAnnotation)
    library(cgdv17)
    file <- system.file("vcf", "NA06985_17.vcf.gz", package = "cgdv17")
    
    ## Explore the file header with scanVcfHeader
    hdr <- scanVcfHeader(file)
    
    info(hdr)

    ## DataFrame with 3 rows and 3 columns
    ##         Number        Type                 Description
    ##    <character> <character>                 <character>
    ## NS           1     Integer Number of Samples With Data
    ## DP           1     Integer                 Total Depth
    ## DB           0        Flag dbSNP membership, build 131

    
    geno(hdr)

    ## DataFrame with 12 rows and 3 columns
    ##             Number        Type                         Description
    ##        <character> <character>                         <character>
    ## GT               1      String                            Genotype
    ## GQ               1     Integer                    Genotype Quality
    ## DP               1     Integer                          Read Depth
    ## HDP              2     Integer                Haplotype Read Depth
    ## HQ               2     Integer                   Haplotype Quality
    ## ...            ...         ...                                 ...
    ## mRNA             .      String                     Overlaping mRNA
    ## rmsk             .      String                  Overlaping Repeats
    ## segDup           .      String Overlaping segmentation duplication
    ## rCov             1       Float                   relative Coverage
    ## cPd              1      String                called Ploidy(level)

     
Convert the gene symbols to gene ids compatible with the 
TxDb.Hsapiens.UCSC.hg19.knownGene annotations. The annotaions 
are used to define the TRPV ranges that will be extracted from
the VCF file.


    ## get entrez ids from gene symbols
    library(org.Hs.eg.db)
    genesym <- c("TRPV1", "TRPV2", "TRPV3")
    geneid <- select(org.Hs.eg.db, keys = genesym, keytype = "SYMBOL", columns = "ENTREZID")

    ## Warning: The 'cols' argument has been deprecated and replaced by 'columns'
    ## for versions of Bioc that are higher than 2.13.  Please use the 'columns'
    ## argument anywhere that you previously used 'cols'

    geneid

    ##   SYMBOL ENTREZID
    ## 1  TRPV1     7442
    ## 2  TRPV2    51393
    ## 3  TRPV3   162514

     
Load the annotation package and keep only the standard chromosomes.


    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    keepStandardChromosomes(txdb)


    ## Modify the seqlevels (chromosomes) in the txdb to match
    those in the VCF file. This step is necessary because we want
    to use ranges from the txdb to extract a subset from the VCF.
    renameSeqlevels(txdb, gsub("chr", "", seqlevels(txdb)))


    ## Create a list of transcripts by gene:
    txbygene = transcriptsBy(txdb, "gene")

    
    ## Create the gene ranges for the TRPV genes
    gnrng <- unlist(range(txbygene[geneid$ENTREZID]), use.names=FALSE)
    names(gnrng) <- geneid$SYMBOL

     
A ScanVcfParam object is used to retrieve data subsets. This object 
can specify genomic coordinates (ranges) or individual VCF elements.
Extractions of ranges (vs fields) requires a tabix index.
See ?indexTabix for details.


    param <- ScanVcfParam(which = gnrng, info = "DP", geno = c("GT", "cPd"))
    param

    ## class: ScanVcfParam 
    ## vcfWhich: 1 elements
    ## vcfFixed: character() [All] 
    ## vcfInfo: DP 
    ## vcfGeno: GT cPd 
    ## vcfSamples:  
    
    ## Extract the TRPV ranges from the VCF file
    vcf <- readVcf(file, "hg19", param)
    ## Inspect the VCF object with the 'fixed', 'info' and 'geno' accessors
    vcf

    ## class: CollapsedVCF 
    ## dim: 405 1 
    ## rowData(vcf):
    ##   GRanges with 5 metadata columns: paramRangeID, REF, ALT, QUAL, FILTER
    ## info(vcf):
    ##   DataFrame with 1 column: DP
    ## info(header(vcf)):
    ##       Number Type    Description
    ##    DP 1      Integer Total Depth
    ## geno(vcf):
    ##   SimpleList of length 2: cPd, GT
    ## geno(header(vcf)):
    ##        Number Type   Description         
    ##    cPd 1      String called Ploidy(level)
    ##    GT  1      String Genotype  
    
    head(fixed(vcf))

    ## DataFrame with 6 rows and 4 columns
    ##              REF                ALT      QUAL      FILTER
    ##   <DNAStringSet> <DNAStringSetList> <numeric> <character>
    ## 1              A                  G       120        PASS
    ## 2              A                            0        PASS
    ## 3          AAAAA                            0        PASS
    ## 4             AA                            0        PASS
    ## 5              C                  T        59        PASS
    ## 6              T                  C       157        PASS

    
    geno(vcf)

    ## List of length 2
    ## names(2): cPd GT


To find the structural location of the variants, use the locateVariants 
function with the TxDb.Hsapiens.UCSC.hg19.knownGene package that was
loaded eariler. 

    
    ## Use the 'region' argument to define the region of interest. See
    ## ?locateVariants for details.
    cds <- locateVariants(vcf, txdb, CodingVariants())
    five <- locateVariants(vcf, txdb, FiveUTRVariants())
    splice <- locateVariants(vcf, txdb, SpliceSiteVariants())
    intron <- locateVariants(vcf, txdb, IntronVariants())
    all <- locateVariants(vcf, txdb, AllVariants())

     
Each row in cds represents a variant-transcript match so multiple rows
per variant are possible. If we are interested in gene-centric questions
the data can be summarized by gene regardless of transcript.


    ## Did any variants match more than one gene?
    table(sapply(split(values(all)[["GENEID"]], 
          values(all)[["QUERYID"]]), function(x) length(unique(x)) > 1))

    ## 
    ## FALSE  TRUE 
    ##   391   11 

    
    ## Summarize the number of variants by gene:
    idx <- sapply(split(values(all)[["QUERYID"]], values(all)[["GENEID"]]), unique)
    sapply(idx, length)

    ## 125144 162514  51393   7442  84690 
    ##      1    172     62    143     35 
    
    ## Summarize variant location by gene:
    sapply(names(idx), function(nm) {
        d <- all[values(all)[["GENEID"]] %in% nm, c("QUERYID", "LOCATION")]
        table(values(d)[["LOCATION"]][duplicated(d) == FALSE])
    })


    ##            125144 162514 51393 7442 84690
    ## spliceSite      0      2     0    1     0
    ## intron          0    153    58  117    19
    ## fiveUTR         0      0     0    0     0
    ## threeUTR        0      0     0    0     0
    ## coding          0      5     3    8     0
    ## intergenic      0      0     0    0     0
    ## promoter        1     12     1   17    16



Amino acid coding for non-synonymous variants can be computed
with the function predictCoding. The BSgenome.Hsapiens.UCSC.hg19 
package is used as the source of the reference alleles. Variant 
alleles are provided by the user.

    ## Convert the VCF and txdb to UCSC seqlevel style to match BSgenome.
    library(BSgenome.Hsapiens.UCSC.hg19)
    seqlevelsStyle(vcf) <- "UCSC"
    seqlevelsStyle(txdb) <- "UCSC"
    aa <- predictCoding(vcf, txdb, Hsapiens)

     
predictCoding returns results for coding variants only. As with 
locateVariants, the output has one row per variant-transcript match
so multiple rows per variant are possible.


    ## Did any variants match more than one gene
    table(sapply(split(values(aa)[["GENEID"]], 
          values(aa)[["QUERYID"]]), function(x) length(unique(x)) > 1))

    ## 
    ## FALSE 
    ##    17

    
    ## Summarize the number of variants by gene
    idx <- sapply(split(values(aa)[["QUERYID"]], values(aa)[["GENEID"]], drop = TRUE), 
        unique)
    sapply(idx, length)

    ## 162514  51393   7442 
    ##      6      3      8

    
    ## Summarize variant consequence by gene
    sapply(names(idx), function(nm) {
        d <- aa[values(aa)[["GENEID"]] %in% nm, c("QUERYID", "CONSEQUENCE")]
        table(values(d)[["CONSEQUENCE"]][duplicated(d) == FALSE])
    })

    ##                162514 51393 7442
    ## nonsynonymous       2     0    2
    ## not translated      1     0    5
    ## synonymous          3     3    1

    
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
    library(PolyPhen.Hsapiens.dbSNP131)
    keys <- keys(PolyPhen.Hsapiens.dbSNP131)
    cols <- columns(PolyPhen.Hsapiens.dbSNP131)

    ## column descriptions are found at ?PolyPhenDbColumns
    cols

    ##  [1] "RSID"        "TRAININGSET" "OSNPID"      "OACC"        "OPOS"       
    ##  [6] "OAA1"        "OAA2"        "SNPID"       "ACC"         "POS"        
    ## [11] "AA1"         "AA2"         "NT1"         "NT2"         "PREDICTION" 
    ## [16] "BASEDON"     "EFFECT"      "PPH2CLASS"   "PPH2PROB"    "PPH2FPR"    
    ## [21] "PPH2TPR"     "PPH2FDR"     "SITE"        "REGION"      "PHAT"       
    ## [26] "DSCORE"      "SCORE1"      "SCORE2"      "NOBS"        "NSTRUCT"    
    ## [31] "NFILT"       "PDBID"       "PDBPOS"      "PDBCH"       "IDENT"      
    ## [36] "LENGTH"      "NORMACC"     "SECSTR"      "MAPREG"      "DVOL"       
    ## [41] "DPROP"       "BFACT"       "HBONDS"      "AVENHET"     "MINDHET"    
    ## [46] "AVENINT"     "MINDINT"     "AVENSIT"     "MINDSIT"     "TRANSV"     
    ## [51] "CODPOS"      "CPG"         "MINDJNC"     "PFAMHIT"     "IDPMAX"     
    ## [56] "IDPSNP"      "IDQMIN"      "COMMENTS"

    
    ## Get the rsids for the non-synonymous variants from the predictCoding
    ## results
    rsid <- unique(names(aa)[values(aa)[["CONSEQUENCE"]] == "nonsynonymous"])
    
    ## Retrieve predictions for non-synonymous variants. Two of the six
    ## variants are found in the PolyPhen database.
    select(PolyPhen.Hsapiens.dbSNP131, keys = rsid, cols = c("AA1", "AA2", "PREDICTION"))

    ##       RSID TRAININGSET   OSNPID   OACC OPOS OAA1 OAA2    SNPID    ACC POS
    ## 1 rs224534      humdiv rs224534 Q8NER1  469    T    I rs224534 Q8NER1 469
    ## 2 rs224534      humvar rs224534   <NA>   NA <NA> <NA> rs224534 Q8NER1 469
    ## 3 rs222747      humdiv rs222747 Q8NER1  315    M    I rs222747 Q8NER1 315
    ## 4 rs222747      humvar rs222747   <NA>   NA <NA> <NA> rs222747 Q8NER1 315
    ## 5 rs322937      humdiv rs322937 Q8NET8  117    R    G rs322937 Q8NET8 117
    ## 6 rs322937      humvar rs322937   <NA>   NA <NA> <NA> rs322937 Q8NET8 117
    ## 7 rs322965      humdiv rs322965 Q8NET8   25    I    V rs322965 Q8NET8  25
    ## 8 rs322965      humvar rs322965   <NA>   NA <NA> <NA> rs322965 Q8NET8  25
    ##   AA1 AA2  NT1  NT2        PREDICTION   BASEDON EFFECT PPH2CLASS PPH2PROB
    ## 1   T   I    C    T            benign alignment   <NA>   neutral    0.006
    ## 2   T   I <NA> <NA>            benign      <NA>   <NA>      <NA>    0.022
    ## 3   M   I    G    T            benign alignment   <NA>   neutral    0.000
    ## 4   M   I <NA> <NA>            benign      <NA>   <NA>      <NA>    0.003
    ## 5   R   G    C    G possibly damaging alignment   <NA>   neutral    0.254
    ## 6   R   G <NA> <NA>            benign      <NA>   <NA>      <NA>    0.178
    ## 7   I   V    A    G            benign alignment   <NA>   neutral    0.000
    ## 8   I   V <NA> <NA>            benign      <NA>   <NA>      <NA>    0.000
    ##   PPH2FPR PPH2TPR PPH2FDR SITE REGION PHAT DSCORE SCORE1 SCORE2 NOBS
    ## 1   0.307   0.967   0.387 <NA>   <NA>   NA  0.565  1.108  0.543   46
    ## 2   0.568   0.962      NA <NA>   <NA>   NA     NA     NA     NA   NA
    ## 3   0.996   1.000   0.665 <NA>   <NA>   NA -0.059  1.009  1.068   48
    ## 4   0.821   0.988      NA <NA>   <NA>   NA     NA     NA     NA   NA
    ## 5   0.152   0.889   0.255 <NA>   <NA>   NA  1.794  1.170 -0.624   12
    ## 6   0.375   0.892      NA <NA>   <NA>   NA     NA     NA     NA   NA
    ## 7   0.996   1.000   0.665 <NA>   <NA>   NA -0.101  0.635  0.736   12
    ## 8   1.000   1.000      NA <NA>   <NA>   NA     NA     NA     NA   NA
    ##   NSTRUCT NFILT PDBID PDBPOS PDBCH IDENT LENGTH NORMACC SECSTR MAPREG DVOL
    ## 1      20    NA  <NA>     NA  <NA>    NA     NA      NA   <NA>   <NA>   NA
    ## 2      NA    NA  <NA>     NA  <NA>    NA     NA      NA   <NA>   <NA>   NA
    ## 3      20    NA  <NA>     NA  <NA>    NA     NA      NA   <NA>   <NA>   NA
    ## 4      NA    NA  <NA>     NA  <NA>    NA     NA      NA   <NA>   <NA>   NA
    ## 5      22    NA  <NA>     NA  <NA>    NA     NA      NA   <NA>   <NA>   NA
    ## 6      NA    NA  <NA>     NA  <NA>    NA     NA      NA   <NA>   <NA>   NA
    ## 7      22    NA  <NA>     NA  <NA>    NA     NA      NA   <NA>   <NA>   NA
    ## 8      NA    NA  <NA>     NA  <NA>    NA     NA      NA   <NA>   <NA>   NA
    ##   DPROP BFACT HBONDS AVENHET MINDHET AVENINT MINDINT AVENSIT MINDSIT
    ## 1    NA    NA   <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>
    ## 2    NA    NA   <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>
    ## 3    NA    NA   <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>
    ## 4    NA    NA   <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>
    ## 5    NA    NA   <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>
    ## 6    NA    NA   <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>
    ## 7    NA    NA   <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>
    ## 8    NA    NA   <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>
    ##   TRANSV CODPOS CPG MINDJNC PFAMHIT IDPMAX IDPSNP IDQMIN         COMMENTS
    ## 1      0      1   0    <NA>      NA 20.732  20.73  88.20 chr17:3486702_GA
    ## 2     NA     NA  NA    <NA>      NA     NA     NA     NA chr17:3486702_GA
    ## 3      1      2   0    <NA>      NA 30.843  30.84  93.80 chr17:3493200_CG
    ## 4     NA     NA  NA    <NA>      NA     NA     NA     NA chr17:3493200_CG
    ## 5      1      0   1    <NA>      NA  4.151     NA  93.42 chr17:3446885_TC
    ## 6     NA     NA  NA    <NA>      NA     NA     NA     NA chr17:3446885_TC
    ## 7      0      0   0    <NA>      NA 93.418  93.42  93.42 chr17:3458072_TC
    ## 8     NA     NA  NA    <NA>      NA     NA     NA     NA chr17:3458072_TC


<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>


<h2 id="install-and-use">Installation and Use</h2>

Follow [installation instructions](/install/) to start using these
packages.  To install VariantAnnotation use 


    library(BiocInstaller)
    biocLite("VariantAnnotation")


Package installation is required only once per R installation. View a
full list of available
[software](/packages/release/bioc/)
and 
[annotation](/packages/release/data/annotation/)
packages.

To use the `VariantAnnotation`, evaluate the commands


    library(VariantAnnotation)


These commands are required once in each R session.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="exploring-package-content">Exploring Package Content</h2>

Packages have extensive help pages, and include vignettes highlighting
common use cases. The help pages and vignettes are available from
within R. After loading a package, use syntax like

    help(package="VariantAnnotation")
    ?predictCoding

to obtain an overview of help on the `VariantAnnotation` package, and 
the `predictCoding` function. View the package vignette with


    browseVignettes(package = "VariantAnnotation")


To view vignettes providing a more comprehensive introduction to
package functionality use


    help.start()



<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>


    sessionInfo()

    ## R version 3.0.0 (2013-04-03)
    ## Platform: x86_64-unknown-linux-gnu (64-bit)
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=C                 LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ##  [1] splines   stats4    parallel  stats     graphics  grDevices utils    
    ##  [8] datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] PolyPhen.Hsapiens.dbSNP131_1.0.2       
    ##  [2] BSgenome.Hsapiens.UCSC.hg19_1.3.19     
    ##  [3] BSgenome_1.29.0                        
    ##  [4] cgdv17_0.0.20                          
    ##  [5] TxDb.Hsapiens.UCSC.hg19.knownGene_2.9.2
    ##  [6] GenomicFeatures_1.13.14                
    ##  [7] GGtools_4.9.19                         
    ##  [8] GGBase_3.23.0                          
    ##  [9] snpStats_1.11.0                        
    ## [10] Matrix_1.0-12                          
    ## [11] lattice_0.20-15                        
    ## [12] survival_2.37-4                        
    ## [13] org.Hs.eg.db_2.9.0                     
    ## [14] RSQLite_0.11.3                         
    ## [15] DBI_0.2-7                              
    ## [16] AnnotationDbi_1.23.15                  
    ## [17] Biobase_2.21.2                         
    ## [18] VariantAnnotation_1.7.25               
    ## [19] Rsamtools_1.13.14                      
    ## [20] Biostrings_2.29.3                      
    ## [21] GenomicRanges_1.13.16                  
    ## [22] XVector_0.1.0                          
    ## [23] IRanges_1.19.8                         
    ## [24] BiocGenerics_0.7.2                     
    ## [25] knitr_1.2                              
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] annotate_1.39.0    biomaRt_2.17.1     bit_1.1-10        
    ##  [4] bitops_1.0-5       digest_0.6.3       evaluate_0.4.3    
    ##  [7] ff_2.2-11          formatR_0.7        genefilter_1.43.0 
    ## [10] grid_3.0.0         RCurl_1.95-4.1     rtracklayer_1.21.5
    ## [13] stringr_0.6.2      tools_3.0.0        XML_3.96-1.1      
    ## [16] xtable_1.7-1       zlibbioc_1.7.0


<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>
