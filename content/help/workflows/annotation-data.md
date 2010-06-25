![](/images/icons/help.gif)Using Bioconductor for Annotation
============================================================

Bioconductor has extensive facilities for mapping between microarray
probe, gene, pathway, gene ontology, homology and other annotations.

Bioconductor has built-in representations of GO, KEGG, vendor, and
other annotations, and can easily access NCBI, Biomart, UCSC, and
other sources.

## Sample Work Flow ##

The following psuedo-code illustrates a typical R / Bioconductor
session. It continues the
[differential expression](/help/workflows/oligo-arrays/) work flow,
taking a 'top table' of differentially expressed probesets and
discovering the genes probed, and the Gene Ontology pathways to which
they belong.

    ## Affymetrix U133 2.0 array IDs of interest; these might be
    ## obtained from
    ##
    ##   tbl <- topTable(efit, coef=2)
    ##   ids <- tbl[["ID"]]
    ##
    ## as part of a more extensive work flow.
    > ids <- c("39730_at", "1635_at", "1674_at", "40504_at", "40202_at")
             
    ## load libraries as sources of annotation
    > library("hgu95av2.db")
    
    ## map probe ids to ENTREZ gene ids...
    > entrez <- hgu95av2ENTREZID[ids]
    > toTable(entrez)
      probe_id gene_id
    1  1635_at      25
    2  1674_at    7525
    3 39730_at      25
    4 40202_at     687
    5 40504_at    5445
    ## ... and to GENENAME
    > genename <- hgu95av2GENENAME[ids]
    ## ... and merge results
    > merge(toTable(entrez), toTable(genename))
      probe_id gene_id                                          gene_name
    1  1635_at      25         c-abl oncogene 1, receptor tyrosine kinase
    2  1674_at    7525 v-yes-1 Yamaguchi sarcoma viral oncogene homolog 1
    3 39730_at      25         c-abl oncogene 1, receptor tyrosine kinase
    4 40202_at     687                              Kruppel-like factor 9
    5 40504_at    5445                                      paraoxonase 2
    
    ## find and extract the GO ids associated with the first id
    > goIds <- mappedRkeys(hgu95av2GO[ids[1]])

    ## use GO.db to find the Terms associated with the goIds, displaying
    ## the head (first six entries) of the result as a data frame
    > library("GO.db")
    > head(as.data.frame(Term(goIds)))
                                                                         Term(goIds)
    GO:0000115 regulation of transcription involved in S-phase of mitotic cell cycle
    GO:0006298                                                       mismatch repair
    GO:0006355                            regulation of transcription, DNA-dependent
    GO:0006464                                          protein modification process
    GO:0007155                                                         cell adhesion
    GO:0007165                                                   signal transduction

## Installation and Use ##

Follow [installation instructions]("/install/"") to start using these
packages.  To install the annotations associated with the Affymetrix
Human Genome U95 V 2.0, and with Gene Ontology, use

    > source("http://bioconductor.org/biocLite.R")
    > biocLite(c("hgu95av2.db", "GO.db"))

Package installation is required only once per R installation. View a
full list of available
[software](http://bioconductor.org/packages/release/Software.html)
and 
[annotation](http://bioconductor.org/packages/release/AnnotationData.html)
packages.

To use the `AnnotationDbi` and `GO.db` package, evaluate the commands

    > library(AnnotationDbi")
    > library("GO.db")

These commands are required once in each R session.

## Exploring Package Content ##

Packages have extensive help pages, and include vignettes highlighting
common use cases; the help pages and vignettes are available from
within R. After loading a package, use syntax like

    > help(package="GO.db")
    > ?GOTERM

to obtain an overview of help on the `GO.db` package, and the `GOTERM`
mapping.  The `AnnotationDbi` package is used by most `.db`
packages. View the vignettes in the `AnnotationDbi` package with

    > browseVignettes(package="AnnotationDbi")

To view vignettes (providing a more comprehensive introduction to
package functionality) in the `AnnotationDbi` package. Use

    > help.start()

To open a web page containing comprehensive help resources.

## Annotation Resources ##

The following guides the user through key annotation packages.  Users
interested in how to create custom chip packages should see the
vignettes in `AnnotationDbi` package. There is additional information
in the annotate package for how to use some of the extra tools
provided.

### Key Packages ###

*
  [AnnotationDbi](http://bioconductor.org/packages/release/bioc/html/AnnotationDbi.html)
  Almost all annotations require the `AnnotationDbi` package. This
  package will be automatically installed for you if you install
  another ".db" annotation package using biocLite(). It contains the code to
  allow annotation mapping objects to be made and manipulated as well
  as code to generate custom chip platforms.
* [Category](http://bioconductor.org/packages/release/bioc/html/Category.html)
  This is the base level package for dealing with annotation questions
  that involve categorical data.
* [GOstats](http://bioconductor.org/packages/release/bioc/html/GOstats.html)
  This builds on what is found in Category so that you can do
  hypergeometric testing using the Gene Ontology found in the GO.db
  package.
* [annotate](http://bioconductor.org/packages/release/bioc/html/annotate.html)
  This package contains many helpful tools for making use of
  annotations.
* [biomaRt](http://bioconductor.org/packages/release/bioc/html/biomaRt.html)
  This package is a great way to pull annotation data directly from
  web based annotation resources. Such data is extremely "current", so
  it is a good idea to save and locally manage the data that you pull
  down from biomaRt so that your code will be reproducible.

### Types of Annotation Packages ###

* Organism annotation packages contain all the gene data for an entire
  organism. All Organism packages are named like this:
  org."Xx"."yy".db. Where "Xx" is the abbreviation for Genus and
  species. And "yy" is the source of the central ID that is used to
  tie all the data together. Some examples are:
  [org.Hs.eg.db](http://www.bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)
  which is for Homo sapiens and is based upon Entrez Gene IDs. And
  [org.At.tair.db](http://www.bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html)
  which is for Arabidopsis thaliana and is based on the tair IDs.
* There are also packages for questions about general systems biology
  data. Some examples of this are:
  [KEGG.db](http://www.bioconductor.org/packages/release/data/annotation/html/KEGG.db.html)
  for accessing data that pertains to Kyoto Encyclopedia of Genes and
  Genomes. [GO.db for](http://www.bioconductor.org/packages/release/data/annotation/html/GO.db.html)
  accessing data that pertains to the Gene
  Ontology. [PFAM.db](http://www.bioconductor.org/packages/release/data/annotation/html/PFAM.db.html)
  for accessing data that pertains to different protein family
  identifiers and how they relate to each other.
* Chip annotation packages are for accessing only the data from one
  specific platform at a time. These packages are named like this:
  "platformName".db.  Where "platformName" is the name of the chip
  platform that these packages refer to. And example would be
  [hgu95av2.db](http://www.bioconductor.org/packages/release/data/annotation/html/hgu95av2.db.html)
  which is for the hgu95av2 platform from Affymetrix.
* Inparanoid homology packages are for accessing inparanoid
  data. hom."Xx".inp.db Where "Xx" is the abbreviation for Genus and
  species. An example is
  [hom.Hs.inp.db](http://www.bioconductor.org/packages/release/data/annotation/html/hom.Hs.inp.db.html)
  which contains inparanoid based mapping data between genes for Homo
  sapiens and 35 other organisms.</li>
* .db0 packages are for making custom platform specific
  packages. These packages are named like this: "name".db0. Where
  "name" is the name of the organism being represented. A list of the
  available .db0 packages can be obtained by calling
  available.db0Pkgs(). There is one of these for each supported
  organism. An example would be
  [human.db0](http://www.bioconductor.org/packages/release/data/annotation/html/Human.db0.html). Users
  should not need these installed unless they plan to make custom chip
  packages according the guidelines in the SQLForge vignette that is
  included with the `AnnotationDbi` package.  These packages must be
  upgraded before you attempt to update your custom chip packages as
  they contain the source databases needed by the SQLForge code.
* For relevant Affymetrix platforms you may also want the cdf and
  probe packages for that platform.  These packages are named using
  the following convention: "platformName"cdf and
  "platformName"probe. Where "platformName" is the name of the chip
  platform that these packages refer to.
