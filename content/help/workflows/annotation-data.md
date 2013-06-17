![](/images/icons/help.gif)Using Bioconductor for Annotation
============================================================

Bioconductor has extensive facilities for mapping between microarray
probe, gene, pathway, gene ontology, homology and other annotations.

Bioconductor has built-in representations of GO, KEGG, vendor, and
other annotations, and can easily access NCBI, Biomart, UCSC, and
other sources.

* [Sample Workflow](#sample-workflow)  
* [Installation and Use](#install-and-use)
* [Exploring Package Content](#exploring-package-content)
* [Annotation Resources](#annotation-resources)

<h2 id="package-type">Package Types</h2> Bioconductor contains many
different types of annotation packages.  You can browse the currently
available types here
[here](http://www.bioconductor.org/packages/release/BiocViews.html#___PackageType)
by simply using the bioconductor web site.

You will see that there are packages that contain annotation data
about a particular microarray platform (ChipDb), there are packages
that contain gene centered data about an organism (OrgDb), and even
packages that contain genome centered data about an organisms
transcriptome (TranscriptDb).  This document will talk about typical
uses for most of these more popular kinds of annotation package.  As
well as describe a newer meta package that wraps access to several
different kinds of packages (OrganismDb).


<h2 id="sample-workflow-ChipDb">Sample ChipDb Workflow</h2>

The following examples illustrates a typical R / Bioconductor
session using a ChipDb style package for information about a specific
type of microarray. It continues the [differential
expression](/help/workflows/oligo-arrays/) workflow, taking a 'top
table' of differentially expressed probesets and discovering the genes
probed, and the Gene Ontology pathways to which they belong.

    ## Affymetrix U133 2.0 array IDs of interest; these might be
    ## obtained from
    ##
    ##   tbl <- topTable(efit, coef=2)
    ##   ids <- tbl[["ID"]]
    ##
    ## as part of a more extensive workflow.
    > ids <- c("39730_at", "1635_at", "1674_at", "40504_at", "40202_at")
             
    ## load libraries as sources of annotation
    > library("hgu95av2.db")
    
    ## To list the kinds of things that can be retrieved, use the cols method.
    > cols(hgu95av2.db)
    
    ## To list the kinds of things that can be used as keys 
    ## use the keytypes method
    > keytypes(hgu95av2.db)

    ## To extract viable keys of a particular kind, use the keys method.
    > head(keys(hgu95av2.db, keytype="ENTREZID"))
    
    ## the select method allows you to mao probe ids to ENTREZ gene ids...
    > select(hgu95av2.db, ids, "ENTREZID", "PROBEID")
       PROBEID ENTREZID
    1 39730_at       25
    2  1635_at       25
    3  1674_at     7525
    4 40504_at     5445
    5 40202_at      687

    ## ... and to GENENAME etc.
    > select(hgu95av2.db, ids, c("ENTREZID","GENENAME"), "PROBEID")
       PROBEID ENTREZID                                           GENENAME
    1 39730_at       25     c-abl oncogene 1, non-receptor tyrosine kinase
    2  1635_at       25     c-abl oncogene 1, non-receptor tyrosine kinase
    3  1674_at     7525 v-yes-1 Yamaguchi sarcoma viral oncogene homolog 1
    4 40504_at     5445                                      paraoxonase 2
    5 40202_at      687                              Kruppel-like factor 9
    
    ## find and extract the GO ids associated with the first id
    > res <- select(hgu95av2.db, ids[1], "GO", "PROBEID")
    > head(res)
       PROBEID         GO EVIDENCE ONTOLOGY
    1 39730_at GO:0000115      TAS       BP
    2 39730_at GO:0000287      IDA       MF
    3 39730_at GO:0003677      NAS       MF
    4 39730_at GO:0003785      TAS       MF
    5 39730_at GO:0004515      TAS       MF
    6 39730_at GO:0004713      IDA       MF

    ## use GO.db to find the Terms associated with those GOIDs
    > library("GO.db")
    > head(select(GO.db, res$GO, "TERM", "GOID"))
            GOID                                                                   TERM
    1 GO:0000115  regulation of transcription involved in S phase of mitotic cell cycle
    2 GO:0000287                                                  magnesium ion binding
    3 GO:0003677                                                            DNA binding
    4 GO:0003785                                                  actin monomer binding
    5 GO:0004515                     nicotinate-nucleotide adenylyltransferase activity
    6 GO:0004713                                       protein tyrosine kinase activity
    
<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>



<h2 id="sample-workflow-OrgDb">Sample OrgDb Workflow</h2>
The organism wide gene centered packages (OrgDb packages) all contain
gene centered data for an organism.  These packages are accessed in
much the same way as the platform based (ChipDb) packages previously
discussed.  But because they are general, they don't contain
infromation like probe IDs that would relate to any specific platform.

But the important thing to understand is that the same methods apply.
So for example you can look up information in this way:

    > library(org.Hs.eg.db)
    > keys <- head(keys(org.Hs.eg.db), n=2)
    > cols <- c("PFAM","GO", "SYMBOL")
    > select(org.Hs.eg.db, keys, cols, keytype="ENTREZID")


<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>



<h2 id="sample-workflow-TranscriptDb">Sample TranscriptDb Workflow</h2>

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>



<h2 id="sample-workflow-OrganismDb">Sample OrganismDb Workflow</h2>

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>



<h2 id="sample-workflow-AnnotationHub">Sample AnnotationHub Workflow</h2>

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>



<h2 id="install-and-use">Installation and Use</h2>

Follow [installation instructions](/install/) to start using these
packages.  To install the annotations associated with the Affymetrix
Human Genome U95 V 2.0, and with Gene Ontology, use

    > source("http://bioconductor.org/biocLite.R")
    > biocLite(c("hgu95av2.db", "GO.db"))

Package installation is required only once per R installation. View a
full list of available
[software](/packages/release/bioc/)
and 
[annotation](/packages/release/data/annotation/)
packages.

To use the `AnnotationDbi` and `GO.db` package, evaluate the commands

    > library(AnnotationDbi")
    > library("GO.db")

These commands are required once in each R session.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="exploring-package-content">Exploring Package Content</h2>

Packages have extensive help pages, and include vignettes highlighting
common use cases. The help pages and vignettes are available from
within R. After loading a package, use syntax like

    > help(package="GO.db")
    > ?select

to obtain an overview of help on the `GO.db` package, and the `select`
method.  The `AnnotationDbi` package is used by most `.db`
packages. View the vignettes in the `AnnotationDbi` package with

    > browseVignettes(package="AnnotationDbi")

To view vignettes (providing a more comprehensive introduction to
package functionality) in the `AnnotationDbi` package. Use

    > help.start()

To open a web page containing comprehensive help resources.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="annotation-resources">Annotation Resources</h2>

The following guides the user through key annotation packages.  Users
interested in how to create custom chip packages should see the
vignettes in the `AnnotationForge` package. There is additional
information in the `AnnotationDbi`, `OrganismDbi` and
`GenomicFeatures` packages for how to use some of the extra tools
provided. You can also refer to the [complete list of annotation
packages](/packages/release/BiocViews.html#___AnnotationData).

### Key Packages ###

* [AnnotationDbi](/packages/release/bioc/html/AnnotationDbi.html)
  Almost all annotations require the `AnnotationDbi` package. This
  package will be automatically installed for you if you install
  another ".db" annotation package using biocLite(). It contains the
  code to allow annotation mapping objects to be made and manipulated
  as well as code to use the select methods etc..
* [OrganismDbi](/packages/release/bioc/html/.html)\\\
  OrganismDbi allows meta packages that enable the user to access
  resources from several different packages as if they were coming
  from one place.  So for example the Homo.sapiens package is enabled
  by OrganismDbi and allows the user to get access to GO.db, the
  associated organism package IDs and the related transcript data for
  the hg19 build of the human genome all as if it were contained in a
  single convenient object.
* [GenomicFeatures](/packages/release/bioc/html/.html)
  GenomicFeatures allows the existance of TranscriptDb objects and
  allows convenient representation of ranges from Transcritomes.
  There are accessors for things like exons, transcripts as well as
  the select method for retrieving data from packages supported by
  GenomicFeatures.
* [AnnotationForge](/packages/release/bioc/html/.html)
  AnnotationForge documents and assists in the creation of some kinds
  of custom annotation packages. 
* [Category](/packages/release/bioc/html/Category.html)
  This is the base level package for dealing with annotation questions
  that involve categorical data.
* [GOstats](/packages/release/bioc/html/GOstats.html)
  This builds on what is found in Category so that you can do
  hypergeometric testing using the Gene Ontology found in the GO.db
  package.
* [annotate](/packages/release/bioc/html/annotate.html)
  This package contains many helpful tools for making use of
  annotations.
* [biomaRt](/packages/release/bioc/html/biomaRt.html)
  This package is a great way to pull annotation data directly from
  web based annotation resources. Such data is extremely "current", so
  it is a good idea to save and locally manage the data that you pull
  down from biomaRt so that your code will be reproducible.

### Types of Annotation Packages ###

* Organism annotation packages contain all the gene based data for an entire
  organism. All Organism packages are named like this:
  org."Xx"."yy".db. Where "Xx" is the abbreviation for Genus and
  species. And "yy" is the source of the central ID that is used to
  tie all the data together. Some examples are:
  [org.Hs.eg.db](/packages/release/data/annotation/html/org.Hs.eg.db.html)
  which is for Homo sapiens and is based upon Entrez Gene IDs. And
  [org.At.tair.db](/packages/release/data/annotation/html/org.At.tair.db.html)
  which is for Arabidopsis thaliana and is based on the tair IDs.
* TransriptDb packages contain range and chromosome information for
  specific transcriptomes.  These are based on a particular genome
  build and are are the place to look for where a
  gene/transcripts/exon coordinate information is relative to a
  genome.  These are also named in a way that tells you about where
  the data came from and can be generated with functions contained in
  the GenomicFeatures package.
* OrganismDb packages are named for the species they represent (such
  as the Homo.sapiens package).  These packages contain references to
  other key annotations packages and can thus represent all the
  underlying data as if it were coming from one place.  So for
  example, the Homo.sapiens package can allow you to retrieve data
  about the ranges of a genes transcripts at the same time that you
  extract it's gene name because it represents both a the
  transcriptome and the relevant org package for Homo sapiens.  These
  can be generated using functions in the OrganismDbi package if you
  have specific packages that you want to link together.
* There are also packages for questions about general systems biology
  data. Some examples of this are:
  [KEGG.db](/packages/release/data/annotation/html/KEGG.db.html)
  for accessing data that pertains to Kyoto Encyclopedia of Genes and
  Genomes. [GO.db](/packages/release/data/annotation/html/GO.db.html) for
  accessing data that pertains to the Gene
  Ontology. [PFAM.db](/packages/release/data/annotation/html/PFAM.db.html)
  for accessing data that pertains to different protein family
  identifiers and how they relate to each other.
* Chip annotation packages are for accessing only the data from one
  specific platform at a time. These packages are named like this:
  "platformName".db.  Where "platformName" is the name of the chip
  platform that these packages refer to. And example would be
  [hgu95av2.db](/packages/release/data/annotation/html/hgu95av2.db.html)
  which is for the hgu95av2 platform from Affymetrix.
* Inparanoid homology packages are for accessing inparanoid
  data. hom."Xx".inp.db Where "Xx" is the abbreviation for Genus and
  species. An example is
  [hom.Hs.inp.db](/packages/release/data/annotation/html/hom.Hs.inp.db.html)
  which contains inparanoid based mapping data between genes for Homo
  sapiens and 35 other organisms.</li>
* .db0 packages are for making custom platform specific
  packages. These packages are named like this: "name".db0. Where
  "name" is the name of the organism being represented. A list of the
  available .db0 packages can be obtained by calling
  available.db0Pkgs(). There is one of these for each supported
  organism. An example would be
  [human.db0](/packages/release/data/annotation/html/human.db0.html). Users
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

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>
