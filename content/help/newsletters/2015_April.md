<!-- %\VignetteEngine{knitr::knitr} %\VignetteIndexEntry{} -->

<!-- To deploy this to the website, do the following: library(knitr)
render_markdown(strict=TRUE) knit("2015_April.Rmd")

This will produce an .md file that can be used  on the website. The {:toc} and
{:no_toc} tags will be rendered by nanoc, the tool that produces the website.
-->

{::options parse_block_html="true" /}

# *Bioconductor* Newsletter
{:.no_toc}

posted by [Valerie Obenchain](mailto:vobencha@fhcrc.org), April 2015

On April 17, the release of `Bioconductor` 3.1 will mark the 24th release of the
software. The project started in 2001 with the first svn commits made in May of
that year:

    r3      rgentlem      2001-05-25 14:08:57 -0700 (Fri, 25 May 2001)
    r2      rgentlem      2001-05-25 08:28:31 -0700 (Fri, 25 May 2001)
    r1      (no author)      2001-05-25 08:28:31 -0700 (Fri, 25 May 2001) 

At the time of the first official
[Bioconductor manuscript](http://genomebiology.com/2004/5/10/r80)
in 2004 the project consisted of 

    "... more than 80 software packages, hundreds of metadata
    packages and a number of experimental data packages ..."

Eleven years later (and after more than 100000 svn commits) `Bioconductor` 
hosts over 990 software, 900 annotation and 230 experimental data packages.

Another quote from the 2004 paper shows that, fortunately, not everything has 
changed,

    "... The group dynamic has also been an important factor in the success of 
    Bioconductor. A willingness to work together, to see that cooperation and 
    coordination in software development yields substantial benefits for the 
    developers and the users and encouraging others to join and contribute to 
    the project are also major factors in our success. ..."

This issue looks at the growing role of proteomics in `Bioconductor` and the use
of web sockets to bridge the gap between workspace data and interactive
visualization. We re-visit Docker with use cases in package development and
managing system administration tasks. We also have a section on new and notable 
functions recently added to base `R` and `Bioconductor`.


## Contents 
{:.no_toc}

* Table of contents will replace this text. 
{:toc}


## Proteomics in Bioconductor 

The diversity of proteomics analysis available in `Bioconductor` continues
to grow steadily and the devel branch now hosts 68 proteomic-based software 
packages. Many individuals have contributed to this area in the form of
packages, web-based workflows and course offerings. One very active member is
Laurent Gatto, head of the Computational Proteomics Unit at the Cambridge Centre
for Proteomics. On a day to day basis he is responsible for developing robust
proteomics technologies applicable to a wide variety of biological questions.

Laurent is the author of many `Bioconductor` packages, including the
new
[ProtGenerics](http://www.bioconductor.org/packages/3.1/bioc/html/ProtGenerics.html). 
Similar in concept to the more general `BiocGenerics` package, 
`ProtGenerics` provides a central location where proteomic-specific S4 generics 
can be defined and reused. He has produced course materials and tools to
help newcomers get started including the detailed
[proteomics workflow](http://www.bioconductor.org/help/workflows/proteomics/)
on the web site and a two publications titled
[Using R and Bioconductor for proteomics analysis](http://www.ncbi.nlm.nih.gov/pubmed/23692960)
and
[Visualisation of proteomics data using R and Bioconductor](http://www.ncbi.nlm.nih.gov/pubmed/25690415).
These publications have a companion experimental data package,
[RforProteomics](http://www.bioconductor.org/packages/release/data/experiment/html/RforProteomics.html)
which illustrates data input/output, data processing, quality control,
visualisation and quantitative proteomics analysis within the
`Bioconductor` framework. Since the first release, `RforProteomics`
has benefited from contributions from additional developers.

The `RforProteomics` data package has 4 vignettes:

* The first vignette offers a 
[current perspective](http://bioconductor.org/packages/devel/data/experiment/vignettes/RforProteomics/inst/doc/HUPO2014poster.pdf) 
on proteomics in `Bioconductor` and is in poster format. It gives an overview of
the `Bioconductor` proteomics infrastructure and mass spectrometry analysis.
Topics covered include raw data manipulation, identification, quantitation, MS
data processing, visualization, statistics and machine learning.

* Also in poster format, the 
[RforProteomics BioC2013](http://bioconductor.org/packages/devel/data/experiment/vignettes/RforProteomics/inst/doc/Bioc2013poster.pdf)
vignette is specific to the `Using R and Bioconductor for proteomics analysis`
publication. Special attention is given to labelled vs label-free quantitation 
and `Bioconductor` packages that offer these methods. 

* [Using R and Bioconductor for Proteomics Data Analysis](http://www.bioconductor.org/packages/release/data/experiment/vignettes/RforProteomics/inst/doc/RforProteomics.pdf)
includes code executed in the `Using R and Bioconductor for proteomics 
analysis` publication.

* [Visualisation of proteomics data using R and Bioconductor](http://www.bioconductor.org/packages/release/data/experiment/vignettes/RforProteomics/inst/doc/RProtVis.html)
includes the code from the `Visualisation of proteomics data using R
and Bioconductor` publication.

Many proteomic packages are worthy of mention in the `Bioconductor` 
repository. Here we highlight a few that have played a primary role in
the growing infrastructure. 

Packages `mzR` and `mzID` read and parse raw and identification MS
data. The former is an R interface to the popular C++
[proteowizard](http://www.ncbi.nlm.nih.gov/pubmed/23051804) toolkit.

Identification methods are offered in `rTANDEM` and `MSGFplus` (and
its `shiny` interface `MSGFgui`) and quantitation methods can be found
in `MSnbase` and `isobar` (isobaric tagging and spectral counting
methods), `synapter`, `xcms` and `MALDIquant` packages (label-free).
Statistical modelling and machine learning are offered in `MSnbase`,
`isobar`, `MSstats` and `msmsTests`.

The `rpx` package provides an interface to the
[ProteomeXchange](http://www.ncbi.nlm.nih.gov/pubmed/24727771)
infrastructure, which coordinates multiple data repositories of
MS-based proteomics data.

The `pRoloc` package contains methods for spatial proteomics analysis,
e.g., machine learning and classification methods for assigning a
protein to an organelle.

The relatively new `Pbase` package contains a `Proteins` class for storing and
manipulating protein sequences and ranges of interest. The package has multiple
vignettes on coordinate mapping. One addresses mapping proteins between
different genome builds and the other mapping from protein to genomic 
coordinates.


## Web Sockets

### Overview 

As our ability to generate volumes of sequencing data grows so does the need for
effective visualization tools. Tools that can summarize large quantities of
information into digestible bits and quickly identify unique features or
outliers are important steps in any analysis pipeline. In the `R` world we
have seen an increase in the use of web sockets to provide an interactive link
between data in the workspace and exploration in the browser.

The analysis capabilities of `R` make it a good fit for the rich and interactive
graphics of HTML5 web browsers. The WebSocket protocol enables more interaction
between a browser and a web site, facilitating live content and the creation of
real-time graphics. 

Web sockets are often described as "a standard for bi-directional, full duplex
communication between servers and clients over a single TCP connection". These
characteristics offer several advantages over HTTP:

- 'bi-directional' means either client or server can send a message to the other
  party. HTTP is uni-directional and the request is always initiated by the
  client.

- 'full duplex' allows client and server to talk independently. In the case of
  HTTP, at any given time, either the client is talking or the server is
  talking.

- Web sockets open a 'single TCP connection' over which the client and server
  communicate for the lifecycle of the web socket connection. In contrast, HTTP
  typically opens a new TCP connection for each round trip; a connection is
  initiated for a request and terminated after the response is received. A new
  TCP connection must be established for each request/response.  The opening and
  closing creates overhead, especially in the case where rapid responses or real
  time interactions are needed. 

For those interested, this 
[blog post](http://blog.arungupta.me/rest-vs-websocket-comparison-benchmarks/)
provides more in-depth details and benchmarking against REST.

### Applications in `Shiny` and `epivizr`

The [Shiny package](http://shiny.rstudio.com/) created by the RStudio team
pioneered the use of web sockets in `R`. Shiny enables the building of
interactive web applications from within an `R` session. Popular applications
are interactive plots and maps that allow real-time manipulation through
widgets. The workhorse behind Shiny is the 
[httpuv package](https://github.com/rstudio/httpuv/), also authored by 
RStudio.  httpuv provides low-level socket and protocol support for handling 
HTTP and WebSocket requests within `R`.

The httpuv infrastructure is also used by the 
[epivizr](http://www.bioconductor.org/packages/3.1/bioc/html/epivizr.html) 
package. In this application, web sockets create a two-way communication 
between the `R` environment and the 
[Epiviz visualization tool](http://epiviz.github.io/).
Objects available in an `R` session can be displayed as tracks or plots on
Epiviz. 

### `BrowserViz`

A slightly different approach is taken in the new 
[BrowserViz](http://bioconductor.org/packages/3.1/bioc/html/BrowserViz.html)
package by Paul Shannon. This application provides access to both the browser 
and an active `R` prompt.

The `BrowserViz` package contains the BrowserViz class whose main purpose is to
provide the necessary `R`, Javascript websocket and JSON infrastructure for
communication.  By loosely coupling `R` and the browser the two environments are
linked but kept maximally ignorant of each other; only simple JSON messages pass
back and forth with no HTML, CSS or Javascript.  The result is access to the
interactive graphics of a web browser in conjunction with an active `R` session.

The companion Javascript library, BrowserViz.js, is also included in the
package. The combination of library and base class provides the infrastructure
necessary for any BrowserViz-style application. The `BrowserVizDemo` and `RCyjs`
packages build on the `BrowserViz` class and will be available in the
`Bioconductor` 3.2 release. `BrowserVizDemo` is a minimal example of interactive
plotting and selection of xy points using the popular d3.js library. The more
full featured RCyjs provides interactive access to the full power of
Cytoscape.js, a richly featured browser-based network visualization library.

More details on the `BrowserViz` class and applications can be found in the 
[package vignette](http://www.bioconductor.org/packages/3.1/bioc/vignettes/BrowserViz/inst/doc/BrowserViz.pdf).


## Infrastructure

### Changes in `AnnotationHub`

This quarter Marc and Sonali continued their work on `AnnotationHub`. Several 
new resources were added and the display method and search navigation were 
substantially reworked.

New resources:

* TF-Target Gene Files from the PAZAR public database 
  (available as `GRanges` objects)
* BioPAX files (Level1 and Level2) from the NCI Pathway Interaction Database
  (available as `biopax` objects)
* Background file for ChEA required for the command line version of ChEA
  (available as `data.frame` object)
* GTF files from Ensembl release 76 to 79
  (available as `GRanges` object)
* Expression Set of raw read counts for the GSE62944 dataset from GEO:
  (available as `ExpressionSet` object)

An improved show method and more flexible data retrieval make interacting with
the 18900+ files straightforward. Sonali has a new [AnnotationHub
video](https://www.youtube.com/watch?v=pFvUOPfR8eA&feature=youtu.be) where she 
gives a tour of the resource with tips and tricks for data access.

Code below was generated with `AnnotationHub` version 1.99.75. The show
method now list fields common for subsetting up front, e.g., providers, 
species and class of `R` object.

    > library(AnnotationHub)
    > hub <- AnnotationHub()
    > hub
    AnnotationHub with 18992 records
    # snapshotDate(): 2015-03-12 
    # $dataprovider: UCSC, Ensembl, BroadInstitute, NCBI, Haemcode, dbSNP, Inpar...
    # $species: Homo sapiens, Mus musculus, Bos taurus, Pan troglodytes, Danio r...
    # $rdataclass: GRanges, FaFile, OrgDb, ChainFile, CollapsedVCF, Inparanoid8D...
    # additional mcols(): taxonomyid, genome, description, tags, sourceurl,
    #   sourcetype 
    # retrieve records with, e.g., 'object[["AH169"]]' 
    
                title                                         
      AH169   | Meleagris_gallopavo.UMD2.69.cdna.all.fa       
      AH170   | Meleagris_gallopavo.UMD2.69.dna.toplevel.fa   
      AH171   | Meleagris_gallopavo.UMD2.69.dna_rm.toplevel.fa
      AH172   | Meleagris_gallopavo.UMD2.69.dna_sm.toplevel.fa
      AH173   | Meleagris_gallopavo.UMD2.69.ncrna.fa          
      ...       ...                                           
      AH28575 | A500002_Erg.csv                               
      AH28576 | A500005_Erg.csv                               
      AH28577 | A500001_IgG.csv                               
      AH28578 | A500004_IgG.csv                               
      AH28579 | GSM730632_Runx1.csv    

Tab completion on a hub object lists all fields available for subsetting:

    > hub$
    hub$ah_id         hub$dataprovider  hub$taxonomyid    
    hub$description   hub$rdataclass    hub$sourcetype    
    hub$title         hub$species       hub$genome        
    hub$tags          hub$sourceurl


Quick discovery of file type and provider:

    > sort(table(hub$sourcetype), decreasing=TRUE)
    
              BED         FASTA    UCSC track           GTF NCBI/blast2GO 
             7855          3876          2208          1606          1145 
            Chain           CSV           VCF        BigWig    Inparanoid 
             1113           406           316           315           268 
           TwoBit  BioPaxLevel2         RData        BioPax         GRASP 
              144             6             4             3             1 
           tar.gz           Zip 
                1             1 

    > sort(table(hub$dataprovider), decreasing=TRUE)
    
                                UCSC                          Ensembl 
                                8746                             4590 
                      BroadInstitute                             NCBI 
                                3146                             1145 
                            Haemcode                            dbSNP 
                                 945                              316 
                         Inparanoid8                            Pazar 
                                 268                               91 
    NIH Pathway Interaction Database                        EncodeDCC 
                                   9                                5 
                              RefNet                             ChEA 
                                   4                                1 
                                 GEO                            NHLBI 
                                   1                                1 

Given the volume and diversity of data available in the hub we encourage
using these files as sample data before creating your own experimental data 
package.

For example, to get an idea of available GRCh37 FASTA from Ensembl:

    >  hub[hub$sourcetype=="FASTA" & hub$dataprovider=="Ensembl" & hub$genome=="GRCh37"] 
    AnnotationHub with 42 records
    # snapshotDate(): 2015-03-26 
    # $dataprovider: Ensembl
    # $species: Homo sapiens
    # $rdataclass: FaFile
    # additional mcols(): taxonomyid, genome, description, tags, sourceurl,
    #   sourcetype 
    # retrieve records with, e.g., 'object[["AH18924"]]' 
    
                title                                    
      AH18924 | Homo_sapiens.GRCh37.73.cdna.all.fa       
      AH18925 | Homo_sapiens.GRCh37.73.dna_rm.toplevel.fa
      AH18926 | Homo_sapiens.GRCh37.73.dna_sm.toplevel.fa
      AH18927 | Homo_sapiens.GRCh37.73.dna.toplevel.fa   
      AH18928 | Homo_sapiens.GRCh37.73.ncrna.fa          
      ...       ...                                      
      AH21181 | Homo_sapiens.GRCh37.72.dna_rm.toplevel.fa
      AH21182 | Homo_sapiens.GRCh37.72.dna_sm.toplevel.fa
      AH21183 | Homo_sapiens.GRCh37.72.dna.toplevel.fa   
      AH21184 | Homo_sapiens.GRCh37.72.ncrna.fa          
      AH21185 | Homo_sapiens.GRCh37.72.pep.all.fa        

Advanced developers may be interested in writing a 'recipe' to add 
additional online resources to `AnnotationHub`. The process involves
writing functions to first parse file metadata and then create `R` objects or 
files from these metadata. Detailed HOWTO steps are in the 
[AnnotationHubRecipes vignette](http://bioconductor.org/packages/3.1/bioc/vignettes/AnnotationHub/inst/doc/AnnotationHubRecipes.html).

### `Rhtslib` package

Nate recently completed work on the
[Rhtslib](http://bioconductor.org/packages/3.1/bioc/html/Rhtslib.html)
package which wraps the [htslib](http://www.htslib.org/) C library from
Samtools. The plan is for `Rhtslib` to replace the Samtools code inside
`Rsamtools`.  `Rhtslib` contains a clean branch of htslib directly from
Samtools, including all unit tests. This approach simplifies maintenance when
new versions or bug fixes become available. The clean API also promises to make
outsourcing to the package more straightforward for both `Rsamtools` and others
wanting access to the native routines.

htslib was developed with a 'linux-centric' approach and getting the library
to build across platforms (specifically Windows) was a challenge. To address this,
Nate chose to use [Gnulib](https://www.gnu.org/software/gnulib/), the GNU
portability library.  Briefly, Gnulib is a collection of modules that package
portability code to enable POSIX-compliance in a transparent manner; the goal
being to supply common infrastructure to enable GNU software to run on a variety
of operating systems.  Modules are incorporated into a project at the source
level rather than as a library that is built, installed and linked against.

Incorporating Gnulib involves (at minimum) the following steps:

* adapt the project to use Autoconf and Automake
* identify and import relevant Gnulib modules using gnulib-tool
* add `#include "config.h"` to source files
* remove (now unnecessary!) preprocessor complier/platform tests from source

For more on specific functions available in `Rhtslib` see the 
[Samtools docs](http://www.htslib.org/doc/) or the API-type headers in the
package, faidx.h, hfile.h, hts.h, sam.h, tbx.h and vcf.h. Headers are located
in Rhtslib/src/htslib/htslib or if the package is installed,

    library(Rhtslib)
    system.file(package="Rhtslib", "include")


## Reproducible Research

### AMI course links

The [course materials](http://www.bioconductor.org/help/course-materials/) web
page has links to several resources including slides, presentations and
packages. Recently Dan started adding an "AMI" link for courses that use them.
The AMI contains the packages, sample data and exact version of
`R`/`Bioconductor` used. This is a convenient, portable way to ensure
reproduciblity.  One can imagine using an AMI or Docker container to
capture the state of a research project or publication which can be easily
shared with colleagues.

### Developing with Docker

Elena Grassi is a Ph.D. student in Biomedical Sciences and Oncology in the
Department of Genetics, Biology and Biochemistry at the University of Torino.
Her research focuses on transcriptional and post transcriptional
regulation with special interest in transcription factors and the alternative 
polyadenylation phenomenon. 

With a background in computer science she is involved in developing
computational pipelines and tools and is the author of `Bioconductor` packages
[roar](http://bioconductor.org/packages/3.1/bioc/html/roar.html) 
(preferential usage of APA sites) and 
[MatrixRider](http://bioconductor.org/packages/3.1/bioc/html/MatrixRider.html)
(propensity of binding protein to interact with a sequence). Elena was one of 
the first to try out the Docker containers and found them useful for both
package development and system administration tasks. I asked a few questions
about her experience and got some interesting answers.

**What motivated you to try Docker when developing `MatrixRider`?**

I heard about docker from some friends last year and I was eager to try it.
During the New Year's Eve holidays I decided to start using it with `R` / 
`Bioconductor` to run different versions on our computational server without
adding burden to the sysadmin work. I started with mere curiosity fiddling with
rocker and, eased by the holiday laziness, I stopped with the idea to begin
working on some ad hoc `R` \ `Bioconductor` containers in January. Imagine my
happiness when I read in the newsletter about the brand new Bioconductor docker
containers: they were ready for me :). I decided to use them to develop
MatrixRider as long as I needed to have working versions both for release and
devel. I work on different computers and using the `Bioconductor`
devel_sequencing container freed me completely from the procedure of getting the
source, building, and installing all needed packages.  Besides this advantage
using docker made me sure that the package I was developing did not have any
dependencies on my local system libraries that would not be available in a
clean installation. This was my first package containing C code and it was nice
to be sure.

**Which image did you use, base, core, sequencing, ... ?**

Mainly devel_sequencing to start with a fully-fledged working environment. I
had to install some other packages (TFBSTools and JASPAR2014) and it worked
flawlessly with `biocLite()`.

**Any unanticipated pros/cons of developing in these containers?**

No. I think that I will continue using the devel containers to develop
and maintain packages.

**Describe how Docker was useful for managing multiple `R` versions
  on your computational server. Was this for multiple users or
  just yourself?**

Right now I'm the only one that needs devel so the version management was for
myself. Eventually I would like to set it up on our server and have it working
for multiple users but it will take a little work to integrate it with our
"pipeline management system".

In the past we have had up to three different `R` versions, one from the package
management system of our distribution, Debian, and two compiled ad hoc. Teaching
new students how to reach them and the related library paths has been hard - I
am pretty sure docker will give a huge hand in these situations, helping
also in tracking which versions of packages were used to perform certain
analyses.

## New and Noteworthy

A number of functions added to `R` (3.2) and `Bioconductor` (3.1) this
quarter have potential for wide-spread use. I thought they were worth a
mention.

*   *base::lengths()*

    Computes the element lengths of a `list` object. In `Bioconductor`,
    S4Vectors::elementLengths performs the same operation on `List` objects.
    (contributed by Michael Lawrence)

*   *base::trimws()*

    Removes leading or trailing whitespace from character strings.
    (contributed by Kurt Hornik)

*   *utils::methods()*

    This function previously worked on S3 generics only and has been enhanced 
    to also handle S4. 
    (enhanced by Martin Morgan)

        > library(Rsamtools)
        > methods("scanBam")
        [1] scanBam,BamFile-method    scanBam,BamSampler-method
        [3] scanBam,BamViews-method   scanBam,character-method 
        see '?methods' for accessing help and source code
        Warning message:
        In findGeneric(generic.function, envir) :
          'scanBam' is a formal generic function; S3 methods will not 
           likely be found

        > methods(class = "BamFile")
         [1] $                   $<-                 asMates            
         [4] asMates<-           close               coerce             
         [7] countBam            filterBam           indexBam           
        [10] initialize          isIncomplete        isOpen             
        [13] obeyQname           obeyQname<-         open               
        [16] path                pileup              qnamePrefixEnd     
        [19] qnamePrefixEnd<-    qnameSuffixStart    qnameSuffixStart<- 
        [22] quickBamFlagSummary scanBam             scanBamHeader      
        [25] seqinfo             show                sortBam            
        [28] testPairedEndBam    updateObject        yieldSize          
        [31] yieldSize<-        
        see '?methods' for accessing help and source code

*   *GenomicFeatures::transcriptLengths()*

    Computes transcripts lengths in a TxDb object with the option to
    include / excluded coding and UTR regions.
    (contributed by Herv&eacute; Pag&egrave;s)

*   *BiocParallel::bpvalidate()*

    Flags undefined symbols in functions intended for parallel,
    distributed memory computations.
    (contributed by Martin Morgan, Valerie Obenchain)

*   *BiocInstaller::biocLite()*

    Now capable of installing git repositories. When the 'pkg' argument 
    contains a forward slash, e.g., "myRepo/myPkg", it is assumed to be a
    repository and is installed with devtools::install_github.
    (contributed and enhanced by Martin Morgan)


## Project Statistics 

### Website traffic

The following compares the number of sessions and new users from the first
quarter of 2015 (January 1 - March 30) with the first quarter of 2014. Sessions
are broken down by new and returning visitors. New visitors correspond to the
total new users.

<table border="0" cellpadding="5" cellspacing="0">
 <caption><b>First Quarter Website Traffic 2015 vs 2014</b></caption>
  <tbody valign="top">
    <tr>
        <td><b>Sessions: Total</b></td>
        <td>24.03%</td>
        <td>(339,283 vs 273,559)</td>
    </tr>
    <tr>
        <td><b>Sessions: Returning Visitor</b></td>
        <td>21.42%</td>
        <td>(213,848 vs 176,128)</td>
    </tr>
    <tr>
        <td><b>Sessions: New Visitor</b></td>
        <td>28.74%</td>
        <td>(125,435 vs 97,431)</td>
    </tr>
    <tr>
        <td><b>New Users</b></td>
        <td>28.74%</td>
        <td>(125,435 vs 97,431)</td>
    </tr>
  </tbody>
</table>

<br/>
Statistics generated with [Google Analytics](http://www.google.com/analytics/).

### Package downloads and new submissions 

The number of unique IP downloads of software packages for January, February and
March of 2015 were 31720, 31956, and 38379, respectively.  For the same time
period in 2014, numbers were 29690, 28993 and 34634. Numbers must be
compared by month (vs sum) because some IPs are the same between months.
See the web site for a full summary of [download
stats](http://bioconductor.org/packages/stats/).

A total of 55 software packages were added in the first quarter of 2015 bringing counts to 991 in devel (`Bioconductor` 3.2) and 936 in release 
(`Bioconductor` 3.1).


## Upcoming Events

See the [events page](http://www.bioconductor.org/help/events/) for a listing
of all courses and conferences.

* [Use R / Bioconductor for Sequence Analysis](https://register.bioconductor.org/Seattle-Apr-2015/):
This intermediate level course is offered 06 - 07 April in Seattle, WA
USA.

* [Advanced RNA-Seq and ChiP-Seq Data Analysis](http://www.ebi.ac.uk/training/course/RNA2015):
Held in Hinxton, UK at EMBL-EBI, 11 - 14 May. 

* [CSAMA 2015 - Statistics and Computing in Genome Data Science](http://www-huber.embl.de/csama/):
Held the 14 - 19 of June in Bressanone-Brixen, Italy.

* [BioC2015](http://www.bioconductor.org/help/course-materials/2015/BioC2015/):
This year the 20 - 22 of July in Seattle, WA, USA.


## Acknowledgements 

Thanks to Laurent Gatto, Elena Grassi and Paul Shannon for contributing to
the Proteomics, Docker and Web Sockets sections. Also thanks to the 
`Bioconductor` team in Seattle for project updates and editorial review.


Send comments or questions to Valerie at 
[vobencha@fredhutch.org](vobencha@fredhutch.org).
