# ![](/images/icons/magnifier.gif)Mentored Projects

A mentored Bioconductor software development project is one in which
experienced programmers work with volunteers to develop new capabilities
needed by the community. <a href=#introduction>More...</a>


## Projects Needing Volunteers
  * [Extending mzR](#extendingMZR)
  * [Expanding MotifDb](#expandingMotifDb)
  * [Galaxy-ification of Useful Scripts](#galaxy)

## Projects in Progress
  * [Add DEXSeq Functionality to easyRNASeq](#easyrnaseq)

## Get Help With Your Own Project
  * [Finish your own package](#your-package)

## Contact Us
  * [Interested in Helping or in Receiving Help?](#contact)


## Completed Projects
  * [Add constructors to the graph package](#graph)
  * [VCF Allele Frequencies](#VCF_alleleFrequency)
  * [VCF Genotypes to Probability-Based SNP Encoding](#VCF_probabilityBasedSnpEncoding)
  * [msGUI](#msGUI)
  * [Create an AnnotationDbi Package for PANTHER](#panther)


<a name="introduction"></a>
## Introduction

A mentored Bioconductor software development project is one in which
experienced programmers work with volunteers to develop new capabilities
needed by the community.

Developers new to Bioconductor may find mentored projects a useful way
to apply, refine and extend their skills.  Projects are identified by
experienced Bioconductor developers.  The projects involve important
but manageable programming tasks.  Experienced developers act as
mentors, providing guidance and oversight.  Successful mentored
projects will be incorporated into the appropriate packages, and
contributors will receive full credit for their work.  Users,
contributors and mentors will all benefit.

We anticipate that mentored projects will usually be run by one or two
experienced Bioconductor-savvy programmers who provide guidance, usually
remotely, to one or more less-experienced programmers.  All the tools
of 'social coding' -- from email and svn to github and skype -- can be
used, at the discretion of the participants.  Except in unusual
circumstances, we expect that participants will have their own
independent funding, most likely as the result of a good fit between the
mentored project and their current employment or academic studies.


Below you will find a list of proposed projects.  We invite your participation.  We welcome your suggestions.

<h2 id='extendingMZR'>Extending mzR</h2>
The `mzR` R/Bioconductor package provides a unified API to the common open and community-driven file formats and parsers available for mass spectrometry data, namely `mzXML`, `mzML` and `mzData` (see [vignette](/packages/devel/bioc/vignettes/mzR/inst/doc/mzR.pdf) for details). It uses `C` and `C++` code from other third party open-source projects and heavily relies on the [`Rcpp`](http://dirk.eddelbuettel.com/code/rcpp.html) package to, notably, provide a direct mapping from `R` to `C++` infrastructure.
Currently, `mzR` provides two actual backends to read Mass Spectrometry raw data:

1. `netCDF` which reads, as the name implies, `netCDF` data
2. `RAMP` to read `mzData` and `mzXML` via the ISB `RAMP` parser. This backend can also read `mzML` through the proteowizard `RAMPadapter` around the proteowizard infrastructure, but this interface is limited to the lowest common denominator between the `mzXML`/`mzData`/`mzML` formats.

This project is intended to add several related backends to `mzR`, by providing a direct wrapper around -- and full access to -- the proteowizard `msdata` object. The candidate will interact closely with [Laurent Gatto](https://github.com/lgatto) and [Steffen Neumann](https://github.com/sneumann), and the [proteowizard](http://www.ncbi.nlm.nih.gov/pubmed/23051804)  and `Rcpp` communities.

Project attributes and estimates:


* Difficulty: medium to difficult, depending on experience and `C++` fluency.
* Length: 3 months.
* Skills needed: intermediate R programming, knowledge of package development helpful, good knowledge of `C` and especially `C++` essential. The candidate will have to familiarise herself with the mass-spectrometry data, the respective data formats and the [proteowizard code base](http://proteowizard.sourceforge.net/dox/index.html).
* Deliverable: pwiz and identificaiton backends to be added to the `mzR` package.
* Mentors: Laurent Gatto and Steffen Neuman, with additional Rcpp support from Dirk Eddelbuettel.
* More information:   <a href="https://github.com/sneumann/mzR/wiki/Extending-mzR">github wiki</a>

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>



<h2 id='expandingMotifDb'>Expanding MotifDb</h2>
The `MotifDb` R/Bioconductor package provides unified access to (currently) seven transcription
factor binding site motif collections, covering 22 organisms.  We wish to expand its holdings by
adding <b>JASPAR 2014</b> and <b>HOCOMOCO</b>, and improving the annotation we offer for 
<b>stamlab</b>.

Project attributes and estimates:

* Difficulty: medium.
* Length: 1-6 weeks if full-time
* Skills needed: intermediate R programming, basics of gene regulation (or willingness to learn)
* Deliverable: new import scripts producing serialized data files for the package
* Mentors: Paul Shannon

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>


<h2 id="galaxy">Galaxy-ification of Useful Scripts</h2>
From the <a href="http://en.wikipedia.org/wiki/Galaxy_%28computational_biology%29">Wikipedia entry for Galaxy</a>:

<blockquote>
Galaxy is a scientific workflow, data integration, and data and
analysis persistence and publishing platform that aims to make
computational biology accessible to research scientists that do not
have computer programming experience. Although it was initially
developed for genomics research, it is largely domain agnostic and is
now used as a general bioinformatics workflow management system.
</blockquote>
The new Bioconductor <a href="http://www.bioconductor.org/packages/devel/bioc/html/RGalaxy.html">RGalaxy</a>
 simplifies the process of exposing an R function in Galaxy so that a user
can run the function using nothing more than a web browser.

This project would involve taking an existing workflow (or conceiving a new workflow)
and exposing it in Galaxy.

Project attributes and estimates:

* Difficulty: Medium
* Length: 3-4 weeks
* Skills needed: R programming
* Mentor: Dan Tenenbaum
* Status: Awaiting volunteers

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>



## <a id="easyrnaseq"></a>Add DEXSeq Functionality to easyRNASeq

The easyRNASeq package facilitates and expedites the processing and
filtering of large RNA-seq datasets for subsequent analysis by
Bioconductor packages edgeR and DESseq, which are concerned with
gene-expression and alternative splicing, respectively.  We propose to
add an output format compatible with DEXSeq, a package for exon-level
differential expression analysis.

Project attributes and estimates:

 * Difficulty: medium/advanced
 * Length: 25 weeks, part-time
 * Skills needed: intermediate R programming; some familiarity with R packages and S4 classes. Familiarity with next generation sequencing and the core RNA-seq codes in Bioconductor
 *  Deliverables: An updated, upgraded package, submitted to Bioconductor.
 *  Funding: Self-funded
 *  Mentee: Vincent Zimmern
 *  Mentor: Nicolas Delhomme
 *  Status: imminent (January 2013)
 *  More information: github repo coming soon

<p class="back_to_top"> [ <a href="#top">Back to top</a> ]</p>






## <a id="your-package"></a>Get Help With Your Package

Package authors sometimes have excellent statistical and bioinformatic
ideas, but are not fully confident in their ability to produce a
robust software package suitable for inclusion in Bioconductor. This
mentored project pairs the package developer with an experienced
programmer to produce quality software. Participants are expected to
have a working version of their package, with the major ideas and
preliminary implementation complete.

Project attributes and estimates:

* Difficulty: medium
* Length: 6-8 weeks
* Skills needed: intermediate R programming; some familiarity with R
  packages and S4 classes.
* Mentor: Various
* Deliverables: A finished package, submitted to Bioconductor.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>


## <a id="orphaned-package"></a>Take Over an 'Orphaned' Package

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

Sometimes the maintainer of an older Bioconductor package is no longer
able to perform that job.  These older packages remain useful but
occasionally need a bug fix or a small change.  We are looking for
volunteers to maintain such packages -- which would otherwise be
abandonded.  Relatively little work is required, the original author
will be available to answer questions, the Bioconductor core team can help,
and the Bioconductor community will benefit.

Current orphans are listed below.
<br>

(No orphans at this time)


<h2 id="contact">Interested?</h2>

Please send mail to pshannon AT fhcrc DOT org if you would like to
help out on any of these projects, or have an idea of your own which
you wish to propose.


<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>






## <a id="completed"></a>Completed projects

<h2 id='panther'> Create an AnnotationDbi Package for PANTHER</h2>

We would like to see PANTHER annotation contained in a Bioconductor AnnotationDbi
package.

PANTHER is found <a href=http://www.pantherdb.org/>here</a>, and summarized:

<blockquote>
The PANTHER (Protein ANalysis THrough Evolutionary Relationships)
Classification System is a unique resource that classifies genes by
their functions, using published scientific experimental evidence and
evolutionary relationships to predict function even in the absence of
direct experimental evidence.  "classifies genes by their function"
</blockquote>

Project attributes and estimates:

* Difficulty: medium
* Length: 6-8 weeks, part-time
* Skills needed:  Familiarity with R, SQL and Panther.
* Mentor: Marc Carlson
* Status: Complete
<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>


### <a id="graph"></a>Add contructors for key classes in the graph package

The [graph package] was developed when users created objects with
calls like `new("graphNEL")`, but there are advantages to hiding this
level of implementation from the user and instead creating a new
instance with `graphNEL()`. The project modernizes this aspects of the
graph package.

* Difficulty: easy
* Length: 1 week
* Skills needed: basic R programming, some familiarity with package
  structure and documentation.
* Mentor: Martin Morgan
* Deliverables: Unit tests, code, and revised documentation for
  constructors `graphNEL()` and `graphAM()`; constructors for
  additional classes may also be provided, e.g., `attrData()`,
  `clusterGraph()`, `distGraph()`, `edgeSet()`, `edgeSetAM()`,
  `edgeSetNEL()`, `renderInfo()`, `simpleEdge()`.
* Mentee: Paul Shannon
* Status: Complete

[graph package]: /packages/devel/bioc/html/graph.html

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

### <a id="VCF_alleleFrequency"></a>VCF Allele Frequencies

The VariantAnnotation package needed a function to compute genotype counts, 
allele frequencies and Hardy-Weinberg estimates from the genotype data in a 
VCF class. 

Project attributes and estimates:

* Difficulty: Easy
* Length: 2 weeks
* Skills needed: R programming, familiarity with S4 classes
* Deliverables: Implement, test and document snpSummary,CollapsedVCF-method
* Mentor: Valerie Obenchain
* Mentee: Chris Wallace
* Status: Complete 
* github: <a
href='https://github.com/Bioconductor/VCF_projects'>VCF_alleleFrequency page</a> and <a
href='https://github.com/Bioconductor/VCF_projects/wiki/VCF_alleleFrequency'>wiki</a>.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

### <a id="VCF_probabilityBasedSnpEncoding"></a>VCF Probability-Based SNP Encoding of
Genotypes

MatrixToSnpMatrix() in the VariantAnnotation package
converts the genotype data in a VCF object into a SnpMatrix 
object. Currently this is done without taking uncertain
uncertain genotype calls into consideration. This project 
involves modifying MatrixToSnpMatrix() to use, when available, 
genotype uncertainty and likelihood information to convert 
genotypes to probability-based SnpMatrix encodings.

Project attributes and estimates:

* Difficulty: Advanced 
* Length: 12 weeks
* Skills needed: R programming, familiarity with S4 classes, statistics
* Deliverables: See wiki page for details 
* Mentors: Valerie Obenchain and Vince Carey
* Status: Complete 
* Mentee: Stephanie Gogarten
* github: <a
href='https://github.com/Bioconductor/VCF_projects'>VCF_probabilityBasedSnpEncoding page</a> and <a
href='https://github.com/Bioconductor/VCF_projects/wiki/VCF_probabilityBasedSnpEncoding'>wiki</a>.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

## <a id="msGUI"></a>msGUI - an interactive mass spectrometry data browser

The aim of this project is to build a simple GUI to navigate raw mass
spectrometry data files. Data input functionality and relevant data
structures are available in the mzR and MSnbase packages. The final
deliverable would be a new R package, that will be submitted to
Bioconductor, implementing the GUI allowing users to directly browse
raw data files as well as MSnExp raw data instances. The overall goal
being to complement programmatic data access with interactive
visualisation.

Project attributes and estimates:

* Difficulty: Easy to medium, depending on experience
* Length: 4 - 6 weeks
* Skills needed: intermediate R programming, experience with GUI  programming is a major advantage, knowledge of package development helpful, but not essential.
* Deliverable: msGUI package
* Mentors: Laurent Gatto and Michael Lawrence
* Status: Complete
* Mentee: Andrius Druzinis
* References: mzR and MSnbase packages, Programming Graphical User Interfaces in R book.
* More information: <a href='https://github.com/lgatto/msGUI'>msGUI github page</a> and <a href='https://github.com/lgatto/msGUI/wiki'>wiki</a>.
<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

