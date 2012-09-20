# ![](/images/icons/magnifier.gif)Mentored Projects

* [Introduction](#introduction)
* [Galaxy-ification of Useful Scripts](#galaxy)
* [New Variant Call Format (VCF) Class Methods](#vcf)
* [Create an AnnotationDbi Package for PANTHER](#panther)
* [Finish your own package](#your-package)
* [Interested?](#contact)
* [Completed projects](#completed)


## <a id="introduction">Introduction</a>

Developers new to Bioconductor may find mentored projects a useful way
to apply, refine and extend their skills.  Projects are identified by
experienced Bioconductor developers.  The projects involve important
but manageable programming tasks.  Experienced developers act as
mentors, providing guidance and oversight.  Successful mentored
projects will be incorporated into the appropriate packages, and
contributors will receive full credit for their work.  Users,
contributors and mentors will all benefit.

Below you will find a list of proposed projects.  We invite your participation.  We welcome your suggestions.

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
The new Bioconductor <a href="packages/devel/bioc/html/RGalaxy.html">RGalaxy</a>
 simplifies the process of exposing an R function in Galaxy so that a user
can run the function using nothing more than a web browser.

This project would involve taking an existing workflow (or conceiving a new workflow)
and exposing it in Galaxy.

Project attributes and estimates:

* Difficulty: Medium
* Length: 1-2 weeks
* Skills needed: R programming
* Mentor: Dan Tenenbaum

<h2 id="vcf">New Variant Call Format (VCF) Class Methods</h2>

The <a href=http://www.bioconductor.org/packages/2.10/bioc/html/VariantAnnotation.html>VariantAnnotation</a> package
provides functions to read data from a VCF file and annotate the variants. 
Several opportunities are available to develop new functionality 
and expand existing methods.

Compute genotype allele frequency :
<blockquote>
Write a function to parse the genotype data in 
a VCF class and compute the allele frequency.
</blockquote>

<blockquote>
Project attributes and estimates:
<ul>
<li>Difficulty: Easy</li>
<li>Length: 2 weeks</li>
<li>Skills needed: R programming, familiarity with S4 classes</li>
<li>Deliverables: Implement, test and document alleleFrequency,VCF</li>
<li>Mentor: Valerie Obenchain</li>
</ul>
</blockquote>

Further development of writeVcf() :
<blockquote>
The current writeVcf() function writes a VCF file from data
stored in a VCF-class object. We would like to expand this 
function to write data from more general structures such 
as matrices, DataFrames or lists.
</blockquote>

<blockquote>
Project attributes and estimates:
<ul>
<li>Difficulty: Medium</li>
<li>Length: 6 weeks</li>
<li>Skills needed: R programming, familiarity with S4 classes</li>
<li>Deliverables: Implement, test and document writeVCF,GRanges</li>
<li>Mentor: Valerie Obenchain</li> 
</ul>
</blockquote>

Convert genotypes to probability-based SnpMatrix encoding :
<blockquote>
MatrixToSnpMatrix() in the VariantAnnotation package
converts the genotype data in a VCF object into a SnpMatrix 
object. Currently this is done without taking uncertain
uncertain genotype calls into consideration. This project 
involves modifying MatrixToSnpMatrix() to use, when available, 
genotype uncertainty and likelihood information to convert 
genotypes to probability-based SnpMatrix encodings.
</blockquote>

<blockquote>
Project attributes and estimates:
<ul>
<li>Difficulty: Medium</li>
<li>Length: 8 weeks</li>
<li>Skills needed: R programming, familiarity with S4 classes, statistics</li>
<li>Mentor: Valerie Obenchain</li>
</ul>
</blockquote>


<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>


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

## <a id="your-package"></a>Finish your own package

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

<h2 id="contact">Interested?</h2>

Please send mail to pshannon AT fhcrc DOT org if you would like to
help out on any of these projects, or have an idea of your own which
you wish to propose.


<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

## <a id="completed"></a>Completed projects

* [Add constructors to the 'graph' package](#graph)

### <a id="graph"></a>Add 'contructors' for key classes in the graph package

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

[graph package]: http://bioconductor.org/packages/devel/bioc/html/graph.html

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>
