![](/images/icons/magnifier.gif)Mentored Projects
==================================================

* [Introduction](#introduction)
* [Galaxy-ification of Useful Scripts](#galaxy)
* [New Variant Call Format (VCF) Class Methods](#vcf)
* [Create an AnnotationDbi Package for PANTHER](#panther)
* [Interested?](#contact)


<h2 id="introduction">Introduction</h2>

The Bioconductor project welcomes code developers who wish to apply, refine and extend their skills, working on
tasks and small projects we have identified as important, but which we have not had the time to tackle ourselves.   We will
offer technical support and advice to those working on any of these projects.   We predict that contributors, users, and established
Bioconductor developers will all benefit.


Below you will find a list of proposed projects.  We welcome other suggestions as well.

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
<li>Length: 1-2 weeks</li>
<li>Skills needed: R programming, familiarity with S4 classes</li>
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
<li>Length: 4-6 weeks</li>
<li>Skills needed: R programming, familiarity with S4 classes</li>
<li>Mentor: Valerie Obenchain</li> 
</ul>
</blockquote>

Convert genotypes to probability-based SnpMatrix encoding :
<blockquote>
The MatrixToSnpMatrix() function in the VariantAnnotation
package converts genotype data in a VCF-class object
into a SnpMatrix object without taking uncertain genotype 
calls into consideration. This project involves modifying
MatrixToSmpMatrix() to use, when available, genotype 
uncertainty and likelihood information to convert 
genotypes to the probability-based SnpMatrix encoding.
</blockquote>

<blockquote>
Project attributes and estimates:
<ul>
<li>Difficulty: Medium</li>
<li>Length: 6-8 weeks</li>
<li>Skills needed: R programming, familiarity with S4 classes, statistics</li>
<li>Mentor: Valerie Obenchain</li>
</ul>
</blockquote>


<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>


<h2 id='panther'> Create an AnnotationDbi Package for PANTHER</h2>

PANTHER is found <a href=http://www.pantherdb.org/>here</a>, and summarized:
<blockquote>

The PANTHER (Protein ANalysis THrough Evolutionary Relationships)
Classification System is a unique resource that classifies genes by
their functions, using published scientific experimental evidence and
evolutionary relationships to predict function even in the absence of
direct experimental evidence.  "classifies genes by their function"
</blockquote>

Project attributes and estimates:

* Difficulty:
* Length
* Skills needed:  
* Mentor:

We would like to see PANTHER annotation contained in a Bioconductor AnnotationDbi
package.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="contact">Interested?</h2>

Please send mail to pshannon AT fhcrc DOT org if you would like to help out on any of these projects, or have an idea of your own which you
wish to propose.


<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>
