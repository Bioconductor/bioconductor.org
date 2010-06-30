![](/images/icons/magnifier.gif)About Bioconductor
==================================================

Bioconductor is an open source and open development software project
to provide tools for the analysis and comprehension of genomic data.

Bioconductor is based primarily on the [R](http://www.r-project.org) programming language, but
does contain contributions in other programming languages. It has two
releases each year that follow the semiannual releases of R. At any
one time there is a [release](http://www.bioconductor.org/packages/release/bioc/) version, which corresponds to the released
version of R, and a [development](http://www.bioconductor.org/packages/devel/bioc) version, which corresponds to the
development version of R. Most users will find the release version
appropriate for their needs. In addition there are a large number of
[meta-data packages](http://www.bioconductor.org/packages/release/data/annotation/) available that are mainly, but not solely, oriented
towards different types of microarrays.

The Bioconductor project was started in the Fall of 2001 and is
overseen by the [Bioconductor core team](/about/core-team/), based
primarily at the [Fred Hutchinson Cancer Research
Center](http://fhcrc.org) with other members coming from various US
and international institutions. It gained widespread exposure in the
groundbreaking Genome Biology 2004 paper [Bioconductor: open software
development for computational biology and bioinformatics](http://genomebiology.com/content/pdf/gb-2004-5-10-r80.pdf). More project
details are available online in the [Bioconductor annual reports](http://merlot2.fhcrc.org:3000/about/annual-reports/).

Bioconductor Packages
---------------------

Most Bioconductor components are distributed as R [packages](http://cran.r-project.org/doc/FAQ/R-FAQ.html#R-Add_002dOn-Packages), which are
add-on modules for R. Initially most of the [Bioconductor software
packages](/packages/release/bioc/) focused primarily on DNA microarray data analysis. As the
project has matured, the functional scope of the software packages
broadened to include the analysis of all types of genomic data, such
as SAGE, sequence, or SNP data.


Goals of the Bioconductor Project
---------------------------------

The broad goals of the Bioconductor project are:

* To provide widespread access to a broad range of powerful statistical and graphical methods for the analysis of genomic data.
* To facilitate the inclusion of biological metadata in the analysis of genomic data, e.g. literature data from PubMed, annotation data from Entrez genes.
* To provide a common software platform that enables the rapid development and deployment of extensible, scalable, and interoperable software.
* To further scientific understanding by producing high-quality [documentation](/help/package-vignettes/) and reproducible research.
* To [train](/help/course-materials/) researchers on computational and statistical methods for the analysis of genomic data.


Main Features of the Bioconductor Project
-----------------------------------------

*  The R Project for Statistical Computing. [R](http://www.r-project.org) and the [R package system](http://cran.r-project.org/doc/FAQ/R-FAQ.html#R-Add_002dOn-Packages) provide a broad range of advantages to the Bioconductor project including:
<ul>
 <li> It contains a high-level interpreted language in which one can easily and quickly prototype new computational methods.</li>
 <li> It includes a well established system for packaging together software components and documentation.</li>
 <li> It can address the diversity and complexity of computational biology and bioinformatics problems in a common object-oriented framework.</li>
 <li> It provides on-line computational biology and bioinformatics data sources.</li>
 <li> It supports a rich set of statistical simulation and modeling activities.</li>
 <li> It contains cutting edge data and model visualization capabilities.</li>
 <li> It has been the basis for pathbreaking research in parallel statistical computing.</li>
 <li> It is under very active development by a dedicated team of researchers with a strong commitment to good documentation and software design.</li>
</ul>
* Documentation and reproducible research. Each [Bioconductor package](/packages/release/bioc/) contains at least one [vignette](/help/package-vignettes/), which is a document that provides a textual, task-oriented description of the package's functionality. These vignettes come in several forms. Many are simple "HowTo"s that are designed to demonstrate how a particular task can be accomplished with that package's software. Others provide a more thorough overview of the package or might even discuss general issues related to the package. In the future, we are looking towards providing vignettes that are not specifically tied to a package, but rather are demonstrating more complex concepts. As with all aspects of the Bioconductor project, users are encouraged to participate in this effort.
* Statistical and graphical methods. The Bioconductor project aims to provide access to a wide range of powerful statistical and graphical methods for the analysis of genomic data. Analysis [packages](/packages/release/bioc/) are available for: pre-processing Affymetrix and cDNA array data; identifying differentially expressed genes; graph theoretical analyses; plotting genomic data. In addition, the [R package system](http://cran.r-project.org/doc/FAQ/R-FAQ.html#R-Add_002dOn-Packages) itself provides implementations for a broad range of state-of-the-art statistical and graphical techniques, including linear and non-linear modeling, cluster analysis, prediction, resampling, survival analysis, and time-series analysis.
* Annotation. The Bioconductor project provides software for associating microarray and other genomic data in real time to biological metadata from web databases such as GenBank, Entrez genes and PubMed ([annotate](http://www.bioconductor.org/packages/release/bioc/html/annotate.html) package). Functions are also provided for incorporating the results of statistical analysis in HTML reports with links to annotation WWW resources.  Software tools are available for assembling and processing genomic annotation data, from databases such as GenBank, the Gene Ontology Consortium, Entrez genes, UniGene, the UCSC Human Genome Project ([AnnotationDbi](http://www.bioconductor.org/packages/release/bioc/html/AnnotationDbi.html) package).  [Data packages](http://www.bioconductor.org/packages/release/bioc/html/AnnotationDbi.html) are distributed to provide mappings between different probe identifiers (e.g. Affy IDs, Entrez genes, PubMed). Customized annotation libraries can also be assembled.
* Bioconductor short courses. The Bioconductor project has developed a program of [short courses](#) on software and statistical methods for the analysis of genomic data. Courses have been given for audiences with backgrounds in either biology or statistics. All [course materials](/help/course-materials/) (lectures and computer labs) are available on this site. [Customized short courses](#) may also be designed for interested parties.
* Open source. The Bioconductor project has a commitment to full open source discipline, with distribution via a SourceForge-like platform. All contributions are expected to exist under an open source license such as Artistic 2.0, GPL2, or BSD. There are many different reasons why open--source software is beneficial to the analysis of microarray data and to computational biology in general. The reasons include:
          o To provide full access to algorithms and their implementation
          o To facilitate software improvements through bug fixing and software extension
          o To encourage good scientific computing and statistical practice by providing appropriate tools and instruction
          o To provide a workbench of tools that allow researchers to explore and expand the methods used to analyze biological data
          o To ensure that the international scientific community is the owner of the software tools needed to carry out research
          o To lead and encourage commercial support and development of those tools that are successful
          o To promote reproducible research by providing open and accessible tools with which to carry out that research (reproducible research is distinct from independent verification)

* Open development. Users are encouraged to become developers, either by contributing [Bioconductor compliant packages](http://wiki.fhcrc.org/bioc/Package_Guidelines) or documentation. Additionally Bioconductor provides a mechanism for linking together different groups with common goals to foster collaboration on software, often at the level of shared development.
* New Users. To see practical application of Bioconductor software, you might consider buying [Bioconductor Case Studies](/help/books/bioconductor-case-studies/).
