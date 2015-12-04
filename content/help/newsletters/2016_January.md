<!-- %\VignetteEngine{knitr::knitr} %\VignetteIndexEntry{} -->

{::options parse_block_html="true" /}

# *Bioconductor* Newsletter
{:.no_toc}

posted by [Valerie Obenchain](mailto:vobencha@roswellpark.org), January 2016


## Contents 
{:.no_toc}

* Table of contents will replace this text. 
{:toc}


## Annotations (need better title)

### `Bioconductor` annotation packages 

`Bioconductor` contains many different types of annotation packages. In this
section we highlight the most heavily used packages and their common uses.

* OrgDb

* ChipDb

* TxDb

* BSgenome

* ...

After individual explanations segue into ... several annotations 
contain similar data or can perform like actions 
(OrgDb and biomaRt; ChipDb and biomaRt ...). 

Explain why some packages are build-specific and some are not.

### Common tasks 

* convert an Affy ID to Entrez Gene ID, 
* find the genomic region for a gene, etc. 
Cover some of the inherent difficulties of annotating things, like 1:many 
mappings, and the tradeoffs you have to make.

### Advanced tasks 

Then if we aren't getting too long, we could move from simple
task-based examples to use cases (I have an RNA-Seq experiment and I
want to generate a list of DE genes with annotation for my PI, etc).



## ExperimentHub

## Infrastructure
- 'generics' package
- SummarizedExperiment ongoing development
- ExperimentHub
- ChiA-PET and Hi-C and InteractionSet class (https://github.com/LTLA/InteractionSet)


## New and Noteworthy

New functions added to `R` (3.3) and `Bioconductor` (3.3) this quarter:

*   *SummarizedExperiment::readKalisto()*

    TODO

*   *GenomicFeatures::mapToIds()* and *GenomicFeatures::mapToRanges()*

    TODO


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


Send comments or questions to Valerie at 
[vobencha@roswellpark.org](vobencha@roswellpark.org).
