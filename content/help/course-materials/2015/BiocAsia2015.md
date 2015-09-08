
## Overview

The _Bioconductor_ project is a widely used collection of over 1000 R
packages for high-throughput genomic analysis and recently surpassed
the milestone of 100,000 commits to its svn repository. Many popular
software tools from the project are developed by researchers based in
the Asia-Pacific region. To enhance collaborations and provide an
avenue for networking, the region's first developers' meeting will be
held directly preceding GIW/InCoB 2015. This event will bring together
both current and prospective package developers to provide a forum for
exchanging ideas and future plans for software in the project.  The
meeting will consist of a number of longer talks selected from
abstracts as well as short 'lightning' presentations to maximize the
opportunities for participants to highlight their work. A Developer
workshop on 'Turning great ideas into _Bioconductor_ packages with
impact' will also be held.

## Organisers

Dr Matt Ritchie, Molecular Medicine Division, The Walter and Eliza
Hall Institute of Medical Research, Melbourne, Australia

Dr. Martin Morgan, Bioconductor / Program in Computational Biology
Fred Hutchinson Cancer Research Center, Seattle, WA USA

## Venue

AIST Tokyo Waterfront Bio-IT Research Building 2-4-7 Aomi, Koto-ku,
Tokyo, 135-0064, Japan

## Key dates

Email proposals for longer talks / lightning talks or posters to
gustin.s near wehi.edu.au (Title and abstract of up to 300 words) by
Friday 31st July 2015.

## Registration

Registration is handled as part of the GIW/InCoB 2015
[registration process](https://perdana.apbionet.org/giw-incob-2015/).
Alternatively email gustin.s near wehi.edu.au if you will be
attending.

## Tentative schedule

|:------------------|:--------------------------------------------------------------------------------|
| 10:00am - 10:50am | Project Updates (Dr Martin Morgan)                                              |
| 10:50am - 11:50am | Getting to know you session: Contributed 'lightning talks'                      |
| 11:50am - 1:00pm  | Lunch                                                                           |
| 1:00pm   - 2:40pm | Contributed talks session I (4 x 25 minute talks)                               |
| 2:40pm   - 3:15pm | Coffee break                                                                    |
| 3:15pm   - 4:15pm | Developer workshop "Turning great ideas into Bioconductor packages with impact" |
| 4:15pm   - 5:15pm | Poster session/social hour                                                      |
| 6:00pm   - 8:00pm | Dinner                                                                          |

## Scholarships

A number of scholarships (up to $400 USD each) are available for students and Bioconductor package developers to help with the cost of travel to attend the meeting.

To apply for a scholarship, please submit a brief statement (200 words of less) describing your interest in the meeting. If you are a package developer or maintainer please tell us the package(s) you work on. The due date for applications is August 21st. Please send applications to gustin.s AT wehi DOT edu DOT au and put SCHOLARSHIP in the subject line.

## Materials

* Introduction [HTML](W1-Introduction.html) | [Rmd](W1-Introduction.Rmd) | [R](W1-Introduction.R)
* Data Representations [HTML](W2-Data-Representations.html) | [Rmd](W2-Data-Representations.html) | [R](W2-Data-Representations.R)
* RNASeq [HTML](W3-RNASeq.html) | [Rmd](W3-RNASeq.Rmd) | [R](W3-RNASeq.R)
* Annotation [HTML](W4-Annotation.html) | [Rmd](W4-Annotation.Rmd) | [R](W4-Annotation.R)

<h2 id="ami">To launch an Amazon Machine Image (AMI) for this course:</h2>

* [Create an Amazon Web Services (AWS) Account](https://aws.amazon.com/) if you
  don't already have one.
* Start the instance <%= ami_url("ami-c5472aa0") %>; See the [documentation for this](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/launching-instance.html). Make sure your [security group](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/using-network-security.html) has
port 80 accessible.
* Paste the Public DNS name of your instance into a web browser.
* Log in to RStudio with username *ubuntu* and password *bioc* .
* Be sure and [terminate your instance](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/terminating-instances.html) when you are done using it, in order to avoid excessive charges.
* Inside RStudio, view the available vignettes with
  `help(package="BiocAsiaPacific2015")` and then click on 
  "User guides, package vignettes and other documentation.".
* For more information, see the
  [Bioconductor AMI page](/help/bioconductor-cloud-ami/).





