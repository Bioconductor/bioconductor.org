##### Provisional agenda

<table>
  <tr><td colspan="2">Monday</td><tr>
  <tr><td>9-11</td><td>Introduction to <em>R</em> and
      <em>Bioconductor</em></td></tr>
  <tr><td>11-12, 1-2</td><td>Working with sequences, ranges and
      annotations</td></tr>
  <tr><td>2-4</td><td>Exploring sequences</td></tr>
  <tr><td colspan="2">Tuesday</td><tr>
  <tr><td>9-11:15</td><td>RNA-seq differential
      representation</td></tr>
  <tr><td>11:15-12, 1-2:30</td><td>Approaches to ChIP-seq</td></tr>
  <tr><td>2:30-4</td><td>Gene-centric annotation</td></tr>
</table>


## [Course Vignette (PDF)](Bioconductor-tutorial.pdf)
## [Annotation Slides (PDF)](AnnotationSlides.pdf)

##### Installing the course package on your computer

We will use R version 2.14 alpha. To install this version of R, visit (<a
href="http://cran.fhcrc.org/bin/windows/base/rtest.html">Windows</a>, <a
href="http://r.research.att.com/">Mac OS</a>, install from <a
href="http://cran.fhcrc.org/src/base-prerelease/">source</a> 
or <a href="http://cran.fhcrc.org/bin/linux/">binary</a> for Linux).  
The course uses additional R and Bioconductor packages. To
install these packages, start R and copy and paste the following
command into your R session:

	source("http://bioconductor.org/scratch-repos/pkgInstall.R")
	pkgInstall(c("SeattleIntro2011", "SeattleIntro2011Data"))


###### Installing R from Source on Ubuntu Linux

    sudo apt-get update
    sudo apt-get install -y make libc6-dev gfortran gfortran-4.3 build-essential \
      libreadline5-dev libx11-dev libxt-dev libcurl4-openssl-dev libxml2-dev \
      texlive-full tcl8.5-dev tk8.5-dev libxss-dev libpng12-dev libjpeg62-dev \
      libcairo2-dev gcj
    curl -O http://cran.fhcrc.org/src/base-prerelease/R-latest.tar.gz
    tar zxf R-latest.tar.gz
    cd R-alpha # or R-beta if this does not work
    ./configure
    make
    sudo make install
