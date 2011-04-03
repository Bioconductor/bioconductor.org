# ![](/images/icons/help.gif)Using Bioconductor for High Throughput Assays

Bioconductor includes packages for analysis of diverse areas of
high-throughput assays such as flow cytometry, quantitative real-time PCR,
mass spectrometry, proteomics and other cell-based data. 

* [Sample Workflow](#sample-workflow)  
* [Installation and Use](#install-and-use)
* [Exploring Package Content](#exploring-package-content)
* [Diverse Assays Resources](#diverse-assays-resources)

<h2 id="sample-workflow"> Sample Workflow</h2>

The following psuedo-code illustrates a typical R / Bioconductor
session. It makes use of the flow cytometry packages to load, transform and
visualize the flow data and gate certain populations in the dataset. 

The workflow loads the `flowCore`, `flowStats` and `flowViz` packages and its
dependencies.  It loads the ITN data with 15 samples, each of which includes,
in addition to FSC and SSC, 5 fluorescence channels: CD3, CD4, CD8, CD69 and
HLADR. 

    ## Load packages
    > library(flowCore)
    > library(flowStats)
    > library(flowViz) # for flow data visualization

    ## Load data
    > data(ITN)
    > ITN
	A flowSet with 15 experiments.

	An object of class "AnnotatedDataFrame"
  	rowNames: sample01, sample02, ..., sample15  (15 total)
  	varLabels and varMetadata description:
    	GroupID:  Group Identifier
    	SiteCode:  Site of examination
    	...: ...
    	name: sample name
    	(7 total)

  	column names:
  	FSC SSC CD8 CD69 CD4 CD3 HLADr Time

First, we need to transform all the fluorescence channels. Using a `workFlow`
object can help to keep track of our progress.

    ## Create a workflow instance and transform data using asinh
    > wf <- workFlow(ITN)
    > asinh <- arcsinhTransform()
    > tl <- transformList(colnames(ITN)[3:7], asinh, 
                          transformationId = "asinh")
    > add(wf, tl)

Next we use the `lymphGate` function to find the T-cells in the CD3/SSC
projection.
    
    ## Identify T-cells population
    > lg <- lymphGate(Data(wf[["asinh"]]), channels=c("SSC", "CD3"),
               preselection="CD4", filterId="TCells", eval=FALSE,
               scale=2.5)
    > add(wf, lg$n2gate, parent="asinh")
    > print(xyplot(SSC ~ CD3| PatientID, wf[["TCells+"]],
                   par.settings=list(gate=list(col="red", 
                   fill="red", alpha=0.3))))
    
![lymphGate](lymphGate.png)
    
A typical workflow for flow cytometry data analysis in Bioconductor flow
packages include data transformation, normalization, filtering, manual gating,
semi-automatic gating and automatic clustering if desired. Details can be
found in [flowWorkFlow.pdf](flowWorkFlow.pdf) or the vignettes of the
[flow cytometry packages](#diverse-assays-resources).

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="install-and-use">Installation and Use</h2>
Follow [installation instructions](/install/) to start using these
packages.  To install the `flowCore` package and all of its
dependencies, evaluate the commands

    > source("http://bioconductor.org/biocLite.R")
    > biocLite("flowCore")

Package installation is required only once per R installation. View a
full list of
[available packages](/packages/release/bioc/).

To use the `flowCore` package, evaluate the command

    > library("flowCore")

This instruction is required once in each R session.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="exploring-package-content">Exploring Package Content</h2>

Packages have extensive help pages, and include vignettes highlighting
common use cases. The help pages and vignettes are available from
within R. After loading a package, use syntax like

    > help(package="flowCore")
    > ?read.FCS

to obtain an overview of help on the `flowCore` package, and the
`read.FCS` function, and

    > browseVignettes(package="flowCore")

to view vignettes (providing a more comprehensive introduction to
package functionality) in the `flowCore` package. Use

    > help.start()

to open a web page containing comprehensive help resources.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<h2 id="diverse-assays-resources">Diverse Assays Resources</h2>

The following provide a brief overview of packages useful for analysis
of high-throughput assays. More comprehensive workflows can be found
in documentation (available from [package
descriptions](/packages/release/bioc/))
and in Bioconductor [publications](/help/publications/).

### Flow Cytometry ###

These packages use standard FCS files, including infrastructure,
utilities, visualization and semi-autogating methods for the
analysis of flow cytometry data.

[flowCore](/packages/release/bioc/html/flowCore.html),
[flowViz](/packages/release/bioc/html/flowViz.html),
[flowQ](/packages/release/bioc/html/flowQ.html),
[flowStats](/packages/release/bioc/html/flowStats.html),
[flowUtils](/packages/release/bioc/html/flowUtils.html),
[flowFP](/packages/release/bioc/html/flowFP.html),
[flowTrans](/packages/release/bioc/html/flowTrans.html),
[iFlow](/packages/release/bioc/html/iFlow.html)

Algorithms for clustering flow cytometry data are found in these packages:

[flowClust](/packages/release/bioc/html/flowClust.html),
[flowMeans](/packages/release/bioc/html/flowMeans.html),
[flowMerge](/packages/release/bioc/html/flowMerge.html),
[SamSPECTRAL](/packages/release/bioc/html/SamSPECTRAL.html)

A typical workflow using the packages `flowCore`, `flowViz`, `flowQ` and
`flowStats` is described in detail in [flowWorkFlow.pdf](flowWorkFlow.pdf).
The data files used in the workflow can be downloaded from
[here](dataFiles.tar).

### Cell-based Assays ###

These packages provide data structures and algorithms for cell-based
high-throughput screens (HTS).

[cellHTS2](/packages/release/bioc/html/cellHTS2.html),
[RNAither](/packages/release/bioc/html/RNAither.html)

This package supports the xCELLigence system which contains a series of
real-time cell analyzer (RTCA).

[RTCA](/packages/release/bioc/html/RTCA.html)

### High-throughput qPCR Assays ###

These package provide algorithm for the analysis of cycle threshold
(Ct) from quantitative real-time PCR data.

[HTqPCR](/packages/release/bioc/html/HTqPCR.html),
[ddCt](/packages/release/bioc/html/ddCt.html),
[qpcrNorm](/packages/release/bioc/html/qpcrNorm.html)

### Mass Spectrometry and Proteomics data ###

These packages provide framework for processing, visualization, and
statistical analysis of mass spectral and proteomics data.

[clippda](/packages/release/bioc/html/clippda.html),
[MassArray](/packages/release/bioc/html/MassArray.html),
[MassSpecWavelet](/packages/release/bioc/html/MassSpecWavelet.html),
[PROcess](/packages/release/bioc/html/PROcess.html),
[flagme](/packages/release/bioc/html/flagme.html),
[xcms](/packages/release/bioc/html/xcms.html)

### Imaging Based Assays ###

These packages provide infrastructure for image-based phenotyping and automation of other image-related tasks:

[EBImage](/packages/release/bioc/html/EBImage.html)

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>
