Flow Cytometry Overview
=======================

Key packages
------------


Open source Bioconductor packages for analysis of flow cytometry data
provides a unified framework for bioinformaticians to develop methods
to analyze and interpret flow cytometry data. A workflow for analysis
of flow cytometry data might involve the following packages.

* flowCore - handles importing, storing, preprocessing and assessment
  of data from flow cytometry experiments
  
* flowViz - provides graphical methods for visualization of flow
  cytometry data.
  
* flowQ - provides quality control and quality assessment tools for
  flow cytometry data
  
* flowStats - provides tools and methods to analyze flow cytometry
  data that are beyond the basic infrastructure provided in the
  flowCore package.
  
* flowUtils - provides utilities, mainly to integrate foreign flow
  cytometry analysis tools.
  
* flowClust - provides tools for robust model based clustering using t
  mixture models with Box-Cox transformation.
  
* flowMerge - provides tools for merging of mixture components for
  model based automated gating of flow cytometry data using the
  flowClust framework.
  
* flowFP - provides tools for fingerprint generation of flow cytometry
  data to facilitate the application of machine learning and data
  mining tools for flow cytometry.

These packages can be installed using biocLite("package name"). Each
package also contains a detailed vignette illustrating typical use
cases as well as man pages. The vignette for a package can be viewed
using openVignette("package name")

A typical workflow using the packages flowCore, flowViz, flowQ and
flowStats is described in detail in flowWorkFlow.pdf. The data files
used in the workflow can be downloaded from here.
