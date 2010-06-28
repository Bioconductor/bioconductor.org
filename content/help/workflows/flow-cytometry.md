![](/images/icons/help.gif)Flow Cytometry Packages
==================================================

Open source Bioconductor packages for analysis of flow cytometry data
provides a unified framework for bioinformaticians to develop methods to
analyze and interpret flow cytometry data. A workflow for analysis of flow
cytometry data might involve the following packages.

* [flowCore](http://bioconductor.org/packages/release/bioc/html/flowCore.html) -
  handles importing, storing, preprocessing and assessment of data from flow
  cytometry experiments.
* [flowViz](http://bioconductor.org/packages/release/bioc/html/flowViz.html) -
  provides graphical methods for visualization of flow cytometry data.
* [flowQ](http://bioconductor.org/packages/release/bioc/html/flowQ.html) -
  provides quality control and quality assessment tools for flow cytometry data.
* [flowStats](http://bioconductor.org/packages/release/bioc/html/flowStats.html) -
  provides tools and methods to analyze flow cytometry data that are beyond
  the basic infrastructure provided in the flowCore package.
* [flowUtils](http://bioconductor.org/packages/release/bioc/html/flowUtils.html) -
  provides utilities, mainly to integrate foreign flow cytometry analysis
  tools.
* [flowClust](http://bioconductor.org/packages/release/bioc/html/flowClust.html) -
  provides tools for robust model based clustering using t mixture models with
  Box-Cox transformation.
* [flowMerge](http://bioconductor.org/packages/release/bioc/html/flowMerge.html) -
  provides tools for merging of mixture components for model based automated
  gating of flow cytometry data using the flowClust framework.
* [flowFP](http://bioconductor.org/packages/release/bioc/html/flowFP.html) -
  provides tools for fingerprint generation of flow cytometry data to facilitate
  the application of machine learning and data mining tools for flow cytometry.
* [iFlow](http://www.bioconductor.org/packages/release/bioc/html/iFlow.html) -
  provides GUI and tools to explore, analyze, and visualize flow cytometry data.
  The tutorial video can be downloaded from
  [here](../flowcytometry/tutorial.mpeg).

These packages can be installed using [biocLite("package name")](/install/)
Each package also contains a detailed vignette illustrating typical use cases
as well as man pages. The vignette for a package can be viewed using
`openVignette("package name")`.

A typical workflow using the packages flowCore, flowViz, flowQ and flowStats
is described in detail in [flowWorkFlow.pdf](../flowcytometry/flowWorkFlow.pdf).
The data files used in the workflow can be downloaded from
[here](../flowcytometry/dataFiles.tar).
