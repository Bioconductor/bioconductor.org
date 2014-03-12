![](/images/icons/magnifier.gif)biocViews
========================================================

Introduction
-------------

Packages being added to the Bioconductor Project require biocViews in their
DESCRIPTION file. 

(Note that the field name "biocViews" is case-sensitive and must begin
with a lower-case 'b'.)

biocViews are "keywords" which are used to describe a given package. They are 
broadly divided into three categories, representing the type of packages present
in the Bioconductor Project -

1. Software

2. Annotation Data

3. Experiment Data

Motivation
------------

One can use biocViews for two broad purposes.

1. A researcher might want to identify all packages in the Bioconductor Project
 which are related to a specific purpose.  For example, one may want to look for
 all packages related to Copy Number Variants.
 
2. During development, a package contributor can "tag" his package with 
   biocViews so that when someone is looking for packages (like in scenario 1),
   they can easily find his package.

Devel and Release branch of biocViews
-------------------------------------

The biocViews for the release branch of Bioconductor can be used to search for 
packages that belong to the current release of Bioconductor Project and can be
viewed at:[Release](http://bioconductor.org/packages/release/BiocViews.html#___Software)


The devel branch of biocViews contains the latest biocViews and can be used
to search for packages that are under the development in Bioconductor Project.
The biocViews and packages present here will become part of the next release.
The biocViews for the devel branch of Bioconductor can be viewed at :[Devel](http://bioconductor.org/packages/devel/BiocViews.html#___Software)



Using biocViews to search for packages 
--------------------------------------

Researchers are often directed to the release and devel page when they are 
looking to find new packages. The tree structure in each category can be 
expanded to explore packages in both release and devel. 

The devel branch has a small check box under the tree structure which states
"Developers: check this box to toggle the visibility of childless biocViews."

When this checkbox is unselected, one can view only the biocViews which have 
packages assigned to them.
 
biocViews during development
----------------------------

### Add BiocViews to your new package
It is recommended to visit the devel branch of biocViews when you're in the 
process of adding biocViews to your new package.

The devel branch has a small check box under the tree structure which states
"Developers: check this box to toggle the visibility of childless biocViews."
When this checkbox is selected, package contributors can view the full 
comprehensive list of biocViews. This list also contains biocViews that are
not tagged with any package. 

We believe that by constantly revising and adding new biocViews to the tree in
the devel branch, we give package contributors the power to contribute richer 
and diverse packages for the analysis and comprehension of high-throughput 
genomic data. 

Please Note: 

1. Since your package will belong to any one part of Bioconductor Project(Software,
Annotation Data or Experiment Data), you are allowed only to add biocViews from 
any one of these categories. 

2. Please do not modify any of these biocViews when adding them to your package
(eg: spelling, lower-case, uppper-case) as these biocViews are checked during
the build process by the Bioconductor Project and the package reviewer. 


### biocViews during build  
When you submit your new package for review , your package is checked and built
by the Bioconductor Project. 

We check the following for biocViews: 

1. Package contributor has added biocViews. 

2. biocViews are valid.

3. Package contributor has added biocViews from only one of the categories. 


If you recieve a "RECOMMENDED" direction for any of these biocViews
after you have submitted your package, you can try correcting them on
your own following the  directions given here or ask your package
reviewer for more information.




