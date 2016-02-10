# ![](/images/icons/magnifier.gif)biocViews

## Introduction

Packages added to the _Bioconductor_ Project require a `biocViews:`
field in their DESCRIPTION file.  The field name "biocViews" is
case-sensitive and must begin with a lower-case 'b'.

biocViews terms are "keywords" used to describe a given package. They
are broadly divided into three categories, representing the type of
packages present in the _Bioconductor_ Project -

1. Software

2. Annotation Data

3. Experiment Data

biocViews are available for the [release][] and [devel][] branches of
_Bioconductor_. The devel branch has a check box under the tree
structure which, when checked, displays biocViews that are defined but
not used by any package, in addition to biocViews that are in use.

## Motivation

One can use biocViews for two broad purposes.

1. A researcher might want to identify all packages in the
   _Bioconductor_ Project which are related to a specific purpose.
   For example, one may want to look for all packages related to "Copy
   Number Variants".

2. During development, a package contributor can "tag" their package
   with biocViews so that when someone looking for packages (like in
   scenario 1) can easily find their package.

## biocViews during new package development

Visit the [devel][] branch of biocViews when you are in the process of
adding biocViews to your new package. Identify as many terms as
appropriate from the hierarchy. Prefer 'leaf' terms at the end of the
hierarchy, over more inclusive terms. Remember to check the box
displaying all available terms.

Please Note:

1. Your package will belong to only one part of _Bioconductor_ Project
   (Software, Annotation Data or Experiment Data), so choose only only
   biocViews from that category.

2. biocViews listed in your package must match exactly (e.g.,
   spelling, capitalization) the terms in the biocViews hierarchy.

When you submit your new package for review , your package is checked
and built by the _Bioconductor_ Project.  We check the following for
biocViews:

1. Package contributor has added biocViews.

2. biocViews are valid.

3. Package contributor has added biocViews from only one of the categories.

If you receive a "RECOMMENDED" direction for any of these biocViews
after you have submitted your package, you can try correcting them on
your own following the directions given here or ask your package
reviewer for more information.

[release]: http://bioconductor.org/packages/release/BiocViews.html
[devel]: http://bioconductor.org/packages/devel/BiocViews.html
