![](/images/icons/magnifier.gif)Version Numbering
=================

All Bioconductor packages should have a version number like x.y.z.
Here are examples of good version numbers:

    1.2.3
    0.0.1
    0.4.5
    2.3.0
    3.12.44
    

For First time developers (submitting a Package to Bioconductor)
------------------------------------------------------------
* Package Developers submit their first version of package in the format 0.99.z
* A detailed package review is returned to the developer, where they are 
 expected to increase the z in "0.99.z" for every new source tarball they submit. 
 They are also expected to build and check the package on their local machine
before re-submitting to Bioconductor. 

Example  

* A package was submitted 'DemoPackage_0.99.0.tar.gz'
* During package review it was recommended that biocViews terms in the DESCRIPTION file, a NEWS file and a vignette be added.
* The package developer made all the recommended changes, bumped the version number to 0.99.1, built and checked the package on his local machine.
* Developer uploaded the tar ball 'DemoPackage_0.99.1.tar.gz' on the tracker.


For packages in their first devel cycle
----------------------------------------------
* Packages accepted on the tracker are then added to the devel branch with the 
current version number of the accepted package. 
* Developers should  bump the `z` portion of their version number every time they commit changes to their package.
* If developers don't bump the version, the changes made to their 
package **will not propagate** to the Bioconductor web site and package 
repository. 
* There is one circumstance in which developers may *want* to avoid bumping
the version number when committing a change. That is when developers want
to see if a given change will build on all platforms on our nightly
build system, and whether dependencies will be affected by the change.
When developers are satisfied that their change is acceptable, then they should
bump the version.

Example  

* The accepted package tarball is added to the devel branch by the Bioconductor   group members. It has its own landing page and Package users download the   "0.99.1" version of "DemoPackage" using biocLite("DemoPackage")  
* The package developer wants to add a new function to his package. The following steps are recommended.  

  * check out the package from the svn repository (using the link, username and password) emailed to the Developer.  
  * add the functionality to the package.  
  * increment the version number from 0.99.1 to 0.99.2  
  * build and check the package  
  * commit back the changes to the svn repository.  
  
* The following day, after the build report is generated - the users can now access the "0.99.2" version from the Bioconductor repository/website( via biocLite("DemoPackage")).

Note -If the developer had not bumped the version from 0.99.1 to 0.99.2, then 
the package users would still get the old version 0.99.1 version of the package
via biocLite.

For packages copied from devel to release cycle
-------------------------------------------------
* Bioconductor has two releases per year - every release, the devel packages are 
copied into the new release branch. 
* The Bioconductor team will bump y to the next even value and set z to zero ,as
part of the release process.
* The Bioconductor team will bump y again for the package in devel. 

Example-

1. Special Case Scenario  
The  "0.99.2" version of our package is copied by the Bioconductor
team to the release branch with version number "1.0.0". The package version
is bumped to "1.1.0" in the devel branch.  

2. Normal Case scenario  
Suppose a package in the devel branch has version number 1.1.25. The new 
release branch now contains a copy of the package with version "1.2.0".  The
devel branch of Bioconductor contains the package whose version number has
been bumped to "1.3.0"

Even Odd Schedule
-----------------
Given a package with version number x.y.z,

* In devel, y should be an odd number.
* In release, y should be an even number.  At the beginning of the
  release, z should be zero.
* As part of the release, packages that will become the new devel line
  will have their y incremented yet again to get back to an odd
  number.

Here is a table showing examples of the version bumping scheme:

<table border="1" cellpadding="5" cellspacing="0">
<thead valign="bottom">
<tr>
  <th class="head">Current Release</th>
  <th class="head">Current Devel</th>
  <th class="head">Just before Next release</th>
  <th class="head">Next release</th>
  <th class="head">Next devel</th>
</tr>
</thead>
<tbody valign="top">
<tr>
  <td>--</td><td>0.99.1</td><td>0.99.2</td><td>1.0.0</td><td>1.1.0</td>
</tr>
<tr>
  <td>1.4.0</td><td>1.5.0</td><td>1.5.4</td><td>1.6.0</td><td>1.7.0</td>
</tr>
<tr>
  <td>1.4.0</td><td>1.5.0</td><td>1.5.1</td><td>1.6.0</td><td>1.7.0</td>
</tr>
<tr>
  <td>1.4.0</td><td>1.5.0</td><td>1.99.3</td><td>2.0.0</td><td>2.1.0</td>
</tr>
<tr>
  <td>1.8.0</td><td>1.9.0</td><td>1.9.1</td><td>1.10.0</td><td>1.11.0</td>
</tr>
</tbody>
</table>



