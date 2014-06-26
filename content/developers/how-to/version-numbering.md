![](/images/icons/magnifier.gif)Version Numbering
=================

All Bioconductor packages should have a version number like x.y.z.
Here are examples of good version numbers:

    1.2.3
    0.0.1
    0.99.5
    2.3.0
    3.12.44
    

First time developers (submitting a package to Bioconductor) should
set Version: 0.99.0 in the DESCRIPTION file; see additional details on
the [Package Submission](/developers/package-submission) page.

Even Odd Schedule
-----------------

Bioconductor has a 'devel' branch where new features are introduced,
and release branches created twice per year.  Given a package with
version number x.y.z,

* y should be an odd number in devel.
* y should be an even number in release.

During regular development of new features

* Authors increment the z version of their package in the 'devel'
  branch by 1 for each svn commit, e.g., from 1.1.0 to 1.1.1, 1.1.2,
  ..., 1.1.10, 1.1.11
* Any bug fixes ported to the release branch are similarly
  incremented, 1.0.1, 1.0.2, ...
* Changes made to SVN without a corresponding version bump **do not
  propagate** to the repository visible to `biocLite()`.

At the time of a release, the Bioconductor team arranges to:

* Create a release branch package with version x.y+1.0
* Increment the devel branch package version to x.y+2.0
* As a special case, any package with version x.99.z is bumped to
  x+1.0.0 in release, and x+1.1.0 in devel. Thus authors wishing to
  signify a major change to their package should set y to 99 in their
  devel package.

Example

1. Normal Case. Suppose a package in the devel branch has version
   number 1.1.25. The new release branch now contains a copy of the
   package with version "1.2.0".  The devel branch of Bioconductor
   contains the package whose version number has been bumped to
   "1.3.0"

2. Special Case.  The "0.99.2" version of our package is copied by the
   Bioconductor team to the release branch with version number
   "1.0.0". The package version is bumped to "1.1.0" in the devel
   branch.

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



