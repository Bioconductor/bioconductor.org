![](/images/icons/magnifier.gif)Version Numbering
=================

Version Format
-----------------
All _Bioconductor_ packages should have a version number in x.y.z format.

Examples of good version numbers:

    1.2.3
    0.99.5
    2.3.0
    3.12.44

Even Odd Schedule
-----------------

_Bioconductor_ has a 'devel' branch where new features are introduced,
and release branches created twice per year.  Given a package with
version number x.y.z,

* y should be an odd number in devel.
* y should be an even number in release.

During regular development of new features

* Authors increment the z version of their package in the 'devel'
  branch by 1 for each GIT commit, e.g., from 1.1.0 to 1.1.1, 1.1.2,
  ..., 1.1.10, 1.1.11
* Any bug fixes ported to the release branch are similarly
  incremented, 1.0.1, 1.0.2, ...
* Changes made to GIT without a corresponding version bump **do not
  propagate** to the repository visible to `BiocManager::install()`.

At the time of a release, the _Bioconductor_ team arranges to:

* Create a release branch package with version: x*.*(y+1)*.*(z-z)
* Increment the devel branch package version to: x*.*(y+2)*.*(z-z)
* As a special case, any package with version x.99.z is bumped to
  (x+1)*.*0*.*(z-z) in release, and (x+1)*.*1*.*0 in devel. Thus authors
  wishing to signify a major change to their package should set y to 99 in
  their devel package.

New packages
-----------------

New packages submitted to _Bioconductor_ should set Version: 0.99.0 in the
DESCRIPTION file. Specifying y=99 triggers a bump in x at the next release
which in this case results in version 1.0.0.

See additional details on 
the [Package Submission](/developers/package-submission) page.

See also the instructions for [Using Bioc Devel][].

[Using Bioc Devel]: /developers/how-to/useDevel/

Examples
-----------------

1. Normal Case. Suppose a package in the devel branch has version
   number 1.1.25. The new release branch now contains a copy of the
   package with version 1.2.0.  The devel branch of _Bioconductor_
   contains the package whose version number has been bumped to
   1.3.0

2. Special Case.  The "0.99.2" version of our package is copied by the
   _Bioconductor_ team to the release branch with version number
   1.0.0. The package version is bumped to 1.1.0 in the devel
   branch.

Examples of the version bumping scheme:

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
<br />

Summary
-----------------

Below is a summary of how version components are bumped and the
key limitations. Bumps "at release time" are done by the _Bioconductor_ 
team and not the package maintainer.

`x` 
    - only modified by the _Bioconductor_ team
    - bumped to x+1 at release time if y=99
`y` 
    - must be even in release and odd in devel
    - must be <=99
    - bumped at release time for all packages to the next
      even number in release and the next odd in devel
`z`
    - should be incremented sequentially during regular package development
    - no limitation on the size of z
    - bumped at release time to 0 for all packages.
