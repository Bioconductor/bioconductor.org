![](/images/icons/magnifier.gif)Version Numbering
=================

All Bioconductor packages should have a version number like x.y.z.
Here are examples of good version numbers:

    1.2.3
    0.0.1
    0.4.5
    2.3.0

The Meaning of x
----------------
Packages in the devel repository with version numbers matching 0.y.z
are assumed to be too young for release and **will not** be included
in any releases. The exception to this rule is the special version
number 0.99.z which we reserve for packages which will become 1.0.0
when bumped into the release.


Even Odd Schedule
-----------------

Given a package with version number x.y.z,

* In devel, y should be an odd number.
* In release, y should be an even number.  At the beginning of the
  release, z should be zero.
* As part of the release, packages that will become the new devel line
  will have their y incremented yet again to get back to an odd
  number.
* If a package has not changed in devel, the version number will not
  be bumped in the release.

The Bioconductor team will bump y as part of the release process.

Example Changes
---------------

Here is a table showing examples of the version bumping scheme:

<table border="1" cellpadding="5" cellspacing="0">
<thead valign="bottom">
<tr>
  <th class="head">Release</th>
  <th class="head">Devel</th>
  <th class="head">Just before release</th>
  <th class="head">Next release</th>
  <th class="head">Next devel</th>
</tr>
</thead>
<tbody valign="top">
<tr>
  <td>--</td><td>0.99.0</td><td>0.99.3</td><td>1.0.0</td><td>1.1.0</td>
</tr>
<tr>
  <td>1.4.0</td><td>1.5.0</td><td>1.5.4</td><td>1.6.0</td><td>1.7.0</td>
</tr>
<tr>
  <td>1.4.0</td><td>1.5.0</td><td>1.5.1</td><td>1.6.0</td><td>1.7.0</td>
</tr>
<tr>
  <td>1.4.0</td><td>1.5.0</td><td>1.5.0</td><td>1.4.0 (no changes)</td><td>1.5.0</td>
</tr>
<tr>
  <td>1.4.0</td><td>1.5.0</td><td>1.99.3</td><td>2.0.0</td><td>2.1.0</td>
</tr>
<tr>
  <td>1.8.0</td><td>1.9.0</td><td>1.9.1</td><td>1.10.0</td><td>1.11.0</td>
</tr>
</tbody>
</table>


The Importance of Incrementing Version Numbers
-----------------------------------------------

Developers should get into the habit of bumping the `z` portion of 
their version number every time they commit changes to their package.

If you don't bump the version, the other changes you made to your
package **will not propagate** to the Bioconductor web site and
package repository. This means, for example, that if you fix a bug
and commit a bug fix, but don't bump the version, users who download
your package will still experience that bug.

There is one circumstance in which you may *want* to avoid bumping
the version number when committing a change. That is when you want 
to see if a given change will build on all platforms on our nightly
build system, and whether dependencies will be affected by the change. 
When you are satisfied that your change is acceptable, you should 
bump the version.



Historical Note
---------------

The even/odd scheme was introduced as part of the BioC 1.7 release.
For this release, all packages going into the release had their y
digit incremented to the next largest even number.  Then in devel,
*all* packages had y incremented again so that all packages begin the
1.8 release cycle with an odd y digit and z equal to zero.
