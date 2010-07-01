Version Numbering
=================

All Bioconductor packages should have a version number like x.y.z.
Here are some examples of **bad** version numbers:
    
    1.2-3
    3.4
    1-2.4

And here are some examples of **good** version numbers:

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
* As part of the release packages that will become the new devel line
  will have their y incremented yet again to get back to an odd
  number.
* If a package has not changed in devel, the version number will not
  be bumped in the release.


Example Changes
---------------

Here is a table showing examples of the version bumping scheme:

<table border="1">
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
</tbody>
</table>

Historical Note
---------------

The even/odd scheme was introduced as part of the BioC 1.7 release.
For this release, all packages going into the release had their y
digit incremented to the next largest even number.  Then in devel,
*all* packages had y incremented again so that all packages begin the
1.8 release cycle with an odd y digit and z equal to zero.
