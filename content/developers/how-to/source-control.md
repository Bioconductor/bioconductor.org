# ![](/images/icons/magnifier.gif)Source Control

The Bioconductor project is maintained in a Subversion source control
system. Maintainers can also access their packages through git and
Github using the [Bioconductor Git-SVN bridge](../git-svn/).

## Subversion Resources

* The [Subversion Book][1] is a key reference; start with the [guided
  tour][2].
* Visit the [Subversion home page][3].
* Read-only access to our svn repository is available with

  * user name: ``readonly``
  * password: ``readonly``

[1]: http://svnbook.red-bean.com/nightly/en/index.html
[2]: http://svnbook.red-bean.com/nightly/en/svn.intro.html
[3]: http://subversion.tigris.org/

## Software Packages

To check out (co) all packages in the software repository (~3 GB) use:

    svn co https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks Rpacks-devel
    
This creates a copy of all packages on your local machine.  Specify a
name other than "Rpacks-devel" if you want a top-level directory with
different name.

To check out the code for the Biobase package use:

    svn co https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/Biobase

The check out command uses the directory name in the repository if the
destination name is not specified.

## Basic svn Operations

Primary commands:

* `svn update` will update your check out from the server to get any new
  changes.
* `svn add` foo will add foo to the repository (note that unlike CVS this is
  a recursive add. Use the -N switch if you don't want this behavior).
* `svn delete` foo will delete foo. If foo is a file it is removed from your
  local copy as well. If it is a directory it is not but is scheduled for
  deletion.
* `svn copy` foo bar will make a copy of foo named bar and copy the history.
* `svn move` foo bar is the same as copy except foo gets deleted.
* `svn commit` commits your changes. Much like CVS you can choose to specify
  a file (or files) or leave it blank and it will commit everything.

Some other commands:

* `svn status` foo will show you information about the file, particularly
  changes that you've made.
* `svn diff` foo will show you the exact diff of your changes to the server
* `svn revert` foo will bring you back to the server copy.
* `svn log` foo will show the log history for that file.

Many of these commands have additional arguments. Get help on diff,
for example, like this:

    svn help diff

The [Subversion Book][1] has more complete documentation and examples
for all the commands and options.

## Where to Commit Changes

Almost all development is on the trunk ('devel') branch of the SVN
repository, as indicated by the 'trunk' part of the URL in the
examples above. All bug fixes, new features and major changes are
introduced in devel. Each commit should include a bump in the `z`
portion of the `x.y.z` package
[versioning scheme](/developers/how-to/version-numbering/).

<!-- UPDATE THIS PARAGRAPH WITH EACH RELEASE (make sure times are correct): -->

If you commit to trunk before 5:20 PM Seattle time, your changes will
build overnight and be reflected in the next day's  [build
report](http://bioconductor.org/checkResults/devel/bioc-LATEST/) which
should appear around 10:00 AM Seattle time.


### Committing Changes to the Release Branch

Only *bug fixes* should be back-ported to the release branch. This is
so that users of the release branch have a stable environment in which
to get their work done.

If you wish to have a bug-fix made in the devel branch also available
in the current release branch, you first need to take note of the
revision number from your commit, for example,

    $ svn commit -m "Sample commit"
    Adding         Rpacks\Biobase\DESCRIPTION
    Sending        Rpacks\Biobase\DESCRIPTION
    Transmitting file data ..
    Committed revision 140.

was revision 140. The changeset you want is `-c140`.

You'll need to have checked out the branches subdirectory, which is
separate from trunk. If you have only checked out the madman
subdirectory previously, you'll need to also check out the appropriate
branch subdirectory:

    svn co https://hedgehog.fhcrc.org/bioconductor/branches/RELEASE_<%=release_branch%>/madman/Rpacks/MYPACKAGENAME

Merge your changes from the trunk to the release branch, check and fix
any conflicts, and commit. So, from your release branch directory
(e.g. RELEASE_<%=release_branch%>/madman/Rpacks/MYPACKAGENAME):

    svn merge -c140 https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/MYPACKAGENAME
    svn status   # Look for C, indicating a conflict
                 # fix conflicts... (remember to use svn resolve for each)
    svn commit -m "merged r140 from trunk"

## Experiment Data Packages

The root of the Bioconductor experiment data svn repository is
[https://hedgehog.fhcrc.org/bioc-data/trunk/experiment](https://hedgehog.fhcrc.org/bioc-data/trunk/experiment).
Experiment data packages are divided into two components: **pkgs**
(containing all but the data files) and **data_store** (containing the
data files). To obtain a specific experiment data package first check
out the package infrastructure then fold the data in using the Python
script
[https://hedgehog.fhcrc.org/bioc-data/trunk/experiment/pkgs/add_data.py](https://hedgehog.fhcrc.org/bioc-data/trunk/experiment/pkgs/add_data.py)
on the exported package. Here are the commands to check out the
affydata package:

    svn export https://hedgehog.fhcrc.org/bioc-data/trunk/experiment/pkgs/add_data.py
    svn co https://hedgehog.fhcrc.org/bioc-data/trunk/experiment/pkgs/affydata
    ./add_data.py affydata

Note that experiment data packages are only built twice a week,
on Wednesdays and Saturdays.

## Having Problems?

Here is a list of possible issues:

* Unrecognized URL scheme:
  If you see "unrecognized URL scheme" when trying to access the
  repository, it may indicate that your svn client does not support
  HTTPS.  You can verify the supported "modes" by examining the
  output of ``svn --version``.  If you do not see support for HTTPS,
  you will need to upgrade your client.

* Username or password not recognized:
  Most usernames we issue are in the form of an email address.  It may
  help to quote the user name: ``svn co
  --username="bob@internet.net"``.  Try double and single quotes.  

* Need help arranging access to your package? Use [contact
  information](/developers/package-submission/#contact-info) for new
  package submission.

* Still unresolved? Ask on the [bioc-devel](/help/mailing-list/)
  mailing list.
