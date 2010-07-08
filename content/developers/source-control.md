![](/images/icons/magnifier.gif)Source Control
===============

The Bioconductor project is maintained in a Subversion source control
system.

Resources for learning about Subversion
---------------------------------------

* The [Subversion Book][1] is *the* reference for using svn. The
  [guided tour][2] is a recommended starting point.
* [Subversion homepage][3]
* Read-only access to our svn repository is available with

  * user name: ``readonly``
  * password: ``readonly``

[1]: http://svnbook.red-bean.com/nightly/en/index.html
[2]: http://svnbook.red-bean.com/nightly/en/svn.intro.html
[3]: http://subversion.tigris.org/

Obtaining a working copy of all BioC software packages
------------------------------------------------------

The root of the Bioconductor svn repository is
[https://hedgehog.fhcrc.org/bioconductor](https://hedgehog.fhcrc.org/bioconductor). 
To checkout (co) all packages in the repository (~3 GB) use:

    svn co https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks Rpacks-devel
    
This will create a copy of all packages on your local machine.
Note that you can specify any destination name if you want a top-level
directory with a name different from "Rpacks-devel".

Here's the command to checkout the code for only the Biobase package:

    svn co https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/Biobase

Note that the checkout command will default to the directory name in
the repository if the destination name is not specified.

Obtaining a working copy of a BioC experimental data package
------------------------------------------------------------

The root of the Bioconductor svn repository is
[https://hedgehog.fhcrc.org/bioc-data/trunk/experiment](https://hedgehog.fhcrc.org/bioc-data/trunk/experiment).
Data experiment packages are typically divided into two components: **pkgs**
(containing all but the data files) and **data_store** (containing the data
files). To obtain a specific experimental data package first check out the
package infrastructure then fold the data in using the Python script
[https://hedgehog.fhcrc.org/bioc-data/trunk/experiment/pkgs/add_data.py](https://hedgehog.fhcrc.org/bioc-data/trunk/experiment/pkgs/add_data.py)
on the exported package. Here's the set of commands to checkout the
affydata package:

    svn export https://hedgehog.fhcrc.org/bioc-data/trunk/experiment/pkgs/add_data.py
    svn co https://hedgehog.fhcrc.org/bioc-data/trunk/experiment/pkgs/affydata
    ./add_data.py affydata

Basic svn operations
--------------------

Primary commands:

* `svn update` will update your checkout from the server to get any new
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

Many of these commands have extra possible arguments. You can get help on
diff, for example, like this:

    svn help diff

The [Subversion Book][1] will have more complete documentation and
examples for all the commands and options.

Committing changes to the release branches
------------------------------------------

If you wish to have a change you made in the main devel branch (aka
the trunk) also be pushed back into the current release branch you
first need to take note of the revision number from your commit, for
example,

    $ svn commit -m"Sample commit"
    Adding         Rpacks\Biobase\DESCRIPTION
    Sending        Rpacks\Biobase\DESCRIPTION
    Transmitting file data ..
    Committed revision 140.

was revision 140. The changeset you want is -r 139:140, i.e. the changes
between r139 and r140.

You'll need to have checked out the branches subdirectory, which is
separate from trunk. If you have only checked out the madman
subdirectory previously, you'll need to also check out the appropriate
branch subdirectory:

    svn co https://hedgehog.fhcrc.org/bioconductor/branches/RELEASE_2_6/madman/Rpacks

Change to the branch directory, merge your changes, check and fix any
conflicts, and commit. So, from your branch directory
(e.g. RELEASE_2_6/madman/Rpacks):

    svn merge -r 139:140 https://hedgehog.fhcrc.org/gentleman/bioconductor/trunk/madman/Rpacks
    svn status   # Look or C, indicating a conflict
                 # fix conflicts... (remember to use svn resolve for each)
    svn commit -m "merged -r139:140 from trunk"

Having problems?
----------------

Here is a list of possible issues:

unrecognized URL scheme:
  If you see "unrecognizedURL scheme" when trying to access the
  repository, it may indicate that your svn client does not support
  HTTPS.  You can verify the supported "modes" by examinging the
  output of ``svn --version``.  If you do not see support for HTTPS,
  you will need to upgrade your client.

Username or password not recognized:
  Most usernames we issue are in the form of an email address.  It may
  help to quote the user name: ``svn co
  --username="bob@internet.net"``.  Try double and single quotes.  
