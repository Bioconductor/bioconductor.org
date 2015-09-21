# Using Git with Bioconductor SVN repositories #

Bioconductor packages are stored in the Bioconductor Subversion (SVN)
repository.  However Bioconductor also maintains Read-only mirrors for each
package located on [GitHub](https://github.com/Bioconductor-mirror).

These instructions detail how to use these mirrors as well as
[git-svn](http://git-scm.com/docs/git-svn) to interact with the Subversion
repositories from git.

*The examples below use HTTPS authentication, however you are free to
substitute using SSH if you prefer, in all cases USER should be replaced by
your GitHub username, and REPO should be replaced by your package name.*

## Install Git-Svn Pre-Requisites ##

To install git with svn support, see the following OS specific
instructions.

### Windows ###

Download and install the [Git windows client](https://www.git-scm.com/download/win).

When using Windows, the below commands should be run using the Git Bash client
which is bundled with git.

### Mac OS X - Homebrew ###

`git-svn` *should* already be installed. Run the command

```
git svn --help
```

If it shows a help page, then git-svn is already installed. Otherwise, 
install [Homebrew](http://brew.sh/) and try:

```
brew update
brew install git
```

### Linux - Debian Based ###

```
sudo apt-get install git-svn
```

## Download the update_remotes.sh script ##

This script will automatically create local branches to track the development
and release branches from svn and the associated git mirrors. You can download this script
to any directory and (optionally) put it in your PATH.

`curl -O https://raw.githubusercontent.com/Bioconductor/mirror/master/update_remotes.sh`

## Setup ##

### Scenario 1: Use Git Locally, No GitHub Repository ###

If you simply want to use git locally on your machine and do not need to have a
publicly accessible git repository on GitHub (or elsewhere) you can simply
clone your package from the mirror directly.

If you are using the Git-Svn bridge, **please** delete your bridge (see
instructions in the next section) before continuing.


  1. `git clone https://github.com/Bioconductor-mirror/REPO` to clone the repository to your machine.
  2. `cd REPO` to switch to the REPO directory.
  3. `bash /path/to/update_remotes.sh` to setup the git remotes.
  4. Commit to git as you normally would.
  5. Each time you want to push git commits to svn run the following commands:
     1. `git svn rebase` to get the latest SVN changes.
     2. `git svn dcommit --add-author-from` to commit your changes to SVN.
     You may be prompted here for your SVN username and password.

### Scenario 2: Set Up Your Own GitHub Repository ###

If you are currently using the Git-Svn Bridge please disable it at
<https://gitsvn.bioconductor.org/>. 
[Log in](https://gitsvn.bioconductor.org/login), 
[list your bridges](https://gitsvn.bioconductor.org/my_bridges),
and delete the one you want to migrate. **If you migrate to
git mirrors without deleting your Git-svn bridge, your
repository could be damaged!**

If you do not already have a public git repository for package REPO the
simplest thing to do is navigate to
`https://github.com/Bioconductor-mirror/REPO` and click the `Fork` button in
the upper right.  This will create a copy of the repository on your personal
account. You may want to re-enable issue tracking in your repository
(it's disabled in the read-only mirrors and forks inherit
this setting). To do this, go to Settings and then click the Issues
checkbox.
Then perform the following steps in your terminal.

  1. `git clone https://github.com/USER/REPO` to clone the repository to your machine.
  2. `cd REPO` to switch to the REPO directory.
  3. `bash /path/to/update_remotes.sh` to setup the git remotes.
  4. Commit to git and push to GitHub as you normally would.
  5. Each time you want to push git commits to svn:
     1. `git checkout devel` to switch to the devel branch. (use release-X.X for release branches)
     2. `git svn rebase` to get the latest SVN changes.
     3. `git merge master --log` to merge your changes from the master branch or skip this step and work directly on the current branch.
     4. `git svn dcommit --add-author-from` to sync and commit your changes to svn.
     You may be prompted here for your SVN username and password.


When you're done, be sure and merge any changes from svn back into the git master branch:

    git checkout master
    git merge devel


 
## FAQs ##

### How do I let users know I am using GitHub for development and contributions?

Add `URL: https://github.com/USER/REPO` and `BugReports:
https://github.com/USER/REPO/issues` to your `DESCRIPTION` file. You can also
mention your repository on the bioc-devel [mailing
list](http://bioconductor.org/help/mailing-list/).

### I don't know my Subversion username and/or password. What do I do? ###

One of the following steps should work:

* Look in your email. Your SVN credentials were originally sent to you
  by a member of the Bioconductor team, probably
  with the subject line "congrats" or "congratulations". The email 
  by a member of the Bioconductor team,
  probably
  with a subject line containing
  "congrats" or "congratulations". The email   should contain the text "Information about your svn account". 
* Go to your `~/.subversion/auth/svn.simple` directory. There should be
  one or more files whose names are long hexadecimal numbers. Use `grep`
  to find out which file contains your username. If you don't know your 
  username,
  it's usually your first initial, a dot, and your last name (all 
  lowercase). So Jill User would be `j.user`. Example:

        $ grep -l j.user *
        81a52e36a28dfd7750bd975f30c7998b

  This indicates that your password can be found in the file called
  `81a52e36a28dfd7750bd975f30c7998b`. Examine that file and you should see 
  something like:

        password
        V 8
        Z7oRUVH6

  In this case, `Z7oRUVH6` is your password.
* If you still can't find your username or password, contact a 
  member of the Bioconductor team at
  `maintainer at bioconductor dot org`. Mention the package(s) that
  you maintain. We cannot send you your password but we can ask for 
  a new one to be generated, and send it to you. It may take 
  a day or two for the request to be processed.

### How do I commit to the release version of my package? ##

If you are cloning the mirrors directly you can switch to the `release-X.X`
branch of the release you would like to commit to, and then proceed as normal.
(Current release version is <%= config[:release_version]%>.)
If you are hosting on GitHub as well, rather than checking out `bioc/master`
checkout `bioc/release-X.X`, then perform the rest of the steps as normal.

### How do I get the SVN revision for a git commit, or the git commit from the SVN revision? ###

`git svn find-rev` can be used for both directions.

```
git svn find-rev r104237
# 3a5e1d5995322fb0569138930c3d2aaa93b1c54d

git svn find-rev 3a5e1d5995322fb0569138930c3d2aaa93b1c54d
# 104237
```


### How can I search all Bioconductor code on GitHub?

[Here](https://github.com/search?utf8=%E2%9C%93&q=user%3Abioconductor-mirror+goana.default&type=Code&ref=searchresults)'s
an example search that looks through all Bioconductor 
software packages for the specified string.

### So can I submit a pull request against any Bioconductor package?

No. The read-only mirrors are set up to automatically
reject any pull request (unfortunately, GitHub does not let
us disable them altogether). If you want to submit a pull request,
you need to submit it to the repository that the maintainer
of the package maintains (if they have chosen to
create one). Check the `BugReports` or
`URL` fields of the `DESCRIPTION` file to find the link
to this repository.


### I have a question or comment.

Please send your feedback to the 
[bioc-devel mailing list](/help/support/#bioc-devel).
 

## Troubleshooting #

### Unable to determine upstream SVN information

The dreadful message indicating that `git` and `svn` got out of sync is `Unable
to determine upstream SVN information from working tree history`. This can
happen, for example if one forgets to `git pull --rebase` before trying to 
`git svn dcommit` and changes were committed to svn independently. Inspect your git 
log with `git log --graph --decorate --pretty=oneline --abbrev-commit --all` to help 
identify such cases.


Useful references to sort such cases out are:

- [http://stackoverflow.com/questions/9805980/unable-to-determine-upstream-svn-information-from-working-tree-history](http://stackoverflow.com/questions/9805980/unable-to-determine-upstream-svn-information-from-working-tree-history)
- [http://eikke.com/importing-a-git-tree-into-a-subversion-repository/](http://eikke.com/importing-a-git-tree-into-a-subversion-repository/)

## Resources #

* [GitHub's Git help](https://help.github.com/)
* [Good Resources For Learning Git and GitHub](https://help.github.com/articles/good-resources-for-learning-git-and-github/)
* [Most common git screwups/questions and solutions](http://41j.com/blog/2015/02/common-git-screwupsquestions-solutions/)
* [Flight rules for git](https://github.com/k88hudson/git-flight-rules)
