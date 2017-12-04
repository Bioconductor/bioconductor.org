# Fix bugs in `devel` and  `release`

__Goal:__ Maintainers will have to fix bugs from time to time, and
make sure the patch is available both in the `master` branch (svn
devel) and the current `release` branch.

## Steps:

1. [Sync existing repositories][].

        git fetch --all
        git checkout master
        git merge upstream/master
        git merge origin/master
        git checkout <RELEASE_X_Y>
        git merge upstream/<RELEASE_X_Y>
        git merge origin/<RELEASE_X_Y>

1. On your local machine, be sure that you are on the `master` branch.

        git checkout master

   Make the changes you need to fix your bug. Add the modified files
   to the commit. Remember to edit the DESCRIPTION file to update the
   version number, and ensure that __only bug-fix changes__ are
   introduced in the commit.

        git add <files changed>
        ## after version bump
        git add DESCRIPTION

   Commit the modified files. It is helpful to tag the commit as bug
   fix.

        git commit -m "bug fix: my bug fix"

1. (Alternative) If the changes are non-trivial, create a new branch
   where you can easily abandon any false starts. Merge the final
   version onto `master`

        git checkout master
        git checkout -b bugfix-my-bug
        ## add and commit to this branch. When the bug fix is complete...
        git checkout master
        git merge bugfix-my-bug

1. Switch to the release branch and cherry-pick the commits from
   master corresponding to the bug fix, e.g., reference the most
   recent commit on the `master` branch with `master`; see examples in
   `git cherry-pick --help`. Remember to edit the DESCRIPTION file to
   update the release version according to _Bioconductor_'s
   [version numbering scheme][].

        git checkout <RELEASE_X_Y>
        git cherry-pick master
        ## Fix version bump by editing 'Version:' field of DESCRIPTION, then
        git add DESCRIPTION
        git commit -m "update RELEASE version number"

    NOTE: If you are patching your release for the first time, you have to make
    a local copy of the RELEASE_X_Y branch, by

        git checkout -b <RELEASE_X_Y> upstream/<RELEASE_X_Y>

    Following this one time local checkout, you may switch between RELEASE_X_Y
    and master with `git checkout <RELEASE_X_Y>`. If you do not use the command
    to get a local copy of the release branch, you will get the message,

        (HEAD detached from origin/RELEASE_X_Y)

1. Push your changes to both the GitHub and _Bioconductor_ `master`
   and `<RELEASE_X_Y>` branches. Make sure you are on the correct
   branch on your local machine.

   For the `master` branch,

        git checkout master
        git push upstream master
        git push origin master

   For the `release` branch,

        git checkout <RELEASE_X_Y>
        git push upstream <RELEASE_X_Y>
        git push origin <RELEASE_X_Y>

[version numbering scheme]: /developers/how-to/version-numbering
[Sync existing repositories]: ../sync-existing-repositories
