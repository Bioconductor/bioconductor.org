# Sync an existing GitHub repository with _Bioconductor_

__Goal:__ Ensure that your local, _Bioconductor_, and GitHub
repositories are all in sync.

## Steps:

1. Clone the GitHub repository to a local machine. Change into the
   directory containing the repository.

1. Configure the "remotes" of the GitHub clone.

        git remote add upstream git@git.bioconductor.org:packages/<YOUR-REPOSITORY>.git

1. Fetch updates from all (_Bioconductor_ and GitHub) remotes. You may
   see "warning: no common commits"; this will be addressed after
   resolving conflicts, below.

        git fetch --all

1. Make sure you are on the master branch.

        git checkout master

1. Merge updates from the GitHub (`origin`) remote

        git merge origin/master

1. Merge updates from the _Bioconductor_ (`upstream`) remote

        git merge upstream/master

   Users of git version >= 2.9 will see an error message ("fatal:
   refusing to merge unrelated histories") and need to use

        git merge --allow-unrelated-histories upstream/master

1. [Resolve merge conflicts][] if necessary.

1. After resolving conflicts and committing changes, look for duplicate
   commits (e.g., `git log --oneline | wc` returns twice as many
   commits as in SVN) and consider following the steps to
   [force Bioconductor `master` to GitHub `master`][force-bioc-to-github].

1. Push to both _Bioconductor_  and GitHub repositories.

        git push upstream master
        git push origin master

1. Repeat for the release branch, replacing `master` with the name of
   the release branch, e.g., `RELEASE_3_6`. It may be necessary to
   create the release branch in the local repository.

        git checkout RELEASE_3_6
        git merge upstream/RELEASE_3_6
        git merge origin/RELEASE_3_6
        git push upstream RELEASE_3_6
        git push origin RELEASE_3_6

    NOTE: If you are syncing your release branch for the first time,
    you have to make a local copy of the `RELEASE_X_Y` branch, by

        git checkout -b <RELEASE_X_Y> upstream/<RELEASE_X_Y>

    Following this one time local checkout, you may switch between
    `RELEASE_X_Y` and `master` with `git checkout <RELEASE_X_Y>`. If you do
    not use the command to get a local copy of the release branch, you
    will get the message,

       (HEAD detached from origin/RELEASE_X_Y)


   Remember that only `master` and the current release branch of
   _Bioconductor_ repositories can be updated.

[Resolve merge conflicts]: ../resolve-conflicts
[force-bioc-to-github]: ../abandon-changes#force-bioconductor--to-github-
