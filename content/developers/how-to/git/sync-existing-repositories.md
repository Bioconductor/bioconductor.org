# Sync existing repositories with _Bioconductor_

__Goal:__ Ensure that your local, _Bioconductor_, and GitHub
repositories are all in sync.

## Steps:

1. Clone the the GitHub repository to a local machine.

1. Configure the "remotes" of the GitHub clone.

        git remote add upstream git@git.bioconductor.org:packages/<YOUR-REPOSITORY-NAME>.git`

1. Fetch updates from all (_Bioconductor_ and GitHub) remotes.

        git fetch --all

1. Make sure you are on the master branch.

        git checkout master

1. Merge updates from the _Bioconductor_ (`upstream`) remote

        git merge upstream/master

   [Resolve merge conflicts][] if necessary.

1. Merge updates from the GitHub (`origin`) remote

        git merge origin/master

1. Push to both _Bioconductor_  and GitHub repositories.

        git push upstream master
        git push origin master

1. Repeat for the release branch, replacing `master` with the name of
   the release branch, e.g., `RELEASE_3_5`. It may be necessary to
   create the release branch in the local repository.
   
        git checkout RELEASE_3_5
        git merge upstream/RELEASE_3_5
        git merge origin/RELEASE_3_5
        git push upstream RELEASE_3_5
        git push origin RELEASE_3_5
   
   Remember that only `master` and the current release branch of
   _Bioconductor_ repositories can be updated.

[Resolve merge conflicts]: ../resolve-conflicts
