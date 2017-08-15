# Pull upstream changes

__Goal:__ Your _Bioconductor_ repository has been updated by the core
team. You want to fetch these commits from _Bioconductor_, merge them
into your local repository, and push them to GitHub.

__NOTE:__ It is always a good idea to fetch updates from
_Bioconductor_ before making more changes. This will help prevent
merge conflicts.

## Steps

These steps update the `master` branch.

1. Make sure you are on the appropriate branch.

        git checkout master

1. Fetch content from _Bioconductor_

        git fetch upstream

1. Merge upstream with the appropriate local branch

        git merge upstream/master

    Get help on [Resolve merge conflicts][] if these occur.

1. Push changes to GitHub's (`origin`) `master` branch

        git push origin master

To pull updates to the current `RELEASE_X_Y` branch, replace `master`
with `RELEASE_X_Y` in the lines above.

See instructions to [Sync existing repositories][] with
changes to both the _Bioconductor_ and GitHub repositories.

[Resolve merge conflicts]: ../resolve-conflicts
[Sync existing repositories]: ../sync-existing-repositories
