# Add collaborators and leverage GitHub features

__Goal:__ You would like to take advantage of the social coding
features provided by GitHub, while continuing to update your
_Bioconductor_ repository.

## Maintaining Collaborators on GitHub

1. [Adding a new collaborator][]

2. [Removing collaborator][]

## Pull requests on GitHub

1. [Merging a pull request][]

## Push GitHub changes to the _Bioconductor_ repository

Once you have accepted pull requests from your package community on
GitHub, you can push these changes to _Bioconductor_.

1. Make sure that you are on the branch to which the changes were
   applied, for example `master`.

        git checkout master

1. Fetch and merge the GitHub changes to your local repository.

        git fetch origin
        git merge

    [Resolve merge conflicts][] if necessary.

1. Push your local repository to the upstream _Bioconductor_ repository.

        git push upstream master
    
    To push GitHub release branch updates to the _Bioconductor_
    release branch, replace `master` with name of the release branch,
    e.g.: `RELEASE_3_6`.

[Adding a new collaborator]: https://help.github.com/articles/inviting-collaborators-to-a-personal-repository/
[Removing collaborator]: https://help.github.com/articles/removing-a-collaborator-from-a-personal-repository/
[Merging a pull request]: https://help.github.com/articles/merging-a-pull-request/
[Resolve merge conflicts]: ../resolve-conflicts
