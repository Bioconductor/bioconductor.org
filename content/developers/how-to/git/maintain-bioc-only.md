# Maintaining a _Bioconductor_-only repository

__Goal:__ Developer wishes to maintain their _Bioconductor_ repository
without using GitHub.

## Steps:

### Clone and setup the package on your local machine.

1. Make sure that you have `SSH` access to the _Bioconductor_
   repository; be sure to [submit your SSH public key][submit-keys] or
   github id to _Bioconductor_.

1. Clone the package to your local machine,

        git clone git@git.bioconductor.org:packages/<ExamplePackage>

    __NOTE:__ If you clone with `https` you will NOT get read+write
    access.

1. View existing remotes

        git remote -v

    which will display

        origin    git@git.bioconductor.org:packages/<ExamplePackage>.git (fetch)
        origin    git@git.bioconductor.org:packages/<ExamplePackage>.git (push)

    This indicates that your git repository has only one remote
    `origin`, which is the _Bioconductor_ repository.

### Commit changes to your local repository

1. Before making changes to your repository, make sure to `pull`
   changes or updates from the _Bioconductor_ repository. This is
   needed to avoid conflicts.

        git pull

1. Make the required changes, then `add` and `commit` your changes to
   your `master` branch.

        git add <files changed>
        git commit -m "My informative commit message"

1. (Alternative) If the changes are non-trivial, create a new branch
   where you can easily abandon any false starts. Merge the final
   version onto `master`

        git checkout -b feature-my-feature
        ## add and commit to this branch. When the change is complete...
        git checkout master
        git merge feature-my-feature

### Push your local commits to the _Bioconductor_ repository

1. Push your commits to the _Bioconductor_ repository to make them
   available to the user community.

   Push changes to the `master` branch using:

        git checkout master
        git push origin master

### (Optional) Merge changes to the current release branch

1. If the changes (bug fixes) need to be available on the
   `<RELEASE_X_Y>` branch. Checkout the branch `<RELEASE_X_Y>` using:

        git checkout -b <RELEASE_X_Y>

1. Pull any changes from the _Bioconductor_ repository to this branch

        git pull

1.  Then, merge `master` into the `<RELEASE_X_Y>` branch

        git merge master

    A new commit message will show up after this step, with the
    message `Merge branch 'master' into <RELEASE_X_Y>`. It will be
    committed to `<RELEASE_X_Y>` branch when you save. Update the
    version number so that it is correct for the release, and commit
    the modified DESCRIPTION file

        git add DESCRIPTION
        git commit -m "Update version bump"

1. Push changes to the `<RELEASE_X_Y>` branch using:

        git push origin <RELEASE_X_Y>

[submit-keys]: https://goo.gl/forms/eg36vcBkIUjfZfLe2
