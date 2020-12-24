# Remote and local branch names differ

__Goal:__ You'd like to working with 'main' as your local repository
branch, while still pushing changes to the 'master' branch of
_Bioconductor_'s git repository.

A common reason for this scenario is because GitHub uses `main` as the
name of the default branch, while the _Bioconductor_ git repository
uses the older term `master` for the default branch. The `master` git
branch of _Bioconductor_ is used for _Bioconductor_'s 'devel' package
nightly builds.

The scenario might also occur because you wish to maintain a GitHub or
other repository with a default branch containing material more
comprehensive than only the package source tree with its `R/`, `man/`,
`vignettes/` and related directories. For instance, the default branch
of your repository might contain material for github.io 'pages', or
with extensive GitHub actions for continuous integration defined on
the branch; you might then maintain on your local and GitHub
repositories a `package-only` branch, periodically merging appropriate
changes from other branches and pushing `package-only` to the `master`
branch of _Bioconductor_.


## Steps:

The following uses symbols in angle brackets `<>` to indicate values
that need to be substituted based on your situation.

- `<package>`: the name of your package.

- `<package-dir>`: the directory on your computer in which your
  package has been cloned.

- `<bioc-branch>`: the name of the branch in the _Bioconductor_ git
  repository. New packages and features are added to the `master`
  branch of _Bioconductor_ git repositories, and the `master` branch
  is used by the nightly builders to generate the 'devel' repository
  of _Bioconductor_ packages. Release branches have the form
  `RELEASE_3_12`, corresponding to the git branch used by the nightly
  builders to produce the (in this case) _Bioconductor_ version 3.12
  packages. The most common scenario is the `<bioc-branch>`
  corresponds to the `master` branch of _Bioconductor_'s git.

- `<local-branch>`: the name of the branch in your local repository
  that you intend to use to correspond to `<bioc-branch>`. The most
  common scenario is that `<local-branch>` corresponds to the `main`
  GitHub branch.


Suppose you already have a local clone of your original GitHub or
other repository.

1.  Add the _Bioconductor_ git as a remote to your local
    repository. By convention we name the remote as `upstream`.

    ```
    cd <package-dir>
    git remote add upstream git@git.bioconductor.org:packages/<package>.git
    ```

1. Make sure that the `<local-branch>` of your repository is in sync with the
   _Bioconductor_ repository.

   ```
   git checkout <local-branch>
   git fetch upstream
   git merge upstream/<bioc-branch>
   ```

   It may be necessary to [resolve merge conflicts][].

1. Commit changes as usual to `<local-branch>`.

1. Push changes from `<local-branch>` to the upstream `<bioc-branch>`;
   remember that changes only propagate with version number
   increments in the DESCRIPTION file.

   ```
   git push upstream <local-branch>:<bioc-branch>
   ```

[resolve merge conflicts]: ../resolve-conflicts
