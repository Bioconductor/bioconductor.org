# New package workflow

__Goal__: You developed a package in GitHub, following the
_Bioconductor_ new package [Contributions README][] guidelines,
[submitted it to _Bioconductor_][], and your package has been
moderated. As part of moderation process, the package to be reviewed
has been added as a repository on the [_Bioconductor_ git server][1].

During and after the review process, package authors must push changes
that include a **version number 'bump'** to the _Bioconductor_ git
repository. This causes the package to be built and checked on Linux,
macOS, and Windows operating systems, and forms the basis for the
review process.

In this document, package authors will learn best practices for
pushing to the _Bioconductor_ git repository.

## Steps:

1. **SSH keys**. As part of the initial moderation step,
   _Bioconductor_ will use SSH 'public key' keys available in
   `https://github.com/<your-github-id>.keys`.

   After the review process is over, additional SSH keys can be added
   and contact information edited using the [BiocCredentials
   application][].

1. **Configure the "remotes" of your local git repository**. You will need
   to push any future changes to your package to the _Bioconductor_ git
   repository to issue a new build of your package. 
   
   Add a remote named `upstream` to your package's local git
   repository using:

        git remote add upstream git@git.bioconductor.org:packages/<YOUR-REPOSITORY-NAME>.git
        
   Check that you have updated the remotes in your repository; you'll
   see an `origin` remote pointing to `github.com`, and an `upstream` remote
   pointing to `bioconductor.org`
   
        $ git remote -v

        origin  <link to your github> (fetch)
        origin  <link to your github> (push)
        upstream git@git.bioconductor.org:packages/<YOUR-REPOSITORY-NAME>.git (fetch)
        upstream git@git.bioconductor.org:packages/<YOUR-REPOSITORY-NAME>.git (push)
        
   NOTE: As a package developer, you must use the SSH protocol (as in
   the above command) to gain read/write access to your package in the
   _Bioconductor_ git repository.

1. **Add and commit changes to your local repository**. During the
   review process you will likely need to update your package. Do this
   in your local repository by first making sure your repository is
   up-to-date with your `github.com` and `git.bioconductor.org`
   repositories.

        git fetch --all

        ## merge changes from git.bioconductor.org
        git merge upstream/master

        ## merge changes from github.com. If your github default branch name
        ## is main, replace origin/master with origin/main
        git merge origin/master

   Make changes to your package master branch and commit them to your
   local repository
   
        git add <files changed>
        git commit -m "<informative commit message>"
        
1. **'Bump' the package version**.  Your package version number is in
   the format 'major.minor.patch'. When the review process starts, the
   version number is `0.99.0`. Increment the `patch` version number by
   1, e.g., to `0.99.1`, `0.99.2`, ..., `0.99.9`, `0.99.10`, ...
   
   Bumping the version number before pushing is essential. It ensures
   that the package is built across platforms.
   
   Remember to add and commit these changes to your local repository.

1. **Push changes to the Bioconductor and github repositories**.  Push
   the changes in your local repository to the _Bioconductor_ and
   github repositories.

        ## push to git@git.bioconductor.org. If your github default
        ## branch name is main, replace master with main:master
        git push upstream master

        ## push to your github repository
        ## If your github default branch name is main, replace master with main
        git push origin master

1. **Check the updated build report**. If your push to
   git.bioconductor.org included a version bump, you'll receive an
   email directing you to visit your issue on github.com,
   `https://github.com/Bioconductor/Contributions/issues/`; a comment
   is also posted on the issue indicating that a build has started.
   
   After several minutes a second email and comment will indicate that
   the build has completed, and that the build report is
   available. The comment includes a link to the build report. Follow
   the link to see whether further changes are necessary.

1. See other scenarios for working with _Bioconductor_ and GitHub repositories, in particular:

    - [Maintain GitHub and _Bioconductor_ repositories][].
    - [Fix bugs in  devel and release][].
    - [Resolve merge conflicts][].

[submit-keys]: https://git.bioconductor.org/BiocCredentials/
[Maintain GitHub and _Bioconductor_ repositories]: ../maintain-github-bioc
[Pull upstream changes]: ../pull-upstream-changes
[Fix bugs in devel and release]: ../bug-fix-in-release-and-devel
[Resolve merge conflicts]: ../resolve-conflicts
[Contributions README]: https://github.com/Bioconductor/Contributions
[_Bioconductor_ git server]: https://git.bioconductor.org
[submitted it to _Bioconductor_]: http://bioconductor.org/developers/package-submission/
[BiocCredentials application]: https://git.bioconductor.org/BiocCredentials/
