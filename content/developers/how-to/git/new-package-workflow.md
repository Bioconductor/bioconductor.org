# New package workflow

__Goal__: You have developed a package in GitHub, following the
_Bioconductor_ new package [Contributions README][] and other
guidelines, and your package has been accepted! The accepted package
has been added to the _Bioconductor_ git repository. Now what?

After your package has been accepted, it is visible on the
[_Bioconductor_ git server][].

- __SSH (developer) read / write access:__ `git@git.bioconductor.org`

- __HTTPS (public) read only access:__ `https://git.bioconductor.org`


## Steps:

1. _Bioconductor_ needs to know your SSH 'public key'. _Bioconductor_
   will use keys in `https://github.com/<your-github-id>.keys`.

   Alternatively, [submit your SSH public key][submit-keys] or github
   id to _Bioconductor_.

1. Configure the "remotes" of your local git repository. You will need
   to push any future changes to your package to the _Bioconductor_
   repository, and pull changes the _Bioconductor_ core team will make
   (e.g., bug fixes or bumping a version number for a new
   release). Add a remote to your machine's local git repository
   using:

        git remote add upstream git@git.bioconductor.org:packages/<YOUR-REPOSITORY-NAME>.git

1. [Pull upstream changes][] made by the _Bioconductor_ core team
   during addition of your repository.

1. See other scenarios for working with _Bioconductor_ and GitHub repositories, in particular:

    - [Maintain GitHub and _Bioconductor_ repositories][].
    - [Fix bugs in  devel and release][].

[submit-keys]: https://git.bioconductor.org/BiocCredentials/
[Maintain GitHub and _Bioconductor_ repositories]: ../maintain-github-bioc
[Pull upstream changes]: ../pull-upstream-changes
[Fix bugs in devel and release]: ../bug-fix-in-release-and-devel
[Contributions README]: https://github.com/Bioconductor/Contributions
[_Bioconductor_ git server]: https://git.bioconductor.org
