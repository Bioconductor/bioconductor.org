# New package workflow

__Goal__: You developed a package in GitHub, following the
_Bioconductor_ new package [Contributions README][] guidelines,
[submitted it to _Bioconductor_][], and your package has been 
moderated. As part of moderation process, the package to be 
reviewed has been added to the _Bioconductor_ git repository. 
Package authors will learn to push to the _Bioconductor_ git
repository, from the time a package has been submitted to 
_Bioconductor_.

After your package has been moderated and assigned a reviewer, 
it is available on the [_Bioconductor_ git server][].

- __SSH (developer) read / write access:__ `git@git.bioconductor.org`

- __HTTPS (public) read only access:__ `https://git.bioconductor.org`

## Steps:

1. _Bioconductor_ needs to know your SSH 'public key'. _Bioconductor_
   will use keys in `https://github.com/<your-github-id>.keys`. This step
   is now mandated by the package submission process, so it should have been
   completed if you have gone through the first step of the submission process.

   Alternatively, [submit your SSH public key][submit-keys] or github
   id to _Bioconductor_ via the [BiocCredentials application][].

1. Configure the "remotes" of your local git repository. You will need
   to push any future changes to your package to the _Bioconductor_ git
   repository to issue a new build of your package. 
   
   Once the package is accepted, you pull changes the _Bioconductor_ core 
   team will make (e.g., bug fixes or bumping a version number for a new
   release). 
   
   Add a remote to your package's local git repository using:

        git remote add upstream git@git.bioconductor.org:packages/<YOUR-REPOSITORY-NAME>.git

1. **Start a new build**: Your package version number is in the format 
   'major.minor.patch'. To start a new build on our server, you must only
   update the 'patch' version. Newly submitted packages will have versions
   similar to `0.99.0`. Then, push to the _Bioconductor_ git repository. 

   To push,
   
        git add <files changed including version bump in DESCRIPTION file>
        
        git commit -m "<informative commit message>"
        
        git push upstream master

1. [Pull upstream changes][] made by the _Bioconductor_ core team
   during addition of your repository. This step is relevant once you
   package is accepted by Bioconductor. You will need to pull the changes
   made by the Bioconductor team.

1. See other scenarios for working with _Bioconductor_ and GitHub repositories, in particular:

    - [Maintain GitHub and _Bioconductor_ repositories][].
    - [Fix bugs in  devel and release][].

[submit-keys]: https://git.bioconductor.org/BiocCredentials/
[Maintain GitHub and _Bioconductor_ repositories]: ../maintain-github-bioc
[Pull upstream changes]: ../pull-upstream-changes
[Fix bugs in devel and release]: ../bug-fix-in-release-and-devel
[Contributions README]: https://github.com/Bioconductor/Contributions
[_Bioconductor_ git server]: https://git.bioconductor.org
[submitted it to _Bioconductor_]: http://bioconductor.org/developers/package-submission/
[BiocCredentials application]: https://git.bioconductor.org/BiocCredentials/
