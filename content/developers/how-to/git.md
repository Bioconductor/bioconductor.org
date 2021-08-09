# ![](/images/icons/magnifier.gif)Source Control

The _Bioconductor_ project is maintained in a Git source control
system. Package maintainers update their packages by pushing changes
to their git repositories.

The essential steps for transitioning from SVN to git are summarized
in

- [New package workflow][]: update your GitHub repository after
  package acceptance.
- [Create a new GitHub repository][] for an existing package.
- [Maintain a _Bioconductor_-only repository][] for an existing
  package.

[New package workflow]: new-package-workflow
[Create a new GitHub repository]: maintain-github-bioc
[Maintain a _Bioconductor_-only repository]: maintain-bioc-only

More scenarios for repository creation:

- [Sync an existing GitHub repository][] with _Bioconductor_.
- [Create a local repository][] for private use.

[Sync an existing GitHub repository]: sync-existing-repositories
[Create a local repository]: create-local-repository

Scenarios for code update:

- [Pull upstream changes][], e.g., introduced by the core team.
- [Push to GitHub and _Bioconductor_][] repositories.
- [Resolve merge conflicts][].
- [Abandon changes][].
- [Fix bugs in  devel and release][].

[Pull upstream changes]: pull-upstream-changes
[Push to GitHub and _Bioconductor_]: push-to-github-bioc
[Resolve merge conflicts]: resolve-conflicts
[Abandon changes]: abandon-changes
[Fix bugs in devel and release]: bug-fix-in-release-and-devel

Github scenarios

- [Add collaborators][] and use Github social coding features.
- [Change Package Maintainer][]
- [Remove Large Data Files][] and clean git tree

See other **[frequently asked questions][]**.

[Add collaborators]: add-collaborators
[Change Package Maintainer]: change-maintainer
[Remove Large Data Files]: remove-large-data
[frequently asked questions]: faq

## Essential work flow

A minimal workflow is to checkout, update, commit, and push changes to
your repository. Using `BiocGenerics` as an example:

    git clone git@git.bioconductor.org:packages/BiocGenerics
    cd BiocGenerics
    ## add a file, e.g., `touch README`
    ## edit file, e.g., `vi DESCRIPTION`
    BiocGenerics$ git commit README DESCRIPTION
    BiocGenerics$ git push

This requires that _Bioconductor_ knows the SSH keys you use to
establish your identity.

Two useful commands are

    BiocGenerics$ git diff     # review changes prior to commit
    BiocGenerics$ git log      # review recent commits

If the repository is already cloned, the work flow is to make sure
that you are on the 'master' branch, pull any changes, then introduce
your edits.

    BiocGenerics$ git checkout master
    BiocGenerics$ git pull
    ## add, edit, commit, and push as above

## Where to Commit Changes

New features and bug fixes are introduced on the master ('devel')
branch of the GIT repository.

    BiocGenerics$ git checkout master
    BiocGenerics$ git pull
    ## edit 'R/foo.R' and commit on master
    BiocGenerics$ git commit R/foo.R  #
    [master c955179] your commit message
    1 file changed, 10 insertions(+), 3 deletions(-)
    BiocGenerics$ git push

To make more extensive changes see [Fix bugs in devel and release][].

Bug fixes can be ported to the current release branch. Use
`cherry-pick` to identify the commmit(s) you would like to port. E.g.,
for release 3.6, porting the most recent commit to master

    BiocGenerics$ git checkout RELEASE_3_6
    BiocGenerics$ git cherry-pick master
    BiocGenerics$ git push

## Checks and version bumps

Each commit pushed to the _Bioconductor_ repository should build and
check without errors or warnings

    BiocGenerics$ cd ..
    R CMD build BiocGenerics
    R CMD check BiocGenerics_1.22.3.tar.gz

Each commit, in either release or devel, should include a bump in the
`z` portion of the `x.y.z` package [versioning scheme][].

Builds occur once per day, and take approximately 24 hours. See the
[build report][] for git commits captured in the most recent build
(upper left corner)

[versioning scheme]: /developers/how-to/version-numbering/
[build report]: https://bioconductor.org/checkResults/devel/bioc-LATEST/

## Annotation packages

Traditional Annotation packages are not stored in GIT due to the size of
annotation files. To update an existing Annotation package please send an email
to maintainer@bioconductor.org. A member of the Bioconductor team will be in
contact to receive the updated package.

Newer annotation packages can be stored in GIT as it is a requirement to use the
[AnnotationHub](https://bioconductor.org/packages/AnnotationHub/) or similar
server hosted data. The larger sized files are not included directly in the package.
To contribute a new Annotation package please contact hubs@bioconductor.org
for guidance and read the documentation on [How to Create A Hub
package](https://bioconductor.org/packages/devel/bioc/vignettes/HubPub/inst/doc/CreateAHubPackage.html).

Currently direct updates to annotation packages, even those stored on git, are
not supported. If you wish to updated an annotation package, make required
changes and push to git.bioconductor.org. Then send an email to
hubs@bioconductor.org or maintainer@bioconductor.org requesting the package be propagated.


## More help

Need more help? Ask on the [bioc-devel](/help/mailing-list/) mailing
list.
