# ![](/images/icons/magnifier.gif)Source Control

The _Bioconductor_ project is maintained in a Git source control
system. Package maintainers update their packages by pushing changes
to their git repositories.

The essential steps for transitioning from SVN to git are summarized
in

- [New package workflow][]: update your GitHub repository after
  package acceptance.
- [Maintain GitHub and Bioconductor repositories][] of an existing
  package.
- [Maintain a _Bioconductor_-only repository][] of an existing
  package.

[New package workflow]: new-package-workflow
[Maintain GitHub and Bioconductor repositories]: maintain-github-bioc
[Maintain a _Bioconductor_-only repository]: maintain-bioc-only

More scenarios for repository creation:

- [Sync existing repositories][] with _Bioconductor_.
- [Create a local repository][] for private use.

[Sync existing repositories]: sync-existing-repositories
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

See other [frequently asked questions][].

[Add collaborators]: add-collaborators
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
for release 3.5, porting the most recent commit to master

    BiocGenerics$ git checkout RELEASE_3_5
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

Due to the size of the Annotation files these packages are not stored in svn.
To update an existing Annotation package please make the new version available 
in dropbox (or similar) and send an email to maintainer@bioconductor.org. A 
member of the Bioconductor team will add the package to the appropriate 
repository.

To contribute a new Annotation package please contact packages@bioconductor.org
for guidance.

## Common problems

### What SSH keys does _Bioconductor_ know about?

### More help

Need more help? Ask on the [bioc-devel](/help/mailing-list/) mailing
list.
