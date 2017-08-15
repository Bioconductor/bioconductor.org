# Abandon changes

__Goal:__ You want to start fresh after failing to resolve conflicts
or some other issue. If you intend to go nuclear, please contact the
bioc-devel@bioconductor.org mailing list.

# Going sub-nuclear: force _Bioconductor_ `master` to GitHub `master`

One way you can ignore your work and make a new branch is by replacing
your local and GitHub repository `master` branch with the
_Bioconductor_ `master` branch.

NOTE: This works only if you haven't pushed the change causing the
issue to the _Bioconductor_ repository.

## Steps:

1. Checkout a new branch, e.g., `master_backup`, with tracking set to
   track the _Bioconductor_ `master` branch `upstream/master`.

        git checkout -b master_backup upstream/master

1. Rename the branches you currently have on your local
   machine. First, rename `master` to `master_deprecated`. Second,
   rename `master_backup` to `master`. This process is called the
   classic __Switcheroo__.

        git branch -m master master_deprecated
        git branch -m master_backup master

1. You will now have to "force push" the changes to your GitHub
   (`origin`) `master` branch.

        git push -f origin master

1. __(Optional)__ If you have commits on your `master_deprecated`
   branch that you would like ported on to your new `master`
   branch. Git has a special feature called `cherry-pick`

   Take a look at which commit you want to cherry-pick on to the new
   master branch, using `git log master_deprecated`, copy the correct
   commit id, and use:

        git cherry-pick <commit id>

   Push these cherry-picked changes to GitHub and _Bioconductor_
   repositories.

# Reset to a previous commit

If you find yourself in a place where you want to abandon changes
__already committed__ to _Bioconductor_ or GitHub, use `reset` to undo
the commits on your local repository and `push -f` to force the
changes to the remotes. Remember that the `HEAD` commit id is the most
recent __parent__ commit of the current state of your local
repository.

    git reset --hard <commit id>

Example:

    git reset --hard e02e4d86812457fd9fdd43adae5761f5946fdfb3                                                        master
    HEAD is now at e02e4d8 version bump by bioc core

To make the changes permanent, you will then need to push the changes
to both GitHub and _Bioconductor_:

    git push -f origin
    git push -f upstream

# Go Nuclear - Delete your local copy and GitHub repo, because nothing is working.

__CAUTION: These instructions come with many disadvantages. You have
been warned.__

## Steps:

1. Delete your local repository, e.g., `rm -rf BiocGenerics`

1. Delete (or rename) your GitHub repository.

1. [Maintain GitHub and _Bioconductor_ repositories][] for an existing
   _Bioconductor_ repository, then [pull upstream changes][].

## Disadvantages of going "nuclear":

1. You will lose all your GitHub issues

1. You will lose your custom collaborator settings in GitHub.

1. You will lose any GitHub-specific changes.

[Maintain GitHub and _Bioconductor_ repositories]: ../maintain-github-bioc
[pull upstream changes]: ../pull-upstream-changes
