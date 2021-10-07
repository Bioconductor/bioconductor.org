# Resolve merge conflicts

__Goal:__ Resolve merge conflicts in branch and push to GitHub and
_Bioconductor_ repositories.

## Steps:

1. You will know you have a merge conflict when you see something like
   this:

        git merge upstream/master
        Auto-merging DESCRIPTION
        CONFLICT (content): Merge conflict in DESCRIPTION
        Automatic merge failed; fix conflicts and then commit the result

    This merge conflict occurs when the package developer makes a
    change, and also a collaborator or a _Bioconductor_ core team
    member makes a change to the same file (in this case the
    `DESCRIPTION` file).

    How can you avoid this? [pull upstream changes][] before
    committing any changes. In other words, `fetch` and `merge` remote
    branches before a `push`.

1. If in spite of this you have conflicts, you need to fix them. See
   which file has the conflict,

        git status

    This will show you something like this:

	    On branch master
	    Your branch is ahead of 'origin/master' by 1 commit.
          (use "git push" to publish your local commits)
	    You have un-merged paths.
	      (fix conflicts and run "git commit")
	      (use "git merge --abort" to abort the merge)

	    Un-merged paths:
	      (use "git add <file>..." to mark resolution)

          both modified:   DESCRIPTION

	    no changes added to commit (use "git add" and/or "git commit -a")

1. Open the file in your favorite editor. Conflicts look like:

	    <<<<<<< HEAD
	    Version: 0.23.2
	    =======
        Version: 0.23.3
	    >>>>>>> upstream/master

	Everything between `<<<<` and `=====` refers to HEAD, i.e your
    current change. And everything between `=====` and `>>>>>` refers
    to the `remote/branch` shown there, i.e `upstream/master`.

	You want to keep the most accurate change, by deleting what is
    necessary. In this case, keep the latest version:

        Version: 0.23.3

1. Add and commit the file as you would any other change.

	    git add DESCRIPTION
	    git commit -m "Fixed conflicts in version change"

1. 	Push to both your GitHub and _Bioconductor_ repositories,

	    git push origin master
	    git push upstream master

## Extra Resources

[Resolving a merge conflict using command line][]

[Based on: Resolving a merge conflict][]

<iframe width="560" height="315" src="https://www.loom.com/embed/03e20ac00992408aadae33a7b7b3384d" frameborder="0" webkitallowfullscreen mozallowfullscreen allowfullscreen></iframe>

[Resolving a merge conflict using command line]: https://help.github.com/articles/resolving-a-merge-conflict-using-the-command-line/
[Based on: Resolving a merge conflict]: https://www.atlassian.com/git/tutorials/using-branches/merge-conflicts
[pull upstream changes]: ../pull-upstream-changes

