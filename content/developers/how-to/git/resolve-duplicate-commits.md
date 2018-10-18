# How to resolve duplicate commits

__Goal:__ You want to get rid of the duplicate (or triplicate) commits
in your git commit history.

Before you begin [Update your package][] to the existing version on
Bioconductor.


## Steps:

1. **IMPORTANT** Make a backup of the branch with duplicate commits,
   call this `master_backup` (or `RELEASE_3_7_backup`). The name of
   the back up branch should be identifiable and specific to the
   branch you are trying to fix (i.e if you want to fix the *master*
   branch or some `<RELEASE_X_Y>` branch).

   On master, (make sure you are on master by `git checkout master`)

		git branch master_backup


2. Identify the commit from which the duplicates have
   originated. These commits are more often than not, `merge` commits.

3. Reset your branch to the commit *before* the merge commit.

		git reset —hard <commit_id>

4. Now cherry pick your commits from the `master_backup` branch.

		git cherry-pick <commit_id>

	a. The commits you cherry-pick should be only 1 version of the
    duplicate commit, i.e don’t cherry-pick the same commit twice.

	b. Include the branching and version bump commits in your
    cherry-pick. Make the package history look as normal as possible.

5. (Optional) In some cases, there are conflicts you need to fix for
   the cherry-pick to succeed. Please read the documentation on how
   to [resolve conflicts][]

6. Finally, you would need to contact the bioc-devel mailing list to
   reach the Bioconductor core team to sync your repository with the
   version on Bioconductor.  This is not possible as `--force` pushes
   which alter the git timeline are not possible for maintainers.


## How to check if your package has duplicate commits

One way is to check the log. You should see the same commit message
with the same changes, but with a different commit ID, if you try

	git log

or

	git log --oneline

## Script to detect duplicate commits

Run this script written in python to [detect duplicate commits][]
which are specific to Bioconductor repositories.

Usage example:

	python detect_duplicate_commits.py /BiocGenerics 1000

	python detect_duplicate_commits.py <path_to_package> <number of commits to check>'



[detect duplicate commits]: https://github.com/Bioconductor/bioc_git_transition/blob/master/misc/detect_duplicate_commits.py
[Update your package]: http://bioconductor.org/developers/how-to/git/resolve-conflicts
[resolve conflicts]: http://bioconductor.org/developers/how-to/git/resolve-conflicts
