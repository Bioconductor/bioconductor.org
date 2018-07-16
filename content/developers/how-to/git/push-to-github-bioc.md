# Push to GitHub and _Bioconductor_ repositories

__Goal:__ During everyday development, you commit changes to your
local repository `master` branch, and wish to push these commits to
both GitHub and _Bioconductor_ repositories.

__NOTE:__ See [Pull upstream changes][] for best practices before
committing local changes.

## Steps:

1. We assume you already have a GitHub repository with the right setup
   to push to _Bioconductor_'s git server
   (git@git.bioconductor.org). If not please see FAQ's on how to get
   access and follow instructions to
   [maintain GitHub and _Bioconductor_ repositories][]. We use a clone
   of the `BiocGenerics` package in the following example.

1.  To check that remotes are set up properly, run the command inside
	your local machine's clone.

		git remote -v

	which should produce the result (where &#60;developer&#62; is your GitHub
	username):

		origin  git@github.com:<developer>/BiocGenerics.git (fetch)
		origin  git@github.com:<developer>/BiocGenerics.git (push)
		upstream    git@git.bioconductor.org:packages/BiocGenerics.git (fetch)
		upstream    git@git.bioconductor.org:packages/BiocGenerics.git (push)

1. Make and commit changes to the `master` branch

		git checkout master
		## edit files, etc.
		git add <name of file changed>
		git commit -m "My informative commit message describing the change"

1. (Alternative) When changes are more elaborate, best practice is to
   use a local branch for development.

		git checkout master
		git checkout -b feature-my-feature
		## multiple rounds of edit, add, commit

	Merge the local branch to master when the feature is 'complete'.

		git checkout master

		# Pull upstream changes before merging
		# http://bioconductor.org/developers/how-to/git/pull-upstream-changes/

		git merge feature-my-feature

1. Push updates to GitHub's (`origin`) `master` branch

		git push origin master

1.  Next, push updates to _Bioconductor_'s (`upstream`) `master` branch

		git push upstream master

1. Confirm changes, e.g., by visiting the GitHub web page for the repository.

[Pull upstream changes]: ../pull-upstream-changes
[maintain GitHub and _Bioconductor_ repositories]: ../maintain-github-bioc
