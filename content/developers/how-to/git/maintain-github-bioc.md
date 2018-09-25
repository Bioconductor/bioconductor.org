# Create a new GitHub repository for an existing _Bioconductor_ package

__Goal:__ As a maintainer, you'd like to create a new GitHub
repository for your existing _Bioconductor_ repository, so that your
user community can engage in the development of your package.

## Steps:

1.  [Create a new GitHub account][] if you don't have one.

1.  Set up remote access to GitHub via SSH or Https.  Please check
    [which-remote-url-should-i-use][] and
    [add your public key to your GitHub account][].

1.  Once you have submitted your keys, you can login to the
    [BiocCredentials application][BiocCredentials] to check if the correct
    keys are on file with _Bioconductor_.

1.  [Create a new GitHub repository][] on your account, with the name
    of the existing _Bioconductor_ package.

    We use "BiocGenerics" as an example for this scenario.

    ![](/images/git/create_repo.png)

    After pressing the 'Create repository' button, __ignore__ the
    instructions that GitHub provides, and follow the rest of this
    document.

1.  On your local machine, clone the empty repository from GitHub.

    Use `https` URL (replace `<developer>` with your GitHub username)

        git clone https://github.com/<developer>/BiocGenerics.git

    or `SSH` URL

        git clone git@github.com:<developer>/BiocGenerics.git

1.  Add a remote to your cloned repository.

    Change the current working directory to your local repository
    cloned in the previous step.

        cd BiocGenerics
        git remote add upstream git@git.bioconductor.org:packages/BiocGenerics.git

1.  Fetch content from remote upstream,

        git fetch upstream

1.  Merge upstream with origin's master branch,

        git merge upstream/master

    __NOTE:__ If you have the error `fatal: refusing to merge
    unrelated histories`, then the repository cloned in step 4 was not
    empty. Either clone an empty repository, or see
    [Sync existing repositories][].

1. Push changes to your origin master,

        git push origin master

    __NOTE:__ Run the command `git config --global push.default
    matching` to always push local branches to the remote branch of
    the same name, allowing use of `git push origin` rather than `git
    push origin master`.

1.  (Optional) Add a branch to GitHub,

        ## Fetch all updates
        git fetch upstream

        ## Checkout new branch RELEASE_3_6, from upstream/RELEASE_3_6
        git checkout -b RELEASE_3_6 upstream/RELEASE_3_6

        ## Push updates to remote origin's new branch RELEASE_3_6
        git push -u origin RELEASE_3_6

1. Check your GitHub repository to confirm that the `master` (and
   optionally `RELEASE_3_6`) branches are present.

1. Once the GitHub repository is established follow
   [Push to GitHub and _Bioconductor_][] to maintain your repository
   on both GitHub and Bioconductor.

[submit-keys]: https://git.bioconductor.org/BiocCredentials/
[BiocCredentials]: https://git.bioconductor.org/BiocCredentials/
[Create a new GitHub account]: https://help.github.com/articles/signing-up-for-a-new-github-account/
[Create a new GitHub repository]: https://help.github.com/articles/create-a-repo/
[Sync existing repositories]: ../sync-existing-repositories
[which-remote-url-should-i-use]: https://help.github.com/articles/which-remote-url-should-i-use/
[add your public key to your GitHub account]: https://help.github.com/articles/connecting-to-github-with-ssh/
[Push to GitHub and _Bioconductor_]: ../push-to-github-bioc
