# Frequently asked questions

1. I can't access my package.

   You will need to log in to the [BiocCredentials][] app. If you have
   not logged in before, you must first **activate** your account.

   There are two steps,

   1. If there is no SSH key registered, you must add one.

   1. If there is already an SSH key registered, check the packages you
      have access to in the 'Profile' interface.

	You can alternatively check if you have access to your package
	using the command line

		  ssh -Tv git@git.bioconductor.org

	If you have access to your package, but cannot git `pull` or `push`,
	please check FAQ #13, #14, and #15.

1. I'm a developer for _Bioconductor_, my package `ExamplePackage` is on
   the new server https://git.bioconductor.org. What do I do next?

    Take a look at [Maintain GitHub and _Bioconductor_ repositories][].
	This will give you the information needed.

	NOTE: This situation is for packages which were previously
	maintained on SVN and have never been accessed through GIT. It is
	not for newly accepted packages through Github.

1. I have a GitHub repository already set up for my _Bioconductor_
   package at `www.github.com/<developer>/<ExamplePackage>` , how do I
   link my repository in GitHub and https://git.bioconductor.org ?

    Take a look at [New package workflow][]. Step 2 gives you
    information on how to add the remote and link both GitHub and
    _Bioconductor_ repositories.

1. I'm unable to `push` or `merge` my updates from my GitHub
   repository to my _Bioconductor_ package on
   `git@git.bioconductor.org`, how do I go about this?

    If you are unable to `push` or `merge` to either your GitHub
    account or _Bioconductor_ repository, it means you do not have the
    correct access rights. If you are a developer for _Bioconductor_,
    you will need to [submit your SSH public key][] to the
    [BiocCredentials][] app.

    You should also make sure to check that your public key is set up
    correctly on GitHub. Follow
    [Adding an SSH key to your GitHub account][].

1. I'm not sure how to fetch the updates from `git.bioconductor.org`
   with regards to my package, how do I do this?

    Take a look at [Sync existing repositories][]. This will give you
    the information needed.

1. I'm just a package user, do I need to do any of this?

    As a package user, you do not need any of these __developer__
    related documentation. Although, it is a good primer if you want
    to be a contributor to _Bioconductor_.

    You can also open [pull requests][] and [issues][] on the
    _Bioconductor_ packages you use, __if__ they have a GitHub
    repository.

1. I'm new to git and GitHub, where should I learn?

    There are many resources where you can learn about git and GitHub.

    * [git-and-github-learning-resources][]
    * [git-scm][]
    * [Guides][]

1. I'm a _Bioconductor_ package maintainer, but I don't have access to
   the _Bioconductor_ server where my packages are being maintained. How
   do I gain access?

    Please submit your SSH public key using the [BiocCredentials][]
    app. Your key will be added to your our server and you will get
    read+ write access to your package.

    All developers of _Bioconductor_ packages are required to do this,
    if they don't already have access. Please identify which packages
    you need read/write access to in the email.

1. What is the relationship between the `origin` and `upstream` remote?

    In `git` lingo __origin__ is just the default name for a remote
    from which a repository was originally cloned. It might equally
    have been called by another name. We recommend that __origin__ be
    set to the developers GitHub repository.

    Similarly, __upstream__ is the name for a remote which is hosted
    on the _Bioconductor_ server.

    It is important that all the changes/updates you have on your
    __origin__ are equal to __upstream__, in other words, you want
    these two remotes to be in sync.

    Follow [Sync existing repositories][] for details on how to
    achieve this.

    Image explaining GitHub and _Bioconductor_ relationship for a
    developer

    ![](/images/git/github-bioc-relationship.png)

1. Can I have more than one upstream remote, if yes, is this recommended?

    You can have as many remotes as you please. But you can have only
    one remote with the name __upstream__. We recommend having the
    remote `origin` set to GitHub, and `upstream` set to the
    _Bioconductor_ git server to avoid confusion.

1. Common names used in the scenario's

    `developer`: This should be your GitHub username, e.g., mine is
    `nturaga`.

    `BiocGenerics`: This is being used as an example to demonstrate
    git commands.

    `ExamplePackage`: This is being used a place holder for a package
    name.

    SVN `trunk` and git `master` branch are now the development
    branches.

1. I'm a _Bioconductor_ developer only on the _Bioconductor_ server. I do
    not have/want a GitHub account. What should I do?

    You do not have to get a Github account if you do not want
    one. But it is a really good idea, to maintain your package
    publicly and interact with the community via the social coding
    features available in Github.

    We highlight this in [Maintain a _Bioconductor_-only repository][]

1. I cannot push to my package. I get the error,

        $ git push -v origin master
        fatal: remote error: FATAL: W any packages/myPackage nobody DENIED by fallthru
        (or you mis-spelled the reponame)

    (you might have renamed the `origin` remote as `upstream`;
    substitute `upstream` for `origin`. Check your remote,

        $ git remote -v
        origin  https://git.bioconductor.org/packages/myPackage.git (fetch)
        origin  https://git.bioconductor.org/packages/myPackage.git (push)

    As a developer you should be using the SSH protocol, but the
    `origin` remote is HTTPS. Use

        git remote add origin git@git.bioconductor.org:packages/myPackage

    to change the remote to the SSH protocol. Note the `:` after the
    host name in the SSH protocol, rather than the `/` in the HTTPS
    protocol. Confirm that the remote has been updated correctly with
    `git remote -v`.

    If your remote is correct and you still see the message, then your
    SSH key is invalid. See the next FAQ.

1. Before sending a question to the Bioc-devel mailing list about
    git, please check the output of the following commands for
    correctness so that we can help you better.

    * As a developer check to make sure, you are using SSH as your
      access protocol.  Check the output of `git remote -v` for
      consistency. Include this in your email to bioc-devel, if you
      are unsure. The remote should look like,

            origin  git@git.bioconductor.org:packages/myPackage.git (fetch)
            origin  git@git.bioconductor.org:packages/myPackage.git (push)

       or

            origin  git@github.com:<github username>/myPackage.git (fetch)
            origin  git@github.com:<github username>/myPackage.git (push)
            upstream  git@git.bioconductor.org:packages/myPackage.git (fetch)
            upstream  git@git.bioconductor.org:packages/myPackage.git (push)

    * Check if you have access to the bioc-git server
      (git@git.bioconductor.org), by using `ssh -T
      git@git.bioconductor.org`.  This will show you a list of
      packages with READ(R) and WRITE(W) permissions.  As a developer
      you should have `R W` next to your package. This is based on the
      SSH public key you are using, the default for ssh authentication
      is `id_rsa`.

            R    	admin/..*
            R    	packages/..*
            R  	admin/manifest
            R  	packages/ABAData
            R  	packages/ABAEnrichment
            R  	packages/ABSSeq
            R W  	packages/ABarray
            R  	packages/ACME
            R  	packages/ADaCGH2
            R  	packages/AGDEX

1. SSH key not being recognized because of different name?

    If you have named your SSH public key differently from `id_rsa` as
    suggested by `ssh-keygen`, you may find it useful to set up a
    `~/.ssh/config` file on your machine.  Simply make a
    `~/.ssh/config` file if it does not exist, and add,

        host git.bioconductor.org
            HostName git.bioconductor.org
            IdentityFile ~/.ssh/id_rsa_bioconductor
            User git

    In this example, my private key is called `id_rsa_bioconductor`
    instead of `id_rsa`.

	You may find it useful to check the [BiocCredentials][] app to see
    what SSH key you have registered.

1. SSH key asking for a password and I don't know it? How do I
    retrieve it?

    There are a few possibilities here,

    * You have set a password. The bioc-devel mailing list cannot help
    you with this.  You have to submit a new key on the
    [BiocCredentials][] app.

    * The permissions on your SSH key are wrong. Verify that the
    permissions on SSH IdentityFile are `400`. SSH will reject, in a
    not clearly explicit manner, SSH keys that are too readable. It
    will just look like a credential rejection.  The solution, in this
    case, is (if your SSH key for bioconductor is called `id_rsa`):

            chmod 400 ~/.ssh/id_rsa

    * You have the wrong remote set up, please check `git remote -v`
    to make sure the SSH access protocol is being used. Your bioc-git
    server remote, should be
    `git@git.bioconductor.org:packages/myPackage`.

1. Can I create and push new branches to my repository on git.bioconductor.org?

    No. Maintainers only have access to `master` and the current `RELEASE_X_Y`.
    New branches cannot be created and pushed to the bioconductor server. We
    recommend maintainers have additional branches on their Github repository
    if they are maintaining one.

1. How can I fix my duplicate commits issue and find the required
    documentation?

	The detailed documentation to [resolve duplicate commits][]
	can be found at the link.

## More questions?


If you have additional questions which are not answered here already,
please send an email to bioc-devel@bioconductor.org.

## Helpful links:


_Bioconductor_ [source control][] overview.

[Adding a new SSH key to your GitHub account](https://help.github.com/articles/adding-a-new-ssh-key-to-your-github-account/)

[Create a pull request on GitHub](https://help.github.com/articles/creating-a-pull-request/)

[Create an issue on GitHub](https://help.github.com/articles/creating-an-issue/)

[Git and GitHub learning resources](https://help.github.com/articles/git-and-github-learning-resources/)

[git-scm manual](https://git-scm.com/)

[GitHub Guides](https://guides.github.com/)

[source control]: ../

[New package workflow]: ../new-package-workflow

[Maintain GitHub and _Bioconductor_ repositories]: ../maintain-github-bioc

[Maintain a _Bioconductor_-only repository]: ../maintain-bioc-only

[Sync existing repositories]: ../sync-existing-repositories

[Adding an SSH key to your GitHub account]: https://help.github.com/articles/adding-a-new-ssh-key-to-your-github-account/

[submit your SSH public key]: https://git.bioconductor.org/BiocCredentials/

[resolve duplicate commits]: ../resolve-duplicate-commits

[Pull requests]: https://help.github.com/articles/creating-a-pull-request/

[issues]: https://help.github.com/articles/creating-an-issue/

[git-and-github-learning-resources]: https://help.github.com/articles/git-and-github-learning-resources/

[git-scm]: https://git-scm.com/

[Guides]: https://guides.github.com/

[BiocCredentials]: https://git.bioconductor.org/BiocCredentials
