 <a name="top"></a>

# Using the Bioconductor Git-SVN bridge


The Git-SVN Bridge allows Github repositories to be in sync with the 
Bioconductor Subversion (SVN) repository.
Once you have created a bridge, you don't need to use Subversion again
if you don't want to. 
Using a bridge also enables 
social coding features of Github such as issue tracking and pull requests.

## How to create a bridge

In order to create a bridge, you must be the maintainer of a Bioconductor 
package, with read/write access to its directory in the
Bioconductor
[Subversion repository](http://bioconductor.org/developers/how-to/source-control/).

You will also need to create a Github repository which will mirror
the Subversion repository. If you already have a Github repository that has 
files in it, that will work too.

Let's assume that your package is called `MyPackage`, your Subversion
username is `j.user`, your Github username is `username`, and your email 
address is `juser@contributor.org`.

Your package will be in Subversion at the URL

    https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/MyPackage

That's the URL for the devel version of the package. You can also create a bridge to the release version of a package (see the
[FAQ](#bridge-to-release-version)).

<a name="step1"></a>
## Step 1: Configure your Github Repository

### Step 1a: Create Github Repository

If you haven't already created a Github repository, please
[do so now](#create-github). Open the repository page in a web browser;
it will have a URL like:

    https://github.com/yourusername/MyPackage

If you are working with a repository that is part of an organization,
see the [FAQ](#org-repos).


### Step 1b: Add Collaborator

Click on the "Settings" link in the right-side nav bar.
It will look like this:

<img src="settings.jpg"/>

Under **Options** in the left-hand nav bar, click "Collaborators".
At this point, you may be asked to enter your Github password. Do so.
Then in the "Add a Friend" box, type 

    bioc-sync

Then click the **Add** button. This allows the Git-SVN bridge to make changes to your github repository in response to Subversion commits.

<a name="step-1c"></a>
### Step 1c: Add Push Hook

Again in the nav bar at left, click on "Webhooks & Services".
Then click on "Add webhook". You may need to confirm your password 
here.


In the Payload URL box, enter:

    http://gitsvn.bioconductor.org/git-push-hook

**Important Note**: This url must start with `http`, **not** `https`.

Leave "Payload version" alone (it should be 
"application/vnd.github-v3+form")

Under "Which events", choose "Just the `push` event". Make sure
the "Active" box is checked. Then click "Update webhook".

This step lets Bioconductor know when there has been a push to 
your Github repository.


**Important Note**: *Both of the above steps **must** be done 
or your Git-SVN bridge will **not** function properly.*

## Step 2: Create the bridge

Open a browser window pointing to 
the [Git-SVN bridge web application](https://gitsvn.bioconductor.org).

In the bridge web app, click 
"[Log In](https://gitsvn.bioconductor.org/login)".

Log in with your SVN Username, SVN password and email address.
See the [FAQ](#svn-password) if you don't remember either of these.

Once you've logged in, click the 
[Create New Github-SVN mapping](https://gitsvn.bioconductor.org/newproject)
link.

Choose the *root* directory path. If you are creating a bridge
for a software package in Bioconductor's `devel` branch, use
the default value of this dropdown 
(`https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/`).

For *Directory Name*, choose the name of your package, e.g.
`MyPackage`.

In the next box, enter the URL for the Github repository you created
in [step 1](#step1), e.g.

    https://github.com/username/MyPackage

<a name="who-wins"></a>
### Deciding which repository takes precedence

<h1 style="color: red;">DANGER ZONE: READ CAREFULLY</h1>

When initially creating the bridge, the process that
takes place is not a "merge" like you may be accustomed
to from `git` or `svn`. It's really more like an `rsync` 
with the `--delete` option. 

What this means is, you decide who is going to "win", (git or svn)
and if you choose git, then all conflicts will be resolved 
in git's favor. Specifically, the following is what will happen:

* If a file has differing contents in svn and in git, 
  the version that's in git will end up in both git and svn
  (bridge creation is a one-way process; changes are only made
  in one repository based on the contents of the other). 
* If a file has differing contents, it doesn't matter which
  one is the most recent or whether the two could easily
  be merged by git or svn. All that matters is whether
  you selected git or svn as the 'winner'.
* If a file exists in svn but not git, it will be deleted from svn.
* Conversely, if a file exists in git but not svn, that 
  file will be added to the svn repository.

If you picked "svn wins", the above is true, but with git and svn
reversed.


* If your git repository is completely empty (i.e., you haven't
  even added a README file, or any files) and you specified
  "git wins," the bridge will override that and declare svn
  the winner. <span style="color: red;">HOWEVER</span>, if you
  have even a single file in your git repos, and you choose
  "git wins," <span style="color: red;">everything in your svn
  repos will be deleted</span>.
* Of course, nothing is ever deleted in git or svn, so you can
  go back in and retrieve files, but it's best to avoid this
  situation in the first place.



You now need to check two boxes: the first confirms that you have 
configured your Github repository as described in [Step 1](#step1),
the second that you will respond to pull requests and issues
filed in your Github repository (see the [FAQ](#responsibilities)).

You may now click the **Create New Project** button.

You should see a message that your bridge was created successfully.
You can click [My Bridges](https://gitsvn.bioconductor.org/my_bridges)
to confirm this.

## Step 3: What now?

Any commits made to your package in Subversion will be mirrored
in the *master* branch of your Github repository.

Any pushes to the *master* branch of your Github repository will
automatically be mirrored in Subversion, and will propagate
to the [Bioconductor Build System](http://bioconductor.org/checkResults/).

The Git-SVN bridge only affects the *master* branch of your Github
repository. It will ignore changes made in any other branch, even
if those branches are pushed to Github.
So you are free to experiment and even break your package
as long as you don't do it in the *master* branch.

## FAQ

<a name="history"></a>
##### Can I see old commit history?

After creating a bridge, you can't see old svn commit information
from prior to bridge creation if you're using git. (You can still see it 
with svn).

Conversely, in svn, you can't see Git commit messages from
before the bridge was created. You can still see them in git.

Once the bridge is created, you'll see subsequent commit messages
from both git and svn, whether you are using git or svn.

This may change in the future.

[[Back To Top]](#top)

<a name="commit-messages"></a>
##### How do I know whether a commit came from Git or SVN?

If a commit was made in svn, it will show up in the output
of `git log` as something like this:

    commit f0c494108cc854a7a7a267c6a40ea8a3bdef2209
    Author: j.user <j.user@bc3139a8-67e5-0310-9ffc-ced21a209358>
    Date:   Tue Dec 24 19:12:36 2013 +0000

        This is my commit message.
        
        git-svn-id: https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/MyPackage0@85102 bc3139a8-67e5-0310-9ffc-ced21a209358

The `git-svn-id` tells you that this commit originated in svn,
and the number after the `@` is the SVN revision number.

If you are working in Subversion, a commit made in git will look 
something like this when you run `svn log -v`:

    ------------------------------------------------------------------------
    r85104 | j.user | 2013-12-24 11:13:59 -0800 (Tue, 24 Dec 2013) | 12 lines
    Changed paths:
       M /trunk/madman/Rpacks/MyPackage/DESCRIPTION

    Commit made by the Bioconductor Git-SVN bridge.
    Consists of 1 commit(s).

    Commit information:

        Commit id: 6132d20eb3615afdeafcb8a086e952e4b9f8977f
        Commit message:
        Bumped the version number
        Committed by Jill User <juser at contributor.org>
        Commit date: 2013-12-24T11:13:49-08:00

Note that the svn user who did the commit will always be
the user you were logged in as when you created the bridge
in [Step 2](#step2). 

The name of the Git user (denoted by the `Committed by` line)
might vary, if you have granted other users "push" access to
your repository, or if you accept a pull request.



[[Back To Top]](#top)


<a name="other-users-push"></a>
##### Other users have push access to my repository. Will it work for them?

Yes. Be sure this is what you want. If you grant another user push access to your repository, they can push to any branch, including `master`, which will then propagate to the Bioconductor build system. If you don't want the user to have that level of access, then don't grant them push access. You can accept pull requests from them instead.

The Git-SVN bridge will correctly record the name of the git user
who made the commits (see above [FAQ](#commit-messages)).


[[Back To Top]](#top)


<a name="howto-list-bridges"></a>
##### How do I know what Git-SVN bridges exist?

Look at the 
[list of bridges](https://gitsvn.bioconductor.org/list_bridges)
maintained by the web app.

[[Back To Top]](#top)

<a name="advertise"></a>
##### How do I advertise my bridge?

Add the Github URL of your repository to the URL: field 
in the DESCRIPTION file of your package. You can also 
mention your bridge on the bioc-devel 
[mailing list](http://bioconductor.org/help/mailing-list/).

[[Back To Top]](#top)

##### I want to contribute to a package using Github, but there is no bridge for it.

You can request that the maintainer of the package create a bridge,
but if they do not wish to do so, you'll need to contribute
via other means. If a maintainer will not review pull requests and issues
filed via Github, then it is pointless to file them.


<a name="responsibilities"></a>
##### What are my responsibilities when I create a bridge?

As implied by the previous question, package maintainers must
respond to pull requests and issues filed in their Github repositories.

[[Back To Top]](#top)



<a name="svn-password"></a>
##### I don't know my Subversion username and/or password. What do I do?

One of the following steps should work:

* Look in your email. Your SVN credentials were originally sent to you
  by a member of the Bioconductor team (probably Marc Carlson), probably
  with the subject line "congrats" or "congratulations". The email 
  should contain the text "Information about your svn account". 
* Go to your `~/.subversion/auth/svn.simple` directory. There should be
  one or more files whose names are long hexadecimal numbers. Use `grep`
  to find out which file contains your username. If you don't know your 
  username,
  it's usually your first initial, a dot, and your last name (all 
  lowercase). So Jill User would be `j.user`. Example:

        $ grep -l j.user *
        81a52e36a28dfd7750bd975f30c7998b

  This indicates that your password can be found in the file called
  `81a52e36a28dfd7750bd975f30c7998b`. Examine that file and you should see 
  something like:

        password
        V 8
        Z7oRUVH6

  In this case, `Z7oRUVH6` is your password.
* If you still can't find your username or password, contact a 
  member of the Bioconductor team at
  `maintainer at bioconductor dot org`. Mention the package(s) that
  you maintain. We cannot send you your password but we can ask for 
  a new one to be generated, and send it to you. It may take 
  a day or two for the request to be processed.


[[Back To Top]](#top)

<a name="create-github"></a>
##### How do I create a Github repository?

* Go to [https://github.com/](https://github.com/).
* Create a Github account if you don't already have one.
* Click the 
  <a href="https://github.com/new"><img src="newrepo.jpg"></a> icon
  next to your login name at the top right corner of the browser window.
* If you are a member of one or more Github organizations, 
  decide if you want this repository to belong to you 
  as an individual or to one of your organizations. Do this by changing
  the value of the `Owner` dropdown.
* Choose a name for your repository. If possible, it should match the
  name of your package (though it may differ if you have bridges
  for the release and devel version of your package, see
  [the next question](#bridge-to-release-version)).
* You can disregard the rest of the fields and click `Create Repository`.


[[Back To Top]](#top)

<a name="bridge-to-release-version"></a>
##### How do I create a bridge to the release version of my package?

Follow the same instructions as given in [Step 1](#step1), but
give your Github repository a name indicating that it's the release
version, i.e. `MyPackage-release`.

In [Step 2](#step2), be sure you choose the release directory
from the `Root Directory` dropdown when creating the bridge.
For the current release, that would be:

    https://hedgehog.fhcrc.org/bioconductor/branches/RELEASE_<%=config[:release_version].sub(".", "_")%>/madman/Rpacks

Your release bridge is completely separate from your devel bridge.
The Github repositories in each are separate from each other, not
branches of each other. 

Shortly before each Bioconductor release (twice a year, usually
in Spring and Fall), we will disable commits to the release branch,
and your release bridge will stop working. You can delete it.

When the new release branch is created, you can create a new 
release bridge pointing to it.


[[Back To Top]](#top)

<a name="org-repos"></a>
##### Working with a Github Organization repository

If you're working with a Github Organization repository, the steps
to set up a collaborator are a little bit different:

* You need to be an administrator of the Github organization.
* Go to the Organization's main page. If the organization is called 
  `myorg`, you'd go to

        https://github.com/organizations/myorg

* Click on "Teams", at the right-hand side of the top navigation bar.
* Click on "New Team." You may be prompted to enter your Github password.
* Choose a descriptive name for the team, such as `bioc-collab`.
* Grant the "Push & Pull" privileges
* In the Members box, add the user `bioc-sync` and click `Add`.
* In the Repositories box, type the organization name, a slash,
  and the repository name, for example: `myorg/MyPackage`, and
  click `Add`.
* Now click `Save Team`.

Now you can go back to the repository "Settings" page and 
[continue](#step-1c) configuring repository settings.


[[Back To Top]](#top)



##### I don't want to use the bridge, I want to keep my repositories in sync manually.

You can do that by following
[these guidelines](https://github.com/Bioconductor/BiocGithubHelp/wiki/Managing-your-Bioc-code-on-hedgehog-and-github). Thanks to Laurent Gatto for providing these instructions
which are also the inspiration for the Git-SVN bridge.

[[Back To Top]](#top)

<a name="experiment-pkgs"></a>
##### Can I create a bridge to an experiment data package?

No.

##### Things aren't working or I have a question.

Contact the bioc-devel 
[mailing list](http://bioconductor.org/help/mailing-list/).

[[Back To Top]](#top)
