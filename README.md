# Site Maintainer README for bioconductor.org

The canonical location for this code is https://github.com/Bioconductor/bioconductor.org

You can setup by cloning this repository and running

```bash
git clone https://github.com/Bioconductor/bioconductor.org
```

Then after committing code locally run the following to commit the changes
and push the commits back to GitHub.

```bash
# commit code to git
git commit -m "My informative commit message"

# push code to github
git push
```

Unix-ish Developer Required Software

## Required software

NOTE: Before reading the following instructions you may want to consider
installing the web site as a Docker container. See the instructions
below.

1.  Make a fork and clone the git repository for bioconductor.org and
    create a new branch to make your changes, (helpful documentation for
    [Creating a pull request from a fork][])

        git clone https://github.com/<your fork>/bioconductor.org

2.  Make your changes on this branch, add content or edit content.

3.  Once the changes are made, you need use the docker image
    `bioconductor/website:latest` and run the
    container. The container has the dependencies installed to `rake`
    the ruby code and host the website on your local machine at
    https://localhost:3000.

        docker run -v /<full_path>/bioconductor.org:/bioconductor.org/ \
            -p 3000:3000 \
                bioconductor/website:latest

    where,

         -p is mapping the container's port 3000 to the host machine's port

         -v mounting a volume, the website (bioconductor.org) directory
            from your local machine is being mounted on the docker container

4.  Then to kill the process, you need to get the CONTAINER ID with,

        docker ps

    and,

        docker kill <CONTAINER ID>

5.  Before you run the docker image again with more changes, make sure
    to clean the artifacts produced by the `rake` command, with

        git clean -xdf

    The output should look like,

        bioconductor.org$ git clean -xdf
        	Removing assets/bioc-devel-version
        	Removing assets/bioc-version
        	Removing assets/config.yaml
        	Removing content/packages/
        	Removing output/
        	Removing tmp/

6.  Once you have reviewed your changes, make a new branch and send a pull
    request to the `devel` branch. The pull request should be made from your
    `my_changes` branch to the [devel branch on GitHub][].

[devel branch on GitHub]: https://github.com/Bioconductor/bioconductor.org
[Creating a pull request from a fork]: https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request-from-a-fork

### Ruby

The site requires ruby 2.2.2 or newer.

There are numerous issues on various platforms when attempting to install
appropriate versions of ruby and the necessary ruby packages. The simplest
way around all of this is to use rbenv, which allows you to switch
between various ruby versions and avoids conflicts between them.
_NOTE_: rbenv works on Unix only; if you are on Windows, skip to
the Windows section.

On ubuntu, before proceeding, make sure the `libsqlite3-dev` package is
installed (`sudo apt-get install libsqlite3-dev`).

The following instructions are adapted from the
[rbenv page](https://github.com/rbenv/rbenv). It's worth reading this
to understand how rbenv works.

_Important note_: Never use `sudo` when working with a ruby that has been
installed by rbenv. rbenv installs everything in your home directory so
you should never need to become root or fiddle with permissions.

0.  Make sure you do not have rvm installed. `which rvm` should not return
    anything. If you do have it installed, refer to
    [this page](http://stackoverflow.com/questions/3950260/howto-uninstall-rvm)
    for instructions on removing it.

1.  Check out rbenv into `~/.rbenv`.

    ```sh
    $ git clone https://github.com/rbenv/rbenv.git ~/.rbenv
    ```

2.  Add `~/.rbenv/bin` to your `$PATH` for access to the `rbenv`
    command-line utility.

    ```sh
    $ echo 'export PATH="$HOME/.rbenv/bin:$PATH"' >> ~/.bash_profile
    ```

    **Ubuntu Desktop note**: Modify your `~/.bashrc` instead of `~/.bash_profile`.

    **Zsh note**: Modify your `~/.zshrc` file instead of `~/.bash_profile`.

3.  Add `rbenv init` to your shell to enable shims and autocompletion.

    ```sh
    $ echo 'eval "$(rbenv init -)"' >> ~/.bash_profile
    ```

    _Same as in previous step, use `~/.bashrc` on Ubuntu, or `~/.zshrc` for Zsh._

4.  Restart your shell so that PATH changes take effect. (Opening a new
    terminal tab will usually do it.) Now check if rbenv was set up:

    ```sh
    $ type rbenv
    #=> "rbenv is a function"
    ```

5.  Install
    [ruby-build](https://github.com/rbenv/ruby-build),
    which provides the `rbenv install` command that simplifies the process of
    installing new Ruby versions:

        ~~~ sh
        git clone https://github.com/rbenv/ruby-build.git ~/.rbenv/plugins/ruby-build
        ~~~

Now you need to install ruby. Go to the
[Ruby Downloads Page](https://www.ruby-lang.org/en/downloads/)
to find out what the current stable version is. As of 3/06/2020 it is
2.7.0 however that particular version still had issues with modules, so I will
use 2.6.5 the current stable version we use for the website in further examples,
but substitute the current stable version for 2.6.5 if you wish to use or test a
different version.

To install this version of ruby in rbenv, type

    rbenv install 2.6.5

Then, to make this the version of ruby that you will use, type:

    rbenv global 2.6.5

If you want to use different versions of ruby in different contexts, read the
[rbenv page](https://github.com/rbenv/rbenv)
for more information.

## Windows Developer Required Software

1. Download and run the one-click ruby installer
   http://rubyinstaller.org/downloads/. Accept all default settings.

2. Also download and install the Development Kit from
   http://rubyinstaller.org/downloads/. Be sure and add the bin dir
   to your path (see devkitvars.bat)

3. If you don't already have it, be sure and install cygwin and
   explicitly install rsync. rsync is required for parts of the web
   site to work.

4. Install git client. https://git-scm.com/downloads

5. Follow the developer setup instructions below.

## Developer Setup

### Checkout the bioconductor.org codebase and set up upstream remote

    git clone https://github.com/Bioconductor/bioconductor.org
    git remote add upstream git@git.bioconductor.org:admin/bioconductor.org

### Installing Necessary Ruby Packages

Ruby packages are called gems and `gem` is the program used to install them.

To save time, ensure your ~/.gemrc file contains the text

    gem: --no-document

This will ensure that gem does not try to install documentation files
that you will not use.

The web site comes with a Gemfile which is similar to an R package DESCRIPTION
file in that it lists all dependencies needed. Gemfiles are read by
the `bundler` gem, so install that as follows, prepending `sudo` if
necessary (remember, _don't_ use `sudo` if your ruby was installed
with `rbenv`:

    cd bioconductor.org
    gem install bundler

Then, assuming you are in the bioconductor.org working copy, issue this command
to install all dependencies, again prepending `sudo` if necessary:

    bundle install

### Build the site

    cd bioconductor.org # if you aren't already in the working copy
    rake

One step in the build process runs 'nanoc', "a Ruby web publishing system  
for building small to medium-sized websites"; it is one of the
gems you installed above. If you ever need to run nanoc explicitly:

    nanoc compile

To run an abbreviated compile, which does not attempt to build all package pages:

    QUICK_NANOC_COMPILE=true nanoc co

Whether run by hand or by rake, the compiled html files are all found in
and below output/, an immediate subdirectory of the bioconductor.org/ directory
you have been working in.

### Start the built-in development server, 'adsf' "A dead-simple file" server:

    cd output
    adsf

### Test in a browser by going to http://localhost:3000/

### Linters

This project includes the following linters

stylelint:

        npx stylelint "**/*.css"

eslint:

        npx eslint <directory/ file>

## Overview of site source code

- README.md :: You are reading this file or a file generated from
  this file.

- Rakefile :: A `Rakefile` is to `rake` as a `Makefile` is to `make`.
  You can see the available targets by running `rake -T`
  in the directory containing `Rakefile`.

- Rules :: This is a Ruby syntax file that describes how site content
  is transformed from its source form into its output form
  (this is called filtering), what layout to use (layouts are
  the shared templates), and where to write the output (this
  is called routing). See the nanoc
  [tutorial](http://nanoc.stoneship.org/tutorial/) and the
  [manual](http://nanoc.stoneship.org/manual/) for details.

- assets :: This directory is not managed by nanoc. It contains files
  that do not undergo any filtering, layout-ing, or routing.
  Contents of the assets directory are copied to the output
  directory using rsync.

- config.yaml :: Nanoc configuration file for the bioconductor.org
  site. This file is written in [YAML](http://www.yaml.org/).

- content :: This is where the bulk of the raw (source form) site
  content lives. Important details:

             - Content always has two related files: a `.yaml` file
               containing item attributes and a `.<extension>` file
               containing the raw source content. You can actually
               use whatever extension you want.

             - The default behavior is that a content file like
               `help.md` is filtered into HTML and then written to
               `output/help/index.html`. This scheme allows for
               clean URLs that avoid having a file extension.

- layouts :: This is where the content templates live.

- lib :: Ruby helper functions and nanoc extensions live here. Files
  in this directory are automatically loaded by nanoc during
  site processing.

- migration :: Documentation and scripts used in the process of
  migrating the bioconductor.org site from Plone to
  nanoc.

- output :: This directory is created when you compile the
  bioconductor.org site using nanoc. It contains the final
  static HTML and other assets. Deploying the site means
  pushing out an update of the contents of output to the
  live server.

- scripts :: Helper scripts for managing the site live here.

## How to add a page

## How to add event

You will use a helper scripts `./scripts/add_event` to add event
to the site using the following steps:

1.  Always run `./scripts/add_event` from the top-level of your
    website working copy
1.  Run `./scripts/add_event EVENT_NAME`
    This will create an EVENT_NAME.yaml file in the
    `./content/help/events/` directory
1.  The default `EVENT_NAME.yaml` file will look like this:

        title: TITLE FOR EVENT_NAME
        location: Seattle, WA, USA
        event_host: FHCRC
        start: 2010-06-29
        end:   2010-06-29
        link:
            text: details and registration
            url: https://secure.bioconductor.org/EVENT_NAME

1.  Edit the `EVENT_NAME.yaml` file
1.  Use git to commit changes and additions by `add_event`

## How to add course material

You will use a helper script `./scripts/course_mgr` to add course
material to the site. PDF files for labs and presentations as well
as course-specific packages and data are _not_ stored in git. The
index pages that describe the course and provide links to the
materials _are_ stored in git. The `course_mgr` script will help
with index file creation and data transfer.

### `course_mgr` workflow and important tips

To add a course, you will typically perform the following steps
(each described in detail below):

0. Always run `./scripts/course_mgr` from the top-level of your
   website working copy.
1. Run `./scripts/course_mgr --create COURSE_NAME`
2. Run `./scripts/course_mgr --index COURSE_NAME`
3. Build and preview site
4. Run `./scripts/course_mgr --push COURSE_NAME`
5. Use git to commit changes and additions made by `course_mgr`

### Using `course_mgr`

1.  Generate a skeleton course directory structure.

         ./scripts/course_mgr --create seattle-intro

    This will create a `seattle-intro/` directory in the top-level
    of your website working copy -- do not add this directory or any
    files within it to git. Inside will be a `course_config.yaml`
    file that will look like this:

          title:
            The title of the course goes here
          start_date: 2010-01-27
          end_date: 2010-01-29
          instructors: ["Someone", "Another"]
          location: "Seattle, USA"
          url: https://secure.bioconductor.org/SeattleJan10/
          tags: ["intro", "seattle", "package"]
          description:
            You can put some description text here.
            Must be indented.

2.  Put course materials as files and directories into the skeleton
    directory. For example, you might end up with a directory like
    that shown below with two subdirectories, `packages` and
    `presentation-slides`, each containing course materials.

          seattle-intro
          |-- course_config.yaml
          |-- packages
          |   |-- day1_0.0.1.tar.gz
          |   |-- day2_0.0.1.tar.gz
          |   `-- day3_0.0.1.tar.gz
          `-- presentation-slides
              |-- First-steps-presentation.pdf
              |-- Microarray-presentation.pdf
              |-- annotation-presentation.pdf
              `-- sequence-presentation.pdf

3.  Now you are ready to create the index files.

           ./scripts/course_mgr --index seattle-intro
           CREATED: content/help/course-notes/2010/01/seattle-intro.(html|yaml)
           COPIED for preview:
             src: ./seattle-intro/*
             dst: output/help/course-notes/2010/01/seattle-intro/
           NEXT STEPS:
           - preview site with 'rake devserver'
             - Use URL: http://localhost:3000/help/course-materials/2010/seattle-intro/
             - edit CREATED files to add descriptions for links
             - if happy, run ./scripts/course_mgr --push 2010/seattle-intro

    This will create a course index content item in content filed
    appropriately based on the metadata provided in
    `course_config.yaml`. It will also copy the files and directories
    you created into the output directory so that you can do a full
    preview after compiling the site.

4.  If everything looks good, you can sync the data files to the web
    server (note that we do not put these files in git because large
    data files are not appropriate for git and they are not likely to
    change):

            ./scripts/course_mgr --push 2010/seattle-intro
            SYNC:
             src: ./seattle-intro
             dst: biocadmin@staging.bioconductor.org:/loc/www/bioconductor-test.fhcrc.org/help/course-materials/2010/
            NEXT STEPS: git add/commit changes in contents

5.  Finally, "git add" the new course index html and yaml files that were generated in the
    content directory and commit.

### Modifying an existing course

You can edit the pages for an existing course by editing the files in
`./content`. If you need to add or modify data files, run:
./scripts/course_mgr --pull 2010/course_to_modify

    This will create a top-level directory called "course_to_modify". You
    can then add or modify course material. When finished, run
      ./scripts/course_mgr --push 2010/course_to_modify

    If you have changed the .md or .yaml files, do the following:
      cp course_to_modify/course_to_modify.* content/help/course_materials/2010
      git commit -m "made changes" content/help/course-materials/2010/course_to_modify

## Adding course material to the spreadsheet

The page
[http://www.bioconductor.org/help/course-materials/](http://www.bioconductor.org/help/course-materials/)
is built from the tab-delimited file `etc/course_descriptions.tsv`.

Add information to this file using a spreadsheet
program (Excel, LibreOffice, etc.). Be sure to save
in the original tsv format. Note that some spreadsheets
insert non-ASCII characters which cause problems.
Before committing your changes, check for this in R
with:

    tools::showNonASCII(readLines("etc/course_descriptions.tsv"))

And if it reports any non-ascii characters (it will
show line numbers) fix these in a text editor before
committing. Usually the culprit is a non-ascii
hyphen that can be replaced with a regular hyphen.

### Staging site scheduled update

The biocadmin user's crontab on staging.bioconductor.org is used to schedule site
updates every twenty minutes. Below are some details on how the test
site is configured.

The site source is located at
`~biocadmin/bioc-test-web/bioconductor.org`. The `deploy_staging`
Rake task deploys site content to the staging server root on staging.

    task :deploy_staging do
      dst = '/loc/www/bioconductor.org'
      site_config = YAML.load_file("./config.yaml")
      output_dir = site_config["output_dir"]
      system "rsync -gvprt --partial --exclude='.git' #{output_dir}/ #{dst}"
    end

An `update_site` shell script updates from git, builds the site,
and deploys it using Rake.

    #!/bin/bash
    cd bioconductor.org && \
      date && \
      git pull && \
      rake real_clean default deploy_staging deploy_production  && \
      date

We keep track of the output of in a local `log/update_site.log` file and handle
log rotation using `logrotate`. For this we need a config file:

    # logrotate.conf
    /home/biocadmin/bioc-test-web/log/update_site.log {
      daily
      copytruncate
      rotate 30
      dateext
      compress
      missingok
      su biocadmin biocadmin
     }

The following crontab entries are used to schedule site update,
deployment, and log rotation (biocadmin user):

    PATH=/usr/bin:/bin:/usr/sbin
    MAILTO=lori.shepherd@roswellpark.org

    # bioconductor-test.fhcrc.org website publishing
    */20 * * * *  cd $HOME/bioc-test-web;./update_site >> log/update_site.log 2>&1
    59 23 * * * /usr/sbin/logrotate -f -s /home/biocadmin/bioc-test-web/logrotateState /home/biocadmin/bioc-test-web/logrotateFiles

### master.bioconductor.org Apache Configuration

A good resource is available [http://en.opensuse.org/Apache_Quickstart_HOWTO](here).
The staging.bioconductor.org builds and stages the website. The website is
hosted on master.bioconductor.org through apache. This discusses some of the
apache set up.

#### Apache module config

#### Apache vhosts config

Edit /etc/apache2/sites-available/000-default.conf

    <VirtualHost *:80>
      ServerAdmin webmaster@localhost
      ServerName master.bioconductor.org

      # Customized ERROR responses
      ErrorDocument 404 /help/404/index.html
      ErrorDocument 403 /help/403/index.html

      # DocumentRoot: The directory out of which you will serve your
      # documents. By default, all requests are taken from this directory, but
      # symbolic links and aliases may be used to point to other locations.
      DocumentRoot /extra/www/bioc

      # if not specified, the global error log is used
      LogFormat "%h %l %u %t \"%r\" %>s %b \"%{Referer}i\" \"%{User-Agent}i\"" awstats
      ErrorLog /var/log/apache2/bioconductor-error.log
      CustomLog /var/log/apache2/bioconductor-access.log awstats
      ScriptAlias /cgi-bin/ "/usr/local/awstats/wwwroot/cgi-bin/"
      ScriptAlias /cgi/ "/usr/local/cgi/"

      # don't loose time with IP address lookups
      HostnameLookups Off

      # needed for named virtual hosts
      UseCanonicalName Off

      ServerSignature On

        # For most configuration files from conf-available/, which are
        # enabled or disabled at a global level, it is possible to
        # include a line for only one particular virtual host. For example the
        # following line enables the CGI configuration for this host only
        # after it has been globally disabled with "a2disconf".
        #Include conf-available/serve-cgi-bin.conf

      # doc root
      <Directory /extra/www/bioc/>
          Options FollowSymLinks
          AllowOverride All
          #Controls who can get stuff from this server
          #Order allow,deny
          #Allow from all
          Require all granted

         AddOutputFilterByType DEFLATE text/html text/css application/javascript text/x-js
         BrowserMatch ^Mozilla/4 gzip-only-text/html
         BrowserMatch ^Mozilla/4\.0[678] no-gzip
         BrowserMatch \bMSIE !no-gzip !gzip-only-text/html

      </Directory>

      <Directory /extra/www/bioc/checkResults/>
          Options Indexes FollowSymLinks
      </Directory>
      <Directory /extra/www/bioc/packages/submitted/>
          Options Indexes
      </Directory>
      <Directory /extra/www/bioc/packages/misc/>
          Options Indexes
      </Directory>
      <Directory /extra/www/bioc/pending/>
          Options Indexes
      </Directory>
      <Directory /extra/www/bioc/data/>
          Options Indexes
      </Directory>
      <Directory /extra/www/bioc/packages/3.6/bioc/src/contrib/Archive/>
          Options Indexes FollowSymLinks MultiViews
          AllowOverride All
      </Directory>
      <Directory /extra/www/bioc/packages/3.7/bioc/src/contrib/Archive/>
          Options Indexes FollowSymLinks MultiViews
          AllowOverride All
      </Directory>
      <Directory /extra/www/bioc/packages/3.8/bioc/src/contrib/Archive/>
          Options Indexes FollowSymLinks MultiViews
          AllowOverride All
      </Directory>
      <Directory /extra/www/bioc/packages/3.9/bioc/src/contrib/Archive/>
          Options Indexes FollowSymLinks MultiViews
          AllowOverride All
      </Directory>
      <Directory /extra/www/bioc/packages/3.10/bioc/src/contrib/Archive/>
          Options Indexes FollowSymLinks MultiViews
          AllowOverride All
      </Directory>
      <Directory /extra/www/bioc/packages/3.11/bioc/src/contrib/Archive/>
          Options Indexes FollowSymLinks MultiViews
          AllowOverride All
      </Directory>

     # configure the Apache web server to work with Solr
     ProxyRequests Off
     <Proxy *>
       Order deny,allow
       Allow from all
     </Proxy>
     ProxyPass /solr/default/select http://localhost:8983/solr/default/select
     ProxyPreserveHost On
     ProxyStatus On

   </VirtualHost>

## How to test for broken links

You can run wget as shown below to get a report on 404s for the site. Note
that this runs against the staging site so will have a lot of false positives.

    wget -r --spider -U "404 check with wget" -o wwwbioc.log http://master.bioconductor.org

## Optimize redirects

Currently the redirects are defined using Apache's mod_rewrite in a
top-level `.htaccess` file. This has the advantage of allowing easy
revision of the rewrite rules via git that are picked up by Apache on
site update. The downside is that using .htaccess files is suboptimal
in terms of performance. So before the site is launched, consider the
following changes:

1. Copy the directives in the top-level .htaccess file to the site's
   vhost config `/etc/apache2/vhosts.d/bioconductor-test.conf`.

2. Remove the .htaccess file

3. Edit the same vhosts.d config file to set Options to None for the
   top-level directory. This should disable .htaccess files as it
   isn't enough just to remove the .htaccess file itself.

### Site Search

The site search contains several moving parts. The search is built on
Apache Solr, which is in turn built on top of Apache Lucene.

#### How to configure Solr

The default SOLR installation works fine, with the exception of the file
example/solr/conf/schema.xml which must be replaced with the version in
this subversion repository at etc/solr.The changes in this file enable
search query highlighting and snippets.

Solr can be started up as follows (SOLR_HOME is assumed to be the location
where the solr tarball has been expanded):

    cd $SOLR_HOME/example; java -jar start.jar

#### How to ensure that Solr is started up at boot time (on master and staging)

On both machines there is an /etc/rc.local and /etc/init.d/rc.local script
which starts Solr as above.

#### How to configure the Apache web server to work with Solr

Using a2enmod, we added support for the "proxy" and "proxy_http" modules
to the Apache web server. Then we added the following if it hasn't already to
/etc/apache2/sites-available/000-default.conf (master):

    ProxyRequests Off
    <Proxy *>
     Order deny,allow
     Allow from all
    </Proxy>
    ProxyPass /solr/default/select http://localhost:8983/solr/default/select
    ProxyPreserveHost On
    ProxyStatus On

This means that all requests starting with "/solr" will go to the
solr server on port 8983. This allows us to make requests to the
search server without violating the "same-origin" policy.

#### How the client-side portion of the search works

The page /search/index.html includes some javascript (in the file
js/search.html) which in turn uses jQuery. The code parses the
arguments in the URL and then makes an AJAX request to the SOLR
server which returns a JSON string. The javascript code converts
that to an object and then renders the search response page.

#### How to rebuild the search index

Note that you typically do not want to do this by hand as it is handled
by cron jobs (see below).

# NOTE: this may need debugging for staging.bioconductor.org

# from transition from merlot2

On staging.bioconductor.org (ssh to staging.bioconductor.org):

    cd ~/biocadmin/bioc-test-web/bioconductor.org
    rake index_production  (see also rake search_index)

What this command does:

- Runs a Ruby class which determines which files need to be (re)indexed.
- This uses a cache file containing the names of each file and their modification times
  as of the last time the script was run. If the cache file does not exist, all files
  are indexed. This class also handles new files and deletions.
- The class actually does not do the indexing itself; it creates another script
  (index.sh -- created by scripts/search_indexer.rb) which does the actual
  indexing, which is accomplished by using curl to post files to the SOLR web app.

To re-index files on master, ssh to staging (not master) and do this:

    cd ~/biocadmin/bioc-test-web/bioconductor.org
    rake index_production

#### Cron jobs for rebuilding the search index/why it is decoupled from site update

Doing "crontab -l" on staging shows how the index us updated on master. Here are the relevant lines:

    # create search index:
    30 */4 * * * cd $HOME/bioc-test-web/bioconductor.org && rake index_production > $HOME/bioc-test-web/production_index.log 2>&1

Notice that the search indexing process is decoupled from the site building process
(which takes place every 30 minutes). Site indexing can be a time-consuming
process (especially on master) and the site rebuilding should be quick. So
the search indexing takes place once a day on staging at 8 pm on
master (where there are many more files to be indexed which originate from the build system).

### BiocViews Pages

The BiocViews pages are generated by a three-step process:

#### Step 0: Obtain manifest git repository

The `manifest` git repository is available from the Bioconductor git server.
It contains a list of all current Bioconductor packages. Make sure that
it is in the same folder as the `bioconductor.org` repository checkout.
To clone it, first ensure appropriate access rights and then run:

    git clone git@git.bioconductor.org:admin/manifest.git

#### Step 1: rake get_json

This is run by a cron job on staging every day at between 3-8 PM EST (presumably
after the build system has finished and copied all its output to master). Here is the cron job:

    */15 15-20 * * * cd $HOME/bioc-test-web; ./get_json.sh > $HOME/bioc-test-web/get_json.log 2>&1

This Rake target runs some R code which calls code in the BiocViews
package, extracting /packagesdata in JSON format and putting it in
assets/packages/json. Then a ruby script is run which processes that
JSON into a format usable by the javascript tree widget.

If you want to run this target on your own machine, you need R
available with the biocViews (Bioconductor) and rjson (CRAN) packages
installed.

#### Step 2: Build package detail pages

This is done by nanoc and handled by the DataSource subclass BiocViews
(found in lib/data_sources/bioc_views.rb). This data source uses the
JSON files generated in the previous step to build a single page for
each page, one for release and one for devel. The pages are rendered
by the partial layouts/\_bioc_views_package_detail.html.

#### Step 3: The BiocViews Hierarchy page

At http://bioconductor.org/packages. This page uses javascript to
build the tree, reading in data generated in step 1. The relevant
Javascript file is assets/js/bioc_views.js. The automatically
generated (by rake) file output/js/versions.js is also sourced.

### Updating the site during a release

Take a look at the config.yaml file in the root of the
bioconductor.org working copy. This should be the only place you need
to make changes.

# Standard Operating Procedures / SOPs / Troubleshooting

## Problem: Web site does not seem to be updating

Symptom: Commits you made are not going through, and/or
the dashboard (http://bioconductor.org/dashboard/) says
that the site has not been updated in over 20 minutes.
It likely means that an error was introduced in a recent
commit. (make sure you haven't forgotten to `git add` any files).

Solution: ssh to biocadmin@staging.bioconductor.org (ask Lori if
you don't have permission to do so). Change directories:

    cd ~/bioc-test-web

Look at the 2015 Oct 29 10:22:19 AM
then its contents are relevant. You can also look at the last
few lines of `./log/update_site.log`.

# Updating Ruby or Gems

See separate README.updatingRubyOrGems
