* Site Maintainer README for bioconductor.org 

** Unix-ish Developer Required Software

*** Required software 

NOTE: Before reading the following instructions you may want to consider 
installing the web site as a Docker container. See the instructions
[here](https://registry.hub.docker.com/u/dtenenba/bioconductor.org/).


**** Ruby

The site requires ruby 1.9.1 or newer.

There are numerous issues on various platforms when attempting to install
appropriate versions of ruby and the necessary ruby packages. The simplest
way around all of this is to use rbenv, which allows you to switch
between various ruby versions and avoids conflicts between them.
*NOTE*: rbenv works on Unix only; if you are on Windows, skip to
the Windows section.

On ubuntu, before proceeding, make sure the `libsqlite3-dev` package is
installed (`sudo apt-get install libsqlite3-dev`).

The following instructions are adapted from the 
[rbenv page](https://github.com/sstephenson/rbenv). It's worth reading this
to understand how rbenv works.

*Important note*: Never use `sudo` when working with a ruby that has been
installed by rbenv. rbenv installs everything in your home directory so
you should never need to become root or fiddle with permissions.

0. Make sure you do not have rvm installed. `which rvm` should not return 
   anything. If you do have it installed, refer to 
   [this page](http://stackoverflow.com/questions/3950260/howto-uninstall-rvm)
   for instructions on removing it.

1. Check out rbenv into `~/.rbenv`.

    ~~~ sh
    $ git clone https://github.com/sstephenson/rbenv.git ~/.rbenv
    ~~~

2. Add `~/.rbenv/bin` to your `$PATH` for access to the `rbenv`
   command-line utility.

    ~~~ sh
    $ echo 'export PATH="$HOME/.rbenv/bin:$PATH"' >> ~/.bash_profile
    ~~~

    **Ubuntu Desktop note**: Modify your `~/.bashrc` instead of `~/.bash_profile`.

    **Zsh note**: Modify your `~/.zshrc` file instead of `~/.bash_profile`.

3. Add `rbenv init` to your shell to enable shims and autocompletion.

    ~~~ sh
    $ echo 'eval "$(rbenv init -)"' >> ~/.bash_profile
    ~~~

    _Same as in previous step, use `~/.bashrc` on Ubuntu, or `~/.zshrc` for Zsh._

4. Restart your shell so that PATH changes take effect. (Opening a new
   terminal tab will usually do it.) Now check if rbenv was set up:

    ~~~ sh
    $ type rbenv
    #=> "rbenv is a function"
    ~~~

5.  Install 
[ruby-build](https://github.com/sstephenson/ruby-build),
which provides the `rbenv install` command that simplifies the process of
installing new Ruby versions:

    ~~~ sh
    git clone https://github.com/sstephenson/ruby-build.git ~/.rbenv/plugins/ruby-build
    ~~~

Now you need to install ruby. Go to the 
[Ruby Downloads Page](https://www.ruby-lang.org/en/downloads/)
to find out what the current stable version is. As of 3/30/2014 it is
2.1.1 so I will use that in further examples, but substitute the current
stable version for 2.1.1 in what follows.

To install this version of ruby in rbenv, type

    rbenv install 2.1.1

Then, to make this the only version of ruby that you will use, type:

    rbenv global 2.1.1

If you want to use different versions of ruby in different contexts, read the
[rbenv page](https://github.com/sstephenson/rbenv)
for more information.



** Windows Developer Required Software

1. Download and run the one-click ruby installer
   http://rubyinstaller.org/downloads/. Accept all default
   settings.

2. Also download and install the Development Kit
   from http://rubyinstaller.org/downloads/.
   Be sure and add the bin dir to your path
   (see devkitvars.bat)

3. If you don't already have it, be sure and install cygwin
   and explicitly install rsync. rsync is required for parts
   of the web site to work.


4. Install subversion client package. Windows packages are listed
   here:

       http://subversion.tigris.org/getting.html#binary-packages

   I installed the collabnet package, but had to go through an annoying
   registration process. The other binaries should be fine and might
   require less hoop jumping. You will need to open a new terminal
   window to pickup the updated config after installing so that you
   will be able to use the svn command.

5. Follow the developer setup instructions below.

** Developer Setup

*** Checkout the bioconductor.org codebase

   svn co https://hedgehog.fhcrc.org/bioconductor/trunk/bioconductor.org

*** Installing Necessary Ruby Packages

Ruby packages are called gems and `gem` is the program used to install them.


To save time, ensure your ~/.gemrc file contains the text

    gem: --no-document

This will ensure that gem does not try to install documentation files
that you will not use.

The web site comes with a Gemfile which is similar to an R package DESCRIPTION
file in that it lists all dependencies needed. Gemfiles are read by
the `bundler` gem, so install that as follows, prepending `sudo` if
necessary (remember, *don't* use `sudo` if your ruby was installed
with `rbenv`:

    cd bioconductor.org
    gem install bundler


Then, assuming you are in the bioconductor.org working copy, issue this command
to install all dependencies, again prepending `sudo` if necessary:

    bundle install


*** Build the site

   cd bioconductor.org # if you aren't already in the working copy
   rake

One step in the build process runs 'nanoc',  "a Ruby web publishing system  
for building small to medium-sized websites"; it is one of the 
gems you installed above.  If you ever need to run nanoc explicitly: 

   nanoc compile

To run an abbreviated compile, which does not attempt to build all package pages:

   QUICK_NANOC_COMPILE=true nanoc co

Whether run by hand or by rake, the compiled html files are all found in 
and below output/, an immediate subdirectory of the bioconductor.org/ directory
you have been working in.  

*** Start the built-in development server, 'adsf' "A dead-simple file" server:

   cd output
   adsf


*** Test in a browser by going to http://localhost:3000/

** Overview of site source code

* README.md :: You are reading this file or a file generated from
                this file.

* Rakefile :: A `Rakefile` is to `rake` as a `Makefile` is to `make`.
              You can see the available targets by running `rake -T`
              in the directory containing `Rakefile`.

* Rules :: This is a Ruby syntax file that describes how site content
           is transformed from its source form into its output form
           (this is called filtering), what layout to use (layouts are
           the shared templates), and where to write the output (this
           is called routing). See the
           [http://nanoc.stoneship.org/tutorial/](nanoc tutorial) and the
           [http://nanoc.stoneship.org/manual/](nanoc manual) for details.

* assets :: This directory is not managed by nanoc. It contains files
            that do not undergo any filtering, layout-ing, or routing.
            Contents of the assets directory are copied to the output
            directory using rsync.

* config.yaml :: Nanoc configuration file for the bioconductor.org
                 site. This file is written in [http://www.yaml.org/](YAML).

* content :: This is where the bulk of the raw (source form) site
             content lives. Important details:

             - Content always has two related files: a `.yaml` file
               containing item attributes and a `.<extension>` file
               containing the raw source content. You can actually
               use whatever extension you want.

             - The default behavior is that a content file like
               `help.md` is filtered into HTML and then written to
               `output/help/index.html`. This scheme allows for
               clean URLs that avoid having a file extension.

* layouts :: This is where the content templates live.

* lib :: Ruby helper functions and nanoc extensions live here. Files
         in this directory are automatically loaded by nanoc during
         site processing.

* migration :: Documentation and scripts used in the process of
               migrating the bioconductor.org site from Plone to
               nanoc.

* output :: This directory is created when you compile the
            bioconductor.org site using nanoc. It contains the final
            static HTML and other assets. Deploying the site means
            pushing out an update of the contents of output to the
            live server.

* scripts :: Helper scripts for managing the site live here.

** How to add a page

** How to add event

You will use a helper scripts `./scripts/add_event` to add event
to the site using the following steps:

0. Always run `./scripts/add_event` from the top-level of your
   website Subversion working copy
1. Run `./scripts/add_event EVENT_NAME`
   This will create an EVENT_NAME.yaml file in the
   `./content/help/events/` directory
2. The default `EVENT_NAME.yaml` file will look like this:

     title: TITLE FOR EVENT_NAME
     location: Seattle, WA, USA
     event_host: FHCRC
     start: 2010-06-29
     end:   2010-06-29
     link:
       text: details and registration
       url: https://secure.bioconductor.org/EVENT_NAME

3. Edit the `EVENT_NAME.yaml` file 
4. Use svn to commit changes and additions by `add_event`
 
** How to add course material

You will use a helper script `./scripts/course_mgr` to add course
material to the site. PDF files for labs and presentations as well
as course-specific packages and data are *not* stored in svn. The
index pages that describe the course and provide links to the
materials *are* stored in svn. The `course_mgr` script will help
with index file creation and data transfer.

*** `course_mgr` workflow and important tips

To add a course, you will typically perform the following steps
(each described in detail below):

0. Always run `./scripts/course_mgr` from the top-level of your
   website Subversion working copy.
1. Run `./scripts/course_mgr --create COURSE_NAME`
2. Run `./scripts/course_mgr --index COURSE_NAME`
3. Build and preview site
4. Run `./scripts/course_mgr --push COURSE_NAME`
5. Use svn to commit changes and additions made by `course_mgr`

*** Using `course_mgr`

1. Generate a skeleton course directory structure.

   ./scripts/course_mgr --create seattle-intro

   This will create a `seattle-intro/` directory in the top-level
   of your website working copy -- do not add this directory or any
   files within it to svn. Inside will be a `course_config.yaml`
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

2. Put course materials as files and directories into the skeleton
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

3. Now you are ready to create the index files.

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

4. If everything looks good, you can sync the data files to the web
   server:

       ./scripts/course_mgr --push 2010/seattle-intro
       SYNC:
        src: ./seattle-intro
        dst: biocadmin@merlot2.fhcrc.org:/loc/www/bioconductor-test.fhcrc.org/help/course-materials/2010/
       NEXT STEPS: svn add/checkin changes in contents

5. Finally, "svn add" the new course index html and yaml files that were generated in the
   content directory and commit.

*** Modifying an existing course

You can edit the pages for an existing course by editing the files in
`./content`. If you need to add or modify data files, run:
       ./scripts/course_mgr --pull 2010/course_to_modify

    This will create a top-level directory called "course_to_modify". You
    can then add or modify course material. When finished, run
      ./scripts/course_mgr --push 2010/course_to_modify

    If you have changed the .md or .yaml files, do the following:
      cp course_to_modify/course_to_modify.* content/help/course_materials/2010
      svn commit -m "made changes" content/help/course-materials/2010/course_to_modify


** http://bioconductor-test.fhcrc.org test site

We run an inside FHCRC only test instance of the Bioconductor website
at the above URL. The site is rebuilt every ten minutes. Here's an
overview of the test site configuration:

* bioconductor-test.fhcrc.org is a DNS CNAME for merlot2.fhcrc.org.

* The site is served by the system installed Apache2 instance on
  merlot2.

* The scheduled svn checkout and rebuild is handled by the biocadmin
  user's crontab.

* biocadmin uses files under ~/bioc-test-web

* Apache serves the site from /loc/www/bioconductor-test.fhcrc.org

*** Staging site scheduled update

The biocadmin user's crontab on merlot2 is used to schedule site
updates every ten minutes. Below are some details on how the test
site is configured.

The site source is located at
`~biocadmin/bioc-test-web/bioconductor.org`. The `deploy_staging`
Rake task deploys site content to the staging server root on merlot2.

    task :deploy_staging do
      dst = '/loc/www/bioconductor.org'
      site_config = YAML.load_file("./config.yaml")
      output_dir = site_config["output_dir"]
      system "rsync -gvprt --partial --exclude='.svn' #{output_dir}/ #{dst}"
    end

An `update_site` shell script updates from svn, builds the site,
and deploys it using Rake.

    #!/bin/bash
    svn update && rake real_clean default deploy_staging

We keep track of the output of in a local `cron.log` file and handle
log rotation using `logrotate`. For this we need a config file:

    # logrotate.conf
    /home/biocadmin/bioc-test-web/cron.log {
      rotate 5
      compress
      daily
    }

The following crontab entries are used to schedule site update,
deployment, and log rotation (biocadmin user):

    PATH=/usr/bin:/bin:/usr/sbin
    MAILTO=devteam-bioc@fhcrc.org

    # bioconductor-test.fhcrc.org website publishing
    ,*/10 * * * *  cd $HOME/bioc-test-web;./update_site >> cron.log 2>&1
    0    0 * * *  logrotate -s $HOME/bioc-test-web/logrotate.state $HOME/bioc-test-web/logrotate.conf

*** Staging site SuSE Apache Configuration

A good resource is available [http://en.opensuse.org/Apache_Quickstart_HOWTO](here).

**** Apache module config

Edit /etc/sysconfig/apache2

Make sure the following modules are listed in the APACHE_MODULES
variable:

* rewrite
* deflate

One way to add them is to do:

   sudo /usr/sbin/a2enmod deflate
   sudo /usr/sbin/a2enmod rewrite

**** Apache vhosts config

Edit /etc/apache2/vhosts.d/bioconductor-test.conf

    <VirtualHost *:80>
      ServerAdmin devteam-bioc@fhcrc.org
      ServerName bioconductor-test.fhcrc.org

      # DocumentRoot: The directory out of which you will serve your
      # documents. By default, all requests are taken from this directory, but
      # symbolic links and aliases may be used to point to other locations.
      DocumentRoot /loc/www/bioconductor-test.fhcrc.org

      # if not specified, the global error log is used
      ErrorLog /var/log/apache2/bioconductor-test.fhcrc.org-error_log
      CustomLog /var/log/apache2/bioconductor-test.fhcrc.org-access_log combined

      # don't loose time with IP address lookups
      HostnameLookups Off

      # needed for named virtual hosts
      UseCanonicalName Off

      ServerSignature On

      # doc root
      <Directory "/loc/www/bioconductor-test.fhcrc.org">
          # The Options directive is both complicated and important. Please see
          # http://httpd.apache.org/docs-2.2/mod/core.html#options
          # for more information.
          Options FollowSymLinks

          # AllowOverride controls what directives may be placed in .htaccess files.
          AllowOverride FileInfo Indexes

          # Controls who can get stuff from this server.
          Order allow,deny
          Allow from all

          # output compression using mod deflate
          AddOutputFilterByType DEFLATE text/html text/css application/javascript text/x-js
          BrowserMatch ^Mozilla/4 gzip-only-text/html
          BrowserMatch ^Mozilla/4\.0[678] no-gzip
          BrowserMatch \bMSIE !no-gzip !gzip-only-text/html
      </Directory>
    </VirtualHost>

*** TODO Add apache2 to rc startup config

*** Start apache using /etc/init.d/apache2 re

*** Staging site nginx installation

We will most likely deploy the test and production sites using
Apache2. A first test setup was configured using nginx. The details
follow.

    ./configure \
      --user=nginx \
      --group=nginx \
      --with-http_ssl_module \
      --with-http_gzip_static_module

    make
    sudo make install

nginx paths:

   path prefix: "/usr/local/nginx"
   binary file: "/usr/local/nginx/sbin/nginx"
   configuration file: "/usr/local/nginx/conf/nginx.conf"
   error log file: "/usr/local/nginx/logs/error.log"
   http access log file: "/usr/local/nginx/logs/access.log"

*** creating an nginx user (SuSE Linux)

     sudo useradd -c "nginx worker" -d /usr/local/nginx -s /bin/false \
                  -g www -G www nginx

*** nginx config

Followed basic config.

    user  nginx www;
    gzip  on;
    gzip_types text/plain text/css text/javascript;

    server {
        listen       80;
        server_name  merlot2.fhcrc.org www.merlot2.fhcrc.org;

        #charset koi8-r;

        #access_log  logs/host.access.log  main;

        location / {
            root   sites/bioconductor.org;
            index  index.html index.htm;
        }

Started nginx as: `sudo /usr/local/nginx/sbin/nginx`

** How to test for broken links

You can run wget as shown below to get a report on 404s for the site. Note
that this runs against the staging site so will have a lot of false positives.

    wget -r --spider -U "404 check with wget" -o wwwbioc.log http://bioconductor-test.fhcrc.org

    *** 404 report for bioconductor-test.fhcrc.org (Tue May 25 09:07:49 2010)

    # TO FIX
    ## depends on packages, etc.
    http://bioconductor-test.fhcrc.org/packages/release/bioc/
    http://bioconductor-test.fhcrc.org/about/publications/compendia/genemetaex/GeneMetaEx_1.0.0.tar.gz
    http://bioconductor-test.fhcrc.org/about/publications/compendia/golubrr/GolubRR_1.3.1.tar.gz
    http://bioconductor-test.fhcrc.org/help/workflows/flowcytometry/flowWorkFlow.pdf
    http://bioconductor-test.fhcrc.org/about/publications/compendia/CompStatViz/CompStatViz_2.0.1.zip
    http://bioconductor-test.fhcrc.org/about/publications/compendia/golubrr/GolubRR_1.3.1.zip
    http://bioconductor-test.fhcrc.org/help/bioconductor-packages/
    http://bioconductor-test.fhcrc.org/about/publications/compendia/CompStatViz/CompStatViz_2.0.1.tar.gz
    http://bioconductor-test.fhcrc.org/help/docs/papers/2003/Compendium/golubEsets_1.0.tar.gz
    http://bioconductor-test.fhcrc.org/help/workflows/flowcytometry/tutorial.mpeg
    http://bioconductor-test.fhcrc.org/help/workflows/flowcytometry/dataFiles.tar

** Note on launching the new site

*** Discuss production setup with Dirk

**** DNS

You want to set things up so that you can move to the new site or
revert to current quickly. Dirk should be able to suggest a way to
achieve this. Ideally, you would not change the bioconductor.org DNS
record as this can take awhile to propagate and doesn't give a quick
way to revert.

**** Site monitoring and alerting

I imagine PHS IT has some monitoring that can be put in place for the
new site. Would also make sense to add an external monitor so that
you will know if the site becomes unreachable from the outside.

**** Squid

I'm not sure what the current status is w.r.t. to Squid proxy/cache.
With the new setup, I would anticipate that a reasonable web server
running Apache will be enough for the load and that if more throughput
or redundancy is desired, setting up a second server and load
balancing would be a good next step.

*** Optimize redirects

Currently the redirects are defined using Apache's mod_rewrite in a
top-level `.htaccess` file. This has the advantage of allowing easy
revision of the rewrite rules via svn that are picked up by Apache on
site update. The downside is that using .htaccess files is suboptimal
in terms of performance. So before the site is launched, consider the
following changes:

1. Copy the directives in the top-level .htaccess file to the site's
   vhost config `/etc/apache2/vhosts.d/bioconductor-test.conf`.

2. Remove the .htaccess file

3. Edit the same vhosts.d config file to set Options to None for the
   top-level directory. This should disable .htaccess files as it
   isn't enough just to remove the .htaccess file itself.

*** Testing the staging site

**** Use a few days worth of access logs

Extract paths from a few days of access logs (make sure to filter for
200 responses) and "replay" these against the staging site. This should
give a good idea of whether the redirects are doing enough and whether
or not basic repository structure has been appropriately mirrored.

**** Use wget to test for broken links on the site

    wget -r -l 20 --spider -U "404 check with wget" -o wwwbioc.log http://bioconductor-test.fhcrc.org

**** Work with Dirk to make the staging site available on the internet

Then you can ask Wolfgang to do some tests to see how the site
performs from Europe. You could also run some site performance
analysis tools like YSlow to get some suggestions for improvements.

**** Staging site performance

You might look into running some simple benchmarks with `ab` or
`httperf`. Might be interesting to compare against the current
Plone-based site.

*** Misc Concerns

In trying to get some test data from the current site using wget,
I've seen a number of cases where a wget request failed, but then
works when I try in a browser. This makes me worried that the
wget-based snapshot may not be as complete as thought. Not sure if
the issue is wget config options or Plone getting overwhelmed with
requests and failing to respond.

*** Site Search

The site search contains several moving parts. The search is built on 
Apache Solr, which is in turn built on top of Apache Lucene. 

**** How to configure Solr
The default SOLR installation works fine, with the exception of the file
example/solr/conf/schema.xml which must be replaced with the version in
this subversion repository at etc/solr.The changes in this file enable
search query highlighting and snippets.

Solr can be started up as follows (SOLR_HOME is assumed to be the location
where the solr tarball has been expanded):

    cd $SOLR_HOME/example; java -jar start.jar


**** How to ensure that Solr is started up at boot time (on merlot2 and krait)
On both machines there is an /etc/rc.d/rc.local script (with symlink at
/etc/rc.d/rclocal) which starts Solr as above. TODO: reboot (at least merlot2)
and make sure this works.

**** How to configure the Apache web server to work with Solr

Using a2enmod, we added support for the "proxy" and "proxy_http" modules
to the Apache web server. Then we added the following to 
/etc/apache2/vhosts.d/bioconductor-test.conf (merlot2) or
bioconductor.conf (krait):

    ProxyRequests Off
    <Proxy *>
     Order deny,allow
     Allow from all
    </Proxy>
    ProxyPass /solr http://localhost:8983/solr
    ProxyPreserveHost On
    ProxyStatus On

This means that all requests starting with "/solr" will go to the
solr server on port 8983. This allows us to make requests to the 
search server without violating the "same-origin" policy.

**** How the client-side portion of the search works

The page /search/index.html includes some javascript (in the file
js/search.html) which in turn uses jQuery. The code parses the 
arguments in the URL and then makes an AJAX request to the SOLR
server which returns a JSON string. The javascript code converts
that to an object and then renders the search response page.

**** How to rebuild the search index (on your own machine, merlot2, or krait)

Note that you typically do not want to do this by hand as it is handled
by cron jobs (see below). 

On merlot2 (ssh to merlot2):

    cd ~/biocadmin/bioc-test-web/bioconductor.org
    rake search_index

What this command does:

* Runs a Ruby class which determines which files need to be (re)indexed.
* This uses a cache file containing the names of each file and their modification times
  as of the last time the script was run. If the cache file does not exist, all files 
  are indexed. This class also handles new files and deletions.
* The class actually does not do the indexing itself; it creates another script
  (index.sh) which does the actual indexing, which is accomplished by using
  curl to post files to the SOLR web app.

To re-index files on krait, ssh to merlot2 (not krait) and do this:
 
    cd ~/biocadmin/bioc-test-web/bioconductor.org
    rake index_production


**** Cron jobs for rebuilding the search index/why it is decoupled from site update

Doing "crontab l" on merlot2 shows how the index us updated 
on both merlot2 and krait. Here are the relevant lines:

    30 */1 * * * cd $HOME/bioc-test-web; ./index_staging.sh > $HOME/bioc-test-web/index_staging.log 2>&1
    30 */4 * * * cd $HOME/bioc-test-web; rake index_production > $HOME/bioc-test-web/production_index.log 2>&1

Notice that the search indexing process is decoupled from the site building process
(which takes place every 10 minutes). Site indexing can be a time-consuming 
process (especially on krait) and the site rebuilding should be quick. So 
the search indexing takes place every hour on merlot2 and every four hours on
krait (where there are many more files to be indexed which originate from the build system).


**** How to get search working on your own development machine

You could set up apache as described above but I think that is overkill.
I use pound (http://www.apsis.ch/pound/) as a simple front end to both
adsf (serving static content built by nanoc on one port) and solr
(java web app running on another). You can use "rake search_index"
to build the search index. You need to define the shell variables
SOLR_HOME and JAVA_HOME. The rake target may require slight modification
to handle the hostname of your local machine.


--todo: make sure people can't do anything bad as solr admin (change password?)
See http://wiki.apache.org/solr/SolrSecurity

*** BiocViews Pages
The BiocViews pages are generated by a three-step process:

**** Step 1: rake get_json

This is run by a cron job on merlot2 every day at 2PM (presumably after the build system
has finished and copied all its output to krait). Here is the cron job:

    0 14 * * * cd $HOME/bioc-test-web; rake get_json > $HOME/bioc-test-web/get_json.log 2>&1

This Rake target runs some R code which calls code in the BiocViews package, extracting 
/packagesdata in JSON format and putting it in assets/packages/json. Then a ruby script
is run which processes that JSON into a format usable by the javascript tree widget.

If you want to run this target on your own machine, you need R available with the biocViews
(Bioconductor) and rjson (CRAN) packages installed.

**** Step 2: Build package detail pages

This is done by nanoc and handled by the DataSource subclass BiocViews (found in
lib/data_sources/bioc_views.rb). This data source uses the JSON files generated in the
previous step to build a single page for each page, one for release and one for devel.
The pages are rendered by the partial layouts/_bioc_views_package_detail.html.

**** Step 3: The BiocViews Hierarchy page

At http://bioconductor.org/packages.
This page uses javascript to build the tree, reading in data generated in step 1.
The relevant Javascript file is assets/js/bioc_views.js. The automatically generated 
(by rake) file output/js/versions.js is also sourced.

*** Updating the site during a release

Take a look at the config.yaml file in the root of the bioconductor.org working copy.
This should be the only place you need to make changes.
