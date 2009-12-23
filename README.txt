Site Maintainer README for bioconductor.org 
===========================================

Developer Setup
---------------

First you need to have a recent version of Ruby installed:

     ruby --version
     ruby 1.8.7 (2008-08-11 patchlevel 72) [x86_64-linux]

Next, you need to have a recent version of rubygems installed:

     gem --version
     1.3.1

If gem is not available and you are on Linux, try:

    sudo zypper install rubygems

Finally, you need to install the following rubygems:

    sudo gem install nanoc3 BlueCloth rack mime-types haml json \
                     rack-cache httparty rake

nginx installation
------------------


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

creating an nginx user (SuSE Linux)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sudo useradd -c "nginx worker" -d /usr/local/nginx -s /bin/false \
             -g www -G www nginx

nginx config
~~~~~~~~~~~~

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

bioconductor.org site code
--------------------------

    svn co  \
    https://hedgehog.fhcrc.org/bioconductor/trunk/bioconductor.org/

scheduled update
----------------

Created a rake task for local deployment:

    task :deploy_merlot2_local do
      dst = '/usr/local/nginx/sites/bioconductor.org'
      site_config = YAML.load_file("./config.yaml")
      output_dir = site_config["output_dir"]
      system "rsync -gvprt --partial --exclude='.svn' #{output_dir}/ #{dst}"
    end

Added a script to an svn checkout working copy:

    #!/bin/bash
    
    svn update && rake default deploy_merlot2_local

Added chrontab entry:

    PATH=/usr/bin:/bin
    MAILTO=sfalcon@fhcrc.org
    
    # test for bioconductor.org nanoc site gen
    */5 * * * *  cd $HOME/src/SVN/bioconductor.org;./update_site >> cron.log 2>&1

Started nginx as: sudo /usr/local/nginx/sbin/nginx

