# Updating Ruby and Gems on Staging!

This can be no trival task!  All ids and dates are accurate as of 02/20/2018.
This assumes you have access to the Bioconductor EC2 instances/snapshots and 
have a locally checkout version of the website code

## TESTING 
It is suggested to make an EC2 testing environment where a test [git] branch can be 
created. Once that is debugged the testing [git] branch can be pulled locally to 
again test and debug if need be. Then it can be deployed. 

### Make a testing EC2 instance 

Make a testing EC2 instance off the latest snapshot of `staging.bioconductor.org`. 
The snapshot is taken every Friday. 

    1. Go to our Bioconductor EC2 instances: https://console.aws.amazon.com/ec2/v2/home?region=us-east-1#Instances:sort=desc:launchTime
    2. Get the `staging.biconductor.org` instance id (i-6552ca8c)
    3. Search for that instance id in the list of create Snapshots: https://console.aws.amazon.com/ec2/v2/home?region=us-east-1#Snapshots:sort=desc:startTime
    4. Choose the most recent and launch
  
This launched EC2 instance is now our testing environment! This instance will have
the same account users and passwords as the original `staging.bioconductor.org`

####  **Immediately:** Disable the crontab for ubuntu and biocadmin users!!!

This disables any tasks or interaction with the main EC2 instances. Login as each
user and disable the crontab immediately!!  Then log on as `biocadmin` for all
updates.

#### Update Ruby

Log in as `biocadmin` to your testing instance. 
The staging environment uses rbenv to install and update ruby. See either https://github.com/Bioconductor/bioconductor.org#ruby or https://github.com/rbenv/rbenv

    1. cd .rbenv
    2. git pull
    3. cd plugins/ruby-build
    4. git pull
    5. cd /home/biocadmin/bioc-test-web/bioconductor.org
    6. rbenv install -l   # this will list all available versions
    7. rbenv install <ruby version>   (ex: rbenv install 2.5.0)
    8. rbenv global <ruby version>  # this will set the default ruby version
    
I ran into ruby installation ERRORS such as the following: 
`rbenv: cannot rehash: /home/biocadmin/.rbenv/shims/.rbenv-shim exists`

If this occurs simply delete this file AND the partially failed download of the 
version and try to download desired version again. 

EX. If this occurred during the 2.5.0 installation, run the following

```
rm /home/biocadmin/.rbenv/shims/.rbenv-shim
rm -rf /home/biocadmin/.rbenv/versions/2.5.0 
rbenv install 2.5.0
rbenv global 2.5.0
```
#### Create a new github branch

Updating gems at the very least will create an updated Gemfile.lock in the
directory `/home/biocadmin/bioc-test-web/bioconductor.org`. Website files or code 
may also need to be updated with newer versions of required gems or debugged if
gems changed internal functionality or object definitions.  It is recommended to start
a new branch early to keep track of all changes and merge once debugged.

`git checkout -b updatingRubyTesting`

#### Update Gems

On the testing EC2 instance, make sure you are using the latest just installed 
version of ruby, remove the current Gemfile.lock else the old 
versions of gems will be installed.

```
ruby --version
rm Gemfile.lock
gem install nanoc
gem install bundler
bundle install
```

#### Rake and debug

You now should be able to rake with the updated version of ruby and gems. Running
rake in `/home/biocadmin/bioc-test-web/bioconductor.org` will start to build the
website locally - as our code is defined in a directory `output`
You can double check the version and installed gem list with the following:

```
ruby --version
gem list
```

**NOTE** 
There is a website script that is used that accesses support.bioconductor.org.
This will fail unless the IP is added to the inbound rules of the EC2 instance 
for support.bioconductor. Since the Testing EC2 instance will have a different 
IP each time in is stopped and restarted it is recommended to disable this 
script in the Rake file for testing and debug locally if necessary. You may 
need your local IP added to the inbound rules; talk to whomever necessary that is
managing these resources for this change. 

``` 
# In Rakefile comment out the line at the top of the file
# require './scripts/get_post_tag_info.rb'
```

For good measure you can delete the output/ directory to make sure it is clean
before doing the rake.

```
cd output
rm -rf *
cd ..
rake
```

Debug and update if necessary.  If it is taking longer than expected (a few days, 
a week or two) you should try and update the gems to make sure all are still current
by doing a gem update. **NOTE** You still MUST delete the Gemfile.lock if you 
want the gems to be updated
```
rm Gemfile.lock
gem update
```

#### Push up git branch

Push your local git testing branch so that you can pull down all the changes 
locally and retest. 

`git push -u origin updatingRubyTesting`

### Testing locally / deploy locally
This also assumes you have a locally checked out version of the website code found at:  https://github.com/Bioconductor/bioconductor.org
Follow the same steps above for updating ruby locally and git pull the testing
branch that was uploaded from the EC2 testing.  You should pull down a new 
Gemfile.lock that will have the updated gem versions. For good measure 
update or install nanoc and bundler 

```
cd bioconductor.org  # your local checkout of website
git checkout -b updateRubyTesting origin/updateRubyTesting   # pull branch
gem install nanoc
gem install bundler
bundle install
```
You now should be able to rake with the updated version of ruby and gems. 
You can double check the version and installed gem list with the following as
we did before:

```
ruby --version
gem list
```

For good measure you can delete the output/ directory to make sure it is clean
before doing the rake.

```
cd output
rm -rf *
cd ..
rake
```

### Merge the git branches and test locally again 

This may seem like overkill but feel like it is a good precautionary step. Once
The testing branch is debugged locally, Merge the local testing branch into master. 

```
git checkout master
git merge updatingRubyTesting
```
Make sure the Gemfile.lock was updated to the most recent versions of everything 
and make sure you can rake again without ERROR. Debug if necessary. If this goes
as planned, push changes to github and git.bioconductor.org repository. 

**Note** Be sure you are ready to deploy if the changes are pushed to the 
`git.bioconductor.org` location.  The repo will start being pulled immediately on
the next scheduled build of the website - This will mean there could be a period
of time where the website is broken if it pulls the changes but the gems and ruby 
haven't been updated yet!


## DEPLOY

Log into the `staging.bioconductor.org` EC2 instance as `biocadmin`. Follow the 
same instructions above for updating Ruby. Navigate to the website code
`/home/biocadmin/bioc-test-web/bioconductor.org` and do a `git pull`. Ensure that
the Gemfile.lock was updated to the newest versions and run the following:


```
gem install nanoc
gem install bundler
bundle install
```
You now should be able to rake with the updated version of ruby and gems. 
You can double check the version and installed gem list with the following as
we did before:

```
ruby --version
gem list
```

To see if the updates are successful, instead of manually running `rake`, it is 
recommended to tail the log of the website build and monitor there if it is 
successful when staging does its scheduled tasks:
`tail -f /home/biocadmin/bioc-test-web/log/update_site.log`

**Word of Caution** Many Bioconductor systems are inter-related. After updating
monitor that the website is building but also keep tabs on other systems that may
connect to staging. A key example: When updating the sequel gem, the SPB was 
affected do to tighter regulations and when connections were tested. 
