# ![](/images/icons/magnifier.gif)How to Create a Bioconductor Mirror Site #

The Bioconductor package repositories may be mirrored with `rsync`.  If
you would like to become a mirror for package and data package
repositories, please use the commands below.

## BioC release repos ##

If you want to mirror the current Bioconductor release version
(currently <%= config[:release_version] %>),
please use the following commands:


### Directory structure

Pick a destination directory where files will be mirrored. Let's say this will be in `/dest`. 
This directory should be served by your web server.
Under that you'll need a directory called `packages`. 
This directory must be present as it is part of the structure
of a Bioconductor repository.
Underneath `packages` should be a directory corresponding to 
the versions of Bioconductor that you will host. The current
release version is <%= config[:release_version] %> and the current
devel version is <%= config[:devel_version] %>. We recommend you
use symlinks called `release` and `devel` that always point to 
the current release and devel versions; this way you will never
have to change your rsync commands. __But__ you should change the
symlink targets with every Bioconductor release (see
the [release schedule](/developers/release-schedule/) for
exact dates).

The following
commands will create the directory structure you'll need (remember
that `/dest` is just an example of the destination directory
you could use; you can put this directory anywhere on your system
where there is enough free space). 

    mkdir -p /dest/packages
    mkdir /dest/packages/<%= config[:release_version] %> # current release
    mkdir /dest/packages/<%= config[:devel_version] %> # current devel
    ln -s /dest/packages/<%= config[:release_version] %> /dest/packages/release # change these links
    ln -s /dest/packages/<%= config[:release_version] %> /dest/packages/devel   # every 6 months (with Bioc release)

### All Bioconductor release repos (RECOMMENDED) ###

    rsync -zrtlv --delete master.bioconductor.org::release /dest/packages/release

### Bioconductor release Software repo ###

    rsync -zrtlv --delete master.bioconductor.org::release/bioc /dest/packages/release /bioc

### Bioconductor release Data repos ###

    rsync -zrtlv --delete master.bioconductor.org::release/data /dest/packages/release/data

### Bioconductor release Extra repo ###

    rsync -zrtlv --delete master.bioconductor.org::release/extra /dest/packages/release/extra


## BioC devel repos ##

If you want to mirror the Bioconductor 
devel repos (currently <%= config[:devel_version] %>),
please use the following commands:

### All Bioconductor <%= config[:devel_version] %> repos (RECOMMENDED) ###

    rsync -zrtlv --delete master.bioconductor.org::devel /dest/packages/devel

### Bioconductor release Software repo ###

    rsync -zrtlv --delete master.bioconductor.org::devel/bioc /dest/packages/devel/bioc

### Bioconductor devel Data repos ###

    rsync -zrtlv --delete master.bioconductor.org::devel/data /dest/packages/devel/data


### Bioconductor devel Extra repo ###

    rsync -zrtlv --delete master.bioconductor.org::devel/extra /dest/packages/devel/extra


## Additional information ##

Make sure the directory above `packages` is served by
a web server. 

Bioconductor is **big** (> 188GB for BioC <%= config[:devel_version] %>). Please check the size of
what will be transferred with e.g. `rsync -avn master.bioconductor.org::release`
and make sure you have enough room on your local disk before you
start.

It is recommended that package repositories be synced once per day,
scheduled with cron.

**Begin** using your new local repository, by making it accessible on
your webserver. See the **"contriburl"** option to
**install.packages()** (utils) for more information.

**Finally**, [contact us](mailto:webmaster@bioconductor.org) if you
would like to have your mirror listed on 
[our mirror page](/about/mirrors/) and in R's
<code>chooseBioCmirror()</code> function.

The [Bioconductor](/) master package repositories reside at [Fred
Hutchinson Cancer Research Center](http://www.fhcrc.org/) in Seattle,
WA, USA.
