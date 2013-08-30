# ![](/images/icons/magnifier.gif)How to Create a Bioconductor Mirror Site #

The Bioconductor package repositories may be mirrored with `rsync`.  If
you would like to become a mirror for package and data package
repositories, please use the commands below.

## BioC <%= config[:release_version] %> repos ##

If you want to mirror the Bioconductor **<%= config[:release_version] %>** repos (the current
release version), please use the following commands:

### All Bioconductor <%= config[:release_version] %> repos ###

    rsync -zrtlv --delete master.bioconductor.org::<%= config[:release_version] %> /dest/bioc_<%= config[:release_version] %>

### Bioconductor <%= config[:release_version] %> Software repo ###

    rsync -zrtlv --delete master.bioconductor.org::<%= config[:release_version] %>/bioc /dest/bioc_<%= config[:release_version] %>/bioc

### Bioconductor <%= config[:release_version] %> Data repos ###

    rsync -zrtlv --delete master.bioconductor.org::<%= config[:release_version] %>/data /dest/bioc_<%= config[:release_version] %>/data

### Bioconductor <%= config[:release_version] %> Extra repo ###

    rsync -zrtlv --delete master.bioconductor.org::<%= config[:release_version] %>/extra /dest/bioc_<%= config[:release_version] %>/extra


## BioC <%= config[:devel_version] %> repos ##

If you want to mirror the Bioconductor **<%= config[:devel_version] %>** repos (the current
devel version), please use the following commands:

### All Bioconductor <%= config[:devel_version] %> repos ###

    rsync -zrtlv --delete master.bioconductor.org::<%= config[:devel_version] %> /dest/bioc_<%= config[:devel_version] %>

### Bioconductor <%= config[:devel_version] %> Software repo ###

    rsync -zrtlv --delete master.bioconductor.org::<%= config[:devel_version] %>/bioc /dest/bioc_<%= config[:devel_version] %>/bioc

### Bioconductor <%= config[:devel_version] %> Data repos ###

    rsync -zrtlv --delete master.bioconductor.org::<%= config[:devel_version] %>/data /dest/bioc_<%= config[:devel_version] %>/data


### Bioconductor <%= config[:devel_version] %> Extra repo ###

    rsync -zrtlv --delete master.bioconductor.org::<%= config[:devel_version] %>/extra /dest/bioc_<%= config[:devel_version] %>/extra


## Additional information ##

Sync the above directories to a directory called `packages` on your
system, and make sure the directory above `packages` is served by
a web server. 

Bioconductor is big (> 64G for BioC <%= config[:release_version] %>). Please check the size of
what will be transferred with e.g. `rsync -avn master.bioconductor.org::<%= config[:release_version] %>`
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
