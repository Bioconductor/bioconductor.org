# How to Create a Bioconductor Mirror Site #

The Bioconductor package repositories may be mirrored with `rsync`.  If
you would like to become a mirror for package and data package
repositories, please use the commands below.

## BioC 2.6 repos ##

If you want to mirror the Bioconductor **2.6** repos (the current
release version), please use the following commands:

### All Bioconductor 2.6 repos ###

    rsync -rtlv --delete bioconductor.org::2.6 /dest

### Bioconductor 2.6 Software repo ###

    rsync -rtlv --delete bioconductor.org::2.6/bioc /dest

### Bioconductor 2.6 Data repos ###

    rsync -rtlv --delete bioconductor.org::2.6/data /dest

### Bioconductor 2.6 Extra repo ###

    rsync -rtlv --delete bioconductor.org::2.6/extra /dest


## BioC 2.7 repos ##

If you want to mirror the Bioconductor **2.7** repos (the current
devel version), please use the following commands:

### All Bioconductor 2.7 repos ### 

    rsync -rtlv --delete bioconductor.org::2.7 /dest

### Bioconductor 2.7 Software repo ### 

    rsync -rtlv --delete bioconductor.org::2.7/bioc /dest

### Bioconductor 2.7 Data repos ### 

    rsync -rtlv --delete bioconductor.org::2.7/data /dest


### Bioconductor 2.7 Extra repo ### 

    rsync -rtlv --delete bioconductor.org::2.7/extra /dest


## Additional information ##

Bioconductor is big (> 50G for BioC 2.4). Please check the size of
what will be transferred with e.g. `rsync -avn bioconductor.org::2.6`
and make sure you have enough room on your local disk before you
start.

It is recommended that package repositories be synced once per day,
scheduled with cron.

**Begin** using your new local repository, by making it accessible on
your webserver. See the **"contriburl"** option to
**install.packages()** (utils) for more information.

**Finally**, [contact us](mailto:webmaster@bioconductor.org) if you
would like to have your mirror listed on this page.

The [Bioconductor](/) master package repositories reside at [Fred
Hutchinson Cancer Research Center](http://www.fhcrc.org/) in Seattle,
WA, USA.

        
