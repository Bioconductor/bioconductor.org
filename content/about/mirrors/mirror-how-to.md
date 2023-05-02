# ![](/images/icons/magnifier.gif)How to a Bioconductor repository #

Bioconductor package repositories may be mirrored with `rsync`.  This
is appropriate if you or your user community requires frequent local
access to many packages, where access to the main repository would be
expensive or impossible (e.g., because users can only access
repositories behind a firewall).

## SSH key and IP address for `rsync`

Using `rsync` requires that you provide Bioconductor with minimal information
that includes (a) with an ssh public key and (b) the IP address(es) from which
you will perform `rsync`. If you would like to request rsync access please fill
out this google form: [private mirror/rsync request form](https://forms.gle/d42JmFCfyJPjpsWT8)

[webmaster@bioconductor.org]: mailto:webmaster@bioconductor.org

## Public mirrors

A mirror is considered "public" if it is an option in R's
`chooseBioCmirror()` function and listed on our [mirrors
page](/about/mirrors/). Public mirrors must support https on their
site. If you are interested in hosting a publicly available mirror
site, please fill out this google form: [public Bioc mirror request form](https://forms.gle/2BREvZfQfwgo2rSR6)

## Structure of the `rsync` command

The overall structure of the `rsync` command is

    rsync SSH_OPTION -e "ssh -i /path/to/ssh" [OPTIONS] SRC DEST

`SSH_OPTION` is required, and tells `rsync` to use SSH during the
transfer. An appropriate command might be

    -e "ssh -i ~/.ssh/"

This can often be abbreviated to `-e "ssh"`, and customized through an
SSH config file.

`[OPTIONS]` determine how files are synchronized with the server. For
a mirror, appropriate values are `-zrtlv --delete`.

`SRC` consists of the account and host for the connection, and the
path on the host to the hierarchy to be synchronized. For instance, to
synchronize the software packages on the current release branch one would use

    bioc-rsync@master.bioconductor.org:<%= config[:release_version] %>/bioc
    
`DEST` is the location on the local file system of the synchronized
repository, e.g.,

    ~/bioconductor_repositories/<%= config[:release_version] %>/bioc

Thus a complete command might be

    rsync -e "ssh -i ~/.ssh" -zrtlv --delete \
        bioc-rsync@master.bioconductor.org:<%= config[:release_version] %>/bioc \
        ~/bioconductor_repositories/<%= config[:release_version] %>/bioc

## Example: `rsync` an entire Bioconductor release

This is appropriate if you are providing your user community with a
version of the entire Bioconductor release, with software, data
annotation, and data experiment packages. The following uses
Bioconductor version <%= config[:release_version] %>, the current
release.

### Directory structure

Pick a destination directory where files will be mirrored. Let's say
this will be in `/dest`.  This directory should be served by your web
server.  Under that you'll need a directory called `packages`.  This
directory must be present as it is part of the structure of a
Bioconductor repository.  Underneath `packages` should be a directory
corresponding to the versions of Bioconductor that you will host. The
current release version is <%= config[:release_version] %> and the
current devel version is <%= config[:devel_version] %>. We recommend
you use symlinks called `release` and `devel` that always point to the
current release and devel versions; this way you will never have to
change your rsync commands. __But__ you should change the symlink
targets with every Bioconductor release (see the
[release schedule](/developers/release-schedule/) for exact dates).

The following commands will create the directory structure you'll need
(remember that `/dest` is just an example of the destination directory
you could use; you can put this directory anywhere on your system
where there is enough free space).

    mkdir -p /dest/packages
    mkdir /dest/packages/<%= config[:release_version] %> # current release
    mkdir /dest/packages/<%= config[:devel_version] %> # current devel
    ln -s /dest/packages/<%= config[:release_version] %> /dest/packages/release # change these links
    ln -s /dest/packages/<%= config[:release_version] %> /dest/packages/devel   # every 6 months (with Bioc release)

### `rsync` the Bioconductor release repository

To synchronize an entire release, use the command

    rsync -e "ssh" -zrtlv --delete bioc-rsync@master.bioconductor.org:release /dest/packages/release

It is also possible to separately synchronize just the software
packages...

    rsync -e "ssh" -zrtlv --delete bioc-rsync@master.bioconductor.org:release/bioc /dest/packages/release/bioc

...or the data annotation and data experiment packages

    rsync -e "ssh" -zrtlv --delete bioc-rsync@master.bioconductor.org:release/data /dest/packages/release/data

### `rsync` the Bioconductor devel repository

To mirror all Bioconductor 'devel' (version <%= config[:devel_version] %>)
repositories:

    rsync -e "ssh" -zrtlv --delete bioc-rsync@master.bioconductor.org:devel /dest/packages/devel

Bioconductor devel software repository:

    rsync -e "ssh" -zrtlv --delete bioc-rsync@master.bioconductor.org:devel/bioc /dest/packages/devel/bioc

Bioconductor devel annotation and experiment data repositories:

    rsync -e "ssh" -zrtlv --delete bioc-rsync@master.bioconductor.org:devel/data /dest/packages/devel/data

## Additional information ##

Make sure the directory above `packages` is served by a web server. 

Bioconductor is **big** (> 188GB for BioC <%= config[:devel_version] %>).
Please check the size of what will be transferred with e.g. `rsync -e "ssh"
-avn bioc-rsync@master.bioconductor.org:release` and make sure you have enough
room on your local disk before you start.

It is recommended that package repositories be synchronized once per
day, scheduled with cron.

**Begin** using your new local repository by making it accessible on your
webserver. See the **"contriburl"** option to **install.packages()** (utils)
for more information.

The [Bioconductor](/) package repositories for the current release and devel
reside at [https://master.bioconductor.org](https://master.bioconductor.org) in
the Amazon cloud.
