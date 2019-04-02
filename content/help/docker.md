# Docker containers for Bioconductor

[Docker](https://www.docker.com) allows software to be packaged into
containers: self-contained environments that contain everything
needed to run the software. Containers can be run anywhere
(containers run in modern Linux kernels, but can be run
on Windows and Mac as well using a virtual machine called
[Docker Toolbox](https://www.docker.com/products/docker-toolbox).
Containers can also be deployed in the cloud using
[Amazon EC2 Container Service](https://aws.amazon.com/ecs/)
or other cloud providers.


<a name="top"></a>

- [Why Use Containers](#intro)
- [Current Containers](#current)
- [Legacy Containers](#legacy)
- [Using Containers](#usage)
  * [Running Containers](#running)
  * [Mounting Additional Volume](#mounting)
- [Modifying Image Container](#modify)
- [Default package with Core2 Container](#core)
- [Acknowledgements](#acknowledgements)

<a name="intro"></a>

## Why use Bioconductor containers

With Bioconductor containers, we hope to enhance

* *reproducibility*: If you run some code in a container today,
  you can run it again in the same container (with the same
  [tag](https://docs.docker.com/userguide/dockerimages/#setting-tags-on-an-image))
  years later and know that nothing in the container has changed.
  You should always take note of the tag you used if you think
  you might want to reproduce some work later.
* *ease of use*: With one command, you can be running the
  latest release or devel Bioconductor. No need to worry
  about whether packages and system dependencies are
  installed.
* *convenience*: Sometimes you just want a fresh R with
  no packages installed, in order to test something; or
  you typically don't have microarray packages installed
  but suddenly you need to do a microarray analysis.
  Containers make this easy.

Our aim is to provide up-to-date containers for the current
release and devel versions of Bioconductor, and some older
versions. Bioconductor's Docker images are stored in [Docker
Hub](https://hub.docker.com/u/bioconductor/); 
the source Dockerfiles are in
[Github](https://github.com/Bioconductor/bioc_docker).

Our release images are based on
[rocker/rstudio](https://github.com/rocker-org/rocker/tree/master/rstudio) and
built when a Biocondcutor Release occurs. Our devel images are based on
[rocker/rstudio-daily](https://github.com/rocker-org/rstudio-daily) and built
weekly with the latest versions of R and Bioconductor packages.

For each supported version of Bioconductor, we provide several
images:

* *base2*: Contains R, RStudio, and a single Bioconductor
  package (`BiocManager`,
  providing the `install()` function for installing additional
  packages).
  Also contains many system dependencies for Bioconductor packages.
  Useful when you want a relatively blank slate for testing purposes.
  R is accessible via the command line or via RStudio Server.
* *core2*: Built on *base2*, so it contains everything in *base2*, plus
  a selection of core Bioconductor packages.
  See [the full list](#the-full-list).
* *protmetcore2*: Built on *core2*, so it contains everything in *core2*, plus
  a selection of core Bioconductor packages recommended for proteomic and
  metabolomics analysis.
* *metabolomics2*: everything in *protmetcore2*, plus select packages from the
  [Metabolomics](/packages/release/BiocViews.html#___Metabolomics) biocView.
* *cytometry2*: Built on *base2*, so it contains everything in *base2*, plus
  [CytoML](/packages/release/bioc/html/CytoML.html) and needed dependencies.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="current"></a>

## Current Containers


###### Maintained by the Bioconductor Core Team: bioc-issue-bot@bioconductor.org
* [bioconductor/devel_base2](https://hub.docker.com/r/bioconductor/devel_base2/)
* [bioconductor/devel_core2](https://hub.docker.com/r/bioconductor/devel_core2/)
* [bioconductor/release_base2](https://hub.docker.com/r/bioconductor/release_base2/)
* [bioconductor/release_core2](https://hub.docker.com/r/bioconductor/release_core2/)

###### Maintained by Steffen Neumann: sneumann@ipb-halle.de
Maintained as part of the "PhenoMeNal, funded by Horizon2020 grant 654241"

* [bioconductor/devel_protmetcore2](https://hub.docker.com/r/bioconductor/devel_protmetcore2/)
* [bioconductor/devel_metabolomics2](https://hub.docker.com/r/bioconductor/devel_metabolomics2/)
* [bioconductor/release_protmetcore2](https://hub.docker.com/r/bioconductor/release_protmetcore2/)
* [bioconductor/release_metabolomics2](https://hub.docker.com/r/bioconductor/release_metabolomics2/)

###### Maintained by Laurent Gatto: lg390@cam.ac.uk
* [bioconductor/devel_mscore2](https://hub.docker.com/r/bioconductor/devel_mscore2/)
* [bioconductor/devel_protcore2](https://hub.docker.com/r/bioconductor/devel_protcore2/)
* [bioconductor/devel_proteomics2](https://hub.docker.com/r/bioconductor/devel_proteomics2/)
* [bioconductor/release_mscore2](https://hub.docker.com/r/bioconductor/release_mscore2/)
* [bioconductor/release_protcore2](https://hub.docker.com/r/bioconductor/release_protcore2/)
* [bioconductor/release_proteomics2](https://hub.docker.com/r/bioconductor/release_proteomics2/)

#### Maintained by RGLab: wjiang2@fredhutch.org
* [bioconductor/devel_cytometry2](https://hub.docker.com/r/bioconductor/devel_cytometry2/)
* [bioconductor/release_cytometry2](https://hub.docker.com/r/bioconductor/release_cytometry2/)


<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="legacy"></a>

## Legacy Containers

The following containers are legacy and no longer updated. They have been kept
to retain previous versions available via tags:

* bioconductor/devel_base
* bioconductor/devel_core
* bioconductor/devel_flow
* bioconductor/devel_microarray
* bioconductor/devel_proteomics
* bioconductor/devel_sequencing
* bioconductor/devel_metabolomics
* bioconductor/release_base
* bioconductor/release_core
* bioconductor/release_flow
* bioconductor/release_microarray
* bioconductor/release_proteomics
* bioconductor/release_sequencing
* bioconductor/release_metabolomics

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="usage"></a>

## Using the containers

A well organized guide to popular docker commands can be found
[here](https://github.com/wsargent/docker-cheat-sheet). For convenience, below
are some commands to get you started. The following examples use the
`bioconductor/devel_base2` image.

**Note:** that you may need to prepend `sudo` to all `docker` commands.


**Prerequisites**: On Linux, you need Docker
[installed](https://docs.docker.com/installation/) and
on [Mac](http://docs.docker.com/installation/mac/)
or [Windows](http://docs.docker.com/installation/windows/)
you need Docker Toolbox installed and running.


##### List which docker machines are available locally
    docker images

##### List running containers
    docker ps

##### List all containers
    docker ps -a

##### Resume a stopped container
    docker start <CONTAINER ID>

##### Shell into a running container
    docker exec -it <CONTAINER ID> /bin/bash

##### Shutdown container
    docker stop <CONTAINER ID>

##### Delete container
    docker rm <CONTAINER ID>

##### Delete image
    docker rmi bioconductor/devel_base2

<a name="running"></a>

#### Running the container

The above commands can be helpful but the real basics of running a Bioconductor
docker involves pulling the public image and running the container. 

##### Get a copy of public docker image 
    docker pull bioconductor/devel_base2

##### To run RStudio Server:

    docker run -e PASSWORD=<pickYourPassword> -p 8787:8787 bioconductor/devel_base2

You can then open a web browser pointing to your docker host on port 8787.
If you're on Linux and using default settings, the docker host is
`127.0.0.1` (or `localhost`, so the full URL to RStudio would be
`http://localhost:8787)`. If you are on Mac or Windows and running
`Docker Toolbox`, you can determine the docker host with the
`docker-machine ip default` command. 

In the above command `-e PASSWORD=` is
setting the rstudio password and is required by the rstudio docker image. It can
be whatever you like except it cannot be `rstudio`.
Log in to RStudio with the username `rstudio` and whatever password was specified.

If you want to run RStudio as a user on your host machine, in order
to read/write files in a host directory, please
[read this](https://github.com/rocker-org/rocker/wiki/Sharing-files-with-host-machine).

##### To run R from the command line:

    docker run -ti --user bioc bioconductor/devel_base2 R

##### To open a Bash shell on the container:

    docker run -ti --user bioc bioconductor/devel_base2 bash


**Note**: The `docker run` command is very powerful and versatile.
For full documentation, type `docker run --help` or visit
the [help page](https://docs.docker.com/reference/run/).

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="mounting"></a>

##### Mounting Additional Volume

One such option for `docker run` is `-v` to mount an
additional volume to the docker image. This might be useful for say mounting a
local R install directory for use on the docker. The path on the docker image
that should be mapped to a local R library directory is
`/usr/local/lib/R/host-site-library`.  The follow example would mount my locally
installed packages to this docker directory. In turn, that path is automatically
loaded in the R `.libPaths` on the docker image and all of my locally installed
package would be available for use.

    docker run \
        -v /home/lori/R/x86_64-pc-linux-gnu-library/3.6-BioC-3.9:/usr/local/lib/R/host-site-library \
        -it --user bioc bioconductor/devel_base2 R

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="modify"></a>


## Modifying the images

There are two ways to modify these images:

1. Making changes in a running container and then committing them
   using the `docker commit` command.
```
      docker commit <CONTAINER ID> <name for new image> 
```
2. Using a Dockerfile to declare the changes you want to make.

The second way is the recommended way. Both ways are
[documented here](https://docs.docker.com/userguide/dockerimages/#creating-our-own-images).

<a name="the-full-list"></a>

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="core"></a>

## List of packages installed on the *core2* container

These packages, plus their dependencies, are installed:

<ul class="inline_list">
<li>BiocManager</li>
<li>OrganismDbi</li>
<li>ExperimentHub</li>
<li>Biobase</li>
<li>BiocParallel</li>
<li>biomaRt</li>
<li>Biostrings</li>
<li>BSgenome</li>
<li>ShortRead</li>
<li>IRanges</li>
<li>GenomicRanges</li>
<li>GenomicAlignment</li>
<li>GenomicFeatures</li>
<li>SummarizedExperiment</li>
<li>VariantAnnotation</li>
<li>DelayedArray</li>
<li>GSEABase</li>
<li>Gviz</li>
<li>graph</li>
<li>RBGL</li>
<li>Rgraphviz</li>
<li>rmarkdown</li>
<li>httr</li>
<li>knitr</li>
<li>BiocStyle</li>
</ul>

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="acknowledgements"></a>

### Acknowledgements

Thanks to the
[rocker](https://github.com/rocker-org/rocker) project for providing the
R/RStudio Server containers upon which ours are based.
