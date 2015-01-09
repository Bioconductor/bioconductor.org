## Docker containers for Bioconductor

[Docker](https://www.docker.com) allows software to be packaged into
containers: self-contained environments that contain everything
needed to run the software. Containers can be run anywhere
(containers run in modern Linux kernels, but can be run
on Windows and Mac as well using a virtual machine called
[boot2docker](http://boot2docker.io/)). Containers can
also be deployed in the cloud using 
[Amazon EC2 Container Service](https://aws.amazon.com/ecs/)
or other cloud providers.

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

### Available Containers

Our aim is to provide up-to-date containers for the current 
release and devel versions of Bioconductor, and (probably, eventually)
some older versions.

Bioconductor's Docker images are stored in 
[Docker Hub](https://registry.hub.docker.com/repos/bioconductor/);
the source Dockerfiles are in 
[Github](https://github.com/Bioconductor/bioc_docker).

For each supported version of Bioconductor, we provide several
images:


* *base*: Based on the latest Ubuntu LTS distribution.
  Contains R and a single Bioconductor package (`BiocInstaller`,
  providing the `biocLite()` function for installing additional
  packages).
  Also contains many system dependencies for Bioconductor packages.
  Useful when you want a relatively blank slate for testing purposes. 
  R is accessible via the command line or via RStudio Server.
  The release and devel versions of these containers (and the
  containers built from them, below) are rebuilt
  daily with the latest versions of R-release or R-devel
  (with previous versions available via tags).
* *core*: Built on *base*, so it contains everything in *base*, plus
  a selection of core Bioconductor packages.
  See [the full list](#the-full-list).
* *flow*: everything in *core*, plus all packages tagged with the
  [FlowCytometry](/packages/release/BiocViews.html#___FlowCytometry) biocView.
* *microarray*: everything in *core*, plus 
  all packages tagged with the 
  [Microarray](/packages/release/BiocViews.html#___Microarray) biocView.
* *proteomics*: everything in *core*, plus all packages tagged with the
  [Proteomics](/packages/release/BiocViews.html#___Proteomics) biocView.
* *sequencing*: everything in *core*, plus all packages tagged with the
  [Sequencing](/packages/release/BiocViews.html#___Sequencing) biocView.

#### List of Containers

At present, the following containers are available:

* bioconductor/devel_base
* bioconductor/devel_core
* bioconductor/devel_flow
* bioconductor/devel_microarray
* bioconductor/devel_proteomics
* bioconductor/devel_sequencing
* bioconductor/release_base
* bioconductor/release_core
* bioconductor/release_flow
* bioconductor/release_microarray
* bioconductor/release_proteomics
* bioconductor/release_sequencing

## Using the containers

The following examples use the `bioconductor/devel_base` container.
Note that you may need to prepend `sudo` to all `docker` commands.

**Prerequisites**: On Linux, you need Docker 
[installed](https://docs.docker.com/installation/) and
on [Mac](http://docs.docker.com/installation/mac/)
or [Windows](http://docs.docker.com/installation/windows/)
you need boot2docker installed and running.

##### To run RStudio Server:

    docker run -p 8787:8787 bioconductor/devel_base

You can then open a web browser pointing to your docker host on port 8787.
If you're on Linux and using default settings, the docker host is
`127.0.0.1` (or `localhost`, so the full URL to RStudio would be
`http://localhost:8787)`. If you are on Mac or Windows and running
`boot2docker`, you can determine the docker host with the
`boot2docker ip` command.

Log in to RStudio with the username `rstudio` and password `rstudio`.

If you want to run RStudio as a user on your host machine, in order
to read/write files in a host directory, please
[read this](https://github.com/rocker-org/rocker/wiki/Sharing-files-with-host-machine#interactive-containers).

##### To run R from the command line:

    docker run -ti bioconductor/devel_base R

##### To open a Bash shell on the container:

    docker run -ti bioconductor/devel_base bash

*Note*: The `docker run` command is very powerful and versatile. 
For full documentation, type `docker run --help` or visit
the [help page](https://docs.docker.com/reference/run/).

## Modifying the images

There are two ways to modify these images:

1. Making changes in a running container and then committing them
   using the `docker commit` command.
2. Using a Dockerfile to declare the changes you want to make.

The second way is the recommended way. Both ways are
[documented here](https://docs.docker.com/userguide/dockerimages/#creating-our-own-images).

<a name="the-full-list"></a>
## List of packages installed on the *core* container

These packages, plus their dependencies, are installed:

<ul class="inline_list">
<li>AnnotationDbi</li>
<li>AnnotationHub</li>
<li>Biobase</li>
<li>BiocParallel</li>
<li>biocViews</li>
<li>biomaRt</li>
<li>Biostrings</li>
<li>BSgenome</li>
<li>epivizr</li>
<li>GenomicFeatures</li>
<li>GenomicRanges</li>
<li>graph</li>
<li>Gviz</li>
<li>IRanges</li>
<li>RBGL</li>
<li>ReportingTools</li>
<li>Rgraphviz</li>
<li>zlibbioc</li>
</ul>

### Acknowledgements

Some code used to create these containers comes from the
[rocker](https://github.com/rocker-org/rocker) project.
