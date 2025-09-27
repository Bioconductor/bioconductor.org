# Docker containers for Bioconductor

[Docker](https:/docs.docker.com/engine/docker-overview/) packages software
into self-contained environments, called containers, that include necessary
dependencies to run. Containers can run on any operating system including
Windows and Mac (using modern Linux kernels) via the
[Docker engine](https://docs.docker.com/engine/) or
[Docker Desktop](https://www.docker.com/products/docker-desktop).

Containers can also be deployed in the cloud using
[Amazon Elastic Container Service](https://aws.amazon.com/ecs/),
[Google Kubernetes Engine](https://cloud.google.com/kubernetes-engine/)
or [Microsoft Azure Container Instances](https://azure.microsoft.com/en-us/services/container-instances/)

## Quick start

### Docker commands

1.  [Install Docker Engine](https://docs.docker.com/engine/install/) or
  [Docker Desktop](https://www.docker.com/products/docker-desktop)

1.  Run container with Bioconductor and RStudio

        docker run \
        	-e PASSWORD=bioc \
        	-p 8787:8787 \
        	bioconductor/bioconductor_docker:devel

This command will run the Docker container
`bioconductor/bioconductor_docker:devel` on your local machine.

RStudio will be available on your web browser at
`http://localhost:8787`. The USER is fixed to always being
`rstudio`. The password in the above command is given as `bioc` but
it can be set to anything. `8787` is the port being mapped between
the Docker container and your host machine. NOTE: password cannot
be `rstudio`.

The user is logged into the `rstudio` user by default.

### The `bioc-run` script

The `bioc-run` script is a convenience script that can be used to run
Bioconductor Docker images. The script is available at
<https://github.com/Bioconductor/bioc-run>. Use `bioc-run -h` or
see the `README.md` for additional details.

Execute the `bioc-run` script to create a container from the `RELEASE_3_20`
Bioconductor Docker image:

	./bioc-run -v RELEASE_3_20

Note that the script also mounts a local directory to the Docker image to
persist installed packages between sessions. It also maps a local directory to
the `/home/rstudio` directory in the Docker image.

<a name="intro"></a>

## Why use Containers

Bioconductor containers enhance:

- **Reproducibility**: Containers run with pre-installed versions of R and
  Bioconductor. These versions do not change and can be run at any time
  in the future as long as the Bioconductor `RELEASE_X_Y` version tag is noted.

- **Ease of use**: With one command can run either release or devel versions
  of Bioconductor with support for nearly all package system dependencies.

- **Convenience**: Run tests with a fresh R session and a minimal set of 
  pre-installed packages. Quickly run analyses from atypical workflows with
  pre-installed system dependencies.

- **Package Installation**: Bioconductor publishes binary packages for
  fast installation of packages within containers (`RELEASE_3_14` and newer).
  Binary packages do not require compilation and install about 7 to 8 times
  faster than source package installations.

Bioconductor provides up-to-date containers for the current release
and devel versions and supports older release versions. Bioconductor's Docker
images are stored in [Docker Hub](https://hub.docker.com/u/bioconductor);
the source `Dockerfile`(s) are on
[Github](https://github.com/Bioconductor/bioconductor_docker).

Our release images and devel images rely on
[Rocker Project](https://www.rocker-project.org/) -
[rocker/rstudio](https://github.com/rocker-org/rocker/tree/master/rstudio)
images and are built after Bioconductor releases.

<a name="goals"></a>

### Goals for new container architecture

A few of our key goals to migrate to a new set of Docker containers are,

- to keep the image size being shipped by the Bioconductor team at a
  manageable size.

- easy to extend, so developers can just use a single image to
  inherit and build their Docker image.

- easy to maintain, by streamlining the Docker inheritance chain.

- Adopt a "best practices" outline so that new community contributed
  Docker images get reviewed and follow standards.

- Adopt a deprecation policy and life cycle for images similar to
  Bioconductor packages.

- Replicate the Linux build machines on the
  `bioconductor/bioconductor_docker:devel` image as closely as
  possible. While this is not fully possible just yet, this image can
  be used by maintainers who wish to reproduce errors seen on the
  Bioconductor Linux build machine and as a helpful debugging tool.

<a name="current"></a>

## Current Containers

For each supported version of Bioconductor, we provide

- **bioconductor/bioconductor_docker:RELEASE_X_Y**

- **bioconductor/bioconductor_docker:devel**

<a name="usage"></a>

## Using the containers

A well organized guide to popular docker commands can be found
[here](https://github.com/wsargent/docker-cheat-sheet). For
convenience, below are some commands to get you started. The following
examples use the `bioconductor/bioconductor_docker:devel` image.

##### List which Docker machines are available locally

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

    docker rmi bioconductor/bioconductor_docker:devel

<a name="running"></a>

### Running a Bioconductor container

The above commands can be helpful but the real basics of running a
Bioconductor Docker involves pulling the public image and running the
container.

##### Download a public Docker image

    docker pull bioconductor/bioconductor_docker:devel

##### To run the RStudio server:

    docker run -e PASSWORD=<password> \
    	-p 8787:8787 \
    	bioconductor/bioconductor_docker:devel

Open a web browser and browse to `http://localhost:8787` or
`http://127.0.0.1:8787` (where `8787` is the port number specified in the
`docker run` command).

Set the RStudio password with `-e PASSWORD=` (required). The password is
arbitrary since the container is running locally but it cannot be
`rstudio`. Log in to RStudio with the username `rstudio` and the password
specified in the `docker run` command.

If you want to run RStudio as a user on your host machine, in order to
read and write files in a host directory, please see the shared volumes
[documentation](https://rocker-project.org/use/shared_volumes.html).

NOTE: If you forget to add the tag `devel` or `RELEASE_X_Y` while
using the `bioconductor/bioconductor_docker` image, it will
automatically use the `latest` tag which points to the latest `RELEASE_X_Y`
version of Bioconductor.

##### To run R from the command line:

    docker run -it --user rstudio bioconductor/bioconductor_docker:devel R

##### To open a Bash shell on the container:

    docker run -it --user rstudio bioconductor/bioconductor_docker:devel bash

**Note**: For full documentation, type `docker run --help` or visit
the [help page](https://docs.docker.com/engine/containers/run/).

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="mounting"></a>

### Mounting Additional Volume

Use the `-v` flag with `docker run` to mount a volume to the Docker image. This
is useful to keep a more permanent R package installation directory to use with
Docker. Always map the local directory to the path on the Docker image that
corresponds to the local R library: `/usr/local/lib/R/host-site-library`.

The following example mounts and maps a designated folder for
container-installed packages to the Docker directory. Using the designated
`host-site-library` folder ensures that the path is automatically picked up by R
and shown via `.libPaths()` in the Docker image.

      docker run -it \
        -e PASSWORD=<password> \
      	-v /home/user/docker-devel-packages:/usr/local/lib/R/host-site-library \
      	-p 8787:8787 \
      	bioconductor/bioconductor_docker:devel bash

  The `-it` flag gives you an interactive tty (shell/terminal) to the
  Docker container.

- Running it with RStudio interface

      docker run \
      	-v /home/user/docker-devel-packages:/usr/local/lib/R/host-site-library \
      	-e PASSWORD=password \
      	-p 8787:8787 \
      	bioconductor/bioconductor_docker:devel

<a name="dockercompose"></a>

### Using docker-compose

**Note**. Ensure that `docker-compose` is installed. Docker recommends
installing `Docker Desktop` to get `docker-compose`.

The `Bioconductor/bioconductor_docker` repository has a `docker-compose.yaml`
file that can be used to run the Bioconductor Docker image with the command:

```
docker-compose up
```

Within the same directory, the user can run `docker-compose` to launch the
Bioconductor Docker image and access it via `http://localhost:8787`.

The `docker-composer.yaml` includes pre-configured settings for the the port
(`8787`), password (default is `bioc`), and the volume for storing container R
packages.

The configuration sets the local library folder to
`$HOME/R/bioconductor_docker/<bioconductor_version>`. For example, if the 
the Bioconductor version is `3.13`, the local library folder will be
`$HOME/R/bioconductor_docker/3.13`. This location is mapped and mounted to
`/usr/local/lib/R/host-site-library` in the container. Ensure that the
`host-site-library` is available by running `.libPaths()` within in the
container.

When the user starts the Docker image using `docker-compose`, it will
recognize previously mounted libraries with the appropriate
Bioconductor version, and save users time re-installing previously
installed packages.

To add another volume, modify the `docker-compose.yml` to include another
volume. For example, to add a volume for the user's home directory, add the 
following line to the `docker-compose.yml` file:

```
volumes:
	- ${HOME}/R/bioconductor_docker/3.13:/usr/local/lib/R/host-site-library
	- ${HOME}/dockerhome:/home/rstudio
```

To run in the background, use the `-d` or `--detach` flag:

```
docker-compose up -d
```

If the image is run in the background, the `container-name` can be used to
`exec` into the container terminal with `root` access. Run the following command
to get the `container-name`:

```
docker ps -a
```

The `docker exec` command allows the user to install additional system-level
dependencies with `root` access:

```
docker exec -it <container-name> bash
```

For more information on how to use `docker-compose`, use the
[official docker-compose reference](https://docs.docker.com/compose/reference/up/).

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="modify"></a>

## Modifying Docker Images

There are two ways to modify Docker images:

1.  Making changes in a running container and then committing them
    using the `docker commit` command.

        docker commit <CONTAINER ID> <name for new image>

2.  Using a `Dockerfile` to formalize the changes you want to make.

Using a `Dockerfile` is the recommended way to modify and even extend
a Docker image. For more details, see the
[documentation](https://docs.docker.com/get-started/docker-concepts/building-images/).

### Use case: Python package installation

An example scenario may be to add the `tensorflow` Python package to the
`bioconductor/bioconductor_docker:devel` image and at the same time
install the `scAlign` Bioconductor package on top of the base Docker
image.

First, the `Dockerfile` could inherit from the
`bioconductor/bioconductor_docker:devel` image. Note that some knowledge of
Linux is required to install the `tensorflow` package on the Ubuntu image.

In the `Dockerfile`, add the following commands:

    # Docker inheritance
    FROM bioconductor/bioconductor_docker:devel

    # Update apt-get
    RUN apt-get update \
    	## Install the python package tensorflow
    	&& pip install tensorflow		\
    	## Remove packages in '/var/cache/' and 'var/lib'
    	## to remove side-effects of apt-get update
    	&& apt-get clean \
    	&& rm -rf /var/lib/apt/lists/*

    # Install required Bioconductor package
    RUN R -e 'BiocManager::install("scAlign")'

Then build the `Dockerfile` with the command:

    docker build -t bioconductor_docker_tensorflow:devel .

Note that the image name is arbitrary (e.g., `bioconductor_docker_tensorflow`).

Once built, the image can be run with the command:

    docker run -p 8787:8787 -e PASSWORD=bioc bioconductor_docker_tensorflow:devel

### Use case: Adding LaTeX to the Docker image

Another potential use case could be to add LaTeX to the Docker image.
LaTeX is required to build vignettes and knit documents into PDF files.
The `Dockerfile` could look like the following:

    # This docker image has LaTeX to build the vignettes
    FROM bioconductor/bioconductor_docker:devel

    # Update apt-get
    RUN apt-get update \
    	&& apt-get install -y --no-install-recommends apt-utils \
    	&& apt-get install -y --no-install-recommends \
    	texlive \
    	texlive-latex-extra \
    	texlive-fonts-extra \
    	texlive-bibtex-extra \
    	texlive-science \
    	texi2html \
    	texinfo \
    	&& apt-get clean \
    	&& rm -rf /var/lib/apt/lists/*

    ## Install BiocStyle
    RUN R -e 'BiocManager::install("BiocStyle")'

Then build the `Dockerfile` with the command:

    docker build -t bioconductor_docker_latex:devel .

Once built, the image can be run with the command:

    docker run -p 8787:8787 -e PASSWORD=bioc bioconductor_docker_latex:devel

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="binarypackages"></a>

### Using Binary Packages

Binary packages are available for Bioconductor containers starting from
`RELEASE_3_14`. This means that, for all `RELEASE_3_14` and newer images,
Bioconductor packages will be pre-compiled and installed via
`BiocManager::install()`, reducing the installation time significantly.
Binary package installations provide a 7 to 8 times speed up compared to
installing from source.

To install binary packages, simply use `BiocManager::install()` within a
Bioconductor container:

    ## Install binary packages on a container
    BiocManager::install(c('Rhtslib','SingleCellExperiment'))

Note that the container needs to be the `bioconductor/bioconductor_docker` or
a derived image. To build custom images, see the
[Modifying Docker Images](#modifying-docker-images) section.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="singularity"></a>

## Singularity

Singularity is an alternative to Docker that is typically used on compute
clusters where admin access may not be available.

To check if `singularity` is installed on your High Performance Cluster (HPC),
run the following:

    module available

If Singularity is available, enable `singularity` with:

    module load singularity

Then, convert a Docker image to a Singularity image with the `singularity pull`
command: 

    singularity pull docker://bioconductor/bioconductor_docker:devel


For specific usage instructions relevant to Singularity containers see
<https://www.rocker-project.org/use/singularity/>.

<a name="msft"></a>

## Microsoft Azure Container Instances

Microsoft Azure users have the option to run containers using images on the
[Microsoft Container Registry](https://github.com/microsoft/ContainerRegistry).

> Microsoft Container Registry (MCR) is the primary Registry for all Microsoft
Published docker images that offers a reliable and trustworthy delivery of
container images with a syndicated catalog

<a name="mcr"></a>

### Using containers hosted on Microsoft Container Registry

Pull the `bioconductor_docker` image from the Microsoft Container
Registry (MCR), specifying a `tag` of choice. Check the
[MCR](https://hub.docker.com/r/microsoft/bioconductor-bioconductor-docker/)
for the list of tags under "Full Tag Listing":

    docker pull mcr.microsoft.com/bioconductor/bioconductor_docker:<tag>

To pull the latest image:

    docker pull mcr.microsoft.com/bioconductor/bioconductor_docker:latest

#### Usage: Run RStudio from the Docker container

To run RStudio server, run the following and open a browser to `127.0.0.1:8787`.
The default user name is `rstudio` and you can specify a password:

    docker run --name bioconductor_docker_rstudio \
    	-v ~/host-site-library:/usr/local/lib/R/host-site-library \
    	-e PASSWORD='bioc'                               \
    	-p 8787:8787                                     \
    	mcr.microsoft.com/bioconductor/bioconductor_docker:latest

To run R in the container terminal:

    docker run --name bioconductor_docker_rstudio \
    	-it                                            \
    	-v ~/host-site-library:/usr/local/lib/R/host-site-library \
    	-e PASSWORD='bioc'                               \
    	-p 8787:8787                                     \
    	mcr.microsoft.com/bioconductor/bioconductor_docker:latest R

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="aci"></a>

### Use Azure Container Instances to run Bioconductor images on-demand on Azure

[Azure Container Instances](https://azure.microsoft.com/en-us/services/container-instances/#features)
(ACI) provide a way to run Docker containers on-demand in a managed,
serverless Azure environment. To learn more, see the
[documentation](https://docs.microsoft.com/en-us/azure/container-instances/container-instances-overview).

### Run Bioconductor images using ACI

**Prerequisites**:

1. [An Azure account and a subscription](https://learn.microsoft.com/en-us/azure/cost-management-billing/manage/create-subscription)

2. [Azure CLI](https://learn.microsoft.com/en-us/cli/azure/install-azure-cli)

3. [Resource group permissions](https://learn.microsoft.com/en-us/azure/azure-resource-manager/management/manage-resource-groups-portal)

You can run Azure CLI or `az cli`
[commands](https://docs.microsoft.com/en-us/cli/azure/?view=azure-cli-latest)
to create, stop, restart or delete container instances running any
official Bioconductor image - either from Bioconductor or  available on the 
[MCR](https://hub.docker.com/r/microsoft/bioconductor-bioconductor-docker).
An Azure account and subscription is required. To create a free account
browse to <https://azure.microsoft.com/en-us/free/>.

To get familiar with Azure Container Instances, follow
[this tutorial](https://learn.microsoft.com/en-us/azure/container-instances/container-instances-quickstart).

To run the Bioconductor image hosted on the MCR, create a new resource group in
the Azure subscription and then run the following command with the Azure CLI:

    az container create \
    	--resource-group resourceGroupName \
    	--name mcr-bioconductor \
    	--image mcr.microsoft.com/bioconductor/bioconductor_docker \
    	--cpu 2 \
    	--memory 4 \
    	--dns-name-label mcr-bioconductor \
    	--ports 8787 \
    	--environment-variables 'PASSWORD'='bioc'

When completed, run this command to get the fully qualified domain name (FQDN):

    az container show \
    	--resource-group resourceGroupName \
    	--name mcr-bioconductor \
    	--query "{FQDN:ipAddress.fqdn,ProvisioningState:provisioningState}" \
    	--out table

Here we expose port `8787` on this publicly accessible FQDN. You may
have to choose a different "dns-name-label" to avoid conflicts. By
default, the username for RStudio is "rstudio" (as in the Bioconductor Docker
image). The password is set to 'bioc' in the environment variable configuration.
The `--cpu` and `--memory` values (in gigabytes; GB) can be configured as necessary.
By default, ACI have 1 CPU and 1.5 GB of memory allocated.

To learn more about configuring and customizing ACI, run:

    az container create --help

#### Mount Azure File Share to persist analysis data between sessions

To persist data between analysis sessions when using ACIs, mount an
[Azure file share](https://learn.microsoft.com/en-us/azure/container-instances/container-instances-volume-azure-files)
to the container. The following are steps to create an ACI that maps the
`/home/rstudio` directory in RStudio from an
[Azure File Share](https://learn.microsoft.com/en-us/azure/storage/files/storage-files-introduction):

1. Create an [Azure Storage](https://learn.microsoft.com/en-us/azure/storage/common/storage-account-create)
   account

2. Create an [Azure file share](https://learn.microsoft.com/en-us/azure/storage/files/storage-how-to-use-files-portal)

3. Get the [storage account key](https://learn.microsoft.com/en-us/cli/azure/storage/account/keys)

Using the Azure CLI, run the following commands to create the storage account
and file share:

```
# Change these four parameters as needed
ACI_PERS_RESOURCE_GROUP=resourceGroupName
ACI_PERS_STORAGE_ACCOUNT_NAME=storageAccountName
ACI_PERS_LOCATION=eastus
ACI_PERS_SHARE_NAME=fileShareName

# Step1: Create the storage account with the parameters
az storage account create \
	--resource-group $ACI_PERS_RESOURCE_GROUP \
	--name $ACI_PERS_STORAGE_ACCOUNT_NAME \
	--location $ACI_PERS_LOCATION \
	--sku Standard_LRS

# Step2: Create the file share
az storage share create \
	--name $ACI_PERS_SHARE_NAME \
	--account-name $ACI_PERS_STORAGE_ACCOUNT_NAME

# Step3: Get the storage account key
STORAGE_KEY=$(az storage account keys list \
	--resource-group $ACI_PERS_RESOURCE_GROUP \
	--account-name $ACI_PERS_STORAGE_ACCOUNT_NAME \
	--query "[0].value" --output tsv)
echo $STORAGE_KEY
```

Mount an Azure file share to an ACI running Bioconductor:

    az container create \
    	--resource-group resourceGroupName \
    	--name mcr-bioconductor-fs \
    	--image mcr.microsoft.com/bioconductor/bioconductor_docker \
    	--dns-name-label mcr-bioconductor-fs \
    	--cpu 2 \
    	--memory 4 \
    	--ports 8787 \
    	--environment-variables 'PASSWORD'='bioc' \
    	--azure-file-volume-account-name storageAccountName \
    	--azure-file-volume-account-key $STORAGE_KEY \
    	--azure-file-volume-share-name fileShareName \
    	--azure-file-volume-mount-path /home/rstudio

When completed, run this command to get the fully qualified domain name or FQDN:

    az container show \
    	--resource-group resourceGroupName \
    	--name mcr-bioconductor-fs \
    	--query "{FQDN:ipAddress.fqdn,ProvisioningState:provisioningState}" \
    	--out table

Note that the `/home/rstudio` directory is mapped to a
persistent Azure file share named "fileShareName" in the storage
account specified. When you stop or restart the ACI, this data will not
be lost.

#### Stop, Start, Restart or Delete containers running on ACI

You can run Azure CLI commands to
[stop, start, restart](https://learn.microsoft.com/en-us/azure/container-instances/container-instances-stop-start)
or
[delete](https://learn.microsoft.com/en-us/azure/container-instances/container-instances-quickstart#clean-up-resources)
container instances on Azure. You can find all the commands and
options at
<https://learn.microsoft.com/en-us/cli/azure/container?view=azure-cli-latest#commands>.

##### Stop the container instance

Note that `containerName` and `resourceGroupName` in the following CLI commands
should be replaced with the actual container name and resource group name.

    az container stop -n containerName -g resourceGroupName

##### Start the container instance

    az container start -n containerName -g resourceGroupName

##### Restart the container instance

    az container restart -n containerName -g resourceGroupName

##### Delete the container instance

    az container delete -n containerName -g resourceGroupName

Use the `-y` flag to avoid the confirmation prompt when deleting the ACI:

    az container delete -n containerName -g resourceGroupName -y

To troubleshoot any issues when using ACIs, see the
[common issues page](https://learn.microsoft.com/en-us/azure/container-instances/container-instances-troubleshooting). For
feedback or further issues, contact
[genomics@microsoft.com](mailto:genomics@microsoft.com).

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>

<a name="contribute"></a>

## How to Contribute

There is a comprehensive list of best practices and standards on how community
members can contribute images:

<https://github.com/Bioconductor/bioconductor_docker/blob/master/best_practices.md>

<a name="deprecation"></a>

## Deprecation Notice

For previous users of Docker containers for Bioconductor, please note
that we are deprecating the following Bioconductor and community maintained
images.

<a name="legacy"></a>

### Legacy Containers

These images are **NO LONGER MAINTAINED**. They will however be available to use
should a user choose.

Bioconductor Core Team: bioc-issue-bot@bioconductor.org

- [bioconductor/devel_base2](https://hub.docker.com/r/bioconductor/devel_base2/)
- [bioconductor/devel_core2](https://hub.docker.com/r/bioconductor/devel_core2/)
- [bioconductor/release_base2](https://hub.docker.com/r/bioconductor/release_base2/)
- [bioconductor/release_core2](https://hub.docker.com/r/bioconductor/release_core2/)

Steffen Neumann: sneumann@ipb-halle.de, Maintained as part of the "PhenoMeNal, funded by Horizon2020 grant 654241"

- [bioconductor/devel_protmetcore2](https://hub.docker.com/r/bioconductor/devel_protmetcore2/)
- [bioconductor/devel_metabolomics2](https://hub.docker.com/r/bioconductor/devel_metabolomics2/)
- [bioconductor/release_protmetcore2](https://hub.docker.com/r/bioconductor/release_protmetcore2/)
- [bioconductor/release_metabolomics2](https://hub.docker.com/r/bioconductor/release_metabolomics2/)

Laurent Gatto: lg390@cam.ac.uk

- [bioconductor/devel_mscore2](https://hub.docker.com/r/bioconductor/devel_mscore2/)
- [bioconductor/devel_protcore2](https://hub.docker.com/r/bioconductor/devel_protcore2/)
- [bioconductor/devel_proteomics2](https://hub.docker.com/r/bioconductor/devel_proteomics2/)
- [bioconductor/release_mscore2](https://hub.docker.com/r/bioconductor/release_mscore2/)
- [bioconductor/release_protcore2](https://hub.docker.com/r/bioconductor/release_protcore2/)
- [bioconductor/release_proteomics2](https://hub.docker.com/r/bioconductor/release_proteomics2/)

RGLab: wjiang2@fredhutch.org

- [bioconductor/devel_cytometry2](https://hub.docker.com/r/bioconductor/devel_cytometry2/)
- [bioconductor/release_cytometry2](https://hub.docker.com/r/bioconductor/release_cytometry2/)

First iteration containers

- bioconductor/devel_base
- bioconductor/devel_core
- bioconductor/devel_flow
- bioconductor/devel_microarray
- bioconductor/devel_proteomics
- bioconductor/devel_sequencing
- bioconductor/devel_metabolomics
- bioconductor/release_base
- bioconductor/release_core
- bioconductor/release_flow
- bioconductor/release_microarray
- bioconductor/release_proteomics
- bioconductor/release_sequencing
- bioconductor/release_metabolomics

<a name="reason"></a>

### Reason for deprecation

The new Bioconductor Docker image `bioconductor/bioconductor_docker` has nearly
all system dependencies pre-installed. This makes it possible to install any
package and avoids the need for multiple images with different system
dependencies. Packages can be installed with `BiocManager::install()`.

Other reasons for deprecation:

- the chain of inheritance of Docker images was too complex and hard
  to maintain

- images were hard to extend because there were multiple flavors 

- naming convention was not consistent

- unused images were not deprecated

<a name="issues"></a>

### Reporting Issues

To report any issues or bugs, go to the GitHub
[issues page](https://github.com/Bioconductor/bioconductor_docker/issues).

Feel free to ask questions related to usage, extension, and enhancement of
Bioconductor images either via GitHub issues or on the
[Bioc-devel](mailto:bioc-devel@r-project.org) mailing list.

<a name="acknowledgements"></a>

## Acknowledgements

Thanks to the [rocker](https://github.com/rocker-org/rocker) project
for providing the R/RStudio Server containers upon which ours are
based.
