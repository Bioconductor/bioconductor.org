# Bioconductor in the cloud
<a name="top"></a>

[Obtain](#first-time-steps) an Amazon Web Services account and 
<b><a target="start_ami"
href="https://console.aws.amazon.com/cloudformation/home?region=us-east-1#cstack=sn~StartBioconductorAMI|turl~https://s3.amazonaws.com/bioc-cloudformation-templates/start_instance.json">start
the AMI</a></b>. Additional instructions below.

## Contents
* <a href="#overview">Overview</a>
* <a href="#preloaded_ami">Preloaded AMI</a>
* <a href="#first-time-steps">First-Time Steps</a>
* <a href="#launching">Launching The AMI</a>
* <a href="#connecting_ssh">Connecting to your AMI using SSH</a>
* <a href="#ami_ids">AMI IDs</a>
* <a href="#scenarios">Scenarios for using your Bioconductor instance</a>
* <a href="#rgraphviz">Using Rgraphviz</a>
* <a href="#parallel">Parallelization using the parallel package</a>
* <a href="#using_cluster">Using the AMI as a cluster</a>
* <a href="#installing_starcluster">Installing StarCluster</a>
* <a href="#configuring_starcluster">Configuring StarCluster</a>
* <a href="#starting_cluster">Starting a Cluster</a>
* <a href="#connecting_cluster">Connecting to a Cluster</a>
* <a href="#terminate_cluster">Terminating the Cluster</a>
* <a href="#cluster_scenarios">Cluster Scenarios</a>
* <a href="#BiocParallel">Using BiocParallel with Sun Grid Engine</a>
* <a href="#ssh_backend">Using SSH as the back end</a>
* <a href="#mpi">Using MPI as the back end</a>
* <a href="#custom">Creating a custom version of the Bioconductor AMI</a>
* <a href="#provisioning">Provisioning a virtual or physical machine for use with Bioconductor</a>
* <a href="#movingdata">Moving data to and from your Bioconductor AMI instance</a>
* <a href="#questions">Questions</a>

<a name="overview"></a>
## Overview

We have developed an Amazon Machine Image (AMI) that is optimized for running Bioconductor in the Amazon Elastic Compute Cloud 
(or EC2) for sequencing tasks.

Here are a few reasons you could use it:

* You do not want to install Bioconductor on your own machine.
* You have a long-running task and you don't want it to tie up the CPU on your own machine.
* You have a parallelizable task and would like to run it (either on multiple CPUs on a single machine, or in a cluster of many machines).
* You want to run R in your web browser (using RStudio Server).
* The AMI contains many packages which can be very
  difficult to install and configure.

See below for more specific scenarios.

<a name="preloaded_ami"></a>
## Preloaded AMI

The AMI comes pre-loaded with the latest release version of R, 
and the following Bioconductor packages (and all their CRAN dependencies):

<ul class="inline_list">
    <li>affxparser</li>
    <li>affy</li>
    <li>affyio</li>
    <li>affylmGUI</li>
    <li>annaffy</li>
    <li>annotate</li>
    <li>AnnotationDbi</li>
    <li>aroma.light</li>
    <li>BayesPeak</li>
    <li>baySeq</li>
    <li>Biobase</li>
    <li>BiocInstaller</li>
    <li>biomaRt</li>
    <li>Biostrings</li>
    <li>BSgenome</li>
    <li>Category</li>
    <li>ChIPpeakAnno</li>
    <li>chipseq</li>
    <li>ChIPseqR</li>
    <li>ChIPsim</li>
    <li>CSAR</li>
    <li>cummeRbund</li>
    <li>DESeq</li>
    <li>DEXSeq</li>
    <li>DiffBind</li>
    <li>DNAcopy</li>
    <li>DynDoc</li>
    <li>EDASeq</li>
    <li>edgeR</li>
    <li>gage</li>
    <li>genefilter</li>
    <li>geneplotter</li>
    <li>GenomeGraphs</li>
    <li>genomeIntervals</li>
    <li>GenomicFeatures</li>
    <li>GenomicRanges</li>
    <li>Genominator</li>
    <li>GEOquery</li>
    <li>GGBase</li>
    <li>GGtools</li>
    <li>girafe</li>
    <li>goseq</li>
    <li>GOstats</li>
    <li>graph</li>
    <li>GSEABase</li>
    <li>HilbertVis</li>
    <li>impute</li>
    <li>IRanges</li>
    <li>limma</li>
    <li>MEDIPS</li>
    <li>multtest</li>
    <li>oneChannelGUI</li>
    <li>PAnnBuilder</li>
    <li>preprocessCore</li>
    <li>qpgraph</li>
    <li>qrqc</li>
    <li>R453Plus1Toolbox</li>
    <li>RBGL</li>
    <li>Repitools</li>
    <li>rGADEM</li>
    <li>Rgraphviz</li>
    <li>Ringo</li>
    <li>Rolexa</li>
    <li>Rsamtools</li>
    <li>Rsubread</li>
    <li>rtracklayer</li>
    <li>segmentSeq</li>
    <li>seqbias</li>
    <li>seqLogo</li>
    <li>ShortRead</li>
    <li>snpStats</li>
    <li>splots</li>
    <li>SRAdb</li>
    <li>tkWidgets</li>
    <li>VariantAnnotation</li>
    <li>vsn</li>
    <li>widgetTools</li>
    <li>zlibbioc</li>
</ul>

Plus the following categories of annotation package:

<ul class="inline_list">
    <li>org.*</li>
    <li>BSgenome.*</li>
    <li>PolyPhen.*</li>
    <li>SIFT.*</li>
    <li>TxDb.*</li>
</ul>


## How To Use It


<a name="first-time-steps"></a>
### First-time steps

First you will need an [Amazon Web Services](http://aws.amazon.com/) (AWS) account if you do not already have one. Sign up for AWS and then click [here](http://aws-portal.amazon.com/gp/aws/developer/subscription/index.html?productCode=AmazonEC2) to sign up for the EC2 service. This will require that you provide credit
card information; however, you will only be charged for services used.
[Some AWS services](http://aws.amazon.com/free/) are free.


That's all that is required if you want to use RStudio Server
to connect to your AMI with a web browser.
If you also want to connect to it with SSH, create a keypair as follows:

#### Creating a Key Pair
<a name="keypairs"></a>
Launch the [AWS Console](https://console.aws.amazon.com/ec2/home?region=us-east-1). 
Click on the [Key Pairs](https://console.aws.amazon.com/ec2/home?region=us-east-1#s=KeyPairs)
link in the lower left-hand corner of the page. Click the "Create Key Pair" button. When prompted, supply a name.
We suggest that the name be a combination of "bioconductor", your first name, and your machine name
(this will avoid conflicts with other people who share your AWS account, or possibly your own account on another machine).
For example, if your name is Bob and your personal computer is named "mylaptop", your key pair name could be "bioconductor-bob-mylaptop". Download the resulting .pem file and keep it in a safe
place.


<a name="launching"></a>
## Launching the AMI

Once you have [created an AWS account](#first-time-steps), you can launch the AMI simply by clicking on this link:

<b><a target="start_ami" href="https://console.aws.amazon.com/cloudformation/home?region=us-east-1#cstack=sn~StartBioconductorAMI|turl~https://s3.amazonaws.com/bioc-cloudformation-templates/start_instance.json">Start AMI</a></b>

Alternative links:

* <a target="start_ami" href="https://console.aws.amazon.com/cloudformation/home?region=us-east-1#cstack=sn~StartBioconductorAMI|turl~https://s3.amazonaws.com/bioc-cloudformation-templates/start_ssh_instance.json">Start AMI with SSH</a>


You'll see the following screen:

<img src="/images/ami/stack1.jpg"/>

Click "continue".

<img src="/images/ami/stack2.jpg"/>


On this screen, you can choose which version of Bioconductor you want to run.
If you are not sure, use the version that is already filled in.
You can also choose an EC2 [instance type](http://aws.amazon.com/ec2/instance-types/).
The default, t1.micro, is free to use under AWS's 
[free usage tier](http://aws.amazon.com/free/) if you use it for less than 750
hours a month.
After choosing Bioconductor version and instance type, click Continue.

<img src="/images/ami/stack3.jpg"/>

Click Continue here to launch the AMI. If you like, you can click Cost
to see how much it will cost to run the AMI with the selected instance type
(if you have selected the t1.micro instance type, be sure and click
the "FREE USAGE TIER" box in the page that comes up).


You'll see the following screen:

<img src="/images/ami/stack4.jpg"/>

Click Close.

In a few moments, the AMI will be ready (when Status changes to
CREATE_COMPLETE). You may need to click the Refresh button in the
upper right-hand corner of the screen (not your browser's refresh button).
You can then click on the Outputs tab to get
the URL and login information for your instance:

<img src="/images/ami/stack5.jpg"/>

Click on the link shown in the Stack Outputs table under URL.
You can then log in to
RStudio server using the username and password shown.

**Important Note**: When you are finished using the AMI, be sure and
shut it down to avoid incurring extra charges. Shut it down by
going to the
[CloudFormation Console](https://console.aws.amazon.com/cloudformation/home?region=us-east-1)
and checking the box next to StartBioconductorAMI.
Then click "Delete Stack" and confirm by clicking "Yes, Delete":

<img src="/images/ami/stack6.jpg">

<a name="connecting_ssh"></a>
## Connecting to your AMI using SSH

Use the following URL to start your AMI:

<b><a target="start_ami_with_ssh" href="https://console.aws.amazon.com/cloudformation/home?region=us-east-1#cstack=sn~StartBioconductorAMIWithSSH|turl~https://s3.amazonaws.com/bioc-cloudformation-templates/start_ssh_instance.json">Start AMI with SSH</a></b>

Follow the same steps as [above](#launching), but give AWS the name
of a key-pair that you created in the [first-time steps](#first-time-steps).
(A list of your keypairs is available 
[here](https://console.aws.amazon.com/ec2/home?region=us-east-1#s=KeyPairs)).

The Outputs tab will display the ssh command you should use to connect
to your instance.

You can vary this command. If you want to use programs or R packages that use X11, be sure and add a -X flag, to make the command something like this:

	ssh -X -i bioconductor-bob-mylaptop.pem ubuntu@ec2-50-16-120-30.compute-1.amazonaws.com

Now you can paste your command line into a terminal or Command Prompt. Make sure you are in the same directory as your
key pair file. 

**Windows Users**: You will need to install a version of the *ssh* and *scp* commands. Graphical programs like
[PuTTY](http://www.chiark.greenend.org.uk/~sgtatham/putty/) and [WinSCP](http://winscp.net/eng/index.php) will work. 
Our examples, however, will use the command-line versions of these programs, which you can obtain by installing 
[Cygwin](http://www.cygwin.com/) (be sure and install the *openssh* package).

Once you have pasted this command into your Terminal or Command Prompt window (and pressed Enter) you should be connected
to your Amazon EC2 instance. 


<a name="ami_ids"></a>
## AMI IDs

Our AMIs have the following IDs.

<table border="1" cellpadding="5" cellspacing="0">
  <thead valign="bottom">
    <tr>
        <th>Bioconductor Version</th>
        <th>R Version</th>
        <th>AMI ID</th>
    </tr>
   </thead>
  <tbody valign="top">
    <tr>
        <td>3.0 (devel)</td>
        <td>3.1.0</td>
        <td><%= config[:ami_ids][:bioc3_0]%></td>
    </tr>
    <tr>
        <td>2.14 (release, <b>recommended</b>)</td>
        <td>3.1.0</td>
        <td><%= config[:ami_ids][:bioc2_14]%></td>
    </tr>
    <tr>
        <td>2.13</td>
        <td>3.0.2</td>
        <td><%= config[:ami_ids][:bioc2_13]%></td>
    </tr>
    <tr>
        <td>2.12</td>
        <td>3.0</td>
        <td><%= config[:ami_ids][:bioc2_12]%></td>
    </tr>
    <tr>
        <td>2.11</td>
        <td>2.15</td>
        <td><%= config[:ami_ids][:bioc2_11]%></td>
    </tr>
    <tr>
        <td>2.10</td>
        <td>2.15</td>
        <td><%= config[:ami_ids][:bioc2_10]%></td>
    </tr>
    <tr>
        <td>2.9</td>
        <td>2.14</td>
        <td><%= config[:ami_ids][:bioc2_9]%></td>
    </tr>
    <tr>
        <td>2.8</td>
        <td>2.13</td>
        <td><%= config[:ami_ids][:bioc2_8]%></td>
    </tr>
  </tbody>
</table>

<br/>
*Please note that AMI IDs may change over time as we update the underlying AMI. Refer to this page for the most
current AMI IDs. These AMIs live in the US-East-1 region.*


<a name="scenarios"></a>
## Scenarios for using your Bioconductor instance

<a name="rgraphviz"></a>
### Using Rgraphviz

Make sure you have connected to your instance either with a web browser, or
using the -X flag of the *ssh* command. Something like:

	ssh -X -i bioconductor-bob-mylaptop.pem ubuntu@ec2-50-16-120-30.compute-1.amazonaws.com

Then, from within R on the remote instance:

	library("Rgraphviz")
	set.seed(123)
	V <- letters[1:10]
	M <- 1:4
	g1 <- randomGraph(V, M, 0.2)
	plot(g1)
	

This should start a graphics device on your local computer displaying a simple graph.

<a name="parallel"></a>
### Paralellization using the parallel package

This works best if you have selected a high-CPU instance type to run. 

This trivial example runs the <code>rnorm()</code> function, but any function would work. Consult the 
parallel documentation for more information.

    library(parallel)
    mclapply(1:30, rnorm)

<a name="using_cluster"></a>
## Using the AMI as a cluster

You can also use the AMI as a cluster, wherein the machines communicate
with each other via one of the following mechanisms:

* SSH
* MPI
* Sun Grid Engine

In order to use the Bioconductor AMI in one of these scenarios,
you need to install [StarCluster](http://star.mit.edu/cluster/),
which is a software package  designed to automate and simplify
the process of building, configuring, and managing clusters of
virtual machines on Amazon’s EC2 cloud.

StarCluster takes care of the details of making sure that the
machines can communicate with each other, including:

* Passwordless SSH
* Shared disk space (using NFS)
* Convenient aliases for host names (such as master and node001)
* Configuration of job scheduler (Sun Grid Engine)


**Note**: Using the Bioconductor AMI for cluster-related tasks
is only supported for the Bioconductor AMI version 2.14 and higher.


<a name="installing_starcluster"></a>
### Installing StarCluster

Install StarCluster by following the
[StarCluster Installation Guide](http://star.mit.edu/cluster/docs/latest/installation.html). This is a simple and fast process.

If you are on Windows, go directly to 
[the Windows section](http://star.mit.edu/cluster/docs/latest/installation.html#installing-on-windows).

Before continuing, it's worth watching the 
[Quick Start Screencast](http://star.mit.edu/cluster/)
or following the 
[Quick-Start Tutorial](http://star.mit.edu/cluster/docs/latest/quickstart.html).

<a name="configuring_starcluster"></a>
### Configuring StarCluster

Before we can use StarCluster with the Bioconductor AMI,
we need to configure it, by editing its `config` file.

You can create the file by issuing the command:

    starcluster help

This will give you three options:

    Options:
    --------
    [1] Show the StarCluster config template
    [2] Write config template to /home/user/.starcluster/config
    [q] Quit

Choose option 2 and note the location of the config file (it will 
be different from what is shown above).

On Unix systems (including Linux and Mac OS X), this file is found
at `~/.starcluster/config`. On Windows systems, the `.starcluster`
folder should be located in your 
[home directory](http://weka.wikispaces.com/Where+is+my+home+directory+located%3F).

Open the `config` file in your favorite text editor, and edit it 
as follows:


#### AWS Credentials and Connection Settings section

You need to fill in values for `AWS_ACCESS_KEY_ID`
and `AWS_SECRET_ACCESS_KEY`. If you don't know these values,
go to the [Security Credentials](https://console.aws.amazon.com/iam/home?#security_credential)
page of the AWS Console and expand the "Access Keys" section.
You can view or create your access keys here. Be sure and
store these credentials in a safe place (in addition to
your StarCluster config file).

The value of `AWS_USER_ID` can also be found on
the [Security Credentials](https://console.aws.amazon.com/iam/home?#security_credential) page, by expanding
the "Account Identifiers" section. Fill in `AWS_USER_ID`
with the number shown as your "AWS Account ID" (this
should be a 12-digit number with hyphens).

#### Defining EC2 Keypairs section

If you haven't already created a keypair in EC2,
please do so now by reading the [keypairs](#keypairs) section.

You can also create a keypair with StarCluster; run the command

    starcluster createkey --help

...for instructions.

Remember the name that you assigned to your keypair.
Change the line 

    [key mykey]

So that `mykey` is replaced by the name you assigned
to your keypair in EC2, and change the following line

    KEY_LOCATION=~/.ssh/mykey.rsa

So that the value of `KEY_LOCATION` is the full path to
the private key downloaded from .ec2 (it probably has a `.pem`
extension).

#### Defining Cluster Templates section

StarCluster allows you to define multiple clusters in the config
file. For now let's just modify the cluster defined as `smallcluster`.

* Change the value of `KEYNAME` to the name of your key pair
  (see keypair section above).
* Optionally change `CLUSTER_SIZE` to the number of machines you
  want to launch. This number includes the master node, so the
  default value of 2 means one master, and one worker. We
  recommend starting with 2 until you have familiarized yourself
  with using StarCluster and Bioconductor.
* Change `CLUSTER_USER` to `ubuntu`.
* Uncomment the line `DNS_PREFIX = True`. This makes your cluster
  instances easier to recognize when using the AWS Console or 
  command line tools.
* Change the `NODE_IMAGE_ID` to the AMI-ID of the AMI you want to use
  This will be listed in the [AMI IDs](#ami_ids) section of this document.
  Note that StarCluster only works with AMIs for Bioconductor
  version 2.14 and higher.
* Optionally change `NODE_INSTANCE_TYPE` to another instance type.
  See the [Instance Types page](https://aws.amazon.com/ec2/instance-types/)
  for more information.
* Under the line reading `#PERMISSIONS = ssh, http`, add the line
  `permissions = http` (note lowercase). This is related to security
  group permissions (more about this below).

You can make additional changes to this section if you
want to further customize your configuration. Refer
to the
[StarCluster documentation](http://star.mit.edu/cluster/docs/latest/manual/configuration.html#creating-the-configuration-file)
for more information.

#### Configuring Security Group Permissions section

Remove the comments (`#` symbol) from the four lines 
starting with `[permission http]`
so that you end up with:

    [permission http]
    IP_PROTOCOL = tcp
    FROM_PORT = 80
    TO_PORT = 80

This allows port 80 on the cluster instances to be open to the
world, allowing us to use Rstudio Server on that port.

<a name="starting_cluster"></a>
### Starting a Cluster

Assuming you have completed the steps above, you can create
a cluster with the command:

    starcluster start smallcluster

After a few moments, the cluster should be available.

<a name="connecting_cluster"></a>
### Connecting to the cluster

There are two ways to connect to the cluster's master node: RStudio Server and
SSH. Unless you have a special need to use SSH, we recommend 
using RStudio Server.

#### Connecting using RStudio Server

First, get the hostname of the master node by issuing the command:

    starcluster listclusters


You can also abbreviate this:

    starcluster lc

This will produce output like the following:

```
-----------------------------------------------
smallcluster (security group: @sc-smallcluster)
-----------------------------------------------
Launch time: 2014-06-16 09:57:54
Uptime: 0 days, 02:19:56
Zone: us-east-1b
Keypair: bioc-default
EBS volumes: N/A
Cluster nodes:
    smallcluster-master running i-46a76c6d ec2-54-91-23-93.compute-1.amazonaws.com
    smallcluster-node001 running i-47a76c6c ec2-54-224-6-153.compute-1.amazonaws.com
Total nodes: 2
```

The line that starts with `smallcluster-master` ends with a hostname
(in this case it's `ec2-54-91-23-93.compute-1.amazonaws.com`; in your
case it will be something different but similar). You can paste
this host name into a web browser (depending on the browser, you may
need to put `http://` in front of the host name).

This should bring you to the RStudio Server login page. You can log 
in with the username `ubuntu` and the password `bioc`.

#### Connecting using SSH

To connect to the master node using ssh, simply issue the command

    starcluster sshmaster --user=ubuntu smallcluster

<a name="terminate_cluster"></a>
### Terminating the Cluster <font color="red">**IMPORTANT!!**</font>

When you are done, you MUST terminate your cluster or you
will continue to be charged money by Amazon Web Services.
To terminate the cluster, do this:

    starcluster terminate smallcluster

This command will prompt you to confirm that you really want to 
terminate the cluster.

<a name="cluster_scenarios"></a>
### Cluster Scenarios

The following scenarios assume that you have started up a cluster 
and that you are connected to the master node.




<a name="BiocParallel"></a>
### Using BiocParallel with Sun Grid Engine

When you start a cluster with StarCluster, it's automatically
configured to use the `BiocParallel` and `BatchJobs` packages
with Sun Grid Engine as the back end. You can demonstrate this
by loading BatchJobs:

    library(BatchJobs)

Among other output, this will say

          cluster functions: SGE

indicating that Sun Grid Engine is the back end.

Here's how to send a simple job to the cluster:

```
library(BatchJobs)
library(BiocParallel)
param <- BatchJobsParam(2, resources=list(ncpus=1))
register(param)
FUN <- function(i) system("hostname", intern=TRUE)
xx <- bplapply(1:100, FUN)
table(unlist(xx))
```
This will produce:

```
 smallcluster-master smallcluster-node001 
                  50                   50 
```

...indicating that SGE ran half the jobs on the master and
the other half on the worker node.

[This presentation](/help/course-materials/2014/ISMB2014/vjcScalable2.pptx)
outlines how to develop a full-scale analysis 
(identification of cis-dsQTL)
using the kind
of cluster we've just created.

<a name="ssh_backend"></a>
### Using SSH as the back end

Here is the same example as above, except using 
SSH instead of Sun Grid Engine as the back end:

```
library(BatchJobs)
library(BiocParallel)
cluster.functions <- makeClusterFunctionsSSH(
    makeSSHWorker(nodename="smallcluster-master"),
    makeSSHWorker(nodename="smallcluster-node001")
)


param2 <- BatchJobsParam(2, resources=list(ncpus=1),
    cluster.functions=cluster.functions)
register(param2)
FUN <- function(i) system("hostname", intern=TRUE)
xx <- bplapply(1:10, FUN)
table(unlist(xx))
```
You should see results like this:

```
 smallcluster-master smallcluster-node001 
                   5                    5 
```


<a name="mpi"></a>
### Using MPI as the back end

When you start a cluster using the above steps, R is automatically
aware of the cluster, as shown in the following example:

    library(Rmpi)
    mpi.universe.size()

With the default cluster configuration, this should return 2, which
make sense, since our cluster consists of two machines (`smallcluster-master`
and `smallcluster-node001`), each of type m1.small, which 
have one core each.

Again using `BiocParallel`, you can run a simple function
on your MPI cluster:

```
FUN <- function(i) system("hostname", intern=TRUE)
```

Create a `SnowParam` instance with the number of nodes equal to
the size of the MPI universe minus 1 (let one node dispatch jobs to
workers), and register this instance as the default:

```
param3 <- SnowParam(mpi.universe.size() - 1, "MPI")
register(param3)

```

Evaluate the work in parallel and process the results:

```
xx <- bplapply(1:10, FUN)
table(unlist(xx))
```

<a name="custom"></a>
### Creating a custom version of the Bioconductor AMI


*Note*: If you make changes to the running Bioconductor AMI, and then terminate the AMI, your changes will be lost. 
Use the steps described here to ensure that your changes are persistent.

If the AMI is missing some packages or features you think it should have, please let us know.

If you want to customize the AMI for your own purposes, it is simple. Just go ahead and customize your 
running instance as you see fit. Typically this will involve installing R packages with <code>biocLite()</code>,
and software packages (at the operating system level) with the Ubuntu package manager <code>apt-get</code>.

When your instance is customized to your liking, issue the following commands:

   
	sudo clean_ami
	exit

This will remove unneeded files from your AMI, clear your shell history, and log you out of your instance.
You may also want to change the password of the "ubuntu" user (because the default password
is publicly known, in order to run RStudio Server) with the command:

	passwd ubuntu
	


Now use the [AWS Console](https://console.aws.amazon.com/ec2/home?region=us-east-1#s=Instances)
to Stop your instance (**important note: do NOT "Terminate" your instance; use the Stop command (under Instance Actions) instead.**)

Then choose "Create Image (EBS AMI)" under the Instance Actions menu. You will be prompted for a name for your AMI. After
entering the name, your AMI will be created and given a unique AMI ID. You can then launch instances of this AMI using the steps above,
being sure to substitute the ID of your own AMI. Your AMI will be private, accessible only to your AWS account, unless you decide
to make it more widely accessible.

Now you should Terminate the Stopped instance of the Bioconductor AMI.

<a name="provisioning"></a>
### Provisioning a virtual or physical machine for use with Bioconductor

The Bioconductor AMI was created using Vagrant and Chef.
The same scripts that were used to create these AMIs
can also be used to provision virtual machines
(Virtualbox or VMWare) or physical machines.

For more information, see the scripts' 
[github repository](https://github.com/Bioconductor/setup-starcluster-image).


<a name="movingdata"></a>
## Moving data to and from your Bioconductor AMI instance

If you are using RStudio Server, you can upload and download files 
from the Files pane in RStudio server.

If you are connected via ssh,
the *scp* command is the most efficient way to move data to and from your EC2 instances.

To copy a file from your computer to a running Bioconductor AMI instance:

* Open a Terminal or Command Prompt window on your computer
* Using "cd", change to the directory where your key pair (.pem) file lives
* Issue a command like the following (your key pair name and the hostname of the AMI instance will be different; you can determine the correct values by clicking on your running instance in the [AWS Console](https://console.aws.amazon.com/ec2/home?region=us-east-1#s=Instances)):

	
	
<span style="display:none">Hidden</span>
	scp -i bioconductor-bob-mylaptop.pem /path/to/myfile ubuntu@ec2-50-16-120-30.compute-1.amazonaws.com:~
	
That will copy the file at "/path/to/myfile" on your local computer to ubuntu's home directory on the remote instance. To copy a file from a 
running instance to your local computer, do something like this (still at your local computer):

	scp -i bioconductor-bob-mylaptop.pem ubuntu@ec2-50-16-120-30.compute-1.amazonaws.com:~/myfile /some/directory

That will copy the file ~/myfile from the running instance to /some/directory on your local machine.

*Reminder*: Files created on a running EC2 instance are **not** persisted unless you do some special steps. So if you are 
generating output files with Bioconductor, you must copy them to your local machine before terminating your instance, or
your files will be lost. If you have a lot of data to move back and forth, you may want to look into 
[Elastic Block Storage](http://aws.amazon.com/ebs/).

<a name="questions"></a>
## Questions


If you have questions about the Bioconductor AMI, please contact us through the [Bioconductor mailing list](/help/mailing-list/#bioconductor).


