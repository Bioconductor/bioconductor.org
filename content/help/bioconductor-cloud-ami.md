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
* <a href="#multiple">Starting Up Multiple Instances</a>
* <a href="#ami_ids">AMI IDs</a>
* <a href="#scenarios">Scenarios for using your Bioconductor instance</a>
* <a href="#rgraphviz">Using Rgraphviz</a>
* <a href="#multicore">Parallelization using multicore</a>
* <a href="#mpi">Using an MPI cluster in the cloud</a>
* <a href="#parallel">Using a parallel cluster in the cloud</a>
* <a href="#custom">Creating a custom version of the Bioconductor AMI</a>
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
* The AMI contains many packages (such as RGraphviz) which can be very
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
* <a target="start_ami" href="https://console.aws.amazon.com/cloudformation/home?region=us-east-1#cstack=sn~StartBioconductorAMI|turl~https://s3.amazonaws.com/bioc-cloudformation-templates/start_multiple_instances.json">Start Multiple Instances</a>
* <a target="start_ami" href="https://console.aws.amazon.com/cloudformation/home?region=us-east-1#cstack=sn~StartBioconductorAMI|turl~https://s3.amazonaws.com/bioc-cloudformation-templates/start_multiple_ssh_instances.json">Start Multiple Instances with SSH</a>


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

	ssh -X -i bioconductor-bob-mylaptop.pem root@ec2-50-16-120-30.compute-1.amazonaws.com

If you do not want to log in as root, you can change "root" to "ubuntu".

Now you can paste your command line into a terminal or Command Prompt. Make sure you are in the same directory as your
key pair file. 

**Windows Users**: You will need to install a version of the *ssh* and *scp* commands. Graphical programs like
[PuTTY](http://www.chiark.greenend.org.uk/~sgtatham/putty/) and [WinSCP](http://winscp.net/eng/index.php) will work. 
Our examples, however, will use the command-line versions of these programs, which you can obtain by installing 
[Cygwin](http://www.cygwin.com/) (be sure and install the *openssh* package).

Once you have pasted this command into your Terminal or Command Prompt window (and pressed Enter) you should be connected
to your Amazon EC2 instance. 

<a name="multiple"></a>
## Starting Up Multiple Instances

Use the following URL to start your AMI:

<b><a target="start_ami" href="https://console.aws.amazon.com/cloudformation/home?region=us-east-1#cstack=sn~StartBioconductorAMI|turl~https://s3.amazonaws.com/bioc-cloudformation-templates/start_multiple_instances.json">Start Multiple Instances</a></b>

To start multiple instances with ssh, use this link:

<b><a target="start_ami" href="https://console.aws.amazon.com/cloudformation/home?region=us-east-1#cstack=sn~StartBioconductorAMI|turl~https://s3.amazonaws.com/bioc-cloudformation-templates/start_multiple_ssh_instances.json">Start Multiple Instances with SSH</a></b>

The procedure is similar to the one described [above](#launching), except
you enter the number of instances you want to start. Also, you will not find the 
instance hostnames in the Outputs panel of the CloudFormation console. You can
find this information in the 
[EC2 Console](https://console.aws.amazon.com/ec2/home?region=us-east-1#s=Instances).


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

	ssh -X -i bioconductor-bob-mylaptop.pem root@ec2-50-16-120-30.compute-1.amazonaws.com

Then, from within R on the remote instance:

	library("Rgraphviz")
	set.seed(123)
	V <- letters[1:10]
	M <- 1:4
	g1 <- randomGraph(V, M, 0.2)
	plot(g1)
	

This should start a graphics device on your local computer displaying a simple graph.

<a name="multicore"></a>
### Paralellization using multicore

This works best if you have selected a high-CPU instance type to run. 

This trivial example runs the <code>rnorm()</code> function, but any function would work. Consult the 
[multicore](http://cran.r-project.org/web/packages/multicore/index.html) documentation for more information.

	library(multicore)
	mclapply(1:30, rnorm)

<a name="mpi"></a>
### Using an MPI cluster in the cloud

You can launch multiple EC2 instances and set them up as an [MPI](http://en.wikipedia.org/wiki/Message_Passing_Interface)
cluster, to parallelize and shorten long-running CPU-intensive jobs. You must explicitly parallelize your
CPU-intensive code using functions in the [Rmpi](http://cran.r-project.org/web/packages/Rmpi/index.html) package. 

The simplest way to start up a cluster is to just click on this URL:

<b><a target="start_ami"
href="https://console.aws.amazon.com/cloudformation/home?region=us-east-1#cstack=sn~StartBioCMPICluster|turl~https://s3.amazonaws.com/bioc-cloudformation-templates/cluster.json">Start MPI Cluster</a></b>

That will start an MPI cluster which you can access via you web browser using RStudio Server.

If you also want to be able to access your cluster via ssh, use this URL instead (you'll
need to provide the name of an ssh keypair that you have <a href="#first-time-steps">
previously set up</a>):

<b><a target="start_ami"
href="https://console.aws.amazon.com/cloudformation/home?region=us-east-1#cstack=sn~StartBioCMPIClusterWithSSH|turl~https://s3.amazonaws.com/bioc-cloudformation-templates/cluster_ssh.json">Start MPI Cluster with ssh access</a></b>


The startup procedure is similar to <a href="#launch">the launch
procedure</a> discussed earlier, except that you are also asked
how many worker instances you want to start. The 
[EC2 instance types](http://aws.amazon.com/ec2/instance-types/) page
tells you how many cores are available with each instance type.
So if you wanted to start a cluster with 40 workers, you could
choose a NumClusterWorkers of 10 and a ClusterInstanceType of m1.xlarge.
Ten machines with four cores each gives you a cluster with 40 
workers. An additional master node will also be started up.

Once you have logged into the RStudio Server using the
link provided by the above step, you'll be able to access
your cluster as follows:

	library(Rmpi)
	mpi.spawn.Rslaves()

You can then run some code on each worker, for example:

	mpi.parLapply(1:mpi.universe.size(), function(x) x+1)

This performs a simple calculation on each worker and returns the 
results as a list. For more complex examples, read on
or consult the 
[Rmpi documentation](http://cran.r-project.org/web/packages/Rmpi/index.html).

Be sure and delete your stack when you are finished using it, 
in order to stop accruing charges.

The above procedure can handle many parallel processing tasks.

If you want your cluster nodes to have a shared disk, follow
the steps below. (This procedure will soon be replaced by 
a simpler one-click procedure like those above.)

This section assumes you have [started up your AMI and connected to it
with ssh](#connecting_ssh).

Here we present a tutorial for using the Bioconductor package
<code>ShortRead</code> with an EC2-based MPI cluster.

To configure an EC2-based MPI cluster, launch single EC2 instance as described 
<a href="#connecting_ssh">above</a>, if you haven't already.
This will be the master node of your cluster. Consider which [instance type](http://aws.amazon.com/ec2/instance-types/)
you want your cluster to consist of. Note that the master and worker nodes will be of the same instance type.

In most cases, you'll want all of your cluster nodes to have access to the same data files. We accomplish this by 
creating an [EBS Volume](http://aws.amazon.com/ebs/) holding our data. This volume is attached to the master node 
and then shared by the workers using NFS. 

We've provided a script that will create an EBS volume already populated with sample data. To use it, run the following
command from your AMI instance:

	create_sample_volume --access-key-id XXX --secret-key yyy

Where "xxx" is your Amazon Access Key ID and "yyy" is your Secret Key. The script will produce output like the following:

	Waiting for volume to be available...
	.
	Volume is available.
	Created volume vol-dec646b6 in availability zone us-east-1c.

Make a note of the volume ID and the availability zone (you can also find this information in the
[Volumes page of the EC2 console](https://console.aws.amazon.com/ec2/home?region=us-east-1#s=Volumes)). This step is not necessary if cluster nodes do not need to share a disk.

Now you're ready to spin up an MPI cluster. Use the <code>mpiutil</code> script. Invoked without arguments, 
<code>mpiutil</code> produces the following:

	Manage MPI clusters
	--access-key-id     -a  Amazon Access Key ID
	--secret-key        -s  Amazon Secret Key ID
	--num-workers       -w  Number of workers to start
	--cluster-name      -n  Name of cluster
	--halt-cluster      -h  Halt cluster
	--instance-type     -t  Instance type
	--volume-id         -v  Volume ID

	Examples:
	  To start a cluster:
	    mpiutil -a XXX -s YYY -w 2 -n "my cluster" -t t1.micro -v vol-9999999

	  To stop a cluster:
	    mpiutil -a XXX -s YYY -n "my cluster" -h -v vol-9999999
	

Let's start a cluster with three workers. Use a command line this:

	mpiutil -a XXX -s YYY -w 3 -n workers -t t1.micro -v vol-9999999

Make the following substitutions:

* For "xxx", substitute your Amazon Access Key ID.
* For "yyy", substitute your Amazon Secret Key.
* For "vol-9999999", substitute the volume ID generated by the <code>create_sample_volume</code> script above. (The -v option is not required if your task does not need access
to shared data on an EBS volume).

This command will mount your EBS volume, launch three worker instances, share the EBS volume with
them using NFS, and do some other housekeeping tasks.

You can verify that your cluster is working with the following command. You may need to wait a few moments 
before all the workers are up and running, but then this command will work:

	mpirun -np 1 --hostfile /usr/local/Rmpi/hostfile R --no-save -f /usr/local/Rmpi/mpiTest.R --args 3
	
This will run an R script which should produce some output like this:

	> library(Rmpi)
	> 
	> mpi.spawn.Rslaves(nsl = nsl)
		3 slaves are spawned successfully. 0 failed.
	master (rank 0, comm 1) of size 4 is running on: ip-10-117-50-18 
	slave1 (rank 1, comm 1) of size 4 is running on: ip-10-117-47-155 
	slave2 (rank 2, comm 1) of size 4 is running on: ip-10-117-45-245 
	slave3 (rank 3, comm 1) of size 4 is running on: ip-10-117-74-57 
	> 
	> mpi.close.Rslaves()
	> mpi.quit()

This output shows that there is a master and three worker nodes running, each with distinct IP addresses.
You may see the same IP address repeated multiple times, once for every processor core available on 
a machine.

Now, to do some actual work with your MPI cluster, run the following command:

	mpirun -np 1 --hostfile /usr/local/Rmpi/hostfile R --no-save -f /usr/local/Rmpi/ShortReadQA.R --args 3

This will run the <code>qa()</code> function from the <code>ShortRead</code> package on three files in parallel, giving each
node one input file to process. Then the script will create a file called "report.tar.gz" in your current directory,
containing a QA report on these files. You can download this file using <code>scp</code> 
(see [Moving data to and from your Bioconductor AMI instance](#movingdata)) and unarchive it with the following
command on your local computer:

	tar zxf report.tar.gz

This will create a "report" directory. Inside it is an "index.html" file that you can open with a web browser.


You can also run R interactively with the MPI cluster connected, by using this command:

	mpirun -np 1 --hostfile /usr/local/Rmpi/hostfile R --no-save


When you are finished with your MPI cluster, you can shut down the worker instances and do other housekeeping 
tasks with the following command, issued on your master node:

	mpiutil -a XXX -s YYY -n workers -h -v vol-9999999

Make the following substitutions:

* For "xxx", substitute your Amazon Access Key ID.
* For "yyy", substitute your Amazon Secret Key.
* For "vol-9999999", substitute the volume ID generated by the <code>create_sample_volume</code> script above. You can omit the -v option if your cluster was not started up with it).


**Note**: As always when working with EC2, be sure to shut down all running instances when you are done with them,
to avoid unnecessary charges. You can quickly check instance status on the 
[Instances Page](https://console.aws.amazon.com/ec2/home?region=us-east-1#s=Instances) of the AWS Console.


<a name="parallel"></a>
### Using a parallel cluster in the cloud

You can create easily create a socket cluster of multiple machines
using R's `parallel` package. 

<b><a target="start_ami"
href="https://console.aws.amazon.com/cloudformation/home?region=us-east-1#cstack=sn~StartBioCParallelCluster|turl~https://s3.amazonaws.com/bioc-cloudformation-templates/parallel_cluster.json">Start parallel Cluster</a></b>

You can also start a parallel cluster with SSH access enabled:

<b><a target="start_ami"
href="https://console.aws.amazon.com/cloudformation/home?region=us-east-1#cstack=sn~StartBioCParallelClusterWithSSH|turl~https://s3.amazonaws.com/bioc-cloudformation-templates/parallel_cluster_ssh.json">Start parallel Cluster with ssh access</a></b>

These links prompt you for instance type and number of workers. For information
about instance types, refer to
[Amazon's Instance Type Page](http://aws.amazon.com/ec2/instance-types/).
"Virtual cores" is the number that tells us how many cores are available
for parallel computation on each instance.

Choosing 2 for "number of workers" results in a cluster of 3 machines: 2
workers plus one master. If the instance type chosen is m1.xlarge, then
you will be starting up a cluster of 12 cores (3 machines times 4 cores
each).

After starting the above stack, you'll see a URL in the Outputs tab
where you can log into R on the master node.

A file `/usr/local/Rmpi/hostfile.plain` describes the machines in this 
cluster and how many cores are on each. It might look like this:

    10.68.155.37 4
    10.50.213.89 4
    10.29.191.43 4

So you can do a parallel operation as follows:

    library(parallel)
    lines <- readLines("/usr/local/Rmpi/hostfile.plain")
    hosts <- character()
    for (line in lines)
    {
        x <- (strsplit(line[[1]], " "))
        hosts <-
            c(hosts, rep.int(x[[1]][1], as.integer(x[[1]][2])))
    }
    cl <- makePSOCKcluster(hosts,
        master=system("hostname -i", intern=TRUE))
    system.time(clusterCall(cl, Sys.sleep, 1))

We know this was done in parallel because it took just over one
second:

    user  system elapsed 
    0.004   0.000   1.005 



<a name="custom"></a>
### Creating a custom version of the Bioconductor AMI


*Note*: If you make changes to the running Bioconductor AMI, and then terminate the AMI, your changes will be lost. 
Use the steps described here to ensure that your changes are persistent.

If the AMI is missing some packages or features you think it should have, please let us know.

If you want to customize the AMI for your own purposes, it is simple. Just go ahead and customize your 
running instance as you see fit. Typically this will involve installing R packages with <code>biocLite()</code>,
and software packages (at the operating system level) with the Ubuntu package manager <code>apt-get</code>.

When your instance is customized to your liking, issue the following commands (as root):

	clean_ami
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
	scp -i bioconductor-bob-mylaptop.pem /path/to/myfile root@ec2-50-16-120-30.compute-1.amazonaws.com:/root
	
That will copy the file at "/path/to/myfile" on your local computer to the /root directory on the remote instance. To copy a file from a 
running instance to your local computer, do something like this (still at your local computer):

	scp -i bioconductor-bob-mylaptop.pem root@ec2-50-16-120-30.compute-1.amazonaws.com:/root/myfile /some/directory

That will copy the file /root/myfile from the running instance to /some/directory on your local machine.

*Reminder*: Files created on a running EC2 instance are **not** persisted unless you do some special steps. So if you are 
generating output files with Bioconductor, you must copy them to your local machine before terminating your instance, or
your files will be lost. If you have a lot of data to move back and forth, you may want to look into 
[Elastic Block Storage](http://aws.amazon.com/ebs/).

<a name="questions"></a>
## Questions


If you have questions about the Bioconductor AMI, please contact us through the [Bioconductor mailing list](/help/mailing-list/#bioconductor).


