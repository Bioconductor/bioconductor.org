# Bioconductor in the cloud

We have developed an Amazon Machine Image (AMI) that is optimized for running Bioconductor in the Amazon Elastic Compute Cloud 
(or EC2) for sequencing tasks.

Here are a few reasons you could use it:

* You do not want to install Bioconductor on your own machine.
* You have a long-running task and you don't want it to tie up the CPU on your own machine.
* You have a parallelizable task and would like to run it (either on multiple CPUs on a single machine, or in a cluster of many machines).

See below for more specific scenarios.

## Preloaded AMI

The AMI comes pre-loaded with R 2.13, and the following Bioconductor (and CRAN) packages (and all their dependencies):

<ul class="inline_list">
	<li>BayesPeak</li>
	<li>baySeq</li>
	<li>ChIPpeakAnno</li>
	<li>chipseq</li>
	<li>ChIPseqR</li>
	<li>ChIPsim</li>
	<li>CSAR</li>
	<li>DEGseq</li>
	<li>DESeq</li>
	<li>edgeR</li>
	<li>gage</li>
	<li>GenomicFeatures</li>
	<li>GenomicRanges</li>
	<li>girafe</li>
	<li>goseq</li>
	<li>MEDIPS</li>
	<li>MotIV</li>
	<li>multicore</li>
	<li>oneChannelGUI</li>
	<li>PICS</li>
	<li>R453Plus1Toolbox</li>
	<li>rGADEM</li>
	<li>Rgraphviz</li>
	<li>Rmpi</li>
	<li>rnaSeqMap</li>
	<li>Rolexa</li>
	<li>Rsamtools</li>
	<li>segmentSeq</li>
	<li>ShortRead</li>
	<li>splots</li>
	<li>SRAdb</li>
</ul>

All the org.\* and BSgenome.\* annotation packages are also loaded.

## How To Use It

### First-time steps

First you will need an [Amazon Web Services](http://aws.amazon.com/) (AWS) account if you do not already have one. Sign up for AWS and then click [here](http://aws-portal.amazon.com/gp/aws/developer/subscription/index.html?productCode=AmazonEC2) to sign up for the EC2 service. (This will require that you provide credit card information).

Now launch the [AWS Console](https://console.aws.amazon.com/ec2/home). 
Click on the [Key Pairs](https://console.aws.amazon.com/ec2/home#s=KeyPairs)
link in the lower left-hand corner of the page. Click the "Create Key Pair" button. When prompted, supply a name.
We suggest that the name be a combination of "bioconductor", your first name, and your machine name
(this will avoid conflicts with other people who share your AWS account, or possibly your own account on another machine).
For example, if your name is Bob and your personal computer is named "mylaptop", your key pair name could be "bioconductor-bob-mylaptop".

Now you will be prompted to download the key pair file you have created. (This is actually an RSA private key).
The suffix ".pem" will be appended to the key pair name you chose earlier.
Be sure and copy this file to a safe place on your hard drive. On Mac and Unix systems at least, you will also need to set
restrictive permissions on the keypair file, as follows:

	chmod 400 bioconductor-bob-mylaptop.pem

## Launching the AMI

Using the [AWS Console](https://console.aws.amazon.com/ec2/home), click the "Launch Instance" button.

Choose the Community AMIs tab. In the text box, paste in the AMI ID of the Bioconductor AMI:

	ami-26b4454f

*Please note that this AMI ID may change over time as we update the underlying AMI. Refer to this page for the most
current AMI ID.*

Click the Select button

The AWS console will now take you through a wizard-like interface: 

<img src="/images/ami/wizard1.jpg" border="0"/>

This screen shows the [different types of instances available](http://aws.amazon.com/ec2/instance-types/). Prices for these
instances vary, so please refer to [current pricing information](http://aws.amazon.com/ec2/pricing/) before launching your 
instance.

Choose an instance type that is appropriate for your task. If you require a lot of memory, choose a high-memory instance type. If your task is CPU-intensive, choose a high-CPU instance type (but be aware that to take advantage of multiple processors, you must explicitly parallelize your CPU-intensive code with the *multicore* or *Rmpi* packages). If you are not concerned about performance, or are just looking around, choose the Micro instance type.

You don't need to choose an Availability Zone unless you are using [EBS (Elastic Block Store) volumes](http://aws.amazon.com/ebs/).

Click the Continue button. On the next screen (Instance Details), you can accept all the defaults and click Continue again. The following 
screen allows you to optionally set some metadata about this instance. For now, just click Continue (though if you share an AWS
account with others, you may want to give your instance a name that will identify it as belonging to you.). On the Create Key pair screen, choose the Key Pair that you created when you first set up your AWS account. Click Continue. On the Configure Firewall screen, choose the "default" security group and click Continue. On the Review screen, make sure all your options look correct, then click Launch.
Then click "View your instances on the Instances page".

In a moment, you will see that your instance is running. You can then put a check box to the left of your running instance, and click on Instance Actions, then Connect. This will show you a dialog box containing a command line you can use to connect to your instance. 
Copy this command line to your clipboard. It will look something like this:

	ssh -i bioconductor-bob-mylaptop.pem root@ec2-50-16-120-30.compute-1.amazonaws.com

*Note* that both the keypair name and machine name here are examples. Your actual command line will look different.

You can vary this command. If you want to use programs or R packages that use X11, be sure and add a -X flag:

	ssh -X -i bioconductor-bob-mylaptop.pem root@ec2-50-16-120-30.compute-1.amazonaws.com

If you do not want to log in as root, you can change "root" to "ubuntu".

Now you can paste your command line into a terminal or Command Prompt. Make sure you are in the same directory as your
key pair file. 

**Windows Users**: You will need to install a version of the *ssh* and *scp* commands. Graphical programs like
[PuTTY](http://www.chiark.greenend.org.uk/~sgtatham/putty/) and [WinSCP](http://winscp.net/eng/index.php) will work. 
Our examples, however, will use the command-line versions of these programs, which you can obtain by installing 
[Cygwin](http://www.cygwin.com/) (be sure and install the *ssh* package).

Once you have pasted this command into your Terminal or Command Prompt window (and pressed Enter) you should be connected
to your Amazon EC2 instance. 

**Important Note**: Be sure and shut down your instance when you no longer need it, or charges will
continue to accrue. You can shut down your instance from the [AWS Console](https://console.aws.amazon.com/ec2/home#s=Instances)
by clicking on Instance Actions, then Terminate.

The following section assumes that you are connected to your EC2 instance. 

## Scenarios for using your Bioconductor instance

### Launching R

Just type R at a command prompt:

	root@ip-10-117-35-246:~# R

	R version 2.13.0 Under development (unstable) (2010-12-17 r53867)
	[...]

	Loading required package: utils
	BioC_mirror = http://www.bioconductor.org
	Change using chooseBioCmirror().
	> 
	
Note that the <code>biocLite()</code> function is automatically available, you do not have to <code>source()</code> a file first
like you would have to on your R installation on your personal computer.

### using Rgraphviz

Make sure you have connected to your instance using the -X flag of the *ssh* command. Something like:

	ssh -X -i bioconductor-bob-mylaptop.pem root@ec2-50-16-120-30.compute-1.amazonaws.com

Then, from within R on the remote instance:

	library("Rgraphviz")
	set.seed(123)
	V <- letters[1:10]
	M <- 1:4
	g1 <- randomGraph(V, M, 0.2)
	plot(g1)
	

This should start a graphics device on your local computer displaying a simple graph.

### Paralellization using multicore

This works best if you have selected a high-CPU instance type to run. 

This trivial example runs the <code>rnorm()</code> function, but any function would work. Consult the 
[multicore](http://cran.r-project.org/web/packages/multicore/index.html) documentation for more information.

	library(multicore)
	mclapply(1:30, rnorm)


### Configuring an MPI cluster in the cloud

You can launch multiple EC2 instances and set them up as an [MPI](http://en.wikipedia.org/wiki/Message_Passing_Interface)
cluster, to parallelize and shorten long-running CPU-intensive jobs. You must explicitly parallelize your
CPU-intensive code using functions in the [Rmpi](http://cran.r-project.org/web/packages/Rmpi/index.html) package.

To configure an EC2-based MPI cluster, do not launch an EC2 instance as described above. Instead, download and build the
<code>RCloudShortRead</code> package. On your own computer (not an EC2 instance), run the following commands:

	svn co https://hedgehog.fhcrc.org/compbio/pilot/CloudShortRead/RCloudShortRead

(When prompted, the username and password are both "readonly".)

	R CMD build RCloudShortRead
	R CMD install RCloudShortRead_0.99.0.tar.gz
	R

then, from within R:

	browseVignettes(package="RCloudShortRead")

Read the <code>RCloudShortRead</code> documentation for further instructions.

### Creating a custom version of the Bioconductor AMI


*Note*: If you make changes to the running Bioconductor AMI, and then terminate the AMI, your changes will be lost. 
Use the steps described here to ensure that your changes are persistent.

If the AMI is missing some packages or features you think it should have, please let us know.

If you want to customize the AMI for your own purposes, it is simple. Just go ahead and customize your 
running instance as you see fit. Typically this will involve installing R packages with <code>biocLite()</code>,
and software packages (at the operating system level) with the Ubuntu package manager <code>apt-get</code>.

When your instance is customized to your liking, use the [AWS Console](https://console.aws.amazon.com/ec2/home#s=Instances)
to Stop your instance (**important note: do NOT "Terminate" your instance; use the Stop command (under Instance Actions) instead.**)

Then choose "Create Image (EBS AMI)" under the Instance Actions menu. You will be prompted for a name for your AMI. After
entering the name, your AMI will be created and given a unique AMI ID. You can then launch this AMI using the steps above, being 
sure to substitute the ID of your own AMI. Your AMI will be private, accessible only to your AWS account, unless you decide
to make it more widely accessible.

Now you should Terminate the Stopped instance of the Bioconductor AMI.

## Moving data to and from your Biocondctor AMI instance

The *scp* command is the most efficient way to move data to and from your EC2 instances.

To copy a file from your computer to a running Bioconductor AMI instance:

* Open a Terminal or Command Prompt window on your computer
* Using "cd", change to the directory where your key pair (.pem) file lives
* Issue a command like the following (your key pair name and the hostname of the AMI instance will be different; you can determine the correct values by clicking on your running instance in the [AWS Console](https://console.aws.amazon.com/ec2/home#s=Instances)):

	
	
<span style="display:none">Hidden</span>
	scp -i bioconductor-bob-mylaptop.pem /path/to/myfile root@ec2-50-16-120-30.compute-1.amazonaws.com/root
	
That will copy the file at "/path/to/myfile" on your local computer to the /root directory on the remote instance. To copy a file from a 
running instance to your local computer, do something like this (still at your local computer):

	scp -i bioconductor-bob-mylaptop.pem root@ec2-50-16-120-30.compute-1.amazonaws.com/root/myfile /some/directory

That will copy the file /root/myfile from the running instance to /some/directory on your local machine.

*Reminder*: Files created on a running EC2 instance are **not** persisted unless you do some special steps. So if you are 
generating output files with Bioconductor, you must copy them to your local machine before terminating your instance, or
your files will be lost. If you have a lot of data to move back and forth, you may want to look into 
[Elastic Block Storage](http://aws.amazon.com/ebs/).

## Questions

If you have questions about the Bioconductor AMI, please (?? contact Dan? or the bioC list??)



