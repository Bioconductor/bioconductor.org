# Bioconductor in the cloud

We have developed an Amazon Machine Image (AMI) that is optimized for running Bioconductor in the cloud for sequencing tasks.

Here are a few reasons you could use it:

* You do not want to install Bioconductor on your own machine.
* You have a long-running task and you don't want it to tie up the CPU on your own machine.
* You have a parallelizable task and would like to run it (either on multiple CPUs on a single machine, in a cluster of many machines).

## Preloaded AMI

The AMI comes pre-loaded with R 2.13, and the following Bioconductor packages (and all their dependencies):

<ul class="inline_list">
	<li>BayesPeak</li>
	<li>CSAR</li>
	<li>ChIPpeakAnno</li>
	<li>ChIPseqR</li>
	<li>ChIPsim</li>
	<li>DEGseq</li>
	<li>DESeq</li>
	<li>GenomicFeatures</li>
	<li>GenomicRanges</li>
	<li>multicore</li>
	<li>MEDIPS</li>
	<li>MotIV</li>
	<li>PICS</li>
	<li>R453Plus1Toolbox</li>
	<li>Rolexa</li>
	<li>Rmpi</li>
	<li>Rsamtools</li>
	<li>SRAdb</li>
	<li>ShortRead</li>
	<li>baySeq</li>
	<li>chipseq</li>
	<li>edgeR</li>
	<li>gage</li>
	<li>girafe</li>
	<li>goseq</li>
	<li>oneChannelGUI</li>
	<li>rGADEM</li>
	<li>Rgraphviz</li>
	<li>rnaSeqMap</li>
	<li>segmentSeq</li>
	<li>splots</li>
</ul>

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


