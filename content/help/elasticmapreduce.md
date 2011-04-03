# Using Bioconductor with Amazon Elastic MapReduce

Bioconductor can be used with [Amazon Elastic MapReduce](http://aws.amazon.com/elasticmapreduce/) to facilitate robust parallelization of 
computational tasks. 

The following tutorial shows how it can be done.

## Prerequisites

You must have an [Amazon Web Services](http://aws.amazon.com/) account and have signed up for the [Elastic MapReduce](http://aws.amazon.com/elasticmapreduce/),
 [Amazon SimpleDB](http://aws.amazon.com/simpledb/) services, and [Simple Storage Service](http://aws.amazon.com/s3/) or S3.

### Key Pair

You need to generate a key pair using the [AWS Console](https://console.aws.amazon.com/ec2/home#s=KeyPairs). Click "Create Key Pair"
and give the key pair a name. Download the resulting file to a safe place on your hard drive.

On Mac and Unix systems you will need to alter the permissions of the key pair file as follows:

	chmod 0400 mykeypair.pem

You only need to do this if you intend to ssh to your Elastic MapReduce instances (recommended for troubleshooting).

## Streaming MapReduce Tutorial

In this tutorial we will use the <code>qa()</code> function from Bioconductor's <code>ShortRead</code> to perform quality
assessment on short read data. (Note that the <code>qa()</code> function also supports other types of parallelization, using the 
R packages <code>multicore</code> and <code>Rmpi</code>. If you would like to use <code>qa()</code> in this way, refer to the
[Bioconductor Cloud AMI Documentation](/help/bioconductor-cloud-ami)).

### Streaming MapReduce

A streaming MapReduce job flow is a task in which the mapper and reducer can be written in any programming language.
Input data is streamed to the standard input of the mapper, which prints its output to standard output. That output
is then streamed to the reducer, which again writes its output to standard output. Input and output data locations are specified
as S3 "buckets". 

This paradigm is ideal for tasks in which the input and output data are textual. In our case, the data we want to examine are
BAM files (though we could use any file type supported by the <code>qa()</code> function--see the [ShortRead 
documentation](/packages/devel/bioc/html/ShortRead.html)),
and the output data is a QA report (a compressed tar file containing html, jpg, and pdf files). 
So we cheat a little bit by specifying lists of BAM files as our input data, and then let our mapper download the actual BAM files
from S3. Our input files are [here](http://bioconductor-mapreduce-example-inputdir.s3.amazonaws.com/file1) and 
[here](http://bioconductor-mapreduce-example-inputdir.s3.amazonaws.com/file2). (You may want to use the <code>curl</code>
utility to view these files, because a web browser will attempt to download them.) 
They are simply text files, each of which contains the name
of a single BAM file, which can be found in an S3 bucket called "bioconductor-bamfiles". (Note that this bucket is a public bucket 
maintained by the Bioconductor team; other parts of this tutorial require that you provide your own buckets, and bucket names must
be unique across all of S3. We'll make it clear when you need to provide your own bucket names).

If you wanted to run this tutorial on your own data, you could put more than one filename in each input file, but 
having a single name in each file ensures that each file will get its own mapper (assuming, of course, that you specify
a matching number of instances when you start your job flow.)

(Note: To better browse S3 buckets, we recommend using the [S3 Console](https://console.aws.amazon.com/s3/home), or
a third party software tool such as [s3cmd](http://s3tools.org/s3cmd).)

Another issue is that Elastic MapReduce uses Amazon Machine Images (AMIs) which are very generic and contain outdated software. Luckily 
Elastic MapReduce provides [Bootstrap Actions](http://aws.typepad.com/aws/2010/04/new-elastic-mapreduce-feature-bootstrap-actions.html),
a feature that allows us to install the latest version of R, the <code>ShortRead</code> package and its dependencies, 
and the s3cmd program for moving files to and from S3. Our bootstrapping script is available
[here](http://bioconductor-emr-bootstrap-scripts.s3.amazonaws.com/bootstrap.sh).

#### Starting the MapReduce Job Flow

There are two ways to launch a MapReduce job flow: through the 
[graphical web console](https://console.aws.amazon.com/elasticmapreduce/home)
or from a command-line utility that can be downloaded [here](http://aws.amazon.com/developertools/2264?_encoding=UTF8&jiveRedirect=1).

First we will discuss the graphical method, and then show how to start the same job flow from the command line.

##### Starting a Job Flow with the Elastic MapReduce Console

Visit the [Elastic MapReduce Console](https://console.aws.amazon.com/elasticmapreduce/home) and click the
"Create New Job Flow" button. Fill out the screen as follows.

<img src="/images/elasticmapreduce/screenshot1.jpg">

Click "Continue". 

Fill out the screen as follows, with the changes noted below:

<img src="/images/elasticmapreduce/screenshot2.jpg">

**Note**: The output location should NOT be what is shown here, but should be the name of a bucket that does not exist yet. Bucket names
must be unique across all of S3, so it can be useful to include a uniquifying string--such as your birthdate--in the bucket name. 
For example, "mybucket-19710223". Bucket names should be all lowercase, contain only letters, numbers and the dash ("-") symbol.
Elastic MapReduce will fail if the output bucket already exists. Use the [S3 Console](https://console.aws.amazon.com/s3/home)
to delete your output bucket if it already exists, before creating your new job flow.

Click "Continue".

Fill out the screen as follows, with the changes noted below:

<img src="/images/elasticmapreduce/screenshot3.jpg">

**Note**: For "Key Pair", use the name of the Key Pair you created in the Prerequisites step. You can see a list of key pairs on the 
[Key Pairs page](https://console.aws.amazon.com/ec2/home#s=KeyPairs).

**Note**: For "Amazon S3 Log Path", choose the name of a bucket you own. It should exist but be empty. See above for 
bucket naming guidelines. You can use the [S3 Console](https://console.aws.amazon.com/s3/home) to create your bucket.
When filling out the form, put "s3n://" in front of the bucket name.

Click "Continue".

Fill out the screen as follows:

<img src="/images/elasticmapreduce/screenshot4.jpg">

Instead of "XXXXXXXXXXXXXXXXXXXX" and "YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY", supply your Amazon 
Access Key ID and Secret Access Key. These are available from the 
[Security Credentials](https://aws-portal.amazon.com/gp/aws/developer/account/index.html?ie=UTF8&action=access-key&openid.assoc_handle=aws&aToken=4|qHyux74J3Pqo5Ml1%2BpIjsf8w8twCjBb4%2BNaiZTTTknBLi8S3ow/WjtPevRc4tiDz3fblfbOshZIXKunrQfXdZJ%2BN/Zba9jpj0PevsDiIj3Oi13k8IHP9YmKnm9jRfy3kRJrWQ2H3OWkTXQz%2BKwtcpqHMvCLPpPtZ/lonhjnU7f5Wyb8TkaQP4aO6FzdQbJ7yt56GzP6NgKqa82QdI%2BCoKxXbfZ/W45ID&openid.claimed_id=https://www.amazon.com/ap/id/amzn1.account.AEM2GUZSOFV2GFYCOAGKQWTMUZ4Q&openid.identity=https://www.amazon.com/ap/id/amzn1.account.AEM2GUZSOFV2GFYCOAGKQWTMUZ4Q&openid.mode=id_res&openid.ns=http://specs.openid.net/auth/2.0&openid.op_endpoint=https://www.amazon.com/ap/signin&openid.response_nonce=2011-01-18T21:28:22Z3420106691568049774&openid.return_to=https://aws-portal.amazon.com/gp/aws/developer/account/index.html%3Fie%3DUTF8%26action%3Daccess-key&openid.signed=assoc_handle,aToken,claimed_id,identity,mode,ns,op_endpoint,response_nonce,return_to,pape.auth_policies,pape.auth_time,ns.pape,signed&openid.ns.pape=http://specs.openid.net/extensions/pape/1.0&openid.pape.auth_policies=http://schemas.openid.net/pape/policies/2007/06/none&openid.pape.auth_time=2011-01-18T21:28:22Z&openid.sig=GFj1mg6Rt5h5mV3DxQzpheGbgJjvCGoE4Wce4Acw350%3D&) page.
These identifiers are necessary so that the mapper and reducer can communicate with S3. 

The last two arguments in the "Optional Arguments" box should be the names of S3 buckets which exist and which you own.
(We earlier told you that the output dir bucket should not exist, but it will have been created once this part of the job flow
is reached). Here we are using the same bucket for both Streaming MapReduce output and the "real" output of our reducer (the report.tar.gz file), but you are not required to do the same thing.

Click "Continue". Your job flow summary should look something like this:

<img src="/images/elasticmapreduce/screenshot5.jpg">

If it all looks correct, click "Create Job Flow".

Then click on "View Job Flows" to see the progress of your job flow.

If all goes well, your job will eventually complete and a file called report.tar.gz will be in your 
output bucket on S3. You can use the [S3 Console](https://console.aws.amazon.com/s3/home) to download
this file. Then you can use the following command to unarchive the file:

	tar zxf report.tar.gz

This will create some files in a "report" directory. You can open "index.html" with your web browser
to see the report generated by ShortRead and Elastic MapReduce.

##### Starting a Job Flow with the command line

Download and install the 
[Elastic MapReduce command-line utility](http://aws.amazon.com/developertools/2264?_encoding=UTF8&jiveRedirect=1).
This utility is written in [Ruby](http://www.ruby-lang.org/en/) and requires that the Ruby Language be installed on your computer.
It is installed by default on Mac OS X machines and on many Linux machines. 

**Note**: The command line utility works best with Ruby 1.8.6 or 1.8.7. It fails to work with newer versions such as 
Ruby 1.9.2.

Once the command line utility is installed, you can start a job flow with a command like the following:

	elastic-mapreduce --create --name "ShortRead QA" --num-instances 2 \
	--slave-instance-type m1.small --master-instance-type m1.small \
	--key-pair gsg-keypair --stream --input s3n://bioconductor-mapreduce-example-inputdir \
	--output s3n://outdir19671025 --mapper s3n://bioconductor-mapreduce-example/mapper-emr.R \
	--reducer s3n://bioconductor-mapreduce-example/reducer-emr.R \
	--jobconf mapred.reduce.tasks=1 \
	--bootstrap-action s3://bioconductor-emr-bootstrap-scripts/bootstrap.sh \
	--bootstrap-name "Custom Boostrapping Action" --arg XXXXXXXXXXXXXXXXXXXX \
	--arg YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY --arg bioconductor-bamfiles \
	--arg tempdir19671025 --arg outdir19671025 --debug -v \
	--log-uri s3n://emr-debugging-19671025 --enable-debugging

It is important to change the values of the following flags:

* <code>--key-pair</code>: Change this to the name of the key pair you create in the Prerequisites step.
* <code>--output</code>: Change this to the name of a bucket that does not yet exist, and put "s3n://" in front of the bucket name.
* <code>--arg XXXXXXXXXXXXXXXXXXXX</code>: Substitute your Amazon Access Key ID.
* <code>--arg YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY</code>: Substitute your Secret Access Key.
* <code>--arg tempdir19671025</code>: Supply the name of a bucket you own which exists
* <code>--arg outdir19671025</code>: Supply the name of a bucket you own which exists, OR the same bucket name you are using 
with the <code>--output</code> flag, which should not exist.
* <code>--log-uri</code>: Change this to the name of a bucket which you own, which exists, and is empty, with "s3n://" in front of the name.

### How it works

The critical part of our [mapper](http://s3.amazonaws.com/bioconductor-mapreduce-example/mapper-emr.R) is ShortRead's <code>qa()</code>
function. While this function can handle multiple 
files, the best parallelization is achieved when each mapper handles just a single file.

Our [reducer](http://s3.amazonaws.com/bioconductor-mapreduce-example/reducer-emr.R) is simply R's <code>rbind()</code> function. This requires that all intermediate results (the qa objects generated
by each of the mappers) be in the same place. We specify a single reducer with the "--jobconf mapred.reduce.tasks=1" flag. 
In the reducer, all 
of the intermediate qa objects are read into a list, <code>rbind()</code> is called on the list, and the report is generated
from the results of the <code>rbind()</code>.

### Support

If you run into any issues, contact us through the [Bioconductor mailing list](/help/mailing-list/#bioconductor).
