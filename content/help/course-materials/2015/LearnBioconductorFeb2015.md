[Course Package](LearnBioconductor_0.1.10.tar.gz) for use with R-3.2.x and Bioconductor 3.1.

<h4 id="ami">To launch an Amazon Machine Image (AMI) for this course:</h4>

* [Create an Amazon Web Services (AWS) Account](https://aws.amazon.com/) if you
  don't already have one.
* Start the instance <%= ami_url("ami-ded79ab6") %>; See the [documentation for this](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/launching-instance.html). Make sure your [security group](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/using-network-security.html) has
port 80 accessible.
* Paste the Public DNS name of your instance into a web browser.
* Log in to RStudio with username *ubuntu* and password *bioc* .
* Be sure and [terminate your instance](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/terminating-instances.html) when you are done using it, in order to avoid excessive charges.

