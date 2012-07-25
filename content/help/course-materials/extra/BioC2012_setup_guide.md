# BioC2012 Participant Setup Guide

Please read this carefully. In order to participate in workshops at the
conference, you'll need to follow some steps which take a few moments.
It's best if you can do this before the conference begins.

## Logistics

You will receive a participant packet when you arrive that contains
complete information about BioC2012.

Limited parking (15 spots) is available in spots marked "Visitor Parking". Check
in with reception when you arrive or you may receive a ticket. 
On-street metered parking is also available.

### Wi-Fi

Wireless access will be available through the network "FHCRC Guest".
The wireless username and password will be written on all conference
whiteboards. When you use your browser after connecting to this 
network you should be prompted to enter this username and password.


### Schedule

A full schedule is in your participant packet which you will
receive upon arrival.

Conference proceedings begin:

* Monday (Developer Day) at 9:00 AM in M1 A303-307 (Arnold Building)
* Tuesday at 8:15 AM in Pelton Auditorium (Weintraub Building)
* Wednesday at 8:30 AM in M1 A303-307 (Arnold Building)

[FHCRC Campus Map](http://www.fhcrc.org/en/contact-us/visit-us.html)

### Links

* [BioC2012 Page](https://secure.bioconductor.org/BioC2012/)
* [BioC2012 Developer Day wiki](http://wiki.fhcrc.org/bioc/Seattle_Dev_Meeting_2012)


## Amazon Web Services Virtual Machine

We will be using Amazon Web Services to provide each participant
with a virtual machine known as an Amazon Machine Image or AMI.
This will allow you to access R through a web browser (using RStudio Server).
This instance of R is pre-configured with necessary packages and workshop materials.

## How much does it cost?

[Amazon Web Services](http://aws.amazon.com/) (AWS) is generously providing
each BioC2012 participant with USD$100 in services. This is more than enough
to pay for the virtual machine you will use during the conference.

*You will be required to provide credit card information* in order to receive
an AWS account. If you use more than $100 worth of services, your credit card
will be charged! This is unlikely to happen, but it is possible if you forget
to turn off any AWS resources (instances, stacks, etc.) that you have turned
on. 

We will discuss how to do that both in this guide, in person at the conference,
and in a followup email to be sent after the conference. 

Provided that you turn off all resources when done with them, and do not
exceed $100 in charges, your AWS usage during the conference will be *free*.

After you have used up the initial $100 from your coupon, you
are free to continue using Amazon Web Services at your own expense.

To give an example of costs, leaving an m1.large instance on for
five days (120 hours) will cost $38.40.

* [AWS EC2 Pricing Guide](http://aws.amazon.com/ec2/pricing/)
* [AWS Pricing Calculator](http://calculator.s3.amazonaws.com/calc5.html)


## Creating an AWS account

You must create an AWS account if you haven't already.

A [useful tutorial](http://www.slideshare.net/simone.brunozzi/amazon-web-services-signup)
is available, though the appearance of some of the pages has changed
since the tutorial was written.

To create your account, visit the [AWS web site](http://aws.amazon.com/).
Click the "Sign Up Now" button.

<img src="signupnow.jpg"/>

If you already have an [amazon.com](http://www.amazon.com/) account,
enter the email address and password associated with that account,
and choose "I am a returning user". Otherwise, enter your email address
and check "I am a new user".

Part of the creation of a new account will involve an identity
verification by telephone, so make sure you have access to a
telephone where you can receive calls.

<img src="callmenow.jpg"/>

You'll receive an email to confirm the creation of your AWS account.

<img src="activating.jpg"/>


## Redeem your coupon code

You received a coupon code via email. Visit
the [coupon redemption page](http://aws.amazon.com/awscredits/)
and enter that code. You should see credits added to your 
account. If you don't do this step, *your credit card will be charged*.

<img src="redeem.jpg"/>

## Start up the BioC2012 AMI

When you are ready to start working with the conference materials,
You can use this link to <a target="BioC2012AMI" href="https://console.aws.amazon.com/cloudformation/home?region=us-east-1#cstack=sn~BioC2012|turl~https://s3.amazonaws.com/bioc-cloudformation-templates/BioC2012.json">
start the BioC2012 AMI</a>.

<img src="createstack.jpg"/>


Click "Continue".

You will be prompted to select the Instance Type of the virtual
machine you will be launching. The m1.large type is pre-selected
and that should be adequate for most needs during the conference. 
Refer to [the complete list of instance types](http://aws.amazon.com/ec2/instance-types/)
to see what other types are available.

Be sure and check the box that says "I acknowledge that this template
may create IAM resources."

<img src="ack_iam.jpg"/>

Click "Continue". 

Review the stack creation information.

<img src="createstack3.jpg"/>


Click "Continue".

Now the stack is being created.

<img src="createstack4.jpg">

Click "Close".

Check the box next to "BioC2012:"

<img src="cf1.jpg"/>

Click the Refresh button at upper right...

<img src="refresh.jpg">

until the stack creation is complete:

<img src="create_complete.jpg">

Then click on the Outputs tab at the bottom of the page:

<img src="outputs.jpg">

You will see a URL for RStudio Server. Click this URL. 


When Rstudio is up, you'll see a login form:

<img src="rstudiosignin.jpg"/>

You can use the
supplied login and password (*workshop* and *bioc*)
to log in to RStudio Server. (If you are a 
workshop presenter, you may also use the login
you received via email).



## RStudio Server tips
* *Browser tip:* RStudio Server does not work well with Internet Explorer.
  We recommend using Firefox, Chrome, or Safari. 
* You can edit a file by clicking on it in the File pane at lower right,
  or by using the `file.edit()` command.
* When a file is in the edit window, you can execute the current line
  or selection by clicking "Run".
* Interrupt a long-running task with the ESC key instead of Control-C.
* `browseVignettes('pkg')` does not work as it should; instead, use 
  `help(package='pkg')` and then click on "Overview of user guides and
  package vignettes" in the Help pane.
* You can upload or download files from the File pane at lower right.
  You should download any files that you wish to keep after the conference.
* Tools/Shell will open a Unix command shell. This is useful for
  running commands like `R CMD check`.
* To navigate the File pane to a specific directory, click the
  triple-dot button in the top right corner of the File Pane:

<img src="threedots.jpg"/>

<a name="deletion"></a>
## Turning off your AWS resources

When you are finished using RStudio, you should delete your 
virtual machine as follows.

From the
[Cloud Formation page](https://console.aws.amazon.com/cloudformation/home),
make sure the "Bioc2012" stack is checked and then click "Delete Stack":

<img src="delete.jpg">

Confirm the deletion by clicking "Yes, Delete".

If you do not do this step, **your credit card will be charged.**

After deletion, all files that were on your virtual machine will
be lost. If there are files you want to keep, use RStudio's
Download button to download them before shutting down
your virtual machine.

## Questions

If you have a question, please talk to a Bioconductor team member
during the course, or email Dan Tenenbaum (dtenenba at fhcrc dot org).
