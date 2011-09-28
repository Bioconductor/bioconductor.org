<?php



session_start();  // Start the session where the code will be stored.

?>
<html>
<head>
  <title>Bioconductor Guest Posting</title>
</head>

<body>
    
    

<?php
if (empty($_POST)) { ?>
        <script language="javascript">
        function valid(item) {
            var x=document.forms["mailform"][item].value;
            if (x == null || x == "") return false;
            return true;
        }
        
        function validateForm()
        {
            var ret = valid("name") && valid("email") && valid("subject") &&
              valid("body") && valid("code") && valid("sessioninfo");
            if (ret == false) alert("Please fill out all form items.");
            return ret;
        }
        </script>

        <h2>Bioconductor Guest Posting</h2>

        <p>You may post to the <a href="http://bioconductor.org/help/mailing-list/#bioconductor">
        Bioconductor mailing list</a> through this form, without being a subscriber.
        You can also 
        <a href="https://stat.ethz.ch/mailman/listinfo/bioconductor">subscribe</a> to the list.
        You can <a href="https://stat.ethz.ch/pipermail/bioconductor/">read the list online</a>.</P>

        <p>Please read our <a href="http://bioconductor.org/help/mailing-list/posting-guide/">
            posting guide</a> before posting. It's especially important
            to include the output of <i>sessionInfo()</i> when reporting a problem.</p>
            
        <p>We do not keep your email address. Your email address will be added to the
            "To:" line in the email that is sent, which means that people who respond
            may "Reply All" and you will get their responses. There is no guarantee 
            they will respond this way, so you should watch 
            <a href="https://stat.ethz.ch/pipermail/bioconductor/">the online list
                archive</a> for responses.</p>
        
        <p><a href="http://bioconductor.org">[Return to Bioconductor web site]</a></p>


<form name="mailform" method="POST" onsubmit="return validateForm()">
Your name:<br/>
<input width="30" type="text" name="name" /><br/>
Your email address:<br/>
<input width="30 type="text" name="email" /><br/><br/>
Subject:<br/>
<input width="50" type="text" name="subject" /><br/><br/>

Email:<br/>
<textarea cols="80" rows="10" name="body"></textarea><br/><br/>

Output of <i>sessionInfo()</i>:<br/>
<textarea cols="80" rows="5" name="sessioninfo"></textarea><br/><br/>


<div style="width: 430px; float: left; height: 90px">
      <img id="siimage" align="left" style="padding-right: 5px; border: 0" src="securimage_show.php?sid=<?php echo md5(time()) ?>" />

        <object classid="clsid:d27cdb6e-ae6d-11cf-96b8-444553540000" codebase="http://download.macromedia.com/pub/shockwave/cabs/flash/swflash.cab#version=9,0,0,0" width="19" height="19" id="SecurImage_as3" align="middle">
			    <param name="allowScriptAccess" value="sameDomain" />
			    <param name="allowFullScreen" value="false" />
			    <param name="movie" value="securimage_play.swf?audio=securimage_play.php&bgColor1=#777&bgColor2=#fff&iconColor=#000&roundedCorner=5" />
			    <param name="quality" value="high" />
			
			    <param name="bgcolor" value="#ffffff" />
			    <embed src="securimage_play.swf?audio=securimage_play.php&bgColor1=#777&bgColor2=#fff&iconColor=#000&roundedCorner=5" quality="high" bgcolor="#ffffff" width="19" height="19" name="SecurImage_as3" align="middle" allowScriptAccess="sameDomain" allowFullScreen="false" type="application/x-shockwave-flash" pluginspage="http://www.macromedia.com/go/getflashplayer" />
			  </object>

        <br/>
        
        <!-- pass a session id to the query string of the script to prevent ie caching -->
        <a tabindex="-1" style="border-style: none" href="#" title="Refresh Image" onclick="document.getElementById('siimage').src = 'securimage_show.php?sid=' + Math.random(); return false"><img src="images/refresh.gif" alt="Reload Image" border="0" onclick="this.blur()" align="bottom" /></a>
</div>
<div style="clear: both"></div>
Code:<br/>

<!-- NOTE: the "name" attribute is "code" so that $img->check($_POST['code']) will check the submitted form field -->
<input type="text" name="code" size="12" /><br/><br/>

<input type="submit" value="Post to List" />
</form>

<?php
} else { //form is posted
  include("securimage.php");
  $img = new Securimage();
  $valid = $img->check($_POST['code']);

  if($valid == true) {
    //echo "<center>Thanks, you entered the correct code.<br/>Click <a href=\"{$_SERVER['PHP_SELF']}\">here</a> to go back.</center>";
    echo "<center>Thank you. Your email will be posted.</center>\n";
    echo "<center><a href='http://bioconductor.org'>Return to Bioconductor Site</a></center>\n";
    $listemail = "bioconductor@r-project.org";
    $guestemail = "guest@bioconductor.org";
    $sender = $_POST['name'] . " [guest] <" . $guestemail . ">";
    $mailheaders = "From: " . $sender . "\n";
    $to = $listemail . ", " . $_POST['email'];
    $subject = $_POST['subject'];
    $body = $_POST['body'] .
     "\n\n -- output of sessionInfo(): \n\n" .
     $_POST['sessioninfo'] .
     "\n\n--\nSent via the guest posting facility at bioconductor.org.";
    $result = mail($to, $subject, $body, $mailheaders);
    
  } else {
    echo "<center>Sorry, the code you entered was invalid.  <a href=\"javascript:history.go(-1)\">Go back</a> to try again.</center>";
  }
}

?>

</body>
</html>
