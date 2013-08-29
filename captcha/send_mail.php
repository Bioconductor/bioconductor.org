<?php
/**
Validate an email address.
Provide email address (raw input)
Returns true if the email address has the email 
address format and the domain exists.
*/
function validEmail($email)
{
   $isValid = true;
   $atIndex = strrpos($email, "@");
   if (is_bool($atIndex) && !$atIndex)
   {
      $isValid = false;
   }
   else
   {
      $domain = substr($email, $atIndex+1);
      $local = substr($email, 0, $atIndex);
      $localLen = strlen($local);
      $domainLen = strlen($domain);
      if ($localLen < 1 || $localLen > 64)
      {
         // local part length exceeded
         $isValid = false;
      }
      else if ($domainLen < 1 || $domainLen > 255)
      {
         // domain part length exceeded
         $isValid = false;
      }
      else if ($local[0] == '.' || $local[$localLen-1] == '.')
      {
         // local part starts or ends with '.'
         $isValid = false;
      }
      else if (preg_match('/\\.\\./', $local))
      {
         // local part has two consecutive dots
         $isValid = false;
      }
      else if (!preg_match('/^[A-Za-z0-9\\-\\.]+$/', $domain))
      {
         // character not valid in domain part
         $isValid = false;
      }
      else if (preg_match('/\\.\\./', $domain))
      {
         // domain part has two consecutive dots
         $isValid = false;
      }
      else if
(!preg_match('/^(\\\\.|[A-Za-z0-9!#%&`_=\\/$\'*+?^{}|~.-])+$/',
                 str_replace("\\\\","",$local)))
      {
         // character not valid in local part unless 
         // local part is quoted
         if (!preg_match('/^"(\\\\"|[^"])+"$/',
             str_replace("\\\\","",$local)))
         {
            $isValid = false;
         }
      }
      if ($isValid && !(checkdnsrr($domain,"MX") || checkdnsrr($domain,"A")))
      {
         // domain not found in DNS
         $isValid = false;
      }
   }
   return $isValid;
}


session_start();  // Start the session where the code will be stored.
if (empty($_POST)) {
    # don't do anything
} else {
    include("securimage.php");
    $img = new Securimage();
    $valid = $img->check($_POST['code']);

    if($valid == true) {
      if ($_POST['subject'] == "" ||
              $_POST['name'] == "" || 
              $_POST['email'] == "" ||
              $_POST['body'] == "" || 
              $_POST['code'] == "" || 
              $_POST['sessioninfo'] == "" ||
              !validEmail($_POST['email'])) {
          header("Location: http://bioconductor.org/help/mailing-list/mailform/missing_items/");
      } else {
          $body = stripslashes($_POST['body'] .
           "\n\n -- output of sessionInfo(): \n\n" .
           $_POST['sessioninfo'] .
           "\n\n--\nSent via the guest posting facility at bioconductor.org.");
          $cc = str_replace(" at ", "@", $_POST['cc']);
          
          if ($_POST['subject'] == "testignore") {
              $listemail =  "dtenenba@fhcrc.org";
              $body = "[CC: " . $cc . "]\n\n" . $body;
              $cc = "";
          } else {
              $listemail =  "bioconductor@r-project.org";
          }
          
          $guestemail = "guest@bioconductor.org";
          $sender = $_POST['name'] . " [guest] <" . $guestemail . ">";
          $mailheaders = "From: " . $sender . "\n";
          if (strlen($cc) > 0) {
              $mailheaders .= "Cc: " . $cc . "\n";
          }
          $to = $listemail . ", " . $_POST['email'];
          $subject = $_POST['subject'];
          $result = mail($to, $subject, $body, $mailheaders);
          // add sender to list of people who can post to list without being subscribed
          $cmd = '/home/webadmin/python/bin/python /extra/www/event_reg/mailform/update_posters_list.py ' . $_POST['email'];
          $last_line = exec($cmd, $output, $result);
          // TODO - we may want to check the value of $result and see if there was an error.
          header("Location: http://bioconductor.org/help/mailing-list/mailform/ok/");
      }
    } else {
        header("Location: http://bioconductor.org/help/mailing-list/mailform/badcaptcha/");
    }
}
?>
