<?php
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
              $_POST['sessioninfo'] == "") {
          header("Location: http://bioconductor.org/help/mailing-list/mailform/missing_items");
      } else {
          $listemail =  "bioconductor@r-project.org";
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
          header("Location: http://bioconductor.org/help/mailing-list/mailform/ok");
      }
    } else {
        header("Location: http://bioconductor.org/help/mailing-list/mailform/badcaptcha");
    }
}
?>
