<?php
session_start();  // Start the session where the code will be stored.
if (empty($_POST)) {
    header("Location: http://bioconductor.org/help/cloud/badcaptcha");
} else {
    include("securimage.php");
    $img = new Securimage();
    $valid = $img->check($_POST['code']);
    $ami_id = $_POST['ami_id'];

    if($valid == true) {
        $output = "";
        $result_code = -1;
        $last_line = exec("/home/webadmin/python/bin/python /extra/www/event_reg/mailform/start_instance.py " + $ami_id, $output, $result_code);
        header("Location: http://bioconductor.org/help/cloud/started?dns=" + $last_line);
    } else {
        header("Location: http://bioconductor.org/help/cloud/badcaptcha");
    }
}
?>
