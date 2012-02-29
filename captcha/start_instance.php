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
        $py = "/home/webadmin/python/bin/python";
        $script = "/extra/www/event_reg/mailform/start_instance.py";
        $cmd = $py . " " . $script . " " . $ami_id;
        $last_line = exec($cmd, $output, $result_code);
        $segs = explode(";", $last_line);
        $dns = $segs[0];
        $key = $segs[1];
        header("Location: http://bioconductor.org/help/cloud/started?dns=" . $dns ."&key=" . $key);
    } else {
        header("Location: http://bioconductor.org/help/cloud/badcaptcha");
    }
}
?>
