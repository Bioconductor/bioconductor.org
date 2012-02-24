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
        $py = "/home/webadmin/python/bin/python";
        $script = "/extra/www/event_reg/mailform/start_instance.py";
        $cmd = $py . " " . $script . " " . $ami_id;
        echo("cmd = " . $cmd . "\n<br/>");
        $last_line = exec($cmd, $output, $result_code);
        echo("last_line = " . $last_line . "\n<br/>");
        echo("output = " . $output . "\n<br/>");
        echo("result_code = " . $result_code . "\n<br/>");
        //header("Location: http://bioconductor.org/help/cloud/started?dns=" . $last_line);
    } else {
        header("Location: http://bioconductor.org/help/cloud/badcaptcha");
    }
}
?>
