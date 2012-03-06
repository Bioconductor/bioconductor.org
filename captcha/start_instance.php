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
        echo("Hello to you!<br/>\n");
        flush();
        $py = "/home/webadmin/python/bin/python";
        $script = "/extra/www/event_reg/mailform/start_instance.py";
        $cmd = $py . " " . $script . " " . $ami_id;
        $last_line = exec($cmd, $output, $result_code);
        $segs = explode(";", $output[0]);
        $dns = $segs[0];
        $key = $segs[1];
        ///*
        echo("size of output: " . count($output) . "<br>\n");
        echo("last_line = " . $last_line . "<br>\n");
        echo("output = " . $output . "<br>\n");
        echo("dns = " . $dns . "<br>\n");
        echo("key = " . $key . "<br>\n");
        foreach ($output as $item) {
            echo("array item: " . $item . "<br>\n");
        }
        //*/
        
        //$refreshtag =  "BLAHmeta http-equiv='refresh' content='0;url=http://bioconductor.org/help/cloud/started?dns='" . $dns . "&key=" . $key . "'>\n"; 
        //echo($refreshtag);
        //header("Location: http://bioconductor.org/help/cloud/started?dns=" . $dns ."&key=" . $key);
        //echo("tryitnowCallbackFunc('" . $dns . "', '" . $key . "');");
    } else {
        header("Location: http://bioconductor.org/help/cloud/badcaptcha");
    }
}
?>
