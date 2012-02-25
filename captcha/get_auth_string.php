<?php
    $url = $_GET["url"];
    $last_line = exec("/usr/bin/curl " . $url,  $output, $result_code);
    $ret = 'processResults("' . $last_line . '")';
    echo $ret;
?>