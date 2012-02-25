<?php
    $url = $_GET["url"];
    $last_line = exec("/usr/bin/curl --connect-timeout 1 --retry 20 --retry-delay 1 " . $url,  $output, $result_code);
    $ret = 'processResults("' . $last_line . '")';
    echo $ret;
?>