<?php
    $url = $GET_["url"];
    echo("url = " . $url . "<br/>\n");
    $last_line = exec("/usr/bin/curl " . $url,  $output, $result_code);
    echo("result code = " . $result_code . "<br/>\n");
    $ret = 'processResults("' . $last_line . '")';
    echo $ret;
?>