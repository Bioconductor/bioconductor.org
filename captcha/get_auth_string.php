<?php
    $url = $GET_["url"];
    $last_line = exec("curl " . $url);
    $ret = 'processResults("' . $last_line . '")';
    echo $ret;
?>