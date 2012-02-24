<?php
    $url = $GET_["auth_url"];
    $last_line = exec("curl " . $url);
    $ret = 'processResults("' . $last_line . '")';
    echo $ret;
?>