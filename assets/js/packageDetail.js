
jQuery(function () {
    
    log("hello");
    jQuery.ajax({
        url: "../packages.json",
        dataType: 'json',
        success: function(data, textStatus, xmlHttpRequest) {
            log("success!");
            packageName = getParameterByName("packageName");
            log("package name = " + packageName);
            var package = data[packageName];
            log("repos root = " + package["reposRoot"]);
            
            jQuery("#packageName").html("Package Detail - " + package["Package"]);
            
            var html = "";
            
            
            
            for (item in package) {
                if (package[item].split) {
                    log(item + " = " + package[item]);
                    
                    html += item + ": " + package[item] + "<br/>\n";
                    
                    
                } else {
                    log(item + ":")
                    html += item + ": ";
                    html += package[item].join(", ");
                    for(var i = 0; i < package[item].length; i++) {
                        log("\t" + package[item][i]);
                    }
                }
                
                jQuery("#packageDetail").html(html);
            }
            
            
            
        },
        error: function(xmlHttpRequest, textStatus, errorThrown) {
            log("failure");
            log("error = " + errorThrown);
            log("text status = " + textStatus);
        }
    });
    
    /*
    jQuery.getJSON("/bioc_views/packageDetail/short.json", function(data, textStatus){
        log("textStatus = " + textStatus);
        packageName = getParameterByName("packageName");
        alert("packageName = " + packageName);  
    });
    */
});
