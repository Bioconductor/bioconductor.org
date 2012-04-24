var packageInfo = {};
var biocVersion;

var displayPackages = function(packageList, nodeName) {
    if (packageList == null) {
        jQuery("#packages").empty();
        return;
    }
    var html = "<h3>Packages</h3>\n";
    
    var parents = findParents(nodeName);

    var category = parents[0];

    var map = {"Software": "bioc", "AnnotationData": "data/annotation", "ExperimentData": "data/experiment"};

    html += "<p>";
    html += parents.join(" &gt ");
    html += "</p>\n";
    
    html += "<table><tr><th>Package</th><th>Maintainer</th><th>Title</th></tr>\n";
    
    
    var tableData = "";
    for (var i = 0; i < packageList.length; i++) {
        var rowClass = (i % 2 == 0) ? "row_odd" : "row_even";
        var pkg = packageList[i];
        var url = getHostUrl() + "/" + biocVersion + "/" + map[category] + "/html/" + pkg + ".html"
        tableData += '<tr class="'+rowClass+'" id="pkg_' + pkg + '">\n';
        tableData += '\t<td><a href="'+url+'">'+pkg+'</a></td>\n';
        var cleanMaintainer = packageInfo[pkg]["Maintainer"].replace(/ +<[^>]*>/g, "");
        tableData += '\t<td>'+cleanMaintainer+'</td>\n';
        tableData += '\t<td>'+packageInfo[pkg]["Title"]+'</td>\n';
        tableData += "</tr>\n";
    }
    html += tableData;
    
    html += "</table>\n";
    
    jQuery("#packages").html(html);
    
}


var jumpToAnchor = function() {
    var tmp = ("" + window.location).split("#");
    document.getElementById('treeTop').scrollIntoView(true);
}

var nodeSelected = function(event, data){
    var nodeName;
    
    // for IE
    nodeName = data['args'][0]['innerText'];
    
    if (data['args'][0]['text'] != undefined) {
        nodeName = data['args'][0]['text'];
    }
    
    
    
    if (nodeName == undefined) {
        nodeName = getNodeName();
    } 

    nodeName = nodeName.trim();
    var bareNodeName = nodeName.split(" ")[0];
    var wl = ("" + window.location.href).split("?")[0];
    wl = wl.split("#")[0];

    var newUrl = "" + wl + "#" + "___" + bareNodeName;
    window.location.href = newUrl;
    
    var tmp = nodeName.split(" ");
    nodeName = tmp[0];
      var packageListStr = jQuery("#" + nodeName).attr("packageList");
      if (packageListStr) {
          var packageList = packageListStr.split(",");
          displayPackages(packageList, nodeName);
      } else {
          displayPackages(null, nodeName);
      }
      //jumpToAnchor();
}


var setBiocVersion = function() {
    
    var href = window.location.href;


    if (href.match(new RegExp("/" + releaseVersion + "/"))) {
        biocVersion = releaseVersion;
    } else if (href.match(new RegExp("/" + develVersion + "/"))) {
        biocVersion = develVersion;
    } else if (href.toLowerCase().match(new RegExp("/release/"))) {
        biocVersion = releaseVersion;
    } else if (href.toLowerCase().match(new RegExp("/devel/"))) {
        biocVersion = develVersion;
    } else if (href.match(/\/packages\/([^/]*)\//)) {
        biocVersion = RegExp.$1;
    }
    
    
    
    var versionText;
    
    if (biocVersion == releaseVersion) {
        versionText = "Release"
    } else if (biocVersion == develVersion) {
        versionText = "Development";
    } else {
      versionText ="";
    }
    
    var versionString = biocVersion;
    if (versionText != "") {
        versionString += " (" + versionText + ")";
    }
    
    
    jQuery("#bioc_version").html(versionString);
}


var findParents = function (nodeId) {
    var parents = jQuery("#" + nodeId).parentsUntil("#tree");
    var ret = [];
    ret.push(nodeId);
    jQuery.each(parents, function(index, value){
        var id = jQuery(value).attr("id");
        if (id.length > 0) {
            ret.push(id);
        }
    });
    ret.reverse();
    return ret;
}

var getNodeName = function() {
    var wlh = window.location.href;
    var segs = [];
    segs = wlh.split("#");
    if (segs.length == 2) {
        return segs[1].replace("___", "");
    } else {
        return "";
    }
}

var init = function() {
    // todo add ajax failure method (possible?)
    
    var initiallySelected = [];
    var nodeName = "";
    nodeName = getNodeName();
    
    if (nodeName != "") {
        initiallySelected.push(nodeName);
    }
    
    jQuery("#tree").jstree({ 
        "core": {
            "animation": 0
        },
        "ui": {
          "initially_select": initiallySelected
        },
            "themes": {
                "theme": "apple",
                "dots": false,
                "icons": false
            },
            "json_data": dataTree,
            "plugins" : [ "themes", "json_data", "ui" ]
        });

        // explicitly add biocViewsTree class because the widget strips it off
        jQuery("#tree").addClass("biocViewsTree");

    /*
    jQuery("#tree").bind("close_node.jstree dblclick.jstree delete_node.jstree deselect_node.jstree destroy.jstree drag_start.vakata drag_stop.vakata get_rollback.jstree init.jstree load_node.jstree loaded.jstree mousedown.jstree move_node.jstree open_node.jstree rename_node.jstree reopen.jstree select_node.jstree set_rollback.jstree ", function (event, data) {
        log("event name: " + event.type);
    })
    */
    
    jQuery("#tree").bind("select_node.jstree", function(event, data){
        nodeSelected(event, data);
    });
    
    jQuery("#tree").bind("loaded.jstree", function(event, data){
        var initiallyOpen = [];
        var openNode = getNodeName();

        if (openNode != "") {
            initiallyOpen = findParents(openNode);
            for(var i = 0; i < initiallyOpen.length; i++) {
                jQuery("#tree").jstree("open_node", "#" + initiallyOpen[i]);
            }
        }
        
        
    });
    
    
}



var getHostUrl = function() {
    var url = window.location.href;
    var segs = url.split("/");
    segs.pop();
    segs.pop();
    url = segs.join("/");
    return(url);
}


var loadPackageData = function() {
    if (typeof bioc_packages != "undefined")
        jQuery.extend(packageInfo, bioc_packages)
    if (typeof data_annotation_packages != "undefined")
        jQuery.extend(packageInfo, data_annotation_packages)
    if (typeof data_experiment_packages != "undefined")
        jQuery.extend(packageInfo, data_experiment_packages)
}


//document ready function
jQuery(function () {
    setBiocVersion();
    loadPackageData();
    init();
});
