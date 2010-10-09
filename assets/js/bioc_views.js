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
    
    
    html += "<ul class='inline_list'>\n";
    
    
    for (var i = 0; i < packageList.length; i++) {

        var url = "/help/bioc-views/" + biocVersion + "/" + map[category] + "/html/" + packageList[i] + ".html"
        
        html += "\t<li>\n"
        html += "\t\t<a class='bioc_package' id='pkg_" +
          packageList[i] +
          "' " +
          "href='" +
          url +
          "'>" +
          packageList[i] +
          "</a>\n";
          html += "\t</li>\n";
    }
    html += "</ul>\n"
    jQuery("#packages").html(html);
}


var jumpToAnchor = function() {
    var tmp = ("" + window.location).split("#");
    document.getElementById('treeTop').scrollIntoView(true);
}

var nodeSelected = function(event, data){
    var nodeName = data['args'][0]['text'];
    
    if (nodeName == undefined) {
        nodeName = getNodeName();
    } 
    
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
    } else if (href.match(/\/bioc-views\/([^/]*)\//)) {
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
            "json_data" : {
                "ajax" : {
                    "url" : "/help/bioc-views/json/" + biocVersion + "/tree.json",
                    "data" : function (n) { 
                        return { id : n.attr ? n.attr("id") : 0 }; 
                    }
                }
            },
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


var loadedPackageData = false;

var loadPackageData = function() {
  var repos = ["bioc", "data/annotation", "data/experiment"];
  var count = 0;
  
  for (var i = 0; i < repos.length; i++) {
      jQuery.getJSON("/help/bioc-views/json/" + biocVersion + "/" + repos[i] +  "/packages.json", function(data){
          jQuery.extend(packageInfo, data);
          if (count == 2) {
              loadedPackageData = true;
          }
          count++;
      });
  }
  
}

//document ready function
jQuery(function () {
    setBiocVersion();
    loadPackageData();
    init();
    jQuery(".bioc_package").live("mouseover", function(){
        var title = jQuery(this).attr("title");
        if (title == "" && loadedPackageData) {
            var tmp = jQuery(this).attr("id");
            var id = tmp.replace("pkg_", "");
            var title = packageInfo[id]["Description"].replace(/\n/g, " ");
            jQuery(this).attr("title", title);
        }
    })
});
