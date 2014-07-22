// bioconductor.js


// global variables
var checkForEncryptInterval;
var gPayload;
var gMod;
var gExp;

// logging functions:
var fb_lite = false;
try {
	if (firebug) {
		fb_lite = true;  
		firebug.d.console.cmd.log("initializing firebug logging");
	}
} catch(e) {
	// do nothing
}



function log(message) {
	if (fb_lite) {  
		console.log(message);
	} else {
		if (window.console) {
			console.log(message);
		} 
	}
	if (window.dump) {
	    dump(message + "\n");
	}
}

// convenience functions
String.prototype.trim = function() {
	return this.replace(/^\s+|\s+$/g,"");
}
String.prototype.ltrim = function() {
	return this.replace(/^\s+/,"");
}
String.prototype.rtrim = function() {
	return this.replace(/\s+$/,"");
}

//utility functions
var getParameterByName = function ( name ) {
  name = name.replace(/[\[]/,"\\\[").replace(/[\]]/,"\\\]");
  var regexS = "[\\?&]"+name+"=([^&#]*)";
  var regex = new RegExp( regexS );
  var results = regex.exec( window.location.href );
  if( results == null )
    return "";
  else
    return decodeURIComponent(results[1].replace(/\+/g, " "));
}


// general-use function to add handlers. use like this:
//    if(document.getElementById('ehs.form')){
//      addEvent(document.getElementById('ehs.form'), 'click', handleRadioClick);
//    }
// JS Bible 6th ed.
function addEvent(elem,evtType,func){
  if(elem.addEventListener){ elem.addEventListener(evtType,func,false); }
  else if(elem.attachEvent){ elem.attachEvent("on"+evtType, func); }
  else { elem["on"+evtType] = func; }
}


// parse the page and pick out div's that have a certain class
// and change those into shaded boxes by adding HTML. this inserts
// table code, but that should be transparent to all users.
function renderShadedBoxes(){

  // prepare the HTML to insert into the divs of target class
  var insert1 = '<table cellspacing="0" cellpadding="0" class="sb"><tr><td class="sb1"></td><td class="sb2"></td><td class="sb3"></td></tr><tr><td class="sb4">&nbsp;</td><td class="sb5">';
  var insert2 = '</td><td class="sb6">&nbsp;</td></tr><tr><td class="sb7"></td><td class="sb8"></td><td class="sb9"></td></tr></table>';

  // obtain all the div's of the target class. note that pre-ie7 doesn't return .getAttribute('class') but does return .getAttribute('className') so we check for that specially
  var oDivs = document.getElementsByTagName('div');
  var className = '';
  for(var i=0;i<oDivs.length;i++){
    className = oDivs.item(i).getAttribute('class')||oDivs.item(i).getAttribute('className'); //alert(className);
    if(className&&className.indexOf('shaded_box')>-1){  //alert(i);        
      oDivs.item(i).innerHTML = insert1 + oDivs.item(i).innerHTML + insert2;    
      oDivs.item(i).className=''; // this removes the shaded_box class from the original div so the styling i just made takes over
    }
  }
  
}

// check each page load to see if there is any shaded_box class
addEvent(window,'load',renderShadedBoxes);



// Masthead site navigation. we have five or more site navigation elements
// appearing at page top, and depending upon the current page url, we want
// the corresponding element to be olive and color unchanged at hover. we do this by pattern matching
// on the page url (client side), and turning the corresponding element olive.
// the position of each of the patterns corresponds to the masthead nav element number,
// e.g., the third element, /help/, which is index 3 (option base 1), matches masthead_nav_element_3 
// we use one Array of matching patterns for each element in case one element needs to match more than one patten. 
// examples are shown below, but adjust for your info architecture.
var masthead_nav_elements = Array( Array(/^\/$/, /^\/index\.html$/),Array(/\/install\//, /install\.html/),
	Array(/\/help\//),Array(/\/developers\//),Array(/\/about\//) );
function checkNav(){
  for(var i=0; i<masthead_nav_elements.length; i++){
  for(var j=0; j<masthead_nav_elements[i].length; j++){
    if( masthead_nav_elements[i]&&masthead_nav_elements[i][j] ){ // skips elements that are blank
	  if (window.location.pathname.match(masthead_nav_elements[i][j])) {
        // match at element i. make it olive
        if( document.getElementById('masthead_nav_element_'+(i+1)) ){
            document.getElementById('masthead_nav_element_'+(i+1)).className='masthead_nav_element masthead_nav_element_selected'; 
            return; // matched, so no need to continue checking.
        }
      }
    }
  }}
}
addEvent(window,'load',checkNav);

Object.size = function(obj) {
    var size = 0, key;
    for (key in obj) {
        if (obj.hasOwnProperty(key)) size++;
    }
    return size;
};


if (!Object.keys) {
  Object.keys = Object.keys || function(o) {
    var result = [];
    for(var name in o) {
        if (o.hasOwnProperty(name))
          result.push(name);
    }
    return result;
  };
}

var tidyWorkflows = function() {
  if (jQuery('#workflows').length > 0) {
    var workflows = [];
    jQuery(".workflow").each(function(index) {
      workflows.push(jQuery(this).html());
    });
    jQuery("#workflows_left").html("");
    jQuery("#workflows_right").html("");
    var rands = {};
    while(Object.size(rands) < 4) {
      var rand = Math.floor(Math.random()* workflows.length);
      rands[rand] = -1;
    }
    var i = 0;
    var keys = Object.keys(rands);
    keys = keys.sort();
    for (var key in keys.sort()) {
      var id = (i < 2) ? "#workflows_left" : "#workflows_right";
      html = jQuery(id).html();
      jQuery(id).html(html + "<li>" + workflows[parseInt(keys[i])] + "</li>");
      i++;
    }
  }
}

var unRebaseMirrors = function() {
    if (!(window.mirror === undefined) && mirror == true) {
        var wlh = window.location.href;
        var segs = wlh.split("/");
        var host = wlh.replace(/^http:\/\//i, "").split("/")[0];
        segs.pop();
        var url = segs.join("/");
        if (segs[3] != "packages") {
            host += "/" + segs[3];
        }
        jQuery.each(jQuery(".do_not_rebase a"), function(index, value){
            var href = jQuery(value).attr("href");
            if (!href.match(/^http:/i)) {
                if (href.match(/^\//)) {
                    jQuery(value).attr("href", "http://" + host + href);
                } else {
                    jQuery(value).attr("href", url + "/" + href);
                }
            }
        });
    } 
}

/*
 * The little file server we use for development does not follow symlinks, so see if we are running 
 * that server (assume we are if we are not on port 80) and change URLs tagged with the "symlink"
 * class (e.g. containing "release" or "devel" to point to the actual file.
 */
var getHrefForSymlinks = function(href) {
  if (window.location.port == "") {
    return href;
  } else {
    var releaseRegex = /\/release\//;
    var develRegex  = /\/devel\//;
    if (href.match(releaseRegex)) {
      return href.replace(releaseRegex, "/" + releaseVersion + "/");
    } else if (href.match(develRegex)) {
      return href.replace(develRegex, "/" + develVersion + "/");
    } else {
      return href;
    }
  }
}


var handleCitations = function()
{
    if (jQuery("#bioc_citation").length ) {
        jQuery("#bioc_citation_outer").hide();
        var url = window.location.href;
        url = url.replace("html","citations");
        var segs = url.split("/");
        var pkg = segs.pop();
        pkg = pkg.replace(".html", "")
        segs.push(pkg);
        segs.push("citation.html");
        url = segs.join("/");
        jQuery.ajax({url: url, dataType: "html", 
            success: function(data,textStatus,jqXHR){
                // working around possible R bug?
                data = data.replace(/}. /g, "");
                data = data.replace(/}.</g, '<');
                data = data.replace(/}."/g, '"'); // ' to pacify my editor

                data = data.replace(' (????)', "");
                jQuery("#bioc_citation").html(data);
                jQuery("#bioc_citation_outer").show();
            }, error: function(data, textStatus, jqXHR){
                console.log("error!");
            }
        });
    }

}

//document ready function                                      
jQuery(function() {
    unRebaseMirrors(); // comment this out if there are issues with rebasing
    tidyWorkflows();
    jQuery.each(jQuery(".symlink"), function(index, value){
      var href = jQuery(value).attr("href");
      jQuery(value).attr("href", getHrefForSymlinks(href));
    });
    jQuery(".rpack").tooltip({tip: "#tooltip"});//{ effect: 'slide'});
    handleCitations();
});

// // another document ready function, for try-it-now
// jQuery(function(){
//     if (jQuery("#tryitnow_instance_started").length > 0) {
//         jQuery("#initially_hidden").hide();
//         var dnsName = getParameterByName("dns");
//         var key = getParameterByName("key");
//         var url = "http://" + dnsName + ":8787";
//         var action = url + "/auth-do-sign-in";
//         var link = "../launch?username=ubuntu&password=bioc&url=" + url;
//         link += "&encrypted=";
//         jQuery("#ami_link").attr("href", link);
//         jQuery("#instance_url").html(url);
        
//         var payload, exp, mod;
//         payload = "ubuntu\nbioc";
//         var chunks = key.split(':', 2);
//         exp = chunks[0];
//         mod = chunks[1];
        
//         var encrypted = encrypt(payload, exp, mod);

//         var link = jQuery("#ami_link").attr("href");
//         jQuery("#ami_link").attr("href", link + encrypted);
//         jQuery("#instance_loading").html("");
//         jQuery("#initially_hidden").show();
//     }
    
    
//     if (jQuery("#launch_tryitnow").length > 0) { // is this launch.md?
//         jQuery("#hide_this_stuff").hide();
//         var username = getParameterByName("username");
//         var password = getParameterByName("password");
//         var encrypted = getParameterByName("encrypted");
//         encrypted = encrypted.replace(/ /g, "+");
//         var url = getParameterByName("url");
//         var action = url + "/auth-do-sign-in";
//         jQuery("form").get(1).setAttribute("action", action);
//         document.getElementById("username").value = username;
//         document.getElementById("password").value = password;
//         //todo change this:
//         document.getElementById('persist').value = document.getElementById('staySignedIn').checked ? "1" : "0";
//         document.getElementById('clientPath').value = window.location.pathname;
//         document.getElementById('package').value = encrypted;
//         document.realform.submit();
//     }
    
//     if (jQuery("#captcha_js").length > 0) {
//         jQuery("#captcha_js").html("<script type='text/javascript' src='http://cloud.bioconductor.org:2112/cgi-bin/get_captcha.js'></script>")
//     }
    
// });

var submit_tryitnow = function() {
    jQuery("#tryitnow_button").attr("disabled", "disabled");
    jQuery("#tryitnow_button").attr("value", "Please wait...");
    return(true);   
}

var processCaptchaResults = function(factoryFilename, captchaKey) {
    var s = "http://cloud.bioconductor.org:2112/cgi-bin/display_captcha.jpg?factoryFilename=";
    s += factoryFilename;
    s += "&captchaKey=";
    s += captchaKey;
    jQuery("#captcha_img").attr("src", s);
    jQuery("#captchaKey").attr("value", captchaKey);
    jQuery("#factoryFilename").attr("value", factoryFilename);
}
