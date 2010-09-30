// bioconductor.js

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
var masthead_nav_elements = Array( Array('home.html'),Array('/install/','install.html'),Array('/help/'),Array('/developers/'),Array('/about/') );
function checkNav(){
  for(var i=0; i<masthead_nav_elements.length; i++){
  for(var j=0; j<masthead_nav_elements[i].length; j++){
    if( masthead_nav_elements[i]&&masthead_nav_elements[i][j] ){ // skips elements that are blank
      if( window.location.href.indexOf(masthead_nav_elements[i][j])>-1 ){
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


var getCorrectUrlForMirrors = function() {
    url = window.location.href.replace(/^http:\/\//i, "");
    segs = url.split("/");
    host = segs[0];
    jQuery(".site_host").html(host);
}

/*
 * The little file server we use for development does not follow symlinks, so see if we are running 
 * that server (assume we are if we are not on port 80) and change URLs tagged with the "symlink"
 * class (e.g. containing "release" or "devel" to point to the actual file.
 */
var getHrefForSymlinks = function(href) {
  if (window.location.port == 80) {
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

//document ready function                                      
jQuery(function() {
    getCorrectUrlForMirrors();
    
    jQuery.each(jQuery(".symlink"), function(index, value){
      var href = jQuery(value).attr("href");
      jQuery(value).attr("href", getHrefForSymlinks(href));
    });
    
    
});
