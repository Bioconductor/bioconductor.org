// bioconductor.js
if (!/\.html$|\/$|#/.test(window.location.href))
  window.location.href = window.location.href + "/";

// global variables
var checkForEncryptInterval;
var gPayload;
var gMod;
var gExp;

// logging functions:
var fb_lite = false;
try {
  // if (firebug) {
  // 	fb_lite = true;
  // 	firebug.d.console.cmd.log("initializing firebug logging");
  // }
} catch (e) {
  // do nothing
}

//Sitehead code that will underline the nav element if the page URL matches the element

const nav_elements = [/^\/about\//, /^\/developers\//, /^\/help\//];

function checkNav() {
  const currentPath = window.location.pathname;
  const navLinks = document.querySelectorAll(".nav-links a");

  navLinks.forEach((link, index) => {
    if (nav_elements[index].test(currentPath)) {
      link.classList.add("active");
    }
  });
}

window.addEventListener("load", checkNav);

//Mobile siteHead js
document.addEventListener("DOMContentLoaded", function () {
  findHeaderTop();
  const hamburger = document.querySelector(".hamburger");
  const navMenu = document.querySelector(".header-nav");

  function mobileMenu() {
    hamburger.classList.toggle("active");
    navMenu.classList.toggle("active");
  }

  hamburger.addEventListener("click", mobileMenu);

  const navLink = document.querySelectorAll(".mobile-link");

  function closeMenu() {
    hamburger.classList.remove("active");
    navMenu.classList.remove("active");
  }

  navLink.forEach((n) => n.addEventListener("click", closeMenu));

  window.addEventListener( "resize", findHeaderTop);

});

function findHeaderTop() {

  const announcementHeight = document.querySelector(".announcement")?.offsetHeight;
  const header = document.querySelector(".site-masthead");

  if(announcementHeight){
    header.style.top = `-${announcementHeight}px`
  }

}

//Changes body and hero background color once clicked on certain links
function changeBackgroundColors() {
  const heroElement = document.querySelector(".hero");
  const installPageRegex = /\/install\//;
  const aboutPageRegex = /\/about\//;

  if (installPageRegex.test(window.location.href)) {
    document.body.style.backgroundColor = "#fff";
    heroElement.style.backgroundColor = "var(--neutral-n50)";
  }

  if (aboutPageRegex.test(window.location.href)) {
    document.body.style.backgroundColor = "#fff";
  }
}

changeBackgroundColors();

window.addEventListener("hashchange", changeBackgroundColors);

//Allow the cursor to drag and scroll the events
function enableDragScroll() {
  const slider = document.querySelector(".events-container");
  let isDown = false;
  let startX;
  let scrollLeft;

  if(slider){

  slider.addEventListener("mousedown", (e) => {
    isDown = true;
    slider.classList.add("active");
    startX = e.pageX - slider.offsetLeft;
    scrollLeft = slider.scrollLeft;
  });

  function handleMouseUpAndLeave() {
    isDown = false;
    slider.classList.remove("active");
  }

  slider.addEventListener("mouseup", handleMouseUpAndLeave);
  slider.addEventListener("mouseleave", handleMouseUpAndLeave);

  slider.addEventListener("mousemove", (e) => {
    if (!isDown) return;
    e.preventDefault();
    const x = e.pageX - slider.offsetLeft;
    const scrollSpeed = x - startX;
    slider.scrollLeft = scrollLeft - scrollSpeed;
  });
}
}

window.addEventListener("load", enableDragScroll);

function log(message) {
  if (fb_lite) {
    //console.log(message);
  } else {
    // if (window.console) {
    // 	console.log(message);
    // }
  }
  if (window.dump) {
    dump(message + "\n");
  }
}

// convenience functions
String.prototype.trim = function () {
  return this.replace(/^\s+|\s+$/g, "");
};
String.prototype.ltrim = function () {
  return this.replace(/^\s+/, "");
};
String.prototype.rtrim = function () {
  return this.replace(/\s+$/, "");
};

//utility functions
var getParameterByName = function (name) {
  name = name.replace(/[\[]/, "\\[").replace(/[\]]/, "\\]");
  var regexS = "[\\?&]" + name + "=([^&#]*)";
  var regex = new RegExp(regexS);
  var results = regex.exec(window.location.href);
  if (results == null) return "";
  else return decodeURIComponent(results[1].replace(/\+/g, " "));
};

// general-use function to add handlers. use like this:
//    if(document.getElementById('ehs.form')){
//      addEvent(document.getElementById('ehs.form'), 'click', handleRadioClick);
//    }
// JS Bible 6th ed.
function addEvent(elem, evtType, func) {
  if (elem.addEventListener) {
    elem.addEventListener(evtType, func, false);
  } else if (elem.attachEvent) {
    elem.attachEvent("on" + evtType, func);
  } else {
    elem["on" + evtType] = func;
  }
}

Object.size = function (obj) {
  var size = 0,
    key;
  for (key in obj) {
    if (obj.hasOwnProperty(key)) size++;
  }
  return size;
};

if (!Object.keys) {
  Object.keys =
    Object.keys ||
    function (o) {
      var result = [];
      for (var name in o) {
        if (o.hasOwnProperty(name)) result.push(name);
      }
      return result;
    };
}

var tidyWorkflows = function () {
  if (jQuery("#workflows").length > 0) {
    var workflows = [];
    jQuery(".workflow").each(function (index) {
      workflows.push(jQuery(this).html());
    });
    jQuery("#workflows_left").html("");
    jQuery("#workflows_right").html("");
    var rands = {};
    while (Object.size(rands) < 4) {
      var rand = Math.floor(Math.random() * workflows.length);
      rands[rand] = -1;
    }
    var i = 0;
    var keys = Object.keys(rands);
    keys = keys.sort();
    for (var key in keys.sort()) {
      var id = i < 2 ? "#workflows_left" : "#workflows_right";
      html = jQuery(id).html();
      jQuery(id).html(html + "<li>" + workflows[parseInt(keys[i])] + "</li>");
      i++;
    }
  }
};

var unRebaseMirrors = function () {
  if (!(window.mirror === undefined) && mirror == true) {
    var wlh = window.location.href;
    var segs = wlh.split("/");
    var host = wlh.replace(/^http:\/\//i, "").split("/")[0];
    segs.pop();
    var url = segs.join("/");
    if (segs[3] != "packages") {
      host += "/" + segs[3];
    }
    jQuery.each(jQuery(".do_not_rebase a"), function (index, value) {
      var href = jQuery(value).attr("href");
      if (!href.match(/^http:/i)) {
        if (href.match(/^\//)) {
          jQuery(value).attr("href", "http://" + host + href);
        } else if (href.match(/^#/)) {
          jQuery(value).attr("href", window.location.href + href);
        } else {
          jQuery(value).attr("href", url + "/" + href);
        }
      }
    });
  }
};

/*
 * The little file server we use for development does not follow symlinks, so see if we are running
 * that server (assume we are if we are not on port 80) and change URLs tagged with the "symlink"
 * class (e.g. containing "release" or "devel" to point to the actual file.
 */
var getHrefForSymlinks = function (href) {
  if (window.location.port == "") {
    return href;
  } else {
    var releaseRegex = /\/release\//;
    var develRegex = /\/devel\//;
    if (href.match(releaseRegex)) {
      return href.replace(releaseRegex, "/" + releaseVersion + "/");
    } else if (href.match(develRegex)) {
      return href.replace(develRegex, "/" + develVersion + "/");
    } else {
      return href;
    }
  }
};

var handleCitations = function () {
  if (jQuery("#bioc_citation").length) {
    jQuery("#bioc_citation_outer").hide();
    var url = window.location.href;
    url = url.replace("html", "citations");
    var segs = url.split("/");
    var pkg = segs.pop();
    pkg = pkg.replace(".html", "");
    segs.push(pkg);
    segs.push("citation.html");
    url = segs.join("/");
    jQuery.ajax({
      url: url,
      dataType: "html",
      success: function (data, textStatus, jqXHR) {
        // working around possible R bug?
        data = data.replace(/}. /g, "");
        data = data.replace(/}.</g, "<");
        data = data.replace(/}."/g, '"'); // ' to pacify my editor

        data = data.replace(" (????)", "");
        jQuery("#bioc_citation").html(data);
        jQuery("#bioc_citation_outer").show();
      },
      error: function (data, textStatus, jqXHR) {
        //console.log("error!");
      },
    });
  }
};

//document ready function
jQuery(function () {
  unRebaseMirrors(); // comment this out if there are issues with rebasing
  tidyWorkflows();
  jQuery.each(jQuery(".symlink"), function (index, value) {
    var href = jQuery(value).attr("href");
    jQuery(value).attr("href", getHrefForSymlinks(href));
  });
  jQuery(".rpack").tooltip({ tip: "#tooltip" }); //{ effect: 'slide'});
  handleCitations();
});

var submit_tryitnow = function () {
  jQuery("#tryitnow_button").attr("disabled", "disabled");
  jQuery("#tryitnow_button").attr("value", "Please wait...");
  return true;
};

var processCaptchaResults = function (factoryFilename, captchaKey) {
  var s =
    "http://cloud.bioconductor.org:2112/cgi-bin/display_captcha.jpg?factoryFilename=";
  s += factoryFilename;
  s += "&captchaKey=";
  s += captchaKey;
  jQuery("#captcha_img").attr("src", s);
  jQuery("#captchaKey").attr("value", captchaKey);
  jQuery("#factoryFilename").attr("value", factoryFilename);
};
