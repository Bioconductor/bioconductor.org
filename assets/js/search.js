// search.js
// don't use jQuery.noConflict(), it conflicts(!) with other
// document.ready functions

// TODO: add previous and next buttons at top (as well as bottom)
// TODO: change page title to include search term

String.prototype.trim = function() {
	return this.replace(/^\s+|\s+$/g,"");
}
String.prototype.ltrim = function() {
	return this.replace(/^\s+/,"");
}
String.prototype.rtrim = function() {
	return this.replace(/\s+$/,"");
}


jQuery(function() {
	initSearch();
});

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

var getSearchUrl = function(query, start) {
	var url = "/solr/select?indent=on&version=2.2&q=" + query + 
	"&fq=&start=" + start +  "&rows=10&fl=id,score,title&qt=standard&wt=json&explainOther=&hl=on&hl.fl=&hl.fragsize=200";
	return url;
}

var initSearch = function() {
	var q = getParameterByName("q");
	
	if (q == "") {
		jQuery("#q").focus();
	} else {
		jQuery("#q").val(q);
	}
	
	jQuery("#if_search_results_present").hide();
	
	var startParam = getParameterByName("start");
	var start = (startParam == "") ? 0 : parseInt(startParam);
	
	
	var url = getSearchUrl(q, start);
	jQuery.getJSON(url, function(data){
		var numFound = data['response']['numFound'];
		jQuery("#numFound").html(numFound);
		jQuery("#search_query").html(q);
		jQuery("#if_search_results_present").show();
		
		if (numFound > 0) {
			var rows = parseInt(data['responseHeader']['params']['rows']);
			var docs = data['response']['docs'];
			var highlighting = data['highlighting'];
			for (var i = 0; i < docs.length; i++) {
				var doc = docs[i];
				var outer = highlighting[doc.id];
				var text = outer['text'];
				var snippet = text[0]; // does this array ever contain more than one element?
				var isHTML = /\/$|\.html/i.test(doc.id); 
				var googleAnalytics = " onClick=\"javascript: pageTracker._trackPageview('" + doc.id + "'); \" ";
            	
				var isR = /\.R$/.test(doc.id);
				var title = (""  + doc.title).trim();
				if (title == "") {
					title = "Untitled";
				}
				var stringToAppend = "";
				stringToAppend += "<p><a ";
				if (!isHTML) {
				    stringToAppend += googleAnalytics;
				}
				stringToAppend += "href='" +
				  doc.id +
				  "'>" +
				  title +
				  "</a> - " +
				  doc.id +
				  "<blockquote>";
				if (isR) {
					stringToAppend += "<pre>"
				}
				stringToAppend += snippet;
				if (isR) {
					stringToAppend += "</pre>";
				}
				stringToAppend += "</blockquote></p>\n";
				
				jQuery("#search_results").append(stringToAppend);
			}
			
			if (start >= rows) {
				var nextStart = start - rows;
				url = "/help/search/index.html?q=" + q + "&start=" + nextStart;
				jQuery("#previous_search_page").html("&lt; <a href='"+url+"'>Previous</a> ");
			}
			
			if (numFound > (start + rows)) {
				var prevStart = start + rows;
				url = "/help/search/index.html?q=" + q + "&start=" + prevStart;
				jQuery("#next_search_page").html(" <a href='"+url+"'>Next </a> &gt;");
			}
		}
	});
	
}

