// search.js
// don't use jQuery.noConflict(), it conflicts(!) with other
// document.ready functions

// TODO: add previous and next buttons at top (as well as bottom)
// TODO: change page title to include search term


var start;
var q;

jQuery(function() {
	initSearch();
});


var getSearchUrl = function(query, start) {
	var url = "//master.bioconductor.org/solr/default/select?indent=on&version=2.2&q=" + query +
	// boost release and devel landing page paths, and fuzzy title matches
	" id:*\\/release\\/bioc\\/html*^1.2 id:*\\/devel\\/bioc\\/html*^1.1 title:" + query + "~" +
	"&fq=&start=" + start +
	"&rows=20&fl=id,score,title&qt=standard&wt=json&explainOther=&hl=on&hl.fl=&hl.fragsize=200";
	return url;
};

var searchResponse = function(data) {
		var numFound = data['response']['numFound'];
		jQuery("#numFound").html(numFound);
		jQuery("#search_query").html(q);
		jQuery("#if_search_results_present").show();
		
		if (numFound == 0)
			$("#search_results").html("");


		if (numFound > 0) {
			var rows = parseInt(data['responseHeader']['params']['rows']);
			var docs = data['response']['docs'];
			var highlighting = data['highlighting'];
			var stringToAppend = "";
			stringToAppend += "<dl>\n";
			for (var i = 0; i < docs.length; i++) {
				var doc = docs[i];
				var outer = highlighting[doc.id];
				var text = outer['text'];
				var snippet = "&nbsp;";
				if (text != undefined) {
    				snippet = text[0]; // does this array ever contain more than one element?
				}
				var isHTML = /\/$|\.html/i.test(doc.id); 
				var googleAnalytics = " onClick=\"javascript: pageTracker._trackPageview('" + doc.id + "'); \" ";
            	
				var isR = /\.R$/.test(doc.id);
				var title = (""  + doc.title).trim();
				if (title == "") {
					title = "Untitled";
				}
				stringToAppend += "<dt><a ";
				if (!isHTML) {
				    stringToAppend += googleAnalytics;
				}
				stringToAppend += "href='" +
				  doc.id +
				  "'>" +
				  title +
				  "</a> - " +
				  doc.id +
				  "</dt><dd>";
				if (isR) {
					stringToAppend += "<pre>";
				}
				stringToAppend += snippet;
				if (isR) {
					stringToAppend += "</pre>";
				}
				stringToAppend += "</dd>\n";
				
				
			}
			stringToAppend += "</dl>\n";
			jQuery("#search_results").html(stringToAppend);
			
			if (start >= rows) {
				var nextStart = start - rows;
				url = "/help/search/index.html?q=" + q + "&start=" + nextStart;
				jQuery(".previous_search_page").html("&lt; <a href='"+url+"'>Previous</a> ");
			}
			
			if (numFound > (start + rows)) {
				var prevStart = start + rows;
				url = "/help/search/index.html?q=" + q + "&start=" + prevStart;
				jQuery(".next_search_page").html(" <a href='"+url+"'>Next </a> &gt;");
			}
		}
};

var ajaxErrHandler = function(jqXHR, textStatus, errorThrown) {
	$("#search_results").html("A timeout or invalid search term resulted in an error.");
};

var initSearch = function() {
	q = getParameterByName("search-bar");
	q = q.replace(/\/$/, "");


	if (q == "") {
		jQuery("#search-bar").focus();
	} else {
		$("#search_results").html("Searching....");
		jQuery("#search-bar").val(q);
	}
	
	jQuery("#if_search_results_present").hide();
	
	var startParam = getParameterByName("start");
	start = (startParam == "") ? 0 : parseInt(startParam);
	
	
	var url = getSearchUrl(q, start);
	jQuery.ajax({'url': url, 'data': null, 'success': searchResponse,
		'timeout': 5000, 'error': ajaxErrHandler,
		'dataType': 'jsonp', 'jsonp': 'json.wrf'});
	
};
