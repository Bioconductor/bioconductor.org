// -------------------- global variables -----------------------

var curTopic="code";
var curChap="";
var curFig="";

// ------------------- end global variables --------------------




// -------------------- helper functions -----------------------

// show/hide arbitrary HTML block elements
function toggleBlockViz(obj){
    var hide = "none";
    var show = "block";
    var state = obj.style.display;
    if(state == hide)
	obj.style.display=show;
    else
	obj.style.display=hide;
}

function revToggleBlockViz(obj){
    var hide = "none";
    var show = "block";
    var state = obj.style.display;
    if(state == show)
	obj.style.display=hide;
    else
	obj.style.display=show;
}

// show/hide table rows
function toggleRowViz(obj){
    var hide = "none";
    var show = "table-row";
    if(navigator.appName == "Microsoft Internet Explorer")
        show = "inline";
    var state = obj.style.display;
    if(state == hide)
	obj.style.display=show;
    else
	obj.style.display=hide;
}

// change src of an iFrame
function linkToFrame(frame, src){
    parent.document.getElementById(frame).src=src;
}

// ------------------- end helper functions --------------------




// ------------------ global event handlers  -------------------

// event for click on 'Figures' button 
function clickFigures(){
    if(curTopic != "figs"){
	var chap;
	if(curChap=="")
	    chap="defaultFig.html";
	else
	    chap="chapter_"+curChap+"/figures/figure_"+curFig;
	document.getElementById("menuFrame").src="allFigs.html";
	document.getElementById("mainFrame").src=chap;
	document.getElementById("main").className="main";
	document.getElementById("codeTop").className="topic";
	document.getElementById("figTop").className="topicSel";
	document.getElementById("solTop").className="topic";
	curTopic="figs";
    }
}


// event for click on 'Code' button 
function clickCode(){
    if(curTopic != "code"){
	var chap;
	if(curChap=="")
	    chap="defaultChap.html";
	else
	    chap="chapter_"+curChap;
	document.getElementById("menuFrame").src="allCode.html";
	document.getElementById("mainFrame").src=chap;
	document.getElementById("main").className="main";
	document.getElementById("codeTop").className="topicSel";
	document.getElementById("figTop").className="topic";
	document.getElementById("solTop").className="topic";
	curTopic="code";
    }
}


// event for click on 'Solutions' button 
function clickSolutions(){
    if(curTopic != "sols"){
	var sol;
	if(curChap=="")
	    sol="defaultSol.html";
	else
	    sol="chapter_"+curChap+"/solution.pdf";
	document.getElementById("menuFrame").src="allSolutions.html";
	document.getElementById("mainFrame").src=sol;
	document.getElementById("main").className="mainNoPadding";
	document.getElementById("codeTop").className="topic";
	document.getElementById("figTop").className="topic";
	document.getElementById("solTop").className="topicSel";
	curTopic="sols";
    }
}

// ---------------- end global event handlers  -----------------


// event for click on a navigation item
function clickNavigation()
{
    toggleBlockViz(document.getElementById("chap0"));
    obj = document.getElementById("chapterList0");
    if(obj.className=="chapterListNavSel")
	obj.className="chapterListNav";
    else
	obj.className="chapterListNavSel";
}

function clickSubNav(i, len, url)
{
    for (var j = 1; j <= len; j++){
     	document.getElementById("nav"+j).className="navList";
      if(navigator.appName == "Microsoft Internet Explorer")
          document.getElementById("nav"+j).style.backgroundColor="#c9d8e5";
    }
    document.getElementById("nav"+i).className="navListSel";
    if(navigator.appName == "Microsoft Internet Explorer")
        document.getElementById("nav"+i).style.backgroundColor="#7896b1";
    linkToFrame("mainFrame", url);
}

// ---------------- end navigation event handlers  -----------------


// -------------- event handlers 'code' mode  ------------------

// event for click on a chapter when in 'code' mode
function clickChapterCode(i, len)
{
    for (var j = 1; j <= len; j++)
	document.getElementById("chapterList"+j).className="chapterList";
    var chap="chapter_"+i;
    linkToFrame('mainFrame', chap+"/code.html");
    curChap=i;
    document.getElementById("chapterList"+i).className="chapterListSel";
}

// ------------ end event handlers 'code' mode  ----------------




// -------------- event handlers 'figure' mode  ------------------

// event for click on a chapter when in 'figure' mode
function clickChapterFig(i, len)
{
    for (var j = 1; j <= len; j++){
	document.getElementById("chapterList"+j).className="chapterList";
	document.getElementById("chap"+j).style.display="none";
	document.getElementById("chapTit"+j).style.display="block";
    }
    var chap="chapter_"+i;
    if(curChap!=i){
	curChap=i;
	document.getElementById("chapterList"+i).className="chapterListSel";
	toggleBlockViz(document.getElementById("chap"+i));
	document.getElementById("chapTit"+i).style.display="none";
    }else{
	document.getElementById("chap"+i).style.display="none";
	document.getElementById("chapTit"+i).style.display="block";
    }
}

function clickFigure(i, len){
    for (var j = 1; j <= len; j++)
	document.getElementById("fig"+curChap+"_"+j).className="figureList";
    var fig="figure_"+i;
    linkToFrame("mainFrame", "chapter_"+curChap+"/figures/"+fig+".html");
    curFig = i;
    document.getElementById("fig"+curChap+"_"+i).className="figureListSel";
}

// ------------ end event handlers 'figure' mode  ----------------



// -------------- event handlers 'solution' mode  ------------------

// event for click on a chapter when in 'solution' mode
function clickChapterSol(i, len)
{
    for (var j = 1; j <= len; j++)
	document.getElementById("chapterList"+j).className="chapterList";
    var chap="chapter_"+i;
    linkToFrame('mainFrame', chap+"/solutions.pdf");
    curChap=i;
    document.getElementById("chapterList"+i).className="chapterListSel";
}

// ------------ end event handlers 'solution' mode  ----------------




// show/hide individual chunks
function toggleChunk(i){
    // show/hide code chunk
    ch = document.getElementById("chunk"+i);
    toggleRowViz(ch);
    // switch control icon
    toggleBlockViz(document.getElementById("hideChunk"+i));
    toggleBlockViz(document.getElementById("showChunk"+i));
    toggleBlockViz(document.getElementById("solChunk"+i));
    ct = document.getElementById("chunkCont"+i);
   
    if(ct.className == "codeContainer" ||
       ct.className == "highlightedCodeContainer"){
	ct.className = "hiddenCodeContainer";
    }else{
	if(ch.className=="chunk")
	    ct.className = "codeContainer";
	else
	    ct.className="highlightedCodeContainer";
    }
}



// highlight individual chunks
function highlightChunk(i){
    ct = document.getElementById("chunkCont"+i);
    if(ct.className != "hiddenCodeContainer"){
	if(ct.className == "codeContainer")
	    ct.className = "highlightedCodeContainer";
	else
	    ct.className = "codeContainer";
    }
}


// enable mouseover hovering in the iframe for IE
function hover(me)
{
    if(navigator.appName == "Microsoft Internet Explorer"){
        me.style.backgroundColor = "#FFC233";
        me.style.cursor = "pointer";
    }
}

function unhover(me)
{
    if(navigator.appName == "Microsoft Internet Explorer"){
        if(me.className.search("Sel") != -1)
            me.style.backgroundColor = "#7896b1";
        else
	       me.style.backgroundColor = "transparent";
    }
}




