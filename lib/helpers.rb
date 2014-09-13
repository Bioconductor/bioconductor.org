include Nanoc3::Helpers::Text
include Nanoc3::Helpers::Rendering
include Nanoc3::Helpers::Breadcrumbs
include Nanoc3::Helpers::XMLSitemap
include Nanoc3::Helpers::HTMLEscape


# coding: utf-8
require 'htmlentities'
require 'time'
require 'date'
require 'rubygems'
require 'httparty'
require 'yaml'
require 'pp'
require 'rexml/document'
require 'json'
require 'open-uri'
require 'nokogiri'
require 'fileutils'
require 'mechanize'
require 'kramdown'

include REXML


class Time
  def to_date
    Date.new(year, month, day)
  rescue NameError
    nil
  end
end


class Date
  def to_time
    Time.local(year, month, day)
  end
end

def myfunc()
  puts "in myfunc"
  "this is the output of myfunc() " + Date.new.to_s
end

class RowIndexer
  @@rownum = 0
  def rowclass()
    @@rownum += 1
    (@@rownum % 2 == 1) ? "row_odd" : "row_even"
  end

end

def get_cran_packages()
  puts "Grabbing list of CRAN packages..."
  ## Is there a non-web-scraping way to get a list of CRAN packages?
  cran_packages = []
  data = begin
           HTTParty.get("http://cran.fhcrc.org/web/packages/available_packages_by_name.html",
                        :timeout => 6)
         rescue Timeout::Error
           puts "Timeout grabbing list of CRAN packages..."
           return []
         rescue TCPSocket::SocketError
           puts "Socket error grabbing list of CRAN packages..."
           return []
         end
  html = data.to_s
  lines = html.split("\n")
  for line in lines
    next unless line =~ /^<td><a href="/
    if line =~ /<a href="[^"]*">([^<]*)<\/a>/
      cran_packages.push $1
    end
  end
  cran_packages
end

$cran_packages = get_cran_packages()

def nav_link_unless_current(text, path)
  if @item_rep && @item_rep.path && ((@item_rep.path == path) ||
    (path.length > 1 && @item_rep.path[1..-1] =~ /^#{path[1..-1]}/))
    %[<a class="current" href="#{path}">#{text}</a>]
  else
    %[<a href="#{path}">#{text}</a>]
  end
end

def verbose_bioc_version(package)
	if (package.has_key? :bioc_version_str and !package[:bioc_version_str].nil?)
	  "#{package[:bioc_version_str]} (#{package[:bioc_version_num]})"
	else
	  package[:bioc_version_num]
	end
end




def base_filename(path)
  #puts path
  return nil if path.nil?
  return nil if path.is_a? Array
  path.split("/").last
end

# This function returns nil if there is no windows binary at all available for the package.
# If there is a windows package available it will return the path to it. 
# The path may have "windows" or "windows64" in it, but you can't draw
# any conclusions from that, because windows64 is a symlink to windows. It does not
# mean that the package is available only for a particular architecture. Use the Archs flag
# to determine that. 
#
# The fields win.binary.ver and win64.binary.ver may have values or not. If neither of them
# have values, the package is not available for Windows. If either of them have a value, the
# package is available. 
#
# 20101215 - I'm changing the behavior of this function to return the "windows" path (if available) instead of
# the "windows64" path. 
# An ordinary 32-bit windows user might wonder why the download path has 64 in it, 
# and there really isn't any good reason. The more generic "windows" is appropriate.
def windows_binary(package)
  return nil unless package.has_key? :"win.binary.ver" or package.has_key? :"win64.binary.ver"
  return nil if (package[:"win.binary.ver"] == "" or package[:"win.binary.ver"] == [])\
   and (package[:"win64.binary.ver"] == "" or package[:"win64.binary.ver"] == [])
  
  win32 = package[:"win.binary.ver"]
  win64 = package[:"win64.binary.ver"]
  
  return win32 unless win32.nil?  or win32.empty?
  win64
end

def win_format(package)
  wb = windows_binary(package)
  return nil if wb.nil?
  
  both = "(32- &amp; 64-bit)"
  _32only = "(32-bit only)"
  _64only = "(64-bit only)"
  ret = ""
  
  #if windows_binary(package) =~ /64/
  #  ret = both
  #else
  #  ret = _32only
  #end
  
  ret = both
  
  if (package.has_key?(:Archs) && !package[:Archs].empty?)
    archs = package[:Archs]
    if (archs =~ /i386/ && archs =~ /x64/)
      ret = both
    elsif (archs =~ /i386/)
      ret = _32only
    elsif (archs =~ /x64/)
      ret = _64only
    end
  end
  
  
  return ret
end

def cjoin(item, sep)
  if (item.respond_to? :join)
    item.join(sep)
  else
    item
  end
end

def to_array(item) # todo - determine if used
  return item if item.respond_to? :join
  return [item]
end

def version_fragment(package)
  return "&version=#{package[:version_num]}&" if package[:bioc_version_str] == "devel"
  ""
end


def munge_email(email)
  @coder = HTMLEntities.new unless defined? @coder
  ret = ""
  email.gsub(/@/, " at ").split("").each do |char|
    if char.ord > 128 # ?
      ret += @coder.encode(char, :hexadecimal)
    else
      ret += "&#x#{char.bytes.first.to_s(16)};"
    end
  end
  ret
end


def filter_emails(str)
  return str if str.nil?
  emails = str.scan( /(<[^>]*>)/).flatten
  for email in emails
    str = str.gsub(email, munge_email(email))
  end
  str
end

def remove_emails(str)
  str.gsub(/<([^>]*)>/,"").gsub("  "," ").gsub(" ,", ",")
end

def linkify(sym, package)
  items = package[sym]
  # the following key gets set in bioc_views.rb#items()
  key = "#{sym.to_s}_repo".to_sym
  repos = package[key]
  output = []
  
  to_array(items).each_with_index do |item, index|
    next if item.nil?
    item=item.gsub("(", " (").gsub("  (", " (")
    linkable, remainder = item.split(" ", 2)
    remainder = "" if remainder.nil?
    remainder = " " + remainder unless remainder.empty?
    repo = repos[index]

    if (repo == false)
      if ($cran_packages.include?(linkable))
        output.push %Q(<a class="cran_package" href="http://cran.fhcrc.org/web/packages/#{linkable}/index.html">#{linkable}</a>#{remainder})
      else
        output.push item.strip
      end
      next
    end
    if package[:repo] == "bioc/"
      jumpup = "../.."
    else
      jumpup = "../../.."
    end
    output.push %Q(<a href="#{jumpup}/#{repo}/html/#{linkable}.html">#{linkable}</a>#{remainder.strip})
  end
  output.join(", ")
end

def doc_object(package)

  # return an array of hashes
  # [{:file => ..., :title => ..., :script => ...}]


  doc_obj = []
  if package.has_key? :vignettes
    package[:vignettes].each_with_index do |vignette, i|
      next if vignette !~ /\.pdf$/i ## FIX this on BiocViews side
      hsh = {}
      hsh[:type] = "PDF"
      hsh[:file] = vignette.split("/").last
      script = vignette.sub(/\.pdf$/, ".R")
      if (package.has_key? :Rfiles and package[:Rfiles].include? script)
        hsh[:script] = script.split("/").last
      end
      if package.has_key? :vignetteTitles
        hsh[:title] = package[:vignetteTitles][i]
      else
        hsh[:title] = hsh[:file]
      end
      doc_obj.push hsh
    end
  end

  if package.has_key? :htmlDocs
    package[:htmlDocs].each_with_index do |htmlDoc, i|
      hsh = {}
      hsh[:type] = "HTML"
      hsh[:file] = htmlDoc.split("/").last
      script = htmlDoc.sub(/\.html$/i, ".R")
      if (package.has_key? :Rfiles and package[:Rfiles].include? script)
        hsh[:script] = script.split("/").last
      end
      if package.has_key? :htmlTitles
        hsh[:title] = package[:htmlTitles][i].gsub(/^"/, "").gsub(/"$/, "").gsub(/""/, '"')
      else
        hsh[:title] = htmlDoc.split("/").last
      end
      doc_obj.push hsh
    end
  end


  doc_obj.sort! do |a,b|
    a[:title] = "" if (a[:title].nil?)
    b[:title] = "" if (b[:title].nil?)
    if a[:type] != b[:type]
      b[:type] <=> a[:type]
    else
      a[:title].downcase <=> b[:title].downcase
    end
  end

  if doc_obj.empty?
    return [:file => "", :title => "", :script => "", :type => ""]
  end

  doc_obj
end

def bioc_views_links(package)
  links = []
  
  
  if package[:repo] == "bioc/"
    jumpup = "../.."
  else
    jumpup = "../../.."
  end
  
  
  bioc_views = to_array(package[:biocViews])
  bioc_views.each do |bioc_view|
    links.push %Q(<a href="#{jumpup}/BiocViews.html#___#{bioc_view}#{version_fragment(package)}">#{bioc_view}</a>)
  end
  
  links.join(", ")
end


def find_item(items, identifier)
  items.find { |i| i.identifier == identifier }
end

def timeago(time, options = {})
  start_date = options.delete(:start_date) || Time.new
  date_format = options.delete(:date_format) || :default
  delta_minutes = (start_date.to_i - time.to_i).floor / 60
  if delta_minutes.abs <= (8724*60)       
    distance = distance_of_time_in_words(delta_minutes)       
    return "#{distance} ago"
  else
    return "on #{DateTime.now.to_formatted_s(date_format)}"
  end
end

def distance_of_time_in_words(minutes)
  case
  when minutes < 90
    "#{minutes} #{pluralize(minutes, "minute")}"
  when minutes < 1080
    "#{(minutes / 60).round} hours"
  else
    days = (minutes / 1440).round
    "#{days} #{pluralize(days, "day")}"
  end
end

def pluralize(count, what)
  case
  when count == 0 || count > 1
    what + "s"
  when count == 1
    what
  end
end

def add_missing_info
  items.each do |item| 
    if item[:file] 
      # nanoc3 >= 3.1 will have this feature, add for older versions 
      item[:extension] ||= item[:file].path.match(/\.(.*)$/)[0]
    end 
  end
end 

def annual_reports
  # FIXME: need a more robust way to obtain assets path
  Dir.chdir("assets/about/annual-reports") do
    Dir.glob("*").map do |f|
      {
        :href => f,
        :name => /[[:alpha:]]+([[:digit:]]+).pdf/.match(f)[1]
      }
    end
  end
end

def doctype
  %[<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">]
end

def course_materials
  top = @items.find { |i| i.identifier == "/help/course-materials/" }
  top.children.sort { |a, b| b[:title] <=> a[:title] }
end

def upcoming_events(events)
  sorted = events.children.sort do |a, b|
      a[:start] <=> b[:start]
  end
  step2 = sorted.select do |e|
    e[:end] >= Time.now.to_date
  end

  ## Make upcoming BioC sticky at top, if sticky == true
  sticky = false
  if sticky
    bioc = step2.find{|i| i[:title] =~ /^BioC2/} 
    unless bioc.nil?
      step2 = step2.reject{|i| i == bioc}
      step2 = step2.unshift(bioc)
    end
  end
  step2
end

def previous_events(events)
  events.children.sort{|a, b| b[:start] <=> a[:start]}.select do |e|
    e[:end] < Time.now.to_date
  end
end

def event_date(e)
  if (e[:start].month == e[:end].month)
    if (e[:start].day == e[:end].day)
      d1 = e[:start].strftime("%d %B %Y")
      d2 = ""
    else
      d1 = e[:start].strftime("%d")
      d2 = e[:end].strftime(" - %d %B %Y")
    end
  else
    d1 = e[:start].strftime("%d %B")
    d2 = e[:end].strftime(" - %d %B %Y")
  end
  d1 + d2
end

def has_subnav?(anItem)
  (anItem[:subnav] || anItem.parent[:subnav]) rescue false
end

def subnav_items(anItem)
  if anItem[:subnav]
    anItem.children
  elsif anItem.parent[:subnav]
    anItem.parent.children
  else
    []
  end
  rescue
    []
end

def get_stats_url(package)
  if (package[:repo] == "data/annotation/")
    repo = "stats/data-annotation/"
  elsif (package[:repo] == "data/experiment/")
    repo = "stats/data-experiment/"
  else
    repo = "stats/#{package[:repo]}"
  end
  "http://bioconductor.org/packages/#{repo}#{package[:Package]}.html"
end

def get_updated_breadcrumbs(old_breadcrumbs, item)
  return old_breadcrumbs unless (old_breadcrumbs.last.identifier =~ /package-pages/)
  index_page = false
  index_page = true if item.identifier =~ /\/package-pages\/all-/
  last_crumb = old_breadcrumbs.last
  home_crumb = old_breadcrumbs.first
  path = item.path
  segs = path.split("/")
  ver = segs[2]
  repo = ["Software", "bioc"] if path =~ /\/bioc\//
  repo = ["Annotation", "data/annotation"] if path =~ /\/data\/annotation\//
  repo = ["Experiment", "data/experiment"] if path =~ /\/data\/experiment\//
  crumbs = []
  crumbs.push home_crumb
  ver_crumb = Nanoc3::Item.new("", {:title => "Bioconductor #{ver}"}, "/packages/#{ver}/BiocViews.html")
  crumbs.push ver_crumb
  repo_crumb = Nanoc3::Item.new("", {:title => "#{repo.first} Packages"}, "/packages/#{ver}/#{repo.last}/")
  crumbs.push repo_crumb unless index_page
  crumbs.push last_crumb
  crumbs
end

def recent_packages()
  begin
    xml = HTTParty.get("http://bioconductor.org/rss/new_packages.rss").body
    doc = Document.new xml
    items = []
    doc.elements.each("rss/channel/item") {|i| items.push i}
    ret = []
    for item in items#.reverse
      h = {}
      title = nil
      item.elements.each("title") {|i|title = i.text}
      description  = nil
      item.elements.each("description") {|i|description = i.text}
      title_segs = title.split(" ")
      title_segs.shift # get rid of url
      pkg = title_segs.shift # get pkg name
      pkgtitle = title_segs.join " "
      desc = description.split("<br/").first
      h[:package] = pkg
      h[:title] = pkgtitle
      h[:description] = desc
      ret.push h
    end
    return ret
  rescue Exception => ex
    return []
  end
end

def get_svn_commits()
  begin
    xml = HTTParty.get("http://bioconductor.org/rss/svnlog.rss").body
    doc = Document.new xml
    items = []
    doc.elements.each("rss/channel/item") {|i| items.push i}
    ret = []
    for item in items 
      next if item.nil?
      next if item.elements.nil?
      next unless item.elements.respond_to? :each
      h = {}
      revision = title = date = author = description = nil
      item.elements.each("title") {|i| title = i.text}
      item.elements.each("pubDate") {|i| date = i.text}
      item.elements.each("author") {|i| author = i.text}
      item.elements.each("description") {|i| description = i.text}
      
      msg = description.gsub(/<div style="white-space:pre">/,"")
      msg = msg.split("</div>").first
      table = "<table>" + description.split("<table>").last
      
      
      rdate = DateTime.strptime(date, "%a, %e %b %Y %H:%M:%S %Z").iso8601
      date = %Q(<abbr class="timeago" title="#{rdate}">#{rdate}</abbr>)

      h[:revision] = title.split(" ")[1]
      h[:date] = date
      h[:author] = author
      h[:msg] = msg
      h[:table] = table
      ret.push h
    end
    return ret
  rescue Exception => ex
    return []
  end
end

def since(package)
  for key in config[:manifest_keys]
    if config[:manifests][key].include? package
      return key
    end
  end
  nil
end


def get_version_from_item_id(item)
  segs = item.identifier.split "/"
  segs.pop
  version = segs.pop
  version
end

def script_tag_for_package_data(item)
  # todo - something sensible if get_json hasn't been run
  segs = item.identifier.split "/"
  segs.pop
  version = segs.pop
  repos = ["bioc", "data/annotation", "data/experiment"]
  if version == config[:devel_version]
    repos = config[:devel_repos]
  end
  s = ""
  for repo in repos
    #<script type="text/javascript" src="/packages/json/<%=get_version_from_item_id(@item)%>/tree.json">
    s += %Q(<script type="text/javascript" src="/packages/json/#{version}/#{repo}/packages.js"></script>\n)
  end
  s
end

def get_tree(item)
  # todo - something sensible if get_json hasn't been run
  segs = item.identifier.split "/"
  segs.pop
  version = segs.pop
  f = File.open "assets/packages/json/#{version}/tree.js"
  json = f.readlines.join("")
  f.close
  return json.strip
end

def r_ver_for_bioc_ver(bioc_ver)
  hsh = config[:r_ver_for_bioc_ver]
  ret = hsh[bioc_ver.to_sym]
  return ret
end

def get_package_maintainers()
  exclude_these_packages = ["limma", "edgeR"]
  ret = []
  json_file \
    = "assets/packages/json/#{config[:release_version]}/bioc/packages.json"
  return ret unless test(?f, json_file)
  f = File.open(json_file)
  json = f.readlines.join
  f.close
  hsh = JSON.parse(json)
  hsh.each_pair do |k, v|
    row = []
    row.push k
    maintainer \
      = v["Maintainer"].gsub(/^[^<]*/,"").gsub(/<|>/, "").gsub("@", " at ")
    ## there should not be more than one maintainer, but if there is,
    ## just pick the first one
    maintainer = maintainer.split(",").first.strip
    
    maintainer = "#{k} Maintainer <#{maintainer}>"
    if exclude_these_packages.include? k
      row.push ""
    else
      row.push maintainer
    end
    ret.push row
  end
  ret.sort! do |a, b|
    a.first.downcase <=> b.first.downcase
  end
  ret
end


def mac_os(pkg)
  if pkg.has_key? :"mac.binary.ver"
    return "Mac OS X 10.6 (Snow Leopard)"
  elsif pkg.has_key? :"mac.binary.leopard.ver"
    return "Mac OS X 10.5 (Leopard)"
  end
  return nil
end


def mac_ver_key(pkg)
  if pkg.has_key? :"mac.binary.ver"
    return pkg[:"mac.binary.ver"]
  elsif pkg.has_key? :"mac.binary.leopard.ver"
    return pkg[:"mac.binary.leopard.ver"]
  end
  return nil
end

def workflow_helper(item)
  w = item.attributes.dup
  id = item.identifier.sub(/\/$/, "")
  pkg = id.split("/").last
  segs = item.identifier.split "/"
  3.times {segs.shift}
  multivig = (segs.length > 1)
  if multivig
    ["source_tarball", "mac_pkg", "win_pkg"].each do |pkgtype|
      w[pkgtype.to_sym] = "../#{w[pkgtype.to_sym]}"
    end
  else
    w["r_source".to_sym] = "#{pkg}.R"
  end
  ["mac_pkg", "win_pkg"].each do |binpkg|
    if w[binpkg.to_sym] == "NOT_SUPPORTED"
      w.delete binpkg.to_sym
    end
  end
  w

end

def release_branch
  config[:release_version].sub(".", "_")
end

# call me with @package[:URL]
def make_package_url_links(url)
    out = ""
    segs = nil

    url = url.gsub /\s+/, "" if url =~ /,\s/

    if url =~ /\s/
        segs = url.split(/\s+/)
    elsif url =~ /,/
        segs = url.split(",")
    else
        segs = [url]
    end
    for seg in segs
        out += %Q(<a href="#{seg}">#{seg}</a> )
    end
    out
end

#FIXME  should gracefully fail (and allow flow to continue) if no internet access
def get_build_summary(version, repo)
    url = "http://bioconductor.org/checkResults/#{version}/#{repo}-LATEST/"
    css_url = "#{url}/report.css"
    html = open(url)
    doc = Nokogiri::HTML(html.read)
    doc.encoding = "ascii"
    dateline = doc.css %Q(p[style="text-align: center;"])
    dateline = dateline.children[1].text
    dateline.sub!(/^This page was generated on /, "")
    dateline = dateline.split("(").first.strip

    rows = doc.css("table.mainrep tr.summary")

    htmlfrag=<<-"EOT"
        <head>
        <base href="#{url}" target="_blank">
        </head
        <body>
        <p><i>Build report generated at #{dateline}</i></p>
        <LINK rel="stylesheet" href="#{css_url}" type="text/css">
        <table>
            #{rows.to_html}
        </table>
        </body>
    EOT
    FileUtils.mkdir_p "output/dashboard"
    File.open("output/dashboard/build_#{version}_#{repo}.html", "w") do |f|
        f.puts htmlfrag
    end
    ret=<<-EOT
    <iframe src="/dashboard/build_#{version}_#{repo}.html" width="80%"></iframe>
    EOT
    ret
end

# FIXME gracefully fail w/o internet access
def get_new_packages_in_tracker()
    return "" unless File.exists?("tracker.yaml")
    url = "http://tracker.fhcrc.org/roundup/bioc_submit/"
    cfg = YAML::load(File.open("tracker.yaml"))
    @agent = Mechanize.new
    page = @agent.post(url, {
        "__login_name" => cfg['username'],
        "__login_password" => cfg['password'], 
        "__came_from" => url,
        "@action" => "login"
    })
    rows = page.search("table.list tr")
    nr = []
    #nr.push rows.first
    header = <<-"EOT"
    <tr>
  <th>ID</th>
   
   <th>Activity</th>
   
   
   <th>Title</th>
   <th>Status</th>
  </tr>
<tr>
EOT
    nr.push header
    for i in 2..12
      t = rows[i].to_s
      t.gsub!(/<td>[^<]+<\/td>\s+<td>[^<]+<\/td>\s+<\/tr>/, "</tr>")

        nr.push t
    end
    s = "<table>\n"
    nr.each do |i|
        #puts i.to_s
        html = i.to_s.sub(%Q(a href="), # i.to_html.sub
                %Q(a href="http://tracker.fhcrc.org/roundup/bioc_submit/))
        s += html
    end
    s += "</table></body>"
    s
end

def get_mailing_list_link(devel=false)
    if devel
        list = "bioc-devel"
    else
        list = "bioconductor"
    end
    d = DateTime.now
    month = d.strftime("%B")
    year = d.strftime("%Y")
    "https://stat.ethz.ch/pipermail/#{list}/#{year}-#{month}/thread.html"
end

def get_search_terms()
    return "" unless File.exists? "analytics_py/client_secrets.json"
    res = nil
    FileUtils.mkdir_p "output/dashboard"

    Dir.chdir("analytics_py") do
        res = `python search_terms.py > ../output/dashboard/search_terms.tsv`
    end
    html=<<-"EOT"
<meta charset="utf-8">
<style>

.bar {
  fill: steelblue;
}

.bar:hover {
  fill: brown;
}

.axis {
  font: 10px sans-serif;
}

.axis path,
.axis line {
  fill: none;
  stroke: #000;
  shape-rendering: crispEdges;
}

.x.axis path {
  display: none;
}

</style>

<script src="http://d3js.org/d3.v3.min.js"></script>
<div id="search_terms_chart"></div>
<script src="/js/search_terms.js"></script>
    EOT
    html
end

def get_hits()
    return "" unless File.exists? "analytics_py/client_secrets.json"
    FileUtils.mkdir_p "output/dashboard"
    Dir.chdir("analytics_py") do
        res = `python hits.py > ../output/dashboard/hits.tsv`
    end
    html=<<-"EOT"
     <script type="text/javascript" src="https://www.google.com/jsapi"></script>
     <script type="text/javascript" src="/js/hits.js"></script>
    <div id="chart_div" style="width: 900px; height: 500px;"></div>

    EOT
    html
end

def get_current_time
  Time.now.utc.iso8601
end

def recent_spb_builds
    begin
      HTTParty.get("http://merlot2.fhcrc.org:8000/recent_builds").body
    rescue Exception => ex
      "Can't connect to merlot2, not building dashboard"
    end
end

def get_last_svn_commit_time()

    xml = `curl -s http://bioconductor.org/rss/svnlog.rss`
    doc = Document.new xml
    items = []
    doc.elements.each("rss/channel/item") {|i| items.push i}
    date = nil
    item = items.first
    item.elements.each("pubDate") {|i| date = i.text}
    rdate = DateTime.strptime(date, "%a, %e %b %Y %H:%M:%S %Z").iso8601
    %Q(<abbr class="timeago" title="#{rdate}">#{rdate}</abbr>)
end

def get_mac_packs(package, item)

    res = []
    os =  ["Mac OS X 10.6 (Snow Leopard)"] 
    osvers = ["mac.binary.ver"] 

    version = item.identifier.split("/")[4].to_f
    if (version > 2.13)
        os <<  "Mac OS X 10.9 (Mavericks)"
        osvers << "mac.binary.mavericks.ver"
    end



    os.each_with_index do |this_os, i|
        h = {}
        h[:os] = this_os
        if package.has_key? osvers[i].to_sym
            h[:url] = "../#{package[osvers[i].to_sym]}"
            h[:basename] = package[osvers[i].to_sym].split("/").last
        end
        res.push h
    end
    res
end

def is_devel(item)
    return true if item.identifier =~ /\/devel\/|\/#{config[:devel_version]}\//
    return false 
end

def is_new_package(package)
    since = since(package[:Package])
    return(since == config[:devel_version])
end

def get_release_url(item_rep)
    item_rep.raw_path.sub(/^output/, "").sub(/\/devel\/|\/#{config[:devel_version]}\//, "/release/")
end


def get_devel_fragment(package, item, item_rep)
    return "" unless is_devel(item)

    str=<<-"EOT"
<p>This is the <b>development</b> version of #{@package[:Package]}
EOT
    unless is_new_package(package)
        str2=<<-"EOT"
; for the stable release version, see 
<a href="#{get_release_url(item_rep)}">#{@package[:Package]}</a>
EOT
        str = str.strip() + str2
    end
    str = str.strip() + ".</p>"
    str
end

def package_has_source_url(item)
    segs = item.identifier.split('/')
    return true if segs[5] == "bioc"
    return true if segs[5] == "data" and segs[6] == "experiment"
    return false
end

def get_source_url(package, item, item_rep)
    url = "https://hedgehog.fhcrc.org/"
    segs = item.identifier.split("/")
    repos = segs[5]
    repos = segs[5] + "/" + segs[6] if segs[5] == "data" \
      and segs[6] == "experiment"
    if repos == "bioc"
        url += "bioconductor/"
    else
        url += "bioc-data/"
    end
    if is_devel(item)
        url += "trunk/"
    else
        pkg_version = segs[4].sub(".", "_")
        url += "branches/RELEASE_#{pkg_version}/"
    end
    if repos == "bioc"
        url += "madman/Rpacks/"
    else
        url += "experiment/pkgs/"
    end
    url += package[:Package] + "/"
    url
end

def get_video_title(video)
   response = HTTParty.get(video, :verify => false)
   doc = Nokogiri::HTML(response.body)
   doc.css("title").text.sub(/ - YouTube$| on Vimeo$/, "")
end

def render_courses()
    lines = File.readlines("etc/course_descriptions.tsv")

    lines = lines.sort do |a,b|
        d1 = a.split("\t").first
        d2 = a.split("\t").first
        b <=> a
    end

    headers = lines.shift.strip.split("\t")
    out=<<-"EOT" # what class?
<table id="course_descriptions">
    <thead>
        <tr>
            <th>Keyword</th>
            <th>Course</th>
            <th>Title</th>
            <th>Materials</th>
            <th>Date</th>
            <th>Bioc/R Version</th>
        </tr>
    </thead>
    <tbody>
    EOT
    lines.each_with_index do |line, idx|
        # break if idx == 1
        lh = {}
        line.strip.split("\t").each_with_index do |seg, i|
            lh[headers[i].strip] = seg
        end
        year = lh["Date"].split(" - ").first.strip.split("-").first
        course_url = "#{year}/#{lh["Course"]}/"
        out += "        <tr>\n"
        out += "<td>" + lh["Keyword"] + "</td>\n"
        out += "<td><a href='#{course_url}'>" + lh["Course"] + "</a></td>\n"
        out += "<td>" + Kramdown::Document.new(lh["Title"].strip + ", "  + lh["Instructor"].strip).to_html +  "</td>\n"
        out += "<td>" + Kramdown::Document.new(lh["Material"]).to_html + "</td>\n"
        out += "<td>" + lh["Date"].split(" - ").first.strip.gsub("-", "&#8209;") + "</td>\n"
        biocver = lh['Bioc version']
        biocver = "3.0" if biocver.strip == "3"
        out += "<td>" + biocver + '/' +  lh["R version"]  + "</td>\n"



        # ...
        out += "        </tr>\n"
    end

    out+=<<-"EOT"
            </tbody>
        </table>
    EOT
    out
end

def ami_url(ami)
    "<a href='https://console.aws.amazon.com/ec2/home?region=us-east-1#launchAmi=#{ami}'>#{ami}</a>"
end