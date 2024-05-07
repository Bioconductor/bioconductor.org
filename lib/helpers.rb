# encoding: utf-8
require 'nanoc'
include Nanoc::Helpers::Text
include Nanoc::Helpers::Rendering
include Nanoc::Helpers::Breadcrumbs
include Nanoc::Helpers::XMLSitemap
include Nanoc::Helpers::HTMLEscape


# coding: utf-8
require 'htmlentities'
require 'time'
require 'date'
require 'rubygems'
require 'httparty'
require 'net/http'
require 'net/https'
require 'yaml'
require 'pp'
require 'rexml/document'
require 'json'
require 'open-uri'
require 'nokogiri'
require 'fileutils'
require 'mechanize'
require 'kramdown'
require 'open3'
require 'open-uri'
require 'socket'
require 'cgi'
require 'octokit'

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
  def initialize()
    @rownum = 0
  end
  def rowclass()
    @rownum += 1
    (@rownum % 2 == 1) ? "row_odd" : "row_even"
  end
end

class TableRower
  def initialize()
    @cellnum = 0
  end
  def betweencells(cells_per_row=4)
    @cellnum += 1
    (@cellnum % cells_per_row == 0) ? "</tr><tr>\n" : ""
  end
end


def get_cran_packages()
  puts "Grabbing list of CRAN packages..."
  ## Is there a non-web-scraping way to get a list of CRAN packages?
  cran_packages = []
  data = begin
           HTTParty.get("http://cran.rstudio.com/web/packages/available_packages_by_name.html",
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

$pkgdata = nil

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

# This function returns nil if there is no windows binary available for the
# package.
# If there is a windows package available it will return the path to it.
# The field win.binary.ver may have a value or not. If it doesn't, then the
# package is not available for Windows.
def windows_binary(package)
  return nil unless package.has_key? :"win.binary.ver"
  wb = package[:"win.binary.ver"]
  if wb.nil? or wb.empty? or wb == "" or wb == []
    return nil
  end
  return wb
end

def win_format(package)
  wb = windows_binary(package)
  return "" if wb.nil?
  if package.has_key? :Archs
    archs = package[:Archs]
    if archs.nil? or archs.empty? or archs == "" or archs == []
      return ""
    end
    i386 = (archs =~ /i386/)
    x64 = (archs =~ /x64/)
    if i386 and x64
      return "(32- &amp; 64-bit)"
    end
    if i386
      return "(32-bit only)"
    end
    if x64
      return "(64-bit only)"
    end
  end
  return ""
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
    if email.include? "orcid"
      email2="<a title='orcid' href='"+email[1...-1]+"'><img src='/images/orcid.png'/></a>"
      email1="("+email+")"
      str = str.gsub(email1, email2)
    else
      str = str.gsub(email, munge_email(email))
    end
  end
  str
end

def remove_emails(str)
  return str if str.nil?
  str.gsub(/<([^>]*)>/,"").gsub("  "," ").gsub(" ,", ",")
end

def linkify(sym, package)
  if defined?($cran_packages).nil?
    $cran_packages = get_cran_packages
  end

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
        output.push %Q(<a class="cran_package" href="http://cran.rstudio.com/web/packages/#{linkable}/index.html">#{linkable}</a>#{remainder})
      else
        output.push item.strip
      end
      next
    end
    if (package[:repo] == "bioc/" || package[:repo] == "workflows/")
      jumpup = "../.."
    else
      jumpup = "../../.."
    end
    output.push %Q(<a href="#{jumpup}/#{repo}/html/#{linkable}.html">#{linkable}</a>#{remainder.strip})
  end
  output.join(", ")
end

def linkify_license(package)
  link = '<a href = "../licenses/%s/LICENSE">LICENSE</a>' % [ package[:Package] ]
  package[:License].sub("LICENSE", link)
end

def doc_object(package)

  # return an array of hashes
  # [{:file => ..., :title => ..., :script => ...}]


  doc_obj = []
  if package.has_key? :vignettes
    package[:vignettes].each_with_index do |vignette, i|
      hsh = {}
      if vignette =~ /\.pdf$/i ## FIX this on BiocViews side
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
      elsif vignette =~ /\.html$/i
          hsh[:type] = "HTML"
          hsh[:file] = vignette.split("/").last
          script = vignette.sub(/\.html$/, ".R")
          if (package.has_key? :Rfiles and package[:Rfiles].include? script)
            hsh[:script] = script.split("/").last
          end
          if package.has_key? :vignetteTitles
            hsh[:title] = package[:vignetteTitles][i].gsub(/^"/, "").gsub(/"$/, "").gsub(/""/, '"')
          else
            hsh[:title] = vignette.split("/").last
          end
          doc_obj.push hsh
      else
        next
      end
    end
  end

  doc_obj.sort! do |a,b|
    a[:title] = "" if (a[:title].nil?)
    b[:title] = "" if (b[:title].nil?)
    a[:title].downcase <=> b[:title].downcase
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

def get_biostar_post_summaries(items)
  # at some point, might need to order these....
  items.find_all{|i| i.identifier =~ /\/biostar_list\//}
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
    Dir.glob("AnnRep*.pdf").map do |f|
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

def top_events(events)
  sorted = events.children.sort do |a, b|
      a[:start] <=> b[:start]
  end
  toplist = sorted[-5..-1]
  toplist.reverse
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
  "http://bioconductor.org/packages/#{repo}#{package[:Package]}/"
end

def get_updated_breadcrumbs(old_breadcrumbs, item)
  return old_breadcrumbs unless (old_breadcrumbs.last.identifier.to_s =~ /package-pages/)
  index_page = false
  index_page = true if item.identifier.to_s =~ /\/package-pages\/all-/

  return old_breadcrumbs if item.identifier.to_s =~ /\/package-pages\/all--/

  last_crumb = old_breadcrumbs.last
  home_crumb = old_breadcrumbs.first
  path = item.path
  segs = path.split("/")
  ver = segs[2]
  repo = ["Software", "bioc"] if path =~ /\/bioc\//
  repo = ["Annotation", "data/annotation"] if path =~ /\/data\/annotation\//
  repo = ["Experiment", "data/experiment"] if path =~ /\/data\/experiment\//
  repo = ["Workflow", "workflows"] if path  =~ /\/workflows\//
  crumbs = []
  crumbs.push home_crumb
  ver_crumb = Nanoc::Int::Item.new("", {:title => "Bioconductor #{ver}"}, "/packages/#{ver}/BiocViews.html")
  crumbs.push ver_crumb
  repo_crumb = Nanoc::Int::Item.new("", {:title => "#{repo.first} Packages"}, "/packages/#{ver}/#{repo.last}")
  crumbs.push repo_crumb unless index_page
  crumbs.push last_crumb
  crumbs
end


def get_git_commits()
   begin
     REXML::Document.entity_expansion_text_limit =
       REXML::Document.entity_expansion_text_limit * 4
     doc = Document.new File.new("assets/developers/rss-feeds/gitlog.xml")
     items = []
     doc.elements.each("rss/channel/item") {|i| items.push i}
     ret = []
     for item in items
       next if item.nil?
       next if item.elements.nil?
       next unless item.elements.respond_to? :each
       h = {}
       title = date = author = description = guid = nil
       item.elements.each("title") {|i| title = i.text}
       item.elements.each("pubDate") {|i| date = i.text}
       item.elements.each("author") {|i| author = i.text}
       item.elements.each("description") {|i| description = i.text}
       item.elements.each("guid") {|i| guid = i.text}

       msg = description.gsub(/<div style="white-space:pre">/,"")
       msg = msg.split("</div>").first

       h[:package] = title
       h[:date] = date
       h[:author] = author
       h[:msg] = msg
       h[:commit] = guid

       ret.push h
     end
     return ret
   end
end


def get_svn_commits()
  begin
    # FIXME - maybe switch to nokogiri to avoid the necessity
    # for this. See https://stackoverflow.com/questions/15593133/rexml-runtimeerror-entity-expansion-has-grown-too-large
    REXML::Document.entity_expansion_text_limit =
      REXML::Document.entity_expansion_text_limit * 4
    xml = HTTParty.get("http://bioconductor.org/rss/gitlog.rss").body
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

def since(package, conf=nil)
  if conf.nil?
    config = @config
  else
    config = conf
  end

  for key in config[:manifest_keys]
    if config[:manifests][key].include? package
      return key
    end
  end
  nil
end

def get_year_shield(package, make_shield=false, conf=nil)
  if conf.nil?
    config = @config
  else
    config = conf
  end
  yib = years_in_bioc(package, config)

  destdir = File.join("assets", "shields", "years-in-bioc")
  FileUtils.mkdir_p destdir
  shield = File.join(destdir, "#{package}.svg")
  now = DateTime.now
  onedayago = now.prev_day
  if ((!File.exists?(shield)) or  DateTime.parse(File.mtime(shield).to_s) < onedayago)
    if is_new_package2(package, config)
      if make_shield
        puts "Downloading years-in-bioc shield for #{package}..."
        FileUtils.cp(File.join('assets', 'images', 'shields',
          'in_bioc', "devel-only.svg"), shield)
      end
    elsif yib.nil?
      return nil
    else
      if make_shield
        puts "Downloading years-in-bioc shield for #{package}..."
        puts "#{yib}"

        yib = yib.gsub("<", "&#60;")
        template = File.read(File.join('assets', 'images',
        'shields', 'in_bioc', 'inbioc-template.svg'))
        newbadge = template.gsub(/9999 years/,  yib)
        newbadge = newbadge.gsub(/textLength=\"(610)\"/, '')
        File.open(shield, "w") { |file| file.write(newbadge) }

      end
    end
  end
  true
end

def years_in_bioc(package, conf=nil)
  if conf.nil?
    config = @config
  else
    config = conf
  end
  since_ver = since(package, config)
  return nil if since_ver.nil?
  key = since_ver.to_sym
  unless config[:release_dates].has_key? key
    return nil
  end
  release_date_str = config[:release_dates][key]
  release_date = Date.strptime(release_date_str, "%m/%d/%Y")
  today = Date.today
  days_since_release = (today - release_date)
  years_since_release = days_since_release / 365.25

  rep = (years_since_release * 2).round / 2.0
  rep = rep.to_i == rep ? rep.to_i : rep
  srep = (rep.is_a? Integer) ? rep.to_i.to_s : rep.to_s
  if srep == "1"
    srep += " year"
  else
    srep += " years"
  end
  srep = "> " + srep if since_ver == "1.6"


  # we probably won't see this:
  if years_since_release <= 0.5
    return "< 6 months"
  end
  srep
end

def get_version_from_item_id(item)
  segs = item.identifier.to_s.split "/"
  segs.pop
  version = segs.pop
  version
end

def script_tag_for_package_data(item)
  # todo - something sensible if get_json hasn't been run
  segs = item.identifier.to_s.split "/"
  segs.pop
  version = segs.pop
  repos = ["bioc", "data/annotation", "data/experiment", "workflows"]
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
  segs = item.identifier.to_s.split "/"
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
    #return "Mac OS X 10.6 (Snow Leopard)"
    return "macOS 10.13 (High Sierra)"
  elsif pkg.has_key? :"mac.binary.mavericks.ver"
    return "Mac OS X 10.9 (Mavericks)"
  elsif pkg.has_key? :"mac.binary.el-capitan.ver"
    return "Mac OS X 10.11 (El Capitan)"
  end
  return nil
end


def mac_ver_key(pkg)
  if pkg.has_key? :"mac.binary.ver"
    return pkg[:"mac.binary.ver"]
  elsif pkg.has_key? :"mac.binary.mavericks.ver"
    return pkg[:"mac.binary.mavericks.ver"]
  elsif pkg.has_key? :"mac.binary.el-capitan.ver"
    return pkg[:"mac.binary.el-capitan.ver"]
  end
  return nil
end

def workflow_helper(item)
  w = item.attributes.dup
  id = item.identifier.to_s.sub(/\/$/, "")
  pkg = id.split("/").last
  segs = item.identifier.to_s.split "/"
  4.times {segs.shift}
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

def get_build_summary(version, repo)
    url = "http://bioconductor.org/checkResults/#{version}/#{repo}-LATEST/"
    if repo == "bioc"
      url_without_protocol = url.sub(/^http:/i, "")
      css_url = "#{url_without_protocol}report.css"
      url = url + "long-report.html"
    else
      url_without_protocol = url.sub(/^http:/i, "")
      css_url = "#{url_without_protocol}report.css" 
    end
    begin
      html = open(url)
    rescue Exception => e
      puts "open(url) failed"
      puts "  url: " + url
      puts "  message: " + e.message
      return ""
    end
    doc = Nokogiri::HTML(html.read)
    doc.encoding = "ascii"
    dateline = doc.css("p.time_stamp")
    return "" if dateline.children[0].nil?
    dateline = dateline.children[0].text
    dateline.sub!(/^This page was generated on /, "")
    #dateline = dateline.split("(").first.strip

    #rows = doc.css("table.mainrep tr.summary")
    rows = doc.css("thead.quickstats")
      
    htmlfrag=<<-"EOT"
        <head>
        <base href="#{url_without_protocol}" target="_blank">
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
    <iframe src="/dashboard/build_#{version}_#{repo}.html" width="80%", height=200></iframe>
    EOT
    ret
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

def get_current_time
  Time.now.utc.iso8601
end

def get_mac_packs(package, item)

    res = []

    os = []
    osvers = []

    version =  Gem::Version.new(package[:bioc_version_num])

    if version < Gem::Version.new('3.3')
      os <<  "Mac OS X 10.6 (Snow Leopard)"
      osvers << "mac.binary.ver"
    end

    if version > Gem::Version.new('2.13') and version < Gem::Version.new('3.5')
        os <<  "Mac OS X 10.9 (Mavericks)"
        osvers << "mac.binary.mavericks.ver"
    end

    if version >= Gem::Version.new('3.5') and version < Gem::Version.new('3.11')
        os <<  "Mac OS X 10.11 (El Capitan)"
        osvers << "mac.binary.el-capitan.ver"
    end

    if version >= Gem::Version.new('3.11') and version < Gem::Version.new('3.15')
        os <<  "macOS 10.13 (High Sierra)"
        osvers << "mac.binary.ver"
    end

    if version >= Gem::Version.new('3.15') and version < Gem::Version.new('3.16')
        os <<  "macOS Binary (x86_64)"
        osvers << "mac.binary.ver"
    end

    if version == Gem::Version.new('3.16')
        os <<  "macOS Binary (x86_64)" << "macOS Binary (arm64)"
        osvers << "mac.binary.ver" << "mac.binary.big-sur-arm64.ver"
    end
    
    if version >= Gem::Version.new('3.17')
        os <<  "macOS Binary (x86_64)" << "macOS Binary (arm64)"
        osvers << "mac.binary.big-sur-x86_64.ver" << "mac.binary.big-sur-arm64.ver"
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

def is_old?(item)
  return false if item.identifier.to_s =~ /\/devel\/|\/release\//
  return false if item.identifier.to_s =~ /\/#{config[:devel_version]}\/|\/#{config[:release_version]}\//
  return true
end


def is_devel?(item)
    return true if item.identifier.to_s =~ /\/devel\/|\/#{config[:devel_version]}\//
    return false
end

def is_removed?(item)
    return true if File.foreach("content/about/removed-packages.md").grep(/\/#{item}\.html/).any?
    return false
end

def isReleaseOrDevel(package)
    config = @config
    if package[:bioc_version_num] == config[:devel_version]
        return "devel"
    else
        return "release"
    end
end

# FIXME eventually replace is_new_package() implementation with this
def is_new_package2(package, conf) # just a string (pkg name)
  keys = conf[:manifest_keys].dup
  keys.pop
  for rel in keys.reverse
    if conf[:manifests][rel].include? package
      return false
    end
  end
  return true
end

def is_new_package(package)
    if $pkgdata.nil? or $pkgdata[package[:repo]].nil?
        if $pkgdata.nil?
          $pkgdata = {}
        end
        if $pkgdata[package[:repo]].nil?
          $pkgdata[package[:repo]] = {}
        end
        dir = File.join("assets", "packages", "json")
        d = Dir.new(dir)
        for entry in d.entries
            next if entry =~ /^\./
            file = File.join(dir, entry,
                package[:repo].sub(/\/$/, "").gsub("/", File::SEPARATOR), "packages.json")
            if (File.exists?(file))
                f = File.open(file)
                $pkgdata[package[:repo]][entry] = JSON.load(f)
                f.close
            end
        end
    end

    obj = $pkgdata[package[:repo]]
    k = obj.keys#.sort {|a,b| b.to_f <=> a.to_f}
    for key in k
        next if key.to_f >= package[:bioc_version_num].to_f
        return false if obj[key].has_key? package[:Package]
    end
    return true
    # remove:
    # since = since(package[:Package])
    # return(since == config[:devel_version])

end

def get_release_url_orig(item_rep)
    item_rep.raw_path.sub(/^output/, "").sub(/\/devel\/|\/#{config[:devel_version]}\//, "/release/")
end

def get_release_url(item)
    item.path.sub(/\/devel\/|\/#{config[:devel_version]}\//, "/release/")
end


def get_fragment(package, item, item_rep)
  return \
    get_removed_fragment(package, item, item_rep) if is_removed? @package[:Package]
  return \
    get_devel_fragment(package, item, item_rep) if is_devel? item
  return \
    get_old_fragment(package, item, item_rep) if is_old? item
  return ""
end

def get_removed_link(item)
  line = File.readlines("content/about/removed-packages.md").select{|l| l.match /#{item}/}.last
  line[/<li>(.*?)<\/li>/,1]
end

def get_removed_fragment(package, item, item_rep)
  segs = item.identifier.to_s.split('/')
  version_seg = segs.find_index "package-pages"
  return "" if version_seg.nil? or version_seg >= (segs.length() -1)
  bioc_version = segs[version_seg+1]
  link = get_removed_link(@package[:Package])
  str=<<-"EOT"
<p>This package is for version #{bioc_version} of Bioconductor.
This package has been removed from Bioconductor.
For the last stable, up-to-date release version, see #{link}.</p>
EOT
end

def get_old_fragment(package, item, item_rep)
  segs = item.identifier.to_s.split('/')
  version_seg = segs.find_index "package-pages"
  return "" if version_seg.nil? or version_seg >= (segs.length() -1)
  bioc_version = segs[version_seg+1]
  str=<<-"EOT"
<p>This package is for version #{bioc_version} of Bioconductor;
for the stable, up-to-date release version, see
<a href="/packages/#{@package[:Package]}/">#{@package[:Package]}</a>.</p>
EOT
  str
end

def get_devel_fragment(package, item, item_rep)
    str=<<-"EOT"
<p>This is the <b>development</b> version of #{@package[:Package]}
EOT
    if is_new_package(package)
        str2=<<-"EOT"
; to use it, please install the <a href="/developers/how-to/useDevel/">devel version</a> of Bioconductor
EOT
    else
        str2=<<-"EOT"
; for the stable release version, see
<a href="#{get_release_url(item)}">#{@package[:Package]}</a>
EOT
    end
    str = str.strip() + str2
    str = str.strip() + ".</p>"
    str
end

def package_has_source_url(item, software_only=false)
    segs = item.identifier.to_s.split('/')
    return true if segs[5] == "bioc"
    unless software_only
      return true if segs[5] == "data" and segs[6] == "experiment"
      return true if segs[5] == "workflows"
    end
    return false
end


def get_source_url(package, item, item_rep, access_type)
    segs = item.identifier.to_s.split('/')
    if access_type == "ssh"
        "git@git.bioconductor.org:packages/" + package[:Package]
    else
        "https://git.bioconductor.org/packages/" + package[:Package]
    end
end

# try to determine if a package is in the code brower based on biocViews
def package_in_code_browser(package, item)
    # Skip anything that's not a software package
    segs = item.identifier.to_s.split('/')
    return false if segs[5] != "bioc"
    # Deprecated package are missing in the code browser
    return false if package[:PackageStatus] == "Deprecated"
    # Skip if we can't determine the git branch
    return false if package[:git_branch].nil?
    # if we get here, it's probably in the code browser
    return true
end

# create URL for package in code browser 
def get_code_browser_url(package, include_branch=true)
    url = "https://code.bioconductor.org/browse/" + package[:Package] + "/"
    if include_branch
      url = url + package[:git_branch] + "/"
    end
    return url
end

def get_video_title(video)
   response = HTTParty.get(video, :verify => false)
   doc = Nokogiri::HTML(response.body)
   doc.css("title").text.sub(/ - YouTube$| on Vimeo$/, "")
end

def render_courses()
    lines = File.readlines("etc/course_descriptions.tsv")


    headers = lines.shift.strip.split("\t")

    lines = lines.sort do |a,b|
        d1 = a.split("\t").first.split(" ").last.strip
        t1 = a.split("\t")[3]
        d2 = b.split("\t").first.split(" ").last.strip
        t2 = b.split("\t")[3]
        [d2, t1] <=> [d1, t2]   # most recent date, then alpha by title
    end


    out=<<-"EOT" # what class?
<table id="course_descriptions">
    <thead>
        <tr>
            <th>Keyword</th>
            <th>Title</th>
            <th>Course</th>
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
        out += "        <tr valign=\"top\">\n"
        out += "<td>" + lh["Keyword"] + "</td>\n"
        out += "<td>" + Kramdown::Document.new(lh["Title"].strip + ", "  + lh["Instructor"].strip).to_html.gsub(/<\/*p>/, "") +  "</td>\n"
        if lh["Course"].start_with? "["
         out += "<td>" + Kramdown::Document.new(lh["Course"]).to_html.gsub(/<\/*p>/, "") + "</td>\n"
        else
          out += "<td><a href='#{course_url}'>" + lh["Course"] + "</a></td>\n"
        end
        out += "<td>" + Kramdown::Document.new(lh["Material"]).to_html.gsub(/<\/*p>/, "") + "</td>\n"
        out += "<td>" + lh["Date"].split(" ")[0].gsub("-", "&#8209;") + "</td>\n"
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

def pad(input)
    input.to_s.rjust(2, '0')
end

def iterate_month_mode(now, code)
    res = []
    month_mode  = true
    (1..12).each do |offset|
        res << code.call(offset, month_mode)
    end
    res
end

def get_pubmed_cache_date
    cachefile = "tmp/pubmed_cache_file.yaml"
    return "" unless File.exists? cachefile
    File.mtime(cachefile).iso8601
end

def render_mirror_urls(mirror)
    out = ""
    if mirror.key?(:https_mirror_url) then
        out += "<#{mirror[:https_mirror_url]}>; "
    end
    out + "<#{mirror[:mirror_url]}>"
end

def render_mirror_contacts(mirror_orig)
    mirror = mirror_orig.dup
    unless mirror[:contact].is_a? Array and mirror[:contact_email].is_a? Array
        mirror[:contact] = [mirror[:contact]]
        mirror[:contact_email] = [mirror[:contact_email]]
    end
    out = ""
    len = mirror[:contact].length
    mirror[:contact].each_with_index do |m, i|
        out += "#{m} &lt;#{mirror[:contact_email][i].sub("@", " at ")}&gt;"
        if len > 1 and i < (len -1)
            out += " or "
        end
    end
    out
end

def mirror_status()
    cachefile = "tmp#{File::SEPARATOR}mirror.cache"
    now = Time.now
    yesterday = now - (60*60*24)
    FileUtils.mkdir_p "tmp"
    if File.exists? (cachefile) and File.mtime(cachefile) > yesterday
        return YAML.load_file(cachefile)
    end
    h = {}
    h[:last_updated] = Time.now.iso8601
    res = []
    for country in config[:mirrors]
        for mirror in country.values.first
            status = {}
            status[:url] = mirror[:https_mirror_url]
            status[:main] = (check_mirror_url(mirror[:https_mirror_url]) == "1") ? "yes" : "no"
            url = mirror[:https_mirror_url]
            url += "/" unless url.end_with? "/"
            ["release", "devel"].each do |version|
                numeric_version = config["#{version}_version".to_sym]
                url_to_check = "#{url}packages/#{numeric_version}/bioc/src/contrib/PACKAGES"
                #puts "URL: " + url_to_check
                begin
                    result = (check_mirror_url(url_to_check) == "1")
                rescue
                    result = false
                end
                status[version.to_sym] = result ? "yes" : "no"
            end

            res << status
        end
    end
    h[:status] = res
    File.open(cachefile, 'w') {|f| f.write h.to_yaml }
    h
end

def get_build_report_link(package)
    repo = package[:repo].sub(/\/$/, "")
    repo = repo.sub "/", "-"
    version = package[:bioc_version_num]
    package_name = package[:Package]
    "https://bioconductor.org/checkResults/#{version}/#{repo}-LATEST/#{package_name}/"
end


def get_build_results(package)
  return nil unless %w(bioc/ data/experiment/ workflows/).include? package[:repo]
  return nil unless current? package
  # return nil unless [config[:release_version],
  #   config[:devel_version]].include? package[:bioc_version_num]
  # build_dbs_dir = File.join(%w(tmp build_dbs))
  repo = package[:repo].sub(/\/$/, "").sub("/", "-")
  h = {config[:release_version] => 'release',
    config[:devel_version] => 'devel'}
  version = h[package[:bioc_version_num]]
  res = {}
  res[:report_url] = "https://bioconductor.org/checkResults/#{version}/#{repo}-LATEST/#{package[:Package]}/"
  res[:repo] = repo
  res [:version] = version
  res
end

# is package in release or devel?
def current? (package)
    [config[:release_version],
    config[:devel_version]].include? package[:bioc_version_num]
end

def on_fhcrc_network?
  io = IO.popen('hostname')
  fh = io.readlines.first.strip
  fh.end_with? "fhcrc.org"
end

def coverage_color(coverage)
  if coverage == "unknown"
    return "AA0088"
  end

  coverage = coverage.to_i
  if coverage >= 90
    return "green"
  elsif coverage >= 75
    return "yellow"
  end
  return "red"
end

def coverage_url(package)
  pkgname = package[:Package]
  branch = "devel"
  if package[:bioc_version_num] == config[:release_version]
    branch = "release-#{config[:release_version]}"
  end
  dirname = "devel"
  dirname = "release" if branch =~ /^release/
  shield = File.join("assets", "shields", "coverage",
    dirname, "#{pkgname}.svg")
  unless File.exists? shield
    return "https://codecov.io/github/Bioconductor-mirror/#{pkgname}/branch/#{branch}"
  end
  content = File.readlines(shield).first
  if (content =~ /unknown/)
    return "/developers/how-to/unitTesting-guidelines/#coverage"
  else
    return "https://codecov.io/github/Bioconductor-mirror/#{pkgname}/branch/#{branch}"
  end
end

def coverage_svg_url(package)
  vers = "devel"
  if package[:bioc_version_num] == config[:release_version]
    vers = "release"
  end
  "/shields/coverage/#{vers}/#{package[:Package]}.svg"
end

def get_short_url(package, full=false)
  url =  "/packages/#{package[:Package]}/"
  if full
    "https://bioconductor.org#{url}"
  else
    url
  end
end

def get_github_url(package)
  if package[:bioc_version_num] == config[:devel_version]
    branch = 'devel'
  else
    branch = "release-#{package[:bioc_version_num]}"
  end
  url = "https://github.com/Bioconductor-mirror/#{package[:Package]}"

  if branch == 'devel'
    url
  else
    "#{url}/tree/#{branch}"
  end
end

def check_mirror_url(url)
  uri = URI(url)
  http = Net::HTTP.new(uri.host, uri.port)
  if uri.port == 443
    http.use_ssl = true
    http.verify_mode = OpenSSL::SSL::VERIFY_NONE
  end
  begin
    response = http.head(uri.path)
    if response.code =~ /^2/
      "1"
    else
      "0"
    end
  rescue
    return "0"
  end
end

def get_url_from_item_identifier(identifier)
  # from /help/bioc-views/package-pages/3.2/bioc/a4/"
  # to   /packages/3.2/bioc/html/a4.html
  out = identifier.to_s.sub("/help/bioc-views/package-pages", "/packages")
  out = out.sub("/bioc/", "/bioc/html/")
  out = out.sub(/\/$/, ".html")
  out
end

def urlescape(url)
  URI.escape(url)
end

def get_socialized_url(identifier)
  url = "https://bioconductor.org#{get_url_from_item_identifier identifier}"
  CGI::escape(url)
end

def get_social_title(item, package)
  urlescape "#{package[:Package]}:#{item[:Title]}"
end

def package_is_release(package)
  if package[:bioc_version_num] == config[:release_version]
    true
  else
    false
  end
end

def package_is_devel(package)
  if package[:bioc_version_num] == config[:devel_version]
    true
  else
    false
  end
end

def package_has_archive(package)
  if (package_is_devel(package))
    return false
  end
  if !(package[:repo] == "bioc/")
    return false
  end
  version = package[:bioc_version_num]
  url = "https://bioconductor.org/packages/#{version}/bioc/src/contrib/Archive/#{package[:Package]}/"
  uri = URI(url)
  http = Net::HTTP.new(uri.host, uri.port)
  if uri.port == 443
    http.use_ssl = true
    http.verify_mode = OpenSSL::SSL::VERIFY_NONE
  end
  response = http.head(uri.path)
  valid_url = response.is_a?(Net::HTTPSuccess) || response.is_a?(Net::HTTPRedirection)
  if valid_url
    return true
  else
    return false
  end
end

def get_archive_url(package, text=false)
  version = package[:bioc_version_num]
  repo = package[:repo]
  url = "/packages/#{version}/#{repo}src/contrib/Archive/#{package[:Package]}/"
  if text
    "Source Archive"
  else
    url
  end
end

def get_last_git_commits(release=true)
  if release
    url = "https://master.bioconductor.org/developers/rss-feeds/gitlog.release.xml"
  else
    url = "https://master.bioconductor.org/developers/rss-feeds/gitlog.xml"
  end
  begin
    xml = HTTParty.get(url).parsed_response["rss"]["channel"]["item"]
    tbl_str = "<table>\n\n"
    uni_pkg = []
    dx = 0
    path = Dir.pwd
    manifest_path = "#{path}/../manifest/software.txt"
    manifest = File.open("#{manifest_path}").read
    lines = manifest.split("\n").drop(1) - [""]
    pkgs = lines.map {|item| item.gsub("Package: ", "")}
    while uni_pkg.length < 10 do
      item = xml[dx]
      if not uni_pkg.include?(item["title"])
        if pkgs.include?(item["title"])
          uni_pkg.push(item["title"])
          line = "<tr><td><a href="+item["link"]+">"+item["title"]+"</a></td><td>"+item["pubDate"]+"</td></tr>"
          tbl_str += line
        end
      end
      dx = dx + 1
    end
    tbl_str += "</table>"
    tbl_str
  rescue Exception => ex
    tbl_str = "<table>\n<tr><td> Can't read / no records in rss feed, not report last git commit time </td></tr>\n"
    i = 0
    while i < 9
      line = "<tr><td> . </td></tr>\n"
      tbl_str += line
      i += 1
    end
    tbl_str += "</table>"
    tbl_str
  end

end

def recent_spb_builds
    begin
      HTTParty.get("http://staging.bioconductor.org:8000/recent_builds").body
    rescue Exception => ex
      "Can't connect to staging.bioconductor.org, not building dashboard"
    end
end

def recent_spb_submissions
  client = Octokit::Client.new(:access_token => ENV["GITHUB_TOKEN"])
  issues = client.issues 'Bioconductor/Contributions'
  issue_names = []
  issue_url = []
  issues.each do |item|
    issue_names.push(item["title"])
    issue_url.push(item["html_url"])
  end
  dx = 0
  max_vl = [issue_names.length, 20].min
  tbl_str = "<table>\n"
  while dx < max_vl do
    line = "<tr><td><a href="+issue_url[dx]+">"+issue_names[dx]+"</a></td></tr>\n"
    tbl_str += line
    dx = dx + 1
  end
  tbl_str += "</table>"
  tbl_str
end

def latest_packages(repo)

  repo_text = repo
  if (repo_text == "bioc")
    repo_text = "software"
  end
  base_url = repo.sub("-", "/")

  path = Dir.pwd
  if (!File.directory?("#{path}/../manifest/"))
    tbl_str = "<table>\n<tr><td> file not found. skipping information </td></tr>\n"
    i = 0 
    while i < 9
      line = "<tr><td> . </td></tr>\n"
      tbl_str += line
      i += 1
    end
     tbl_str += "</table>"   
  else
    manifest_path = "#{path}/../manifest/#{repo_text}.txt"
    manifest = File.open("#{manifest_path}").read

    lines = manifest.split("\n").drop(1) - [""]
    pkgs = lines.map {|item| item.gsub("Package: ", "")}
    pkgs.reverse!

    trunc_pkgs = pkgs[0..9]
  
    tbl_str = "<table>\n"
    trunc_pkgs.each do |pkg|
      short_url = "https://bioconductor.org/packages/#{pkg}/"
      long_url = "http://bioconductor.org/packages/devel/#{base_url}/html/#{pkg}.html"

      url = URI.parse(long_url)
      req = Net::HTTP.new(url.host, url.port)
      res = req.request_head(url.path)
      if res.code == "200"
          line = "<tr><td><a href= #{short_url} > #{pkg} </a></td></tr>\n"
      else
          line = "<tr><td> #{pkg} </td></tr>\n"
      end
      tbl_str += line
    end
    tbl_str += "</table>"
    tbl_str
  end
  tbl_str
end

# this is an older function that is broken
# but needed currently to rake website
def recent_packages()
  begin
    xml = HTTParty.get("http://master.bioconductor.org/rss/new_packages.rss").body
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
