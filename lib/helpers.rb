include Nanoc3::Helpers::Text
include Nanoc3::Helpers::Rendering
include Nanoc3::Helpers::Breadcrumbs
include Nanoc3::Helpers::XMLSitemap
include Nanoc3::Helpers::HTMLEscape


require 'time'
require 'rubygems'
require 'httparty'
require 'yaml'
require 'pp'
require 'rexml/document'

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

def get_cran_packages()
  puts "Grabbing list of CRAN packages..."
  ## Is there a non-web-scraping way to get a list of CRAN packages?
  cran_packages = []
  data = begin
           HTTParty.get("http://cran.fhcrc.org/web/packages/",
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
  #return nil if path.nil?
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

def to_a(item)
  return item if item.respond_to? :join
  return [item]
end

def version_fragment(package)
  return "&version=#{package[:version_num]}&" if package[:bioc_version_str] == "devel"
  ""
end


## Should we really filter out maintainer email addresses from the package home page?
## These emails are available publically in svn (although this requires a little more
## work on the part of a spam-harvester, who would have to supply the appropriate
## login credentials).
def filter_emails(str)
  str.gsub(/<[^>]*>/,"").gsub("  "," ").gsub(" ,", ",")
end

def linkify(sym, package)
  items = package[sym]
  # the following key gets set in bioc_views.rb#items()
  key = "#{sym.to_s}_repo".to_sym
  repos = package[key]
  output = []
  
  items.each_with_index do |item, index|
    repo = repos[index]



    if (repo == false)
      if ($cran_packages.include?(item))
        output.push %Q(<a class="cran_package" href="http://cran.fhcrc.org/web/packages/#{item}/index.html">#{item}</a>)
      else
        output.push item
      end
      next
    end

    output.push %Q(<a href="/packages/#{package[:bioc_version_num]}/#{repo}/html/#{item}.html">#{item}</a>)
  end
  output.join(", ")
end

def doc_object(package)
  # return an array of hashes
  # [{:file => ..., :title => ..., :script => ...}]
  
  
  
  
  doc_obj = []
  
  
  if (package[:vignetteTitles].respond_to? :join)
    
    package[:vignetteTitles].each_with_index do |title, i|
      item = {}
      if (package.has_key? :vignetteFiles)
        item[:file] = (package[:vignetteFiles].respond_to? :join)  ? package[:vignetteFiles][i] : package[:vignetteFiles]
        item[:title] = title
      else
        item[:file] = item[:title] = title
      end
      
      if (package.has_key? :vignetteScripts and package[:vignetteScripts][i] != "")
        item[:script] = package[:vignetteScripts][i]
      end

      doc_obj.push item
    end
    
    doc_obj.sort!{|a,b|a[:title].downcase <=> b[:title].downcase}
    
    return doc_obj
  elsif package[:vignetteTitles].is_a? String and package[:vignetteTitles] != "NA"
    item = {}
    item[:file] = item[:title] = package[:vignetteTitles]
    item[:script] = ""
    doc_obj.push item
  else
    
    return [:file => "", :title => "", :script => ""]
  end

  doc_obj.sort!{|a,b|a[:title].downcase <=> b[:title].downcase}
  doc_obj
end

def bioc_views_links(package)
  links = []
  
  bioc_views = to_a(package[:biocViews])
  bioc_views.each do |bioc_view|
    links.push %Q(<a href="/packages/#{package[:bioc_version_num]}/BiocViews.html#___#{bioc_view}#{version_fragment(package)}">#{bioc_view}</a>)
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
  events.children.sort{|a, b| a[:start] <=> b[:start]}.select do |e|
    e[:end] >= Time.now.to_date
  end
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
    repo = "dataann-stats/data-annotation/"
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
  ver_crumb = Nanoc3::Item.new(nil, {:title => "Bioconductor #{ver}"}, "/packages/#{ver}/BiocViews.html")
  crumbs.push ver_crumb
  repo_crumb = Nanoc3::Item.new(nil, {:title => "#{repo.first} Packages"}, "/packages/#{ver}/#{repo.last}/")
  crumbs.push repo_crumb unless index_page
  crumbs.push last_crumb
  crumbs
end

def recent_packages()
  config = YAML.load_file("./config.yaml")
  devel_version = config["devel_version"]
  manifest_url = \
    "https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/bioc_#{devel_version}.manifest"
  begin
    raw_manifest = `curl -s -u readonly:readonly #{manifest_url}`
    manifest = raw_manifest.split("\n").reject{|i| i =~ /#/ or i.empty?}.reverse
    manifest = manifest.map{|i| i.gsub(/^Package: /, "")}
    new_packages = manifest[0..5]
    # todo  - make sure pkg home pages actually exist
    #puts;pp new_packages;puts
    return new_packages
  rescue Exception => ex
    return []
  end
end

def get_svn_commits()
  begin
    cmd = "svn --xml --username readonly --password readonly log -v --limit 10 " +
      "https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/"
    raw_xml = `#{cmd}`
    doc = REXML::Document.new(raw_xml)
    entries = []
    doc.elements.each("log/logentry") {|i| entries.push i}
    ret = []
    for entry in entries
      h = {}
      entry.elements.each('author'){|i| h[:author] = i.text}
      h[:revision] = entry.attributes["revision"]
      entry.elements.each('date'){|i| h[:date] = i.text}
      # todo convert date...
      entry.elements.each('msg'){|i| h[:msg] = i.text}
      paths = []
      entry.elements.each('paths/path') do |item|
        path = {}
        path[:action] = item.attributes['action']
        path[:path] = item.text
        path[:copyfrom_path] = item.attributes['copyfrom-path']
        path[:copyfrom_rev] = item.attributes['copyfrom-rev']
        paths.push path
      end
      h[:paths] = paths
      ret.push h
    end
    return ret
  rescue Exception => ex
    puts "caught exception trying to get svn log"
    pp ex
    return []
  end
end