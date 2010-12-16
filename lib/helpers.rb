include Nanoc3::Helpers::Text
include Nanoc3::Helpers::Rendering
include Nanoc3::Helpers::Breadcrumbs
include Nanoc3::Helpers::XMLSitemap



require 'time'

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

def filter_emails(str)
  str.gsub(/<[^>]*>/,"").gsub("  "," ").gsub(" ,", ",")
end

#todo - show links for stuff that is not necessarily a bioconductor package 
#(e.g. CRAN or R)
def linkify(sym, package)
  items = package[sym]
  # the following key gets set in bioc_views.rb#items()
  key = "#{sym.to_s}_repo".to_sym
  repos = package[key]
  output = []
  
  items.each_with_index do |item, index|
    repo = repos[index]
    if (repo == false)
      output.push item
      next
    end
    output.push %Q(<a href="/help/bioc-views/#{package[:bioc_version_num]}/#{repo}/html/#{item}.html">#{item}</a>)
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
    links.push %Q(<a href="/help/bioc-views/#{package[:bioc_version_num]}/BiocViews.html#___#{bioc_view}#{version_fragment(package)}">#{bioc_view}</a>)
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
