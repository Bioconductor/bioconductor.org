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

def base_filename(path)
  #return nil if path.nil?
  path.split("/").last
end

# show only one windows build - the 64-bit (if available) or the 32 (because we have fat (dual-arch) packages now)
def windows_binary(package)
  win32 = package[:"win.binary.ver"]
  win64 = package[:"win64.binary.ver"]
  # assuming that package has both keys
  return win64 unless win64.nil?  or win64.empty?
  win32
end

def win_format(package)
  if windows_binary(package) =~ /64/
    "(32- &amp; 64-bit)"
  else
    "(32-bit only)"
  end
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


def linkify(sym, package)
  items = package[sym]
  key = "#{sym.to_s}_repo".to_sym
  repos = package[key]
  output = []
  
  items.each_with_index do |item, index|
    repo = repos[index]
    if (repo == false)
      output.push item
      next
    end
    output.push %Q(<a href="/help/bioc-views/#{package[:bioc_version_str]}/#{item}">#{item}</a>)
  end
  output.join(", ")
end

def doc_object(package)
  # return an array of hashes
  # [{:file => ..., :title => ..., :script => ...}]
  doc_obj = []
  if (@package[:vignetteTitles].respond_to? :join)
    @package[:vignetteTitles].each_with_index do |title, i|
      item = {}
      if (@package.has_key? :vignetteFiles)
        item[:file] = @package[:vignetteFiles][i]
        item[:title] = title
      else
        item[:file] = item[:title] = title
      end
      
      if (@package.has_key? :vignetteScripts and @package[:vignetteScripts][i] != "")
        item[:script] = @package[:vignetteScripts][i]
      end
      doc_obj.push item
    end
    return doc_obj
  elsif @package[:vignetteTitles].is_a? String and @package[:vignetteTitles] != "NA"
    item = {}
    item[:file] = item[:title] = @package[:vignetteTitles]
    item[:script] = ""
    doc_obj.push item
  else
#    item = {}
#    item[:file] = 
    return [:file => "", :title => "", :script => ""]
  end
  doc_obj
end

def bioc_views_links(package)
  links = []
  
  bioc_views = to_a(package[:biocViews])
  bioc_views.each do |bioc_view|
    links.push %Q(<a href="/help/bioc-views/?openNode=#{bioc_view}#{version_fragment(package)}">#{bioc_view}</a>)
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
