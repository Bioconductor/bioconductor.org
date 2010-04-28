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
        :name => f
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
