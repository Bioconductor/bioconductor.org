include Nanoc3::Helpers::Text
include Nanoc3::Helpers::Rendering
include Nanoc3::Helpers::Breadcrumbs

require 'time'


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
