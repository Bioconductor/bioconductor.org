#!/usr/bin/env ruby


# if the fsevent gem is properly installed, you should not need the following three lines:

#fsevent_dir = "~/src/ruby-fsevent"
#$LOAD_PATH.unshift(fsevent_dir + '/lib')
#$LOAD_PATH.unshift(fsevent_dir + '/ext')  

require 'rubygems'
require 'fsevent'

class PrintChange < FSEvent
  def on_change(directories)
    puts "Detected change in: #{directories.inspect}"
    system "nanoc compile"
  end

  def start
    puts "watching #{registered_directories.join(", ")} for changes"
    super
  end
end

printer = PrintChange.new
printer.latency = 0.2
printer.watch_directories %W(#{Dir.pwd}/content #{Dir.pwd}/layouts)
printer.start
