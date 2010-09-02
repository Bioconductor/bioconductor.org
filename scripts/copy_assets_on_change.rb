#!/usr/bin/env ruby


require 'rubygems'
require 'fsevent'

class CopyAssetsOnChange < FSEvent
  def on_change(directories)
    puts "Detected change in: #{directories.inspect}"
    system "rake copy_assets"
  end

  def start
    puts "watching #{registered_directories.join(", ")} for changes"
    super
  end
end

printer = CopyAssetsOnChange.new
printer.latency = 0.2
printer.watch_directories %W(#{Dir.pwd}/assets)
printer.start
