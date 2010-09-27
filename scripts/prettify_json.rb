#!/usr/bin/env ruby

unless ARGV.size == 1
  STDERR.puts "usage: #{$0} file"
  exit
end

require 'rubygems'
require 'json'
require 'fileutils'
require 'tempfile'

f = File.open(ARGV.first)
json = f.readlines.join("\n")
obj = JSON.parse(json)

pretty = JSON.pretty_generate(obj)

tmpfile = Tempfile.new("prettify_json_in_place")
tmpfile.puts pretty
tmpfile.close

FileUtils.rm ARGV.first
FileUtils.mv tmpfile.path, ARGV.first
