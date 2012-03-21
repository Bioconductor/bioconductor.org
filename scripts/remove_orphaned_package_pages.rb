#!/usr/bin/env ruby

if ARGV.size != 1
  puts "supply web root"
  exit
end

webroot = ARGV.first

require 'rubygems'
require 'json'
require 'pp'
require 'fileutils'

json_root = "#{webroot}/packages/json"

file_list = `find #{json_root}|grep "packages\\.json$"`
files = file_list.split("\n")

deletion_list = []

for file in files
  
  path = file.gsub("#{json_root}/", "")
  segs = path.split("/")
  version = segs.shift
  segs.pop
  repo = segs.join("/")
  
  #puts "\nversion = #{version}, repo = #{repo}\n"
  
  f = File.open(file)
  json = f.readlines.join("\n")
  begin
    obj = JSON.parse(json)
  rescue JSON::ParserError => ex
    puts("failed to parse JSON in #{file}, skipping...")
    next
  end
  packages = obj.keys
  next if packages.empty? # for safety
  
  dir = "#{webroot}/packages/#{version}/#{repo}/html"
  
  next unless (File.exists? dir and File.directory? dir)
  
  existing_files = `ls -1 #{dir}`.split("\n")
  
  existing_packages = []
  for pkg in existing_files
    existing_packages.push pkg.sub(/\.html$/, "")
  end
  
  files_to_delete = existing_packages - packages
  
  next if files_to_delete.empty?

  for item in files_to_delete
    next if item == "package-detail.css"
    next if item == "PAN"
    deletion_list.push "#{dir}/#{item}.html"
  end
  
end

for item in deletion_list
  puts "deleting #{item}"
  FileUtils.rm item
end