#!/usr/bin/env ruby

# this file exists because I can't figure out how biocViews works. I will figure it out and make this file go away.

unless ARGV.size == 1
  puts "usage: #{$0} repospath"
  exit
end

require 'pp'
require 'rubygems'
require 'json'

file_str = `find #{ARGV.first.gsub(/\/$/,"")} |grep -i "\.pdf$"`

files = file_str.split("\n")

hsh = {}

for file in files
  segs = file.split("/")
  
  rnw_file = file.sub(/\.pdf$/i, ".Rnw")
  
  
  base_file = segs.last#.gsub(/rnw$/i, "pdf")
  rfile = base_file.sub(/\.pdf$/, ".R")
  pkg_name = segs[segs.length - 4]
  segs.pop
  dir = segs.join("/")

  
  hsh[pkg_name] = {:vignetteTitles => [], :vignetteScripts => [], :vignetteFiles => []} unless hsh.has_key?(pkg_name)
  
  
  if (test(?f, rnw_file))
    cmd = %Q(grep VignetteIndexEntry #{rnw_file})
    result = `#{cmd}`
    result =~ /\{([^}]*)\}/
    title = $1

    title = base_file if title.nil?
    #next if title.nil?
    hsh[pkg_name][:vignetteTitles].push title
  else
    hsh[pkg_name][:vignetteTitles].push base_file
  end
  
  
  hsh[pkg_name][:vignetteFiles].push base_file
  
  
  #puts "looking for #{dir}/#{rfile}"
  if (test(?f, "#{dir}/#{rfile}"))
    hsh[pkg_name][:vignetteScripts].push rfile
  else
    hsh[pkg_name][:vignetteScripts].push ""
  end
end


json = JSON.pretty_generate(hsh)
puts json