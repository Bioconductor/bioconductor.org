#!/usr/bin/env ruby

# this file exists because I can't figure out how biocViews works. I will figure it out and make this file go away.

unless ARGV.size == 1
  puts "usage: #{$0} repospath"
  exit
end

require 'pp'
require 'rubygems'
require 'json'

file_str = `find #{ARGV.first.gsub(/\/$/,"")} |grep -i "\.Rnw$"`

files = file_str.split("\n")

hsh = {}

for file in files
  segs = file.split("/")
  base_file = segs.last.gsub(/rnw$/i, "pdf")
  pkg_name = segs[segs.length - 4]
  segs.pop
  dir = segs.join("/")

  
  hsh[pkg_name] = {:vignetteTitles => [], :vignetteScripts => [], :vignetteFiles => []} unless hsh.has_key?(pkg_name)
  cmd = %Q(grep VignetteIndexEntry #{file})
  result = `#{cmd}`
  result =~ /\{([^}]*)\}/
  title = $1
  
  
  next if title.nil?
  hsh[pkg_name][:vignetteTitles].push title
  rfile = base_file.sub(/\.pdf$/, ".R")
  
  if (test(?f, "#{dir}/#{base_file}"))
      hsh[pkg_name][:vignetteFiles].push base_file
    else
      hsh[pkg_name][:vignetteFiles].push ""
  end
  
  
  #puts "looking for #{dir}/#{rfile}"
  if (test(?f, "#{dir}/#{rfile}"))
    hsh[pkg_name][:vignetteScripts].push rfile
  else
    hsh[pkg_name][:vignetteScripts].push ""
  end
end


json = JSON.pretty_generate(hsh)
puts json