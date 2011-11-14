#!/usr/bin/env ruby

require 'rubygems'
require 'json'
require 'pp'

def add_extra_files(json_dir, version, repo, type)
  working_dir = "#{json_dir}/#{version}/#{repo}"
  if type == :readme
    what = "READMEs"
    dir = "readmes"
    file = "readmes.txt"
    jsonkey = 'hasREADME'
    pattern = /README/
  elsif type == :news
    what = "NEWS"
    dir = "news"
    file = "news.txt"
    jsonkey = 'hasNEWS'
    pattern = /NEWS/
  end
  puts "getting #{what} for #{repo} version #{version}..."
  system %Q(ssh webadmin@bioconductor.org "find /extra/www/bioc/packages/#{version}/bioc/#{dir}" > #{working_dir}/#{file})
  return unless (test(?f, "#{working_dir}/#{file}"))
  f = File.open("#{working_dir}/#{file}")
  lines = f.readlines()
  f.close()
  pkgs_with_goodies = []
  for line in lines
    next unless line =~ pattern
    segs = line.split("/")
    segs.pop
    pkgs_with_goodies.push segs.last
  end
  
  pkg_file = "#{working_dir}/packages.json"
  f = File.open(pkg_file)
  lines = f.readlines
  f.close()
  json = lines.join("\n")
  json_obj = JSON.parse(json)
  json_obj.each_pair do |key, value|
    if  pkgs_with_goodies.include? key
      json_obj[key][jsonkey] = true
    end
  end
  
  outfile = File.open("#{working_dir}/packages.json", "w")
  
  outfile.puts JSON.pretty_generate(json_obj)
  outfile.close
  
end

def add_readmes(json_dir, version, repo)
  add_extra_files(json_dir, version, repo, :readme)
end

def add_news(json_dir, version, repo)
  add_extra_files(json_dir, version, repo, :news)
end