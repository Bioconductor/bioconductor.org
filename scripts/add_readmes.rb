#!/usr/bin/env ruby

require 'rubygems'
require 'json'
require 'pp'

def add_readmes(json_dir, version, repo)
  working_dir = "#{json_dir}/#{version}/#{repo}"
  puts "getting READMEs for #{repo} version #{version}..."
  system %Q(ssh webadmin@bioconductor.org "find /extra/www/bioc/packages/#{version}/bioc/readmes" > #{working_dir}/readmes.txt)
  return unless (test(?f, "#{working_dir}/readmes.txt"))
  f = File.open("#{working_dir}/readmes.txt")
  lines = f.readlines()
  f.close()
  pkgs_with_readmes = []
  for line in lines
    next unless line =~ /README/
    segs = line.split("/")
    segs.pop
    pkgs_with_readmes.push segs.last
  end
  
  pkg_file = "#{working_dir}/packages.json"
  f = File.open(pkg_file)
  lines = f.readlines
  f.close()
  json = lines.join("\n")
  json_obj = JSON.parse(json)
  json_obj.each_pair do |key, value|
    if  pkgs_with_readmes.include? key
      puts "adding readme for #{key}..."
      json_obj[key]['hasREADME'] = true
    end
  end
  
  outfile = File.open("#{working_dir}/packages.json", "w")
  
  outfile.puts JSON.pretty_generate(json_obj)
  outfile.close
  
end