#!/usr/bin/env ruby

require 'pp'
require 'fileutils'
require 'yaml'
require 'open3'

include FileUtils
include Open3

class SearchIndexer

  def initialize()
    cache = {}
    cache_exists = false

    curl_path = `which curl`.chomp

    pwd = FileUtils.pwd()


    cachefilename = "#{pwd}/search_indexer_cache.yaml"

    if File.exists?(cachefilename) #todo also make sure file is not empty
      cache = YAML.load_file(cachefilename)
      cache_exists = true
    end


    #puts "cache_exists= #{cache_exists}"; exit if true

    dirtoindex = "./output"

    FileUtils.cd(dirtoindex)


    allfilestr = `find .`


    extensions_to_index = %w{html doc pdf R}
    regex = Regexp.new("\\." + extensions_to_index.join("$|\\.") + "$")


    allfiles = allfilestr.split("\n")


    goodfiles = allfiles.grep(regex)

    url = "http://localhost:8983/solr/update"

    puts "#/bin/sh\n"
    
    goodfiles.each do |file|
      cleanfile = file.gsub(/^\./,"")
      mtime = File.stat(file).mtime().to_i()
      if ( \
        (!cache_exists) \
        or \
        (cache_exists and cache.has_key?(file) and mtime > cache[file]) \
        or \
        (cache_exists and !cache.has_key?(file)) \
      )
        nice_name = cleanfile.gsub(/index\.html$/,"")

        puts %Q(echo "indexing #{nice_name}")
        cmd = %Q(#{curl_path} -s "#{url}/extract?literal.id=#{nice_name}&commit=false" -F "myfile=@#{dirtoindex}#{cleanfile}")
        puts cmd
        #result = system(cmd)
        #{}`#{cmd}`
        #puts "#{result}\t#{nice_name}"
        
        #puts "#{nice_name}"
        #Open3.popen3(cmd) do |stdin, stdout, stderr|
        #  puts "stderr:"
        #  puts stderr.readlines
        #  puts "stdout:"
        #  puts stdout.readlines
        #  puts
        #end
        
      end
      cache[file] = mtime
      # todo - deal with files that have been removed since last run

    end



    cmd = "curl -s #{url} --data-binary '<commit/>' -H 'Content-type:text/xml; charset=utf-8'"
    #system(cmd)
    puts cmd
    #`#{cmd}`

    FileUtils.cd(pwd)
    
    
    File.open(cachefilename, "w") do |cachefile|
      YAML.dump(cache, cachefile)
    end
    
  end
  
  def index(file)
  end


  # url = 
  # http://localhost:8983/solr/select?indent=on&version=2.2&q=oligonucleotide&fq=&start=0&rows=10&fl=id,score,title&qt=standard&wt=json&explainOther=&hl=on&hl.fl=






end  

si = SearchIndexer.new
