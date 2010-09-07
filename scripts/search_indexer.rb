#!/usr/bin/env ruby

require 'pp'
require 'fileutils'
require 'yaml'

include FileUtils

class SearchIndexer

  def initialize()
    cache = {}
    cache_exists = false

    pwd = FileUtils.pwd()


    cachefilename = "#{pwd}/search_indexer_cache.yaml"

    if File.exists?(cachefilename) #todo also make sure file is not empty
      cache = YAML.load_file(cachefilename)
      cache_exists = true
    end

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
      if ((!cache_exists) or (cache_exists && cache.has_key?(file) && mtime > cache[file]))
        nice_name = file.gsub(/index\.html$/,"")

        puts %Q(echo "indexing #{nice_name}")
        cmd = %Q(curl -s "#{url}/extract?literal.id=#{nice_name}&commit=false" -F "myfile=@#{dirtoindex}#{cleanfile}")
        puts cmd
        #  result = `#{cmd}`
        #  puts result
        
      end
      cache[file] = mtime


    end



    cmd = "curl -s #{url} --data-binary '<commit/>' -H 'Content-type:text/xml; charset=utf-8'"
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
