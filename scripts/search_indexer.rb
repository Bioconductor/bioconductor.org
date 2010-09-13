#!/usr/bin/env ruby

require 'pp'
require 'fileutils'
require 'yaml'
require 'open3'

include FileUtils
include Open3

class SearchIndexer

  def self.is_solr_running?()
    cmd = "curl http://localhost:8983"
    system(cmd)
  end

  def initialize(args)
    
    directory_to_index, path_to_cache_file, path_to_output_script = args
    
    
    script_file_name = "#{path_to_output_script}/index.sh"
    FileUtils.rm_f script_file_name
    
    script_file = File.open(script_file_name, "w")
    
    java_home = ENV["JAVA_HOME"]
    post_jar_home = "/Users/dante/apache-solr-1.4.1/example/exampledocs" #todo - customize for staging & production  
    
    cache = {}
    cachecopy = {}
    cache_exists = false

    curl_path = `which curl`.chomp

    pwd = FileUtils.pwd()


    cachefilename = "#{path_to_cache_file}/search_indexer_cache.yaml"

    if File.exists?(cachefilename) #todo also make sure file is not empty
      cache = YAML.load_file(cachefilename)
      cache_exists = true
    end


    FileUtils.cd(directory_to_index)


    allfilestr = `find .`


    extensions_to_index = %w{html doc pdf R}
    regex = Regexp.new("\\." + extensions_to_index.join("$|\\.") + "$")


    allfiles = allfilestr.split("\n")


    goodfiles = allfiles.grep(regex)

    url = "http://localhost:8983/solr/update"

    script_file.puts "#/bin/sh\n"
    
    goodfiles.each do |file|
      cachecopy[file] = 1
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

        script_file.puts %Q(echo "indexing #{nice_name}")
        cmd = %Q(#{curl_path} -s "#{url}/extract?literal.id=#{nice_name}&commit=false" -F "myfile=@#{directory_to_index}#{cleanfile}")
        script_file.puts cmd
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

    end

    # deal with files that have been removed since last run
    if (cache_exists)
      a = cache.keys.sort
      b = cachecopy.keys.sort
      to_be_deleted = a - b 
      return if to_be_deleted.empty?
      script_file.puts "cd #{post_jar_home}" 
      #puts "to be deleted:"
      #puts to_be_deleted.join("\n")
      for item in to_be_deleted
        script_file.puts %Q(#{java_home}/bin/java -Ddata=args -Dcommit=no -jar ./post.jar "<delete><id>#{nice_name(item)}</id></delete>")
        cache.delete(item)
      end
      script_file.puts %Q(#{java_home}/bin/java -jar ./post.jar) #commit
      script_file.puts "cd #{pwd}"
    end
    


    cmd = "curl -s #{url} --data-binary '<commit/>' -H 'Content-type:text/xml; charset=utf-8'"
    #system(cmd)
    script_file.puts cmd
    #`#{cmd}`

    FileUtils.cd(pwd)
    
    
    File.open(cachefilename, "w") do |cachefile|
      YAML.dump(cache, cachefile)
    end
    
  end
  
  def index(file)
  end

  def nice_name(name)
    name.gsub(/^\./,"").gsub(/index\.html$/,"")
  end

  # url = 
  # http://localhost:8983/solr/select?indent=on&version=2.2&q=oligonucleotide&fq=&start=0&rows=10&fl=id,score,title&qt=standard&wt=json&explainOther=&hl=on&hl.fl=






end  

#unless (ARGV.length == 3)
#  puts "usage: #{$0} [directory-to-index] [path-to-cache-file] [path-to-output-script]"
#  exit 1
#end

#si = SearchIndexer.new(ARGV)
