#!/usr/bin/env ruby

require 'pp'
require 'fileutils'
require 'yaml'
require 'open3'
require 'tempfile'

include FileUtils
include Open3

class CranSearchIndexer


  def initialize(args)
    
    unless args.size == 4
      puts "invalid arguments!"
      exit 1
    end
    directory_to_index, path_to_cache_file, path_to_output_script, site_url = args
    
    
    script_file_name = "#{path_to_output_script}/index_cran.sh"
    FileUtils.rm_f script_file_name
    
    script_file = File.open(script_file_name, "w")
    
    java_home = ENV["JAVA_HOME"]
    post_jar_home = "#{ENV["SOLR_HOME"]}/example/exampledocs"
    
    cache = {}
    cachecopy = {}
    cache_exists = false

    curl_path = `which curl`.chomp



    cachefilename = "#{path_to_cache_file}/cran_search_indexer_cache.yaml"

    if File.exists?(cachefilename) #todo also make sure file is not empty
      cache = YAML.load_file(cachefilename)
      cache_exists = true
    end

    allurls, allfiles = get_list_of_files_to_index()

    url = "http://localhost:8983/solr/update"

    script_file.puts "#!/bin/sh\n"

    pwd = FileUtils.pwd()
    FileUtils.cd(directory_to_index) # necessary?
    
    
    allfiles.each_with_index do |file, i|
      cachecopy[file] = 1
      mtime = File.lstat(file).mtime().to_i()
      if ( \
        (!cache_exists) \
        or \
        (cache_exists and cache.has_key?(file) and mtime > cache[file]) \
        or \
        (cache_exists and !cache.has_key?(file)) \
      )
        puts "adding #{file} to indexing script"
        cmd = %Q(#{curl_path} -s "#{url}/extract?literal.id=#{allurls[i]}&commit=false&boost.text=2" -F "myfile=@#{directory_to_index}#{file}")
        script_file.puts "echo '#{file}'"
        script_file.puts cmd
      end
      cache[file] = mtime

    end

    FileUtils.cd pwd # necessary?

    # deal with files that have been removed since last run
    if (cache_exists)
      a = cache.keys.sort
      b = cachecopy.keys.sort
      to_be_deleted = a - b 
      unless to_be_deleted.empty?
        script_file.puts "cd #{post_jar_home}" 
        #puts "to be deleted:"
        #puts to_be_deleted.join("\n")
        for item in to_be_deleted
          script_file.puts "echo 'deleting #{item}...'"
          script_file.puts %Q(#{java_home}/bin/java -Ddata=args -Dcommit=no -jar ./post.jar "<delete><id>#{item}</id></delete>")
          cache.delete(item)
        end
        script_file.puts "echo 'committing deletions...'"
        script_file.puts %Q(#{java_home}/bin/java -jar ./post.jar) #commit
        script_file.puts "cd #{pwd}"
      end
    end
    


    cmd = "curl -s #{url} --data-binary '<commit/>' -H 'Content-type:text/xml; charset=utf-8'"
    script_file.puts "echo 'committing changes...'"
    script_file.puts cmd
    script_file.close()
    File.chmod(0777, script_file_name)
    
    
    File.open(cachefilename, "w") do |cachefile|
      YAML.dump(cache, cachefile)
    end
    
  end
  

  
  def get_list_of_files_to_index()
    urls = []
    files = []
    index_file = File.open("/extra/www/cran-mirror/web/packages/available_packages_by_name.html")
    while (line = index_file.gets)
      next unless line.start_with? %Q(<td><a href="../../web/packages/)
      pkg_name = line.gsub(%Q(<td><a href="../../web/packages/),"").split("/").first
      urls.push  "http://cran.r-project.org/web/packages/#{pkg_name}/"
      files.push("/extra/www/cran-mirror/web/packages/#{pkg_name}/index.html")
    end
    return urls, files
  end


end  

