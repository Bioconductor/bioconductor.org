#!/usr/bin/env ruby

require 'pp'
require 'fileutils'
require 'yaml'
require 'open3'
require 'tempfile'

include FileUtils
include Open3

class SearchIndexer

  def self.is_solr_running?()
    orig_stdout = STDOUT.clone
    orig_stderr = STDERR.clone
    
    $stdout.reopen("/dev/null", "w")
    $stderr.reopen("/dev/null", "w")
    
    
    cmd = "curl http://localhost:8983"
    result = system(cmd)
    
    $stdout.reopen orig_stdout
    $stderr.reopen orig_stderr
    result
  end

  def get_boost(url)
    release_regex =  /\/packages\/(release|#{@release_version})\/(bioc|data\/annotation|data\/experiment)\/html/
    devel_regex =  /\/packages\/(devel|#{@devel_version})\/(bioc|data\/annotation|data\/experiment)\/html/
    return 100 if url =~ release_regex
    return 50 if url =~ devel_regex
    #return 3 if url =~ /\/packages\/release/ or url =~ /\/packages\/#{@release_version}/
    #return 2 if url =~ /\/packages\/devel/ or url =~ /\/packages\/#{@devel_version}/
    return 1
  end
  
  
  def throw_out_bad_files(file_list)
    ret = []
    for file in file_list
       #e.g. /packages/release/bioc/html or /packages/2.6/Software.html
      if (file =~ /^\.\/packages/ && file =~ /\/html\/|\.html$/)
        puts "KEEPING #{file}"
        ret.push file
      else
        puts "THROWING OUT #{file}"
      end
    end
    ret
  end

  def initialize(args)
    
    site_config = YAML.load_file("./config.yaml")
    @release_version = site_config["release_version"]
    @devel_version = site_config["devel_version"]
    
    unless args.size == 4
      puts "invalid arguments!"
      exit 1
    end
    directory_to_index, path_to_cache_file, path_to_output_script, site_url = args
    
    
    script_file_name = "#{path_to_output_script}/index.sh"
    FileUtils.rm_f script_file_name
    
    script_file = File.open(script_file_name, "w")
    
    java_home = ENV["JAVA_HOME"]
    post_jar_home = "#{ENV["SOLR_HOME"]}/example/exampledocs"
    
    cache = {}
    cachecopy = {}
    cache_exists = false

    curl_path = `which curl`.chomp



    cachefilename = "#{path_to_cache_file}/search_indexer_cache.yaml"

    if File.exists?(cachefilename) #todo also make sure file is not empty
      cache = YAML.load_file(cachefilename)
      cache_exists = true
    end



    #allfilestr = `find .`


    extensions_to_index = %w{html doc pdf R}
    regex = Regexp.new("\\." + extensions_to_index.join("$|\\.") + "$")


    #allfiles = allfilestr.split("\n")

    allfiles = get_list_of_files_to_index(directory_to_index, site_url)

    goodfiles = allfiles.grep(regex)
    
    ###goodfiles = throw_out_bad_files(goodfiles)
    #exit if true

    url = "http://localhost:8983/solr/update"

    script_file.puts "#!/bin/sh\n"

    pwd = FileUtils.pwd()
    FileUtils.cd(directory_to_index)
    
    
    goodfiles.each do |file|
      
#      puts file
      
      cachecopy[file] = 1
      cleanfile = file.gsub(/^\./,"")
      bail = false
      begin
        bail = false
        mtime = File.lstat(file).mtime().to_i()
      rescue Exception => ex # TODO, make a note of files that don't exist
        bail = true
      end
      next if bail
      if ( \
        (!cache_exists) \
        or \
        (cache_exists and cache.has_key?(file) and mtime > cache[file]) \
        or \
        (cache_exists and !cache.has_key?(file)) \
      )
        nice_name = cleanfile.gsub(/index\.html$/,"")
        puts "adding #{nice_name} to indexing script"
        boost = get_boost(nice_name)
        boost_frag = (boost==1) ? "" : "&boost.text=#{boost}"
        script_file.puts %Q(echo "indexing #{nice_name}")
        cmd = %Q(#{curl_path} -s "#{url}/extract?literal.id=#{nice_name}&commit=false#{boost_frag}" -F "myfile=@#{directory_to_index}#{cleanfile}")
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

    FileUtils.cd pwd

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
          script_file.puts %Q(#{java_home}/bin/java -Ddata=args -Dcommit=no -jar ./post.jar "<delete><id>#{nice_name(item)}</id></delete>")
          cache.delete(item)
        end
        script_file.puts "echo 'committing deletions...'"
        script_file.puts %Q(#{java_home}/bin/java -jar ./post.jar) #commit
        script_file.puts "cd #{pwd}"
      end
    end
    


    cmd = "curl -s #{url} --data-binary '<commit/>' -H 'Content-type:text/xml; charset=utf-8'"
    #system(cmd)
    script_file.puts "echo 'committing changes...'"
    script_file.puts cmd
    script_file.close()
    File.chmod(0777, script_file_name)
    
    #`#{cmd}`

    #FileUtils.cd(pwd)
    
    
    File.open(cachefilename, "w") do |cachefile|
      YAML.dump(cache, cachefile)
    end
    
  end
  
  def index(file)
  end

  def nice_name(name)
    name.gsub(/^\./,"").gsub(/index\.html$/,"")
  end
  
  
  def get_list_of_files_to_index(directory_to_index, site_url)
    f = File.open("links.txt")
    lines = f.readlines
    for line in lines
      line.chomp!
    end
    lines
  end
  
  def get_list_of_files_to_index_old(directory_to_index, site_url)
    spider_tmpfile = Tempfile.new("spider_tmpfile")
    spider_output_file = spider_tmpfile.path()
    spider_tmpfile.close()
    
    puts "spidering..."
    system(%Q(wget -r --spider -U "link check with wget" -o #{spider_output_file} #{site_url}))
    
    pwd = FileUtils.pwd()
    
    f = File.open(spider_output_file)


    url_hash = {}

    urls = []
    broken_links = []

    broken_link_mode = false

    while (line = f.gets)
      line.chomp!
      if line =~ /^--/
        url_hash[line.split(/\s/).last()] = 1
      end


      if line =~ / broken links\.$/
        broken_link_mode = true
      end

      if broken_link_mode and line =~ /^http:|^https:/i
        broken_links.push(line)
      end


    end

    urls = url_hash.keys.sort

    clean_urls = []

    urls.each do |url|
      clean_urls.push url.gsub(/^#{site_url}/, ".").gsub(/\/$/, "/index.html").gsub("%20"," ")
    end

    filtered = (urls - broken_links)


    file_list_file = Tempfile.new("file_list")
    file_list_filename = file_list_file.path()
    file_list_file.close()
    
    FileUtils.cd directory_to_index
    system("find -L . -type f > #{file_list_filename}")
    FileUtils.cd(pwd)
    

    f = File.open(file_list_filename)
    all = f.readlines.map{|i|i.chomp}
    #puts "!#{all[1]}!"

    orphans = all - clean_urls
    non_orphans = all - orphans

    #puts non_orphans.join("\n")
        

    #    puts "#{filtered.size} total links found, #{broken_links.size} broken links found, #{filtered.size} filtered links found"
    #    puts "#{all.size} files on filesystem, #{orphans.size} orphans, non-orphans = #{non_orphans.size}"
    
    
    
    
    non_orphans
    
    
  end

  # url = 
  # http://localhost:8983/solr/select?indent=on&version=2.2&q=oligonucleotide&fq=&start=0&rows=10&fl=id,score,title&qt=standard&wt=json&explainOther=&hl=on&hl.fl=






end  

#unless (ARGV.length == 3)
#  puts "usage: #{$0} [directory-to-index] [path-to-cache-file] [path-to-output-script] [site-url]"
#  exit 1
#end

#si = SearchIndexer.new(ARGV)
 