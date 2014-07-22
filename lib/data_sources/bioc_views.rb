class BiocViews < Nanoc3::DataSource

  require 'pp'
  require 'rubygems'
  require 'json'
  
  identifier :bioc_views
  
  def hard_coded_repos()
    {"bioc/" => "Software", "data/annotation/" => "AnnotationData", "data/experiment/" => "ExperimentData"}
  end

  def get_repos(version, repos)
    return hard_coded_repos() unless version == @site_config["devel_version"]
    h = {}
    
    for repo in @site_config["devel_repos"]
      key = "#{repo}/"
      h[key] = repos[key]
    end
    h
  end
  
  
  # todo - find out if there is a way to skip items() altogether if things are not found in up()
  def up
    @bad_packages = ["snpMatrix2"] # don't process these

    #@repos = {"bioc/" => "Software", "data/annotation/" => "AnnotationData", "data/experiment/" => "ExperimentData"}
    @repos = hard_coded_repos()
    
    @good_to_go = true
    
    @site_config = YAML.load_file("./config.yaml")
    
    
    for version in @site_config["versions"]
      @repos = get_repos(version, @repos)
      @repos.each_pair do |k,v|
        
        
        dir = "#{config[:json_dir]}/#{version}/#{k}"
        
        @good_to_go = false unless test(?f, "#{dir}/packages.json") 

        #todo remove
        #@good_to_go = false unless test(?f, "#{dir}/vignette_titles.json") 


      end
    end
    
    
    @all_packages = {}



    
    if @good_to_go
      for version in @site_config["versions"]
        hsh = {}
        
        @repos = get_repos(version, @repos)
        
        
        @repos.each_pair do |k,v|
          # todo remove this when all 2.9 repos are present
          # next if version == "2.9" and k != "bioc/"
          
          key = k.gsub(/\/$/, "")
          dir = "#{config[:json_dir]}/#{version}/#{k}"
          json_file = File.open("#{dir}/packages.json")
          obj = JSON.parse(json_file.readlines.join("\n"))
          #@all_packages[key] = obj
          hsh[key] = obj
          for bad in @bad_packages
            hsh.delete(bad)
          end
        end
        @all_packages[version] = hsh
      end
    else
      puts "BiocViews data source: json files not present, skipping initialization"
    end
    
      
  end

  
  
  def get_index_page(packages, repo, version)
    item = Nanoc3::Item.new("", {}, "all-#{repo}-#{version}")
    item[:rebase] = true
    rep = Nanoc3::ItemRep.new(item, :package_index_page)
    
    
    item[:package_index_page] = true
    
    info = []
    
    packages.each_pair do |k,v|
      hsh = {}
      hsh[:name] = k
      hsh[:Maintainer] = v["Maintainer"]
      hsh[:Title] = v["Title"]
      info.push hsh
    end

    
    item[:info] = info.sort{|a,b| a[:name].downcase <=> b[:name].downcase}
    
    
    item[:bioc_version_num] = version
    if (version == @site_config["release_version"])
      item[:bioc_version_str] = "Release"
    elsif (version == @site_config["devel_version"])
      item[:bioc_version_str] = "Development"
    else
      item[:bioc_version_str] = nil
    end
    
    
    item[:repo] = repo
    item[:title] = "#{version} #{repo} Packages"
    
    item[:subnav] = []
    item[:subnav].push({:include => "/_bioc_release_packages/"})
    item[:subnav].push({:include => "/_bioc_devel_packages/"})
    item[:subnav].push({:include => "/_bioc_older_packages/"})
    
    item
  end
  
  def items
    
    unless @good_to_go
      puts "BiocViews_DataSource: no JSON file(s) found. Package detail pages will not be built"
      return []
    end

    items = []
    
    link_list = [:Depends, :Imports, :Suggests, :dependsOnMe, :importsMe, :suggestsMe]
    
    for version in @site_config["versions"]
      @repos = get_repos(version, @repos)
      @repos.each_pair do |k,v|
        dir = "#{config[:json_dir]}/#{version}/#{k}"
        
        
        json_file = File.open("#{dir}/packages.json")



        packages = JSON.parse(json_file.readlines.join("\n"))
        for bad in @bad_packages
          packages.delete(bad)
        end

        items.push(get_index_page(packages, v, version))


        for package in packages.keys
          repo = k
          id = "/#{version}/#{repo}#{package}/"
          item = Nanoc3::Item.new("", packages[package], id)
          
          item[:rebase] = true
          item[:subnav] = []
          item[:subnav].push({:include => "/_workflows/"})
          item[:subnav].push({:include => "/_mailing_list/"})

          item[:title] = "#{item[:Package]}"
          if  version == @site_config["devel_version"]
            item[:title] += " (development version)"
          end
          item[:repo] = repo
          item[:bioc_version_num] = version
          if (version == @site_config["release_version"])
            item[:bioc_version_str] = "Release"
          elsif (version == @site_config["devel_version"])
            item[:bioc_version_str] = "Development"
          else
            item[:bioc_version_str] = nil
          end
          rep = Nanoc3::ItemRep.new(item, :unique_name)

          for sym in link_list
            new_sym = "#{sym.to_s}_repo".to_sym
            item[new_sym] = []
            for x in to_array(item[sym])
              item[new_sym].push(is_bioc_package?(x, version))
            end
          end


          items.push item
        end
      end
    end
    
    
    
    items
  end
  
  # if "packageName" is a bioC package, return its repo directory, otherwise return false
  def is_bioc_package?(packageName, version)
    #puts "packageName = #{packageName}" 
    return false if packageName.nil?
    packageName = packageName.split(/[ |(]/).first
    @all_packages[version].each_pair do |k,v|
      return k if v.has_key?(packageName)
    end
    false
  end
  
  
end
