class BiocViews < Nanoc3::DataSource

  require 'pp'
  require 'rubygems'
  require 'json'
  
  identifier :bioc_views
  
  # todo - write an index page for each repo containing links to all packages
  
  
  
  # todo - find out if there is a way to skip items() altogether if things are not found in up()
  def up
    @bad_packages = ["snpMatrix2"] # don't process these

    @repos = {"bioc/" => "Software", "data/annotation/" => "AnnotationData", "data/experiment/" => "ExperimentData"}
    
    @good_to_go = true
    
    @site_config = YAML.load_file("./config.yaml")
    
    
    for version in @site_config["versions"]
      @repos.each_pair do |k,v|
        
        
        dir = "#{config[:json_dir]}/#{version}/#{k}"
        @good_to_go = false unless test(?f, "#{dir}/packages.json") 

        #todo remove
        @good_to_go = false unless test(?f, "#{dir}/vignette_titles.json") 


      end
    end
    
    
    
    
    
    @all_packages = {}
    
    if @good_to_go
      for version in @site_config["versions"]
        hsh = {}
        @repos.each_pair do |k,v|
          
          
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

  
  #todo remove
  def do_bad_stuff(obj, dir)
    
    tf = File.open("#{dir}/vignette_titles.json")
    json = tf.readlines.join("\n")
    vt = JSON.parse(json)
    for bad in @bad_packages
      vt.delete(bad)
    end
    
    obj.each_pair do |k,v|

      if (vt.has_key? k)
        
          
        if (vt[k].respond_to? :has_key?)
          
          
          v["vignetteTitles"] = vt[k]["vignetteTitles"] unless vt[k]["vignetteTitles"].nil? or vt[k]["vignetteTitles"].empty?
          v["vignetteScripts"] = vt[k]["vignetteScripts"] unless vt[k]["vignetteScripts"].nil? or vt[k]["vignetteScripts"].empty?
          v["vignetteFiles"] = vt[k]["vignetteFiles"] unless vt[k]["vignetteFiles"].nil? or vt[k]["vignetteFiles"].empty?

          
          v[:vignetteTitles] = v["vignetteTitles"]
          v[:vignetteScripts] = v["vignetteScripts"]
          v[:vignetteFiles] = v["vignetteFiles"]
          
        end
        
        
      end
      
    end
    obj
  end
  # end remove
  
  
  def get_index_page(packages, repo, version)
    item = Nanoc3::Item.new(nil, {}, "all-#{repo}-#{version}")
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
      @repos.each_pair do |k,v|
        dir = "#{config[:json_dir]}/#{version}/#{k}"
        
        
        json_file = File.open("#{dir}/packages.json")



        packages = JSON.parse(json_file.readlines.join("\n"))
        for bad in @bad_packages
          packages.delete(bad)
        end

        items.push(get_index_page(packages, v, version))


        #todo remove
        packages = do_bad_stuff(packages,dir)
        # end remove




        for package in packages.keys
          repo = k
          item = Nanoc3::Item.new(nil, packages[package], package)
          
          item[:rebase] = true
          item[:subnav] = []
          item[:subnav].push({:include => "/_workflows/"})
          item[:subnav].push({:include => "/_mailing_list/"})

          item[:title] = "#{item[:Package]}"
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
            for x in item[sym]
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
    @all_packages[version].each_pair do |k,v|
      return k if v.has_key?(packageName)
    end
    false
  end
  
  
end
