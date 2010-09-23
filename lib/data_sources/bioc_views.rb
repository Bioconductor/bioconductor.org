class BiocViews < Nanoc3::DataSource

  require 'pp'
  require 'rubygems'
  require 'json'
  
  identifier :bioc_views
  
  # todo - write an index page for each repo containing links to all packages
  
  # todo - find out if there is a way to skip items() altogether if things are not found in up()
  def up
    
    @repos = {"bioc/" => "Software", "data/annotation/" => "AnnotationData", "data/experiment/" => "ExperimentData"}
    
    @good_to_go = true
    
    @repos.each_pair do |k,v|
      dir = "#{config[:json_dir]}/#{k}"
      @good_to_go = false unless test(?f, "#{dir}/packages.json") 
    end
    
    unless @good_to_go
      puts "BiocViews data source: json files not present, skipping initialization"
    end
    
      
  end

  
  def items
    
    unless @good_to_go
      puts "BiocViews_DataSource: no JSON file(s) found. Package detail pages will not be built"
      return []
    end

    items = []
    
    @repos.each_pair do |k,v|
      dir = "#{config[:json_dir]}/#{k}"
      json_file = File.open("#{dir}/packages.json")
      
      packages = JSON.parse(json_file.readlines.join("\n"))
      for package in packages.keys
        repo = v
        packages[package][:repo] = repo
        packages[package][:bioc_version_str] = config[:version_str]
        packages[package][:bioc_version_num] = config[:version_num]
        packages[package][:title] = packages[package][:Package]
        items.push Nanoc3::Item.new(nil, packages[package], package)
      end
    end
    items
  end
  
  
end
