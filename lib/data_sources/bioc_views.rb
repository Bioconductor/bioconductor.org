class BiocViews < Nanoc3::DataSource

  require 'pp'
  require 'rubygems'
  require 'json'
  
  identifier :bioc_views
  
  # todo - find out if there is a way to skip items() altogether if things are not found in up()
  def up
    bioc_views_file = config[:json_file].gsub("packages.json", "biocViews.json")
    if (test(?f, config[:json_file]) and test(?f, bioc_views_file))
      f = File.open(bioc_views_file)
      json = f.readlines.join("\n")
      @bioc_views = JSON.parse(json)
      
      @immediate_parents = {}
      
      @bioc_views.each_pair do |k,v|
        hsh = v.first
        next unless hsh['packageList'].respond_to? :keys
        hsh['packageList'].each_pair do |pk, pv|
          @immediate_parents[pk] = hsh["parentViews"]
        end
      end
      
    else
      puts "BiocViews data source: json files not present, skipping initialization"
    end
  end

  
  def items
    
    bioc_views_file = config[:json_file].gsub("packages.json", "biocViews.json")
    
    if (test(?f, config[:json_file]) and test(?f, bioc_views_file))
      json_file = File.open(config[:json_file])
      packages = JSON.parse(json_file.readlines.join("\n"))
      items = []
      for package in packages.keys
        
        repo = find_repo(package)
        packages[package][:repo] = repo
        packages[package][:bioc_version_str] = config[:version_str]
        packages[package][:bioc_version_num] = config[:version_num]
        packages[package][:title] = packages[package][:Package]
        items.push Nanoc3::Item.new(nil, packages[package], package)
      end
      items
    else
      puts "BiocViews_DataSource: no JSON file(s) found. Package detail pages will not be built"
      []
    end
  end
  
  
  private
  
  def find_repo(pkg_name)
    ancestor = @immediate_parents[pkg_name]
    return "bioc" if ancestor == "BiocViews"
    loop do
      ary = @bioc_views[ancestor]
      break if ary.nil?
      hash = ary.first
      parent = hash["parentViews"]
      if (parent == "BiocViews")
        break
      end
      ancestor = @bioc_views[parent].first["name"]
    end
    choices = {"Software" => "bioc", "AnnotationData" => "data/annotation", "ExperimentData" => "data/experiment",
      "BiocViews" => "bioc"}
    #puts "ancestor = #{ancestor}"
    choices[ancestor]
  end
  
  
end
