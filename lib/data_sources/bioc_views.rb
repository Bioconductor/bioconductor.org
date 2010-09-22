class BiocViews < Nanoc3::DataSource

  require 'pp'
  require 'rubygems'
  require 'json'
  
  identifier :bioc_views
  
  def up
#    puts "hi from bv.up"
#    puts "config:"
#    pp config
  end
  
  def items
    if (test(?f, config[:json_file]))
      json_file = File.open(config[:json_file])
      packages = JSON.parse(json_file.readlines.join("\n"))
      items = []
      for package in packages.keys
        packages[package][:bioc_version_str] = config[:version_str]
        packages[package][:bioc_version_num] = config[:version_num]
        packages[package][:title] = packages[package][:Package]
        items.push Nanoc3::Item.new(nil, packages[package], package)
      end
      items
    else
      puts "BiocViews_DataSource: no JSON file found. Package detail pages will not be built"
      []
    end
  end
  
  
  
  def layouts
    layouts = []
    f = File.open("layouts/_bioc_views_package_detail.html")
    content = f.readlines.join("\n")


    layouts.push Nanoc3::Layout.new(content, {}, "/_bioc_views_package_detail/")
    #layouts.push Nanoc3::Layout.new(nil, nil, "/default/")
    layouts
  end
  
#  def code_snippets
#  end
  
#  def down
#    puts "hi from bv.down"
#    puts "layouts:"
#    i  = @site.layouts.first.inspect
#    puts i
#    exit if true
#  end
end
