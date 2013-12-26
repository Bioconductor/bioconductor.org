#!/usr/bin/env ruby

require 'rubygems'
require 'json'
require 'pp'

BAD_PACKAGES = ["snpMatrix2"]


class ParseBiocViews
  
  def get_subtree(obj, repo)
    clean_hash = {}
    obj.each_pair do |k,v|
      clean_hash[k] = clean(v)
    end

 
    keys_to_delete = []
 
    clean_hash.each_pair do |k,v|
      
      
      
      if (v['subViews'] != "")
        v['children'] = []

        #puts v['subViews'] if v['subViews'] == "qPCR"
        if (v['subViews'].respond_to? :each)
          for subview in v['subViews']
            v['children'].push clean_hash[subview]
          end
        end
        v['children'].sort!{|a,b|a['data'].downcase<=>b['data'].downcase}
      end
      v.delete("subViews")
      
      unless k == "BiocViews"
        if (!v.has_key?("children") or v['children'].empty? or v['children'].nil?)
          keys_to_delete.push k
        end
        unless (v.has_key?("attr") )
          keys_to_delete.push
        end
        if (v.keys.length == 1)
          keys_to_delete.push k
        end
      end
      
    end
 
    for key_to_delete in keys_to_delete
      clean_hash.delete(key_to_delete)
    end

    
    root = clean_hash['BiocViews']['children']

    repos = {"bioc" => "Software", "data/annotation" => "AnnotationData", "data/experiment" => "ExperimentData"}
    
    key = root.find{|i|i['data'] =~ /^#{repos[repo]}/}
    
    
    key
 
  end
  
  
  
  def initialize(args)
    reposlist, inputdata, outputfile = args
    #inputfiles, outputfile = args
    
    
    obj = {}
    
    
    top_level_tree = []

    inputdata.each_with_index do |this_one, i|
      repo = reposlist[i]
      
      subtree = get_subtree(this_one, repo)
      top_level_tree.push subtree
    end
    

    
    top_level_tree.sort! do |a, b|
      if a['data'] =~ /^Software/
        -1
      elsif b['data'] =~ /^Software/
        1
      else
        a['data'] <=> b['data']
      end
    end
    
    outfile = File.open(outputfile, "w")
    
    outfile.puts JSON.pretty_generate(top_level_tree)
    
    
  end
end




def clean(arg)
  item = arg
#  item = arg.first
  item['data'] = item['name']
  item.delete("name")
  item.delete("parentViews")
  
  if (item.has_key? "childnum" and (!item['childnum'].nil?) and item['childnum'] > 0)
    item["data"] += " (#{item['childnum']})"
    item.delete('childnum')
  end


  
  if (item.has_key?("packageList") && !item['packageList'].nil? &&item["packageList"] != "" && !item["packageList"].empty?)
    for bad in BAD_PACKAGES
      item['packageList'].delete(bad)
    end
    
    item["attr"] = {"packageList" => item["packageList"].sort{|a,b|a.downcase<=>b.downcase}.join(","),
      "id" => item['data'].split(" ").first}
    
  end
  item.delete("packageList")

  item
end

if __FILE__ == $0
#  puts "We are running directly as a script."
else
#  puts "We are not running directly as a script."
end


__END__
unless ARGV.size == 2
  puts "usage: #{$0} inputfile outputfile"
  exit
end

pbv = ParseBiocViews.new(args)



This is the target data format we want:

[
   {
      "data":"dan's root node",
      "children":[
         {
            "data":"First child"
         },
         {
            "data":"Second child"
         }
      ]
   },
   {
      "data": "Ajax node"
   }
]