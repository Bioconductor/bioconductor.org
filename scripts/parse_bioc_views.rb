#!/usr/bin/env ruby

require 'rubygems'
require 'json'
require 'pp'



class ParseBiocViews
  def initialize(args)
    inputfile, outputfile = args
    f = File.open(inputfile)
    
    outfile = File.open(outputfile, "w")
    
    lines = f.readlines
    json = lines.join("\n")
    obj = JSON.parse(json)
    
    clean_hash = {}
    obj.each_pair do |k,v|
      clean_hash[k] = clean(v)
    end
 
    keys_to_delete = []
 
    clean_hash.each_pair do |k,v|
      
      
      
      if (v['subViews'] != "")
        v['children'] = []
        
        for subview in v['subViews']
          v['children'].push clean_hash[subview]
        end
        v['children'].sort!{|a,b|a['data'].downcase<=>b['data'].downcase}
      end
      v.delete("subViews")
      
      unless k == "BiocViews"
        if (!v.has_key?("children") or v['children'].empty? or v['children'].nil?)
          keys_to_delete.push k
        else
          keep = []
          for child in v['children']
            keep.push child if (child['data'] =~ /\(/)
          end
          v['children'] = keep
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
    
    outfile.puts root.to_json
    
    
  end
end




def clean(arg)
  item = arg.first
  item['data'] = item['name']
  item.delete("name")
  item.delete("parentViews")
  
  
  if (item.has_key?("packageList") && item["packageList"] != "" && !item["packageList"].empty?)
    item["attr"] = {"packageList" => item["packageList"].keys.sort{|a,b|a.downcase<=>b.downcase}.join(","), "id" => item['data']}
    item['data'] += " (#{item["packageList"].length})"
  end
  item.delete("packageList")

  item
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