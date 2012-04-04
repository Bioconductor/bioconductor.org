#!/usr/bin/env ruby

require 'rubygems'
require 'rgl/adjacency'
require 'sqlite3'
require 'json'
require 'dcf'
require 'net/http'
require 'uri'
require 'pp'

class RGL::DirectedAdjacencyGraph
  def children(parent)
    children = []
    for edge in self.edges
      children.push edge.target if edge.source == parent
    end
    children
  end
  
  def parents(child)
    parents = []
    for edge in self.edges
      parents.push edge.source if edge.target == child
    end
    parents
  end
end


class GetJson
  
  
  def get_dcfs(repo, version)
    url = URI.parse("http://bioconductor.org/packages/#{version}/#{repo}/VIEWS")
    req = Net::HTTP::Get.new(url.path)
        res = Net::HTTP.start(url.host, url.port) {|http|
          http.request(req)
        }
    views = res.body
    view_dcfs = []
    view_lines = views.split("\n")
    dcf = ""
    view_lines.each_with_index do |line, i| 
      if i == (view_lines.length() -1)
        if !line.empty?
          dcf = dcf + "\n" + line
          line = ""
        end
      end
      if line.empty?
        pdcf = Dcf.parse(dcf)
        view_dcfs.push pdcf
        dcf = ""
      else
        if dcf == ""
          dcf = line
        else
          dcf = dcf + "\n" + line
        end
      end
    end
    return clean_dcfs(view_dcfs)
  end
  
  
  def clean_dcfs(dcfs)
    ret = {}
    plural_fields = ["Depends", "Suggests", "Imports", "Enhances", "biocViews"]
    for dcf in dcfs
      obj = dcf.first
      for key in obj.keys
        if plural_fields.include? key
          ary = obj[key].split(",")
          ary = ary.map{|i| i.gsub(/^\s*|\s$/, "")}
          obj[key] = ary
        end
      end
      key = obj["Package"]
      ret[key] = obj
    end
    ret
  end
  
  
  
  
  
  def get_bioc_views(dcfs, repo)
    default_view = "Software" if repo == "bioc"
    default_view = "AnnotationData" if repo == "data/annotation"
    default_view = "ExperimentData" if repo == "data/experiment"
    
    dbfile = `R --vanilla --slave -e "cat(system.file('extdata','biocViewsVocab.sqlite',package='biocViews'))"`
    db = SQLite3::Database.new(dbfile)
    rows = db.execute("select * from biocViews")
    g = RGL::DirectedAdjacencyGraph.new()
    sort_order = []
    for row in rows
      g.add_edge(row.first, row.last)
      sort_order.push row.first
    end
    dcfs.each_pair do |key, value|
      if value.has_key? "biocViews"
        for view in value["biocViews"]
          if g.has_vertex? (view)
            g.add_edge(view, )
          else
            print "[WARNING] non-existent biocView #{view} in package #{key}."
          end
        end
      else
        print "[WARNING] Package #{key} has no biocViews."
        # deal with default_view
      end
    end
  end
  
  def initialize(repo, version, outdir)
    dcfs = get_dcfs(repo, version)
    #pp dcfs
    json = JSON.pretty_generate(dcfs)
    pkgs_file = File.open("#{outdir}/packages.json", "w")
    pkgs_file.print(json)
    pkgs_file.close
    bioc_views = get_bioc_views(dcfs, repo)
  end
end


j = GetJson.new("bioc", "2.11", "/Users/dtenenba/tmp")

