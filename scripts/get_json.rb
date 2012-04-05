#!/usr/bin/env ruby

require 'rubygems'
require 'rgl/adjacency'
require 'sqlite3'
require 'json'
## this DCF parser is slow. Here's a faster one:
## http://gist.github.com/117293
## but it does not like unescaped colons in values,
## and also requires spaces after field names (and colons, presumably)
## so we write our own below, based on python code from BBS.
#require 'dcf'
require 'net/http'
require 'uri'
require 'pp'


module Dcf
  def self.parse(input)
    ret = {}
    lines = input.split("\n")
    key = nil
    for line in lines
      if !(line =~ /^\s+/)
        key, val = line.split(":", 2)
        ret[key] = val.strip
      else
        ret[key] += " " + line.strip
      end
    end
    ret
  end
end

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
      line = line.force_encoding("UTF-8")
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
      for key in dcf.keys
        if plural_fields.include? key
          ary = dcf[key].split(",")
          ary = ary.map{|i| i.gsub(/^\s*|\s$/, "")}
          dcf[key] = ary
        end
      end
      key = dcf["Package"]
      ret[key] = dcf
    end
    ret
  end
  
  
  def get_bioc_views(dcfs, repo, version)
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
    node_attrs = {}
    dcfs.each_pair do |key, value|
      if value.has_key? "biocViews"
        if value["biocViews"].nil? or value["biocViews"].empty?
          bioc_views = []
        else
          bioc_views = value["biocViews"] << default_view
        end
      else
        bioc_views = []
      end
      bioc_views.uniq!
      bioc_views.sort!
      value['biocViews'] = bioc_views
      for bioc_view in bioc_views
        if g.has_vertex? bioc_view
          node_attrs[bioc_view] = {} unless node_attrs.has_key? bioc_view
          node_attrs[bioc_view]['packageList'] = [] unless node_attrs[bioc_view].has_key? 'packageList'
          node_attrs[bioc_view]['packageList'].push key
        else
          puts "[Warning] #{repo} package #{key} has invalid view #{bioc_view}. (BioC #{version})"
        end
      end
    end
    ret = {}
    for node in g.vertices
      h = {}
      h["name"] = node
      h["parentViews"] = g.parents(node)
      h["subViews"] = g.children(node)
      if (node_attrs.has_key? node)
        h["packageList"] = node_attrs[node]['packageList']
      else
        h['packageList'] = []
      end
    
      ret[node] = h
    end
    ret
  end
  
  def initialize(repo, version, outdir)
    FileUtils.mkdir_p outdir unless Kernel.test(?d, outdir)
    dcfs = get_dcfs(repo, version)
    json = JSON.pretty_generate(dcfs)
    pkgs_file = File.open("#{outdir}/packages.json", "w")
    pkgs_file.print(json)
    pkgs_file.close
    bioc_views = get_bioc_views(dcfs, repo, version)
    json = JSON.pretty_generate(bioc_views)
    views_file = File.open("#{outdir}/biocViews.json", "w")
    views_file.print(json)
    views_file.close
  end
end



