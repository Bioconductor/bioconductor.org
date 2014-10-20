#!/usr/bin/env ruby

# coding: utf-8
require 'htmlentities'
#require 'unicode'
#$KCODE = 'UTF-8'

require 'rubygems'
require 'tmpdir'
require 'rgl/adjacency'
require 'rgl/traversal'
require 'sqlite3'
require 'json'
require 'net/http'
require 'uri'
require 'pp'
require 'yaml'

class String
  def to_boolean()
    return true if self.downcase == "true"
    return false if self.downcase =="false"
    return self # maybe return nil?
  end
end

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
    ret.each_pair do |k,v|
      ret[k] = v.strip
      if v == "TRUE" or v == "FALSE"
        ret[k] = v.to_boolean
      end
    end
    ret
  end
  
  def self.get_value_as_array(value, remove_version_specifiers=false)
    return [] if value.empty? or value.nil? 
    value.gsub!(",,", "%%ESCAPED_COMMA%%")
    value.gsub!(", ", ",")
    segs = value.split(",")
    segs = segs.map{|i| i.gsub("%%ESCAPED_COMMA%%", ",")}
    segs = segs.map{|i| i.gsub(/^\s*|\s*$/, "")}
    if remove_version_specifiers
      segs = segs.map{|i| i.split(" ").first}
    end
    segs
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
    url = URI.parse("http://master.bioconductor.org/packages/#{version}/#{repo}/VIEWS")
    req = Net::HTTP::Get.new(url.path)
    res = Net::HTTP.start(url.host, url.port) {|http|
      http.request(req)
    }
    views = res.body
    view_dcfs = []
    view_lines = views.split("\n")
    dcf = ""
    view_lines.each_with_index do |line, i| 
      line = line.force_encoding("UTF-8") if line.respond_to? :force_encoding
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
    plural_fields = ["Depends", "Suggests", "Imports", "Enhances", "biocViews",
      "vignettes", "vignetteTitles", "Rfiles", "htmlDocs", "htmlTitles"]
    for dcf in dcfs
      for key in dcf.keys
        if plural_fields.include? key
          dcf[key] = Dcf.get_value_as_array(dcf[key])
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
    rows = nil
    Dir.mktmpdir do |dir|
      if version == @config["devel_version"]
        branch = "trunk"
      else
        branch = "branches/RELEASE_#{version.gsub(".", "_")}"
      end
      url = "https://hedgehog.fhcrc.org/bioconductor/#{branch}/madman/Rpacks/biocViews/inst/extdata/biocViewsVocab.sqlite"
      auth = {:username => "readonly", :password => "readonly"}
      File.open("#{dir}/biocViewsVocab.sqlite", "wb") do |f|
          resp = HTTParty.get(url, :verify => false, :basic_auth => auth)
          # require 'pry';binding.pry
          f.write resp
      end
      #`curl -s -u  readonly:readonly #{url} > #{dir}/biocViewsVocab.sqlite`
      dbfile = "#{dir}/biocViewsVocab.sqlite"      
      db = SQLite3::Database.new(dbfile)
      rows = db.execute("select * from biocViews")
    end
    #dbfile = `R --vanilla --slave -e "cat(system.file('extdata','biocViewsVocab.sqlite',package='biocViews'))"`
    #db = SQLite3::Database.new(dbfile)
    #rows = db.execute("select * from biocViews")
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
          bioc_views = [default_view]
        else
          bioc_views = value["biocViews"] << default_view
        end
      else
        bioc_views = [default_view]
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

    nodes_to_delete = []
    for node in g.vertices
      delete_parent = true
      g.depth_first_visit(node) do |n|
        if node_attrs.has_key? n
          delete_parent = false
          break
        end
      end
      nodes_to_delete.push node if delete_parent
      
    end
    
    
    # todo - better way to delete nodes-- build up a nodeList from the 
    # biocViews of packages, uniq it, then delete all nodes not in that list.
    for node_to_delete in nodes_to_delete.uniq
      # disabling deletion, because we want to include orphan nodes
      # let javascript handle filtering them out
  #####    g.remove_vertex node_to_delete
    end
    
    
    
    for node in g.vertices
      pkgs = []
      g.depth_first_visit(node) do |n|
        pkgs += node_attrs[n]['packageList'] if node_attrs.has_key? n and \
          node_attrs[n].has_key? 'packageList'
      end
      
      node_attrs[node] = {} unless node_attrs.has_key? node
      node_attrs[node]['packageList'] = \
        [] unless node_attrs[node].has_key? "packageList"
      if node_attrs[node]['packageList'].empty?
        node_attrs[node]['packageList'] = pkgs
      end
      node_attrs[node]['childnum'] = pkgs.uniq.length
    end
    
    
    for node in g.vertices
      h = {}
      h["name"] = node
      h["parentViews"] = g.parents(node)
      h["subViews"] = g.children(node)
      if (node_attrs.has_key? node)
        h["packageList"] = node_attrs[node]['packageList']
        h["childnum"] = node_attrs[node]['childnum']
      else
        h['packageList'] = []
      end
    
      ret[node] = h
    end
    ret
  end
  
  
  def initialize(version, outdir)
    @coder = HTMLEntities.new
    nanoc_dir = File.expand_path $0
    segs = nanoc_dir.split "/"
    2.times {segs.pop}
    nanoc_dir = segs.join("/")
    config = YAML.load_file("./config.yaml")
    @config = config
    repos = ["bioc", "data/experiment", "data/annotation"]
    repos = config["devel_repos"] if version == config["devel_version"]
    packages_data = []
    biocviews_data = []
    for repo in repos
      p, b = handle_repo(repo, version, outdir)
      packages_data.push p
      biocviews_data.push b
    end
    
    packages_data = get_reverse_dependencies(packages_data)
    
    repos.each_with_index do |repo, i|
      write_packages_file(packages_data[i], version, repo, outdir)
    end
    
    write_tree_file(repos, biocviews_data, version, outdir)
    
  end
  
  def write_tree_file(repos, data, version, outdir)
    args = [repos, data, "#{outdir}/#{version}/tree.json"]
    p = ParseBiocViews.new(args)
  end
  
  
  def write_packages_file(data, version, repo, outdir)
    fulloutdir = "#{outdir}/#{version}/#{repo}"
    json = JSON.pretty_generate(data)
    # out = []
    # for i in json.split ""
    #   # if i.ord > 128 # ??
    #   #   out << @coder.encode(i, :named)
    #   #   # puts "converting #{i} to #{out.last}"
    #   # else
    #     # out << i
    #   # end
    # end
    # json = out.join
    pkgs_file = File.open("#{fulloutdir}/packages.json", "w")
    pkgs_file.print(json)
    pkgs_file.close
  end
  
  
  # the biocViews packages builds reverse dependencies, but
  # not across repositories, so we handle that here
  def get_reverse_dependencies(data)
    fields = ["Suggests", "Depends", "Imports"]
    reverse_deps = {"Suggests" => {}, "Depends" => {}, "Imports" => {}}
    for datum in data
      datum.each_pair do |k,v|
        for field in fields
          next unless v.has_key? field and !v[field].empty?
          items = v[field].map{|i| i.split(/ |\(/).first}
          for item in items
            reverse_deps[field][item] = [] unless reverse_deps[field].has_key? item
            reverse_deps[field][item].push k
          end
        end
      end
    end
    
    
    fields = {"Suggests" => "suggestsMe", "Depends" => "dependsOnMe", "Imports" => "importsMe"}
    ret = []
    for datum in data
      datum.each_pair do |k, v|
        fields.each_pair do |fk, fv|
          if reverse_deps[fk].has_key? k
            rdeps = reverse_deps[fk][k]
            rdeps.sort! do |a, b|
              a.downcase <=> b.downcase
            end
            datum[k][fv] = rdeps
          end
        end
      end
      ret.push datum
    end
    ret
  end
  
  
  def handle_repo(repo, version, outdir)
    fulloutdir = "#{outdir}/#{version}/#{repo}"
    FileUtils.mkdir_p fulloutdir unless Kernel.test(?d, fulloutdir)
    dcfs = get_dcfs(repo, version)
    bioc_views = get_bioc_views(dcfs, repo, version)
    
    ## todo - remove this
    json = JSON.pretty_generate(bioc_views)
    views_file = File.open("#{fulloutdir}/biocViews.json", "w")
    views_file.print(json)
    views_file.close
    ## end of todo: remove this
    [dcfs, bioc_views]
  end
end



