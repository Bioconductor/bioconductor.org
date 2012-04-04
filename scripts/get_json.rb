#!/usr/bin/env ruby

require 'rubygems'
require 'rgl/adjacency'
require 'sqlite3'
require 'json'
require 'dcf'
require 'net/http'
require 'uri'
require 'pp'

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
    ret = []
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
      ret.push obj
    end
    ret
  end
  
  def initialize(repo, version, outdir)
    dbfile = `R --vanilla --slave -e "cat(system.file('extdata','biocViewsVocab.sqlite',package='biocViews'))"`
    db = SQLite3::Database.new(dbfile)
    rows = db.execute("select * from biocViews")
    @g = RGL::DirectedAdjacencyGraph.new()
    for row in rows
      @g.add_edge(row.first, row.last)
    end
    dcfs = get_dcfs(repo, version)
    #pp dcfs
    json = JSON.pretty_generate(dcfs.first)
    pkgs_file = File.open("#{outdir}/packages.json", "w")
    pkgs_file.print(json)
    pkgs_file.close
  end
end


j = GetJson.new("bioc", "2.11", "/Users/dtenenba/tmp")

