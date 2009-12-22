# encoding: utf-8

require 'rubygems'
require 'httparty'
require 'time'
require 'yaml'

class GmaneList < Nanoc3::DataSource
  identifier :gmane_list
  def fetch
    # FIXME: make this URL part of site config
    bioc_list_url = ('http://rss.gmane.org/topics/excerpts/' +
                     'gmane.science.biology.informatics.conductor')
    data = HTTParty.get(bioc_list_url)
    data["rdf:RDF"]["item"].map do |item|
      attributes = {
        :title => item["title"],
        :date => (Time.parse(item["dc:date"]) rescue Time.now),
        :link => item["link"],
        :author => item["dc:creator"]
      }
      content = item["description"]
      mtime = nil
      identifier = "/#{File.basename(item["link"])}/"
      Nanoc3::Item.new(content, attributes, identifier, mtime)
    end.sort { |a, b| b[:date] <=> a[:date] }
  end
  
  def items
    @items ||= fetch()
  end
end

