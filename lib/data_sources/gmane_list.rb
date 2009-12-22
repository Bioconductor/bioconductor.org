# encoding: utf-8

require 'rubygems'
require 'httparty'
require 'time'
require 'yaml'

class GmaneList < Nanoc3::DataSource
  identifier :gmane_list
  def fetch
    bioc_list_url = self.config[:gmane_rss_url]
    data = HTTParty.get(bioc_list_url)
    data["rdf:RDF"]["item"].map do |item|
      attributes = {
        :title => item["title"],
        :date => fix_date(item),
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

  private

  def fix_date(rss_item)
    Time.parse(rss_item["dc:date"] + " GMT").utc
  rescue
    Time.new.utc
  end
end

