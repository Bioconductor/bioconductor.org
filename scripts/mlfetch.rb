#!/usr/bin/env ruby

require 'rubygems'
require 'httparty'
require 'yaml'

BIOC_LIST_URL = ('http://rss.gmane.org/topics/excerpts/' +
                 'gmane.science.biology.informatics.conductor')

def fetch
  data = HTTParty.get(BIOC_LIST_URL)
  data["rdf:RDF"]["item"].map do |item|
    {
      :title => item["title"],
      :date => item["dc:date"],
      :link => item["link"],
      :author => item["dc:creator"]
    }
  end
end

zz = fetch()
zz[0..4].each do |m|
  puts m[:title]
end

xx = {
  'title' => 'Bioconductor Mailing Recent Activity',
  'bioc_list_main_url' => 'https://stat.ethz.ch/mailman/listinfo/bioconductor',
  'bioc_list_search' => 'http://dir.gmane.org/gmane.science.biology.informatics.conductor',
  'posts' => zz[0..5]
}

puts xx.to_yaml

