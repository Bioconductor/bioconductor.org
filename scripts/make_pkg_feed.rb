#!/usr/bin/env ruby

require "pp"
require "rubygems"
require "dcf"
require "rexml/document"
require "yaml"
require 'nokogiri'
require 'httparty'
require 'net/http'

include REXML


nanoc_path = File.expand_path($0)
segs = nanoc_path.split("/")
segs.pop
segs.pop
nanoc_path = segs.join("/")

pkg_limit = 100

rpacks_url = 'https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/'



site_config = YAML.load_file("#{nanoc_path}/config.yaml")

devel_version = site_config["devel_version"]


manifest = HTTParty.get("#{rpacks_url}bioc_#{devel_version}.manifest",
  :verify => false, :basic_auth => {:username => 'readonly', :password => 'readonly'}).body


lines = manifest.split("\n")

pkgs = []
for line in lines
  next if line =~ /^#/
  next if line.empty?
  pkgs.push line.gsub(/^Package: /, "")
end

pkgs.reverse!

descs = []


for pkg in pkgs
  # HTTParty.head doesn't work for some reason
  url = URI("http://bioconductor.org/packages/devel/bioc/html/#{pkg}.html")
  resp = nil
  Net::HTTP.start(url.host, url.port){|http|
   response = http.head('/')
   resp = response
  }

  if resp.code == "200" # is there a pkg homepage?
    response = HTTParty.get("#{rpacks_url}#{pkg}/DESCRIPTION",
        :verify => false, :basic_auth => {:username => 'readonly', :password => 'readonly'})
    next unless response.code == 200
    description = response.body
#    3.times {|i| description.gsub!(/\r\n/, "")}
    description.gsub!(/\r\n?/, "\n")
    description.gsub!(/\n{2,}/, "\n")
    description = description.gsub(/^\n+/, "")
    parse_result = Dcf.parse(description)
    dcs = parse_result.first
    descs.push dcs
    
    break if descs.size == pkg_limit
  end
end



#descs.reverse!

startxml = <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<rss version="2.0">
  <channel>
    <title>New Packages: Bioconductor</title>
    <link>http://www.bioconductor.org/</link>
    <description>New Bioconductor Packages (devel branch)</description>
    <language>en-us</language>
    <copyright></copyright>
    <generator>ruby</generator>
  </channel>
</rss>
EOF

doc = REXML::Document.new startxml

elem = nil
doc.elements.each("rss/channel") {|i| elem = i}



for desc in descs
  url = "http://bioconductor.org/packages/devel/bioc/html/#{desc["Package"]}.html"

  item = Element.new "item"
  title = Element.new "title"
  title.text = "#{url} #{desc["Package"]} #{desc["Title"]}"
  pubdate = Element.new "pubdate"
  xml = `svn log --xml --username readonly --password readonly --no-auth-cache --non-interactive -v  --limit 1 -r HEAD:1 https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/#{desc['Package']}`

  xdoc = Document.new xml
  pubdate.text = xdoc.elements["log/logentry/date"].text

  #puts "Package: #{desc['Package']}, pubdate.text: #{pubdate.text}"
  #pubdate.text =  "Sat, 1 Jan 2012 00:00:00 GMT"  #"2011-01-01"
  author = Element.new "author"
  author.text = desc["Author"]
  description  = Element.new "description"
  description.attributes["disable-output-escaping"] = "yes"
  
  # attr?
  
  link = %Q(\n<br/><a href="#{url}">link</a>)
  
  description.text = "#{desc["Description"]}#{link}"
  
  
  
  item.add_element title
  item.add_element pubdate
  item.add_element author
  item.add_element description
  elem.add_element item
end

# for some reason, indentation breaks on merlot2
#doc.write($stdout, 2)
doc.write($stdout)

