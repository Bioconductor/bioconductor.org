#!/usr/bin/env ruby

require "pp"
require "rubygems"
require "dcf"
require "rexml/document"

include REXML

pkg_limit = 100

rpacks_url = 'https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/'

list =  `svn --non-interactive --no-auth-cache --username readonly --password readonly list #{rpacks_url}`

lines = list.split("\n")
mfs = lines.find_all{|i| i =~ /^bioc/ && i =~ /\.manifest$/}



sorted_manifests = mfs.sort do |a, b|
  a = a.gsub(/bioc_/, "").gsub(/\.manifest/, "")
  b = b.gsub(/bioc_/, "").gsub(/\.manifest/, "")
  major_a, minor_a = a.split(".")
  major_b, minor_b = b.split(".")
  if major_a == major_b
    Integer(minor_a) <=> Integer(minor_b)
  else
    Integer(major_a) <=> Integer(major_b)
  end
    
end

mf = sorted_manifests.last

manifest = `curl -s -u readonly:readonly #{rpacks_url}#{mf}`

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
  raw_result = `curl -s -I http://bioconductor.org/packages/devel/bioc/html/#{pkg}.html`
  first_line  = raw_result.split("\n").first
  if (first_line =~ /200/) # is there a pkg homepage?
    description = `curl -s -u readonly:readonly #{rpacks_url}#{pkg}/DESCRIPTION`
#    3.times {|i| description.gsub!(/\r\n/, "")}
    description.gsub!(/\r\n?/, "\n")
    parse_result = Dcf.parse(description)
    dcs = parse_result.first
    descs.push dcs
    
    break if descs.size == pkg_limit
  end
end



descs.reverse!

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

doc = Document.new startxml

elem = nil
doc.elements.each("rss/channel") {|i| elem = i}



for desc in descs
  item = Element.new "item"
  title = Element.new "title"
  title.text = "#{desc["Package"]} - #{desc["Title"]}"
  pubdate = Element.new "pubdate"
  pubdate.text = "2011-01-01"
  author = Element.new "author"
  author.text = desc["Author"]
  description  = Element.new "description"
  description.attributes["disable-output-escaping"] = "yes"
  
  # attr?
  
  link = %Q(\n<br/><a href="http://bioconductor.org/packages/devel/bioc/html/#{desc["Package"]}.html">link</a>)
  
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