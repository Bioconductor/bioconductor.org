#!/usr/bin/env ruby


require "rubygems"
require "hpricot"
require "pp"



if (ARGV.size != 1)
  puts "supply the root directory"
  exit 1
end


@rootdir = ARGV.first.dup
@rootdir.gsub!(/\/$/, "")
startfile = "index.html"

@link_map = {}

@fnf = {}

def get_links(filename)
  #puts "in get_links with filename #{filename}"
  filename.sub!(/^\//, "")
  begin
    fh = File.new("#{@rootdir}/#{filename}")
  rescue Exception => ex
    @fnf[filename] = 1
    return
  end
  filename = "/#{filename}" unless filename =~ /^\//
  @link_map[filename] = 1
  lines = fh.readlines
  html = lines.join
  links  = []
  doc = Hpricot(html)
  doc.search("a").each {|i| links.push i}
  for link in links
    next if link.nil? or link.empty?
    #next unless link.attributes.has_key? "href"
    href = link.attributes["href"]
    next if href.nil? or href.empty?
    next if href =~ /:\/\//
    next if href =~ /^#/
    if href =~ /#/
      href = href.split("#").first
    end
    if href =~ /\/$/
      href = href + "index.html"
    end
    href = "/#{href}" unless href =~ /^\//
#    puts "href = #{href}, text = #{link.inner_text}"
    if (href =~ /\.html$/)
      get_links(href) unless @link_map.has_key? href
    elsif (href.downcase() =~ /\.r$|\.pdf$/)
      @link_map[href] = 1
    end
  end
  
end


get_links(startfile)

for item in @link_map.keys.sort
  puts item
end

$stderr.puts "# files not found:"
for item in @fnf.keys.sort
  $stderr.puts item
end

