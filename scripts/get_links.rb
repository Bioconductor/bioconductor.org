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
    href = link.attributes["HREF"] if href.nil?
    next if href.nil? or href.empty?
    href.sub!(/^http:\/\/bioconductor\.org\//i, "")
    href.sub!(/^http:\/\/www\.bioconductor\.org\//i, "")
    next if href =~ /:\/\//
    next if href =~ /^mailto:/i
    next if href =~ /\.\./ # TODO: deal with this later
    next if href =~ /^#/
    if href =~ /#/
      href = href.split("#").first
    end
    if href =~ /\/$/
      href = href + "index.html"
    end
    #href = "/#{href}" unless href =~ /^\//
    unless href =~ /^\//
      # take the filename,
      # remove the root from it,
      # remove the last segment
      # add href
      temp = filename.sub(@rootdir, "")
      segs = temp.split "/"
      segs.pop
      segs.push href
      href = segs.join("/")
      #puts "RELATIVE LINK = #{href}"
    end
#    puts "href = #{href}, text = #{link.inner_text}"
    if (href =~ /\.html$/)
      get_links(href) unless @link_map.has_key? href
    elsif (href =~ /\.r$|\.pdf$/i)
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

