#!/usr/bin/env ruby


require "rubygems"
require "hpricot"
require "pp"
require "yaml"
require 'find'


if (ARGV.size != 1)
  puts "supply the root directory"
  exit 1
end

@site_config = YAML.load_file("./config.yaml")
@release_version = @site_config["release_version"]
@devel_version = @site_config["devel_version"]


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
    unless href =~ /^\//
      #puts "got relative link: #{href}"
      # take the filename,
      # remove the root from it,
      # remove the last segment
      # add href
      href = href.split("/").last
      temp = filename.sub(@rootdir, "")
      segs = temp.split "/"
      segs.pop 
      if href =~ /\.html$/
        unless segs.last == "html"
          segs.push "html"
        end
      end
      segs.push href
      href = segs.join("/")
      #puts "RELATIVE LINK = #{href}"
    end
#    puts "href = #{href}, text = #{link.inner_text}"
    if (href =~ /\.html$/)
      get_links(href) unless @link_map.has_key? href
    elsif (href =~ /\.r$|\.pdf$|\.doc$/i)
      @link_map[href] = 1
    end
  end
  
end


get_links(startfile)


# TODO - add something that knows that if devel == 2.9,
# then /packages/devel/... and /packages/2.9/... are the same,
# and remove duplicates. (same with release). Remove numbers instead
# of words, so users are taken to latest devel (or release).

def cleanup(release)

  if (release)
    num = @release_version
    str = "release"
  else
    num = @devel_version
    str = "devel"
  end

  keys_to_delete = []
  @link_map.keys.each do |k|
    if k =~ /packages\/#{num}/
      newkey = k.sub "packages/#{num}", "packages/#{str}"
      @link_map[newkey] = 1
      keys_to_delete.push k
    end
  end
  for key in keys_to_delete
    @link_map.delete key
  end
end

cleanup(true)
cleanup(false)

ok = Regexp::new("^\.\/packages\/#{@release_version}\/|^\.\/packages\/#{@devel_version}\/")
exts = Regexp::new(/\.html$|\.R$|\.pdf$|\.doc$/i)

Dir.chdir(@rootdir) do
  Find.find(".") do |path|
    next unless path =~ ok
    next unless path =~ exts
    path = path.gsub(/^\.\/packages\/#{@release_version}\//, "./packages/release/")
    path = path.gsub(/^\.\/packages\/#{@devel_version}\//, "./packages/devel/")
    key = path.gsub(/^\./, "")
    @link_map[key] = 1
  end
end


for item in @link_map.keys.sort
  puts ".#{item}"
end

$stderr.puts "# files not found:"
for item in @fnf.keys.sort
  $stderr.puts item
end

