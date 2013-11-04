#!/usr/bin/env ruby

require 'find'
require 'rss'
require 'pp'
require 'yaml'
require 'fileutils'
require 'rubygems'
require 'uuid'
require 'rexml/document'
include REXML



if ARGV.empty?
    $repo = "bioc"
else
    unless ["bioc", "data-experiment"].include? ARGV.first
        puts "argument must be 'bioc' or 'data-experiment'"
        exit 1
    else
        $repo = ARGV.first
    end
end

$uuid = UUID.new
BASEURL="http://bioconductor.org/checkResults"

if $repo == "bioc"
    DCFDIR="tmp/build_dcfs"
    OUTDIR="assets/rss/build"
else
    DCFDIR="tmp/data_build_dcfs"
    OUTDIR="assets/rss/build/data"
end




$urls = []

def get_status(path)
    f = File.open(path)
    lines = f.readlines
    statusline = lines.find {|i| i =~ /^Status: /}
    statusline = statusline.chomp.sub(/^Status: /, "")
    statusline
end

def tweak(rss, outfile)
    #return rss.to_s if true
    outfile.gsub! /#{OUTDIR}/, ""
    outfile.gsub! /^\//, ""
    spec = ($repo=="bioc") ? "" : "/data"
    url = "http://bioconductor.org/rss/build#{spec}/#{outfile}"
    xml = Document.new rss.to_s
    hub_link =  Element.new "link" #e.add_element("link")
    hub_link.attributes["rel"] = "hub"
    hub_link.attributes["href"] = $hub_url
    self_link = Element.new "link" #e.add_element("link")
    self_link.attributes["rel"] = "self"
    self_link.attributes["href"] = url
    self_link.attributes["type"] = "application/atom+xml"
    $urls.push url

    xml.root.insert_after "//feed/updated", self_link
    xml.root.insert_after "//feed/updated", hub_link
    #xml.root.elements["entry"].add Text.new("Bioconductor Build Information")
    tmp = XPath.match xml, "//entry"
    for thing in tmp
        date = XPath.first(thing, "dc:date")
        thing.delete date
    end

    return xml.to_s
end


def make_problem_feed(pkglist, config, problems, outfile)
    rss = RSS::Maker.make("atom") do |maker|
        maker.channel.author = "Bioconductor Build System"
        maker.channel.updated = Time.now.to_s
        ## FIXME: add content here:
        maker.channel.about = "http://bioconductor.org/developers/rss-feeds/"
        if problems.include? "WARNINGS"
            maker.channel.title = "Build Problems (including warnings)"
        else
            maker.channel.title = "Build Problems (excluding warnings)"
        end


        for key in pkglist.keys.sort  {|a,b| b.downcase <=> a.downcase}
            bad = pkglist[key].find_all {|i| problems.include? i[:status] }
            for b in bad
                maker.items.new_item do |item|
                    item.link = "#{BASEURL}/#{b[:version]}/#{$repo}-LATEST/#{key}/#{b[:node]}-#{b[:phase]}.html"
                    item.title = "#{b[:status]} in #{b[:version]} version of #{key} on node #{b[:node]}"
                    item.summary = item.title
                    item.updated = Time.now.to_s
                    item.id = "#{item.link}?id=#{Time.now.to_i}_#{$uuid.generate}"
                end
            end
        end
    end
    FileUtils.mkdir_p OUTDIR
    f = File.open("#{OUTDIR}/#{outfile}", "w")
    tweaked = tweak(rss, outfile)
    f.puts tweaked 
    f.close
end

def make_individual_feed(pkglist, config)
    rootdir =  "#{OUTDIR}/packages" 
    FileUtils.mkdir_p rootdir
    for key in pkglist.keys
        filename = "#{rootdir}/#{key}.rss"
        file_exists = File.exist?(filename) and File.file?(filename)
        bad = pkglist[key].find_all {|i| i[:status] != "OK"}
        rss = RSS::Maker.make("atom") do |maker|
            maker.channel.author = "Bioconductor Build System"
            maker.channel.title = "#{key} Build Problems"
            maker.channel.updated = Time.now.to_s
            ## FIXME: add content here:
            maker.channel.about = "http://bioconductor.org/developers/rss-feeds/"

            if bad.empty? and not file_exists
                maker.items.new_item do |item|
                    if pkglist[key].find {|i| i[:version] == "release"}
                        version = "release"
                    else
                        version = "devel"
                    end
                    item.link = "#{BASEURL}/#{version}/#{$repo}-LATEST/#{key}/"
                    item.updated = Time.now.to_s
                    item.title = "No build problems for #{key}."
                    item.summary = item.title
                    item.id = "#{item.link}?id=#{Time.now.to_i}_#{$uuid.generate}"
                end
            else
                relprobs = bad.find_all {|i| i[:version] == "release"}
                devprobs = bad.find_all {|i| i[:version] == "devel"}
                os = {"linux" => 1, "windows" => 2, "mac" => 3}
                for ary in [relprobs, devprobs]
                    machines = ary == relprobs ? config["active_release_builders"] : config["active_devel_builders"]
                    ary = ary.find_all{|i| machines.values.include? i[:node]}
                    ary.sort! do |a, b|
                        nodea = a[:node]
                        nodeb = b[:node]
                        osa = machines.find{|k,v| v == nodea}.first
                        osb = machines.find{|k,v| v == nodeb}.first
                        if (os[osa] > os[osb])
                            1
                        elsif os[osa] < os[osb]
                            -1
                        else
                            0
                        end
                    end                    
                    next if ary.empty?
                    version = ary.first[:version]
                    probs = ary.collect{|i| i[:status]}
                    nodes = ary.collect{|i| i[:node]}
                    maker.items.new_item do |item|
                        item.link = "#{BASEURL}/#{version}/#{$repo}-LATEST/#{key}/"
                        nword = (nodes.length > 1) ? "nodes" : "node"
                        item.title = "#{key} #{probs.join "/"} in #{version} on #{nword} #{nodes.join "/"}"
                        item.summary = item.title
                        item.updated = Time.now.to_s
                        item.id = "#{item.link}?id=#{Time.now.to_i}_#{$uuid.generate}"

                    end
                end
            end

        end
        if (not bad.empty?) or (not file_exists)
            f = File.open(filename, "w")
            #puts filename
            tweaked = tweak(rss, filename)
            f.puts tweaked #rss 
            f.close
        end
    end
end

def runit()
    pkglist = {}

    config = YAML.load_file("./config.yaml")
    $hub_url = config["rss_hub_url"]
    Dir.chdir DCFDIR do
        Find.find "." do |path|
            next unless path =~ /-summary\.dcf$/
            segs = path.split("/")
            segs2 = segs.last.split(".")
            segs2.pop
            segs2.pop
            pkg = segs2.join(".")
            status = get_status(path)
            pkglist[pkg] = [] unless pkglist.has_key? pkg
            segs = path.split("/")
            node =  segs[3]
            version = segs[1]
            phase = segs[4]
            pkglist[pkg].push({:version => version, :node => node,
                :phase => phase, :status => status})
        end
    end

    make_problem_feed(pkglist, config, ["ERROR", "WARNINGS", "TIMEOUT"],
        "problems.rss")
    make_problem_feed(pkglist, config, ["ERROR", "TIMEOUT"],
        "errors.rss")
    make_individual_feed(pkglist, config)
    puts "Done at #{Time.now.to_s}"
    rssfile = ($repo == "bioc") ? "tmp/rss_urls.txt" : "tmp/data_rss_urls.txt"
    FileUtils.rm_f rssfile
    urlfile = File.open(rssfile, "w")
    for url in $urls
        urlfile.puts url
    end
    urlfile.close
    pkglist
end

#runit()
if $0 == __FILE__
    runit()
end
