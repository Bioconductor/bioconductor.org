#!/usr/bin/env ruby

require 'find'
require 'rss'
require 'pp'
require 'yaml'
require 'fileutils'
require 'rexml/document'
include REXML

DCFDIR="tmp/build_dcfs"
OUTDIR="assets/rss/build"
BASEURL="http://bioconductor.org/checkResults"

$urls = []

def get_status(path)
    f = File.open(path)
    lines = f.readlines
    statusline = lines.find {|i| i =~ /^Status: /}
    statusline = statusline.chomp.sub(/^Status: /, "")
    statusline
end

def tweak(rss, outfile)
    outfile.gsub! /#{OUTDIR}/, ""
    outfile.gsub! /^\//, ""
    url = "http://bioconductor.org/rss/build/#{outfile}"
    xml = Document.new rss.to_s
    xml.elements.each("feed") do |e|
        hub_link = e.add_element("link")
        hub_link.attributes["rel"] = "hub"
        hub_link.attributes["href"] = $hub_url
        self_link = e.add_element("link")
        self_link.attributes["rel"] = "self"
        self_link.attributes["href"] = url
        self_link.attributes["type"] = "application/atom+xml"
        $urls.push url
        return xml.to_s
    end
end


def make_problem_feed(pkglist, config, problems, outfile)
    rss = RSS::Maker.make("atom") do |maker|
        maker.channel.author = "Bioconductor Build System"
        maker.channel.updated = Time.now.to_s
        ## FIXME: add content here:
        maker.channel.about = "http://bioconductor.org/rss/build"
        if problems.include? "WARNINGS"
            maker.channel.title = "Build Problems (including warnings)"
        else
            maker.channel.title = "Build Problems (excluding warnings)"
        end


        for key in pkglist.keys.sort  {|a,b| b.downcase <=> a.downcase}
            bad = pkglist[key].find_all {|i| problems.include? i[:status] }
            for b in bad
                maker.items.new_item do |item|
                    item.link = "#{BASEURL}/#{b[:version]}/bioc-LATEST/#{key}/#{b[:node]}-#{b[:phase]}.html"
                    item.title = "#{b[:status]} in #{b[:version]} version of #{key} on node #{b[:node]}"
                    item.updated = Time.now.to_s
                end
            end
        end
    end
    FileUtils.mkdir_p OUTDIR
    f = File.open("#{OUTDIR}/#{outfile}", "w")
    tweaked = tweak(rss, outfile)
    f.puts tweaked # rss
    f.close
end

def make_individual_feed(pkglist, config)
    rootdir = "#{OUTDIR}/packages"
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
            maker.channel.about = "http://bioconductor.org/rss/build"

            if bad.empty? and not file_exists
                maker.items.new_item do |item|
                    if pkglist[key].find {|i| i[:version] == "release"}
                        version = "release"
                    else
                        version = "devel"
                    end
                    item.link = "#{BASEURL}/#{version}/bioc-LATEST/#{key}/"
                    item.updated = Time.now.to_s
                    item.title = "No build problems for #{key}."
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
                        item.link = "#{BASEURL}/#{version}/bioc-LATEST/#{key}/"
                        item.title = "#{key} #{probs.join "/"} in #{version} on nodes #{nodes.join "/"}"
                        item.updated = Time.now.to_s
                    end
                end
                # for b in bad
                #     maker.items.new_item do |item|
                #         #item.link = "#{BASEURL}/#{b[:version]}/bioc-LATEST/#{key}/#{b[:node]}-#{b[:phase]}.html"
                #         item.link = "#{BASEURL}/#{version}/bioc-LATEST/#{key}/"
                #         item.title = "#{b[:status]} in #{b[:version]} version of #{key} on node #{b[:node]}"
                #         item.updated = Time.now.to_s
                #     end
                # end
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
            next unless path =~ /\.dcf$/
            segs = path.split("/")
            pkg = segs.last.split(".").first
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
    make_problem_feed(pkglist, config, ["ERROR", "TIMEOUT"], "errors.rss")
    make_individual_feed(pkglist, config)
    puts "Done at #{Time.now.to_s}"
    FileUtils.rm_f "tmp/rss_urls.txt"
    urlfile = File.open("tmp/rss_urls.txt", "w")
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
