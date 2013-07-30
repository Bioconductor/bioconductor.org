#!/usr/bin/env ruby

require 'find'
require 'rss'
require 'pp'
require 'yaml'
require 'fileutils'

DCFDIR="tmp/build_dcfs"
OUTDIR="assets/rss/build"
BASEURL="http://bioconductor.org/checkResults"

def get_status(path)
    f = File.open(path)
    lines = f.readlines
    statusline = lines.find {|i| i =~ /^Status: /}
    statusline = statusline.chomp.sub(/^Status: /, "")
    statusline
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
    f.puts rss
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
            maker.channel.title = "Build Problems for #{key}"
            maker.channel.updated = Time.now.to_s
            ## FIXME: add content here:
            maker.channel.about = "http://bioconductor.org/rss/build"

            if bad.empty? and not file_exists
                maker.items.new_item do |item|
                    item.link = "#{BASEURL}/release/bioc-LATEST/#{key}/#{config["master_release_builder"]}-buildsrc.html"
                    item.updated = Time.now.to_s
                    item.title = "No build problems for #{key}."
                end
            else
                for b in bad
                    maker.items.new_item do |item|
                        item.link = "#{BASEURL}/#{b[:version]}/bioc-LATEST/#{key}/#{b[:node]}-#{b[:phase]}.html"
                        item.title = "#{b[:status]} in #{b[:version]} version of #{key} on node #{b[:node]}"
                        item.updated = Time.now.to_s
                    end
                end
            end

        end
        if (not bad.empty?) or (not file_exists)
            f = File.open(filename, "w")
            f.puts rss
            f.close
        end
    end
end

def runit()
    pkglist = {}

    config = YAML.load_file("./config.yaml")

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
    pkglist
end

#runit()
if $0 == __FILE__
    runit()
end
