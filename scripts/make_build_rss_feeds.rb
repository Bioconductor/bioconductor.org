#!/usr/bin/env ruby

require 'find'
require 'rss'
require 'pp'
require 'yaml'

DCFDIR="tmp/build_dcfs"

def get_status(path)
    f = File.open(path)
    lines = f.readlines
    statusline = lines.find {|i| i =~ /^Status: /}
    statusline = statusline.chomp.sub(/^Status: /, "")
    statusline
end


def runit()
    pkglist = {}

    config = YAML.load_file("./config.yaml")

    Dir.chdir DCFDIR do
        Find.find "." do |path|
            next unless path =~ /\.dcf$/
            segs = path.split("/")
            pkg = segs.last.split(".").first
    #        puts path
            status = get_status(path)
    #        puts status unless status == "OK"
            pkglist[pkg] = [] unless pkglist.has_key? pkg
            segs = path.split("/")
            node =  segs[3]
            version = segs[1]
            phase = segs[4]
            pkglist[pkg].push({:version => version, :node => node,
                :phase => phase, :status => status})
        end
    end

    #puts pkglist.keys
    #pp pkglist
    pkglist
end

#runit()
if $0 == __FILE__
    runit()
end
