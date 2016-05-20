# encoding: utf-8


require "net/https"
require 'uri'
require 'rss'



class BiostarList < Nanoc3::DataSource
    identifier :biostar_list

    def fetch
        ret = []
        # which rss feed to use? latest posts or latest questions?
        # https://support.bioconductor.org/info/rss/
        # for now let's go with latest posts
        url = "https://support.bioconductor.org/feeds/latest/"
        body = HTTParty.get(url, :verify => false).body

        feed = RSS::Parser.parse(body)
        num_wanted = 5
        feed.items.each_with_index do |item, i|
            attributes = {
                :title => item.title,
                :date => item.pubDate.utc,
                :link => item.link,
                :author => "unused"
            }
            content = "unused"
            identifier = item.link.sub("https://", "/biostar_list/#{i}/")
            mtime = nil
            ret.push Nanoc3::Item.new(content, attributes, identifier, mtime)
            break if i == (num_wanted - 1)
        end
        ret
    end


    def items
        fetch()
    end

end
