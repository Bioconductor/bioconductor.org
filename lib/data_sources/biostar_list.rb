# encoding: utf-8


require "net/https"
require 'uri'
require 'rss'
require 'nanoc'


class BiostarList < Nanoc::DataSource
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
            identifier = Nanoc::Identifier.new(item.link.sub("https://", "/biostar_list/#{i}"), type: :legacy)
            mtime = nil
            ret.push new_item(content, attributes, identifier)
            break if i == (num_wanted - 1)
        end
        ret
    end


    def items
        fetch()
    end

end
