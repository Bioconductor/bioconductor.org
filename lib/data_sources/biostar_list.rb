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

        uri = URI.parse(url)
        http = Net::HTTP.new(uri.host, uri.port)
        http.use_ssl = true
        http.verify_mode = OpenSSL::SSL::VERIFY_NONE

        request = Net::HTTP::Get.new(uri.request_uri)

        response = http.request(request)
        body = response.body
        feed = RSS::Parser.parse(response.body)
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
