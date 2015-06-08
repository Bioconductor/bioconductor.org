#!/usr/bin/env ruby
require 'twitter'
require 'yaml'
require 'pp'
require 'nokogiri'
require 'open-uri'



config_dir = File.expand_path(File.dirname($0))

twitter_config = YAML::load(File.open("#{config_dir}/twitter.yaml"))

@user = "bioconductor" # bioconductor


@client = Twitter::REST::Client.new do |config|
  config.consumer_key    =  twitter_config["consumer_key"]
  config.consumer_secret =  twitter_config["consumer_secret"]
  config.access_token    =  twitter_config["access_token"]
  config.access_token_secret = twitter_config["access_token_secret"]
end



def collect_with_max_id(collection=[], max_id=nil, &block)
  response = yield max_id
  collection += response
  response.empty? ? collection.flatten : collect_with_max_id(collection, response.last.id - 1, &block)
end

def get_all_tweets(user)
  collect_with_max_id do |max_id|
    options = {:count => 3200, :include_rts => true}
    options[:max_id] = max_id unless max_id.nil?
    @client.user_timeline(user, options)
  end
end

tutorials_in_rss = []
tweet_hash = {}
doc = Nokogiri::XML(open("https://support.bioconductor.org/feeds/tag/tutorial/").read)
doc.xpath("//item").each do |item|
    title = item.children.find{|i| i.name == 'title'}.text
    title = title.sub(/^Tutorial/, "tutorial")
    link = item.children.find{|i| i.name == 'link'}.text
    tutorials_in_rss.push "#rstats / #Bioconductor new #{title} #{link}"
end

tweets = get_all_tweets(@user)

# one was already tweeted manually
tweeted_tutorials = ["#rstats / #Bioconductor new tutorial: Epigenomics RoadMap Project files now accessible via AnnotationHub https://support.bioconductor.org/p/68039/"]

for tweet in tweets
    text = tweet.text
    if text =~ /#rstats \/ #Bioconductor new tutorial/#%r(^http://t.co)
        tweeted_tutorials.push text
    end
end


def truncate(input)
    preamble = "#rstats / #Bioconductor new tutorial: "
    url = input.split(" ").last
    title = input.sub(preamble, "").sub(url, "")
    tweet_length =  preamble.length + title.length + 1 + 22
    return input if tweet_length <= 140
    diff = tweet_length - 140
    ok = title.length - diff
    preamble + title[0..ok] + " " + url
end




tutorials_in_rss = tutorials_in_rss.map{|i| truncate(i)}

tutorials_in_rss_without_urls = tutorials_in_rss.map{|i| segs = i.split(" "); segs.pop; segs.join(" ")}
tweeted_tutorials_without_urls = tweeted_tutorials.map{|i| segs = i.split(" "); segs.pop; segs.join(" ")}

tutorials_to_tweet = []


tutorials_in_rss_without_urls.each_with_index do |item, i|
  unless tweeted_tutorials_without_urls.include? item
    tutorials_to_tweet.push tutorials_in_rss[i]
  end

end


#tutorials_to_tweet = tutorials_in_rss - tweeted_tutorials


if tutorials_to_tweet.length >= 100
    puts "Can't send more than 100 tweets at once!"
    exit 1
end

for item in tutorials_to_tweet
    begin
        @client.update(item)
        puts "successfully tweeted #{item}..."
    rescue Exception => ex
        puts "Error tweeting #{item}: #{ex.message}"
    end
end

