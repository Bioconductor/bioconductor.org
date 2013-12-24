#!/usr/bin/env ruby
require 'twitter'
require 'yaml'
require 'pp'
require "rexml/document"
require 'open-uri'

include REXML


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

pkgs_in_rss = []
tweet_hash = {}
doc = REXML::Document.new(open("http://www.bioconductor.org/rss/new_packages.rss").read)
doc.elements.each("rss/channel/item/title") do |element|
    key = element.text.split(" ")[1]
    pkgs_in_rss.push key
    tweet_hash[key] = element.text
end

tweets = get_all_tweets(@user)

tweeted_pkgs = []

for tweet in tweets
    text = tweet.text
    if text =~ %r(^http://t.co)
        tweeted_pkgs.push(text.split(" ")[1])
    end
end

pkgs_to_tweet = pkgs_in_rss - tweeted_pkgs

# fixme? if there is more than one url in the tweet this will be thrown off.
def calculate_length(tweet)
    url = tweet.split(" ").first
    txt = tweet.sub(%r(^#{url}), "")
    22 + txt.length
end


def truncate(input)
    url = input.split(" ").first
    txt = input.sub(%r(^#{url}), "")
    tweet_length = (22 + txt.length)
    return input if tweet_length <= 140
    ok = 140 - 22
    url + txt[0...ok]
end


if pkgs_to_tweet.length >= 100
    puts "Can't send more than 100 tweets at once!"
    exit 1
end

for item in pkgs_to_tweet
    tweet = truncate(tweet_hash[item])
    begin
        @client.update(tweet)
        puts "successfully tweeted #{item}..."
    rescue Exception => ex
        puts "Error tweeting #{item}: #{ex.message}"
    end
end

