# encoding: utf-8

require 'rubygems'
require 'httparty'
require 'time'
require 'yaml'
require 'fileutils'

class GmaneList < Nanoc3::DataSource
  identifier :gmane_list
  def fetch
    fetch_entries().map do |e|
      Nanoc3::Item.new(e[:content], e[:attributes], e[:identifier], e[:mtime])
    end.sort { |a, b| b[:date] <=> a[:date] }
  end

  def items
    @items ||= fetch()
  end

  def cache_file
    self.config[:cache_file]
  end

  # I would have thought that these setup actions would
  # be placed in the 'setup' method, but it doesn't seem to get
  # called.  up and down methods do get called.
  def up
    if cache_file().nil?
      raise "please add a 'cache_file' entry to config.yaml"
    end
    FileUtils.mkdir_p(File.dirname(cache_file))
  end

  private

  def fetch_entries
    cache_data = read_cache()
    entries = cache_data[:entries]
    if expired?(cache_data)
      puts "GMANE cache expired"
      res = fetch_entries_conditional(cache_data[:last_modified])
      if res
        write_cache(res)
        entries = res[:entries]
        puts "GMANE 200 UPDATE saved"
      else
        write_cache(cache_data)
      end
    else
      puts "using cached GMANE data"
    end
    entries
  end

  def fetch_entries_conditional(have_last_mod)
    bioc_list_url = self.config[:gmane_rss_url]
    data = begin
             HTTParty.get(bioc_list_url,
                          :headers => if_modified_since(have_last_mod),
                          :timeout => 6)
           rescue Timeout::Error
             puts "GMANE timeout error"
             return nil
           rescue TCPSocket::SocketError
             puts "Socket error"
             return nil
           end
    if data.code == 304
      puts "GMANE 304 content not modified"
      return nil
    end
    if data.code >= 400
      puts "GMANE error: HTTP #{data.code}"
      return nil
    end
    if (data["rdf:RDF"].nil?)
      parent = "RDF"
    else
      parent = "rdf:RDF"
    end
    entries = data[parent]["item"].map do |item|
      attributes = {
        :title => item["title"],
        :date => fix_date(item),
        :link => item["link"],
        :author => item["dc:creator"]
      }
      content = item["description"]
      mtime = nil
      identifier = "/#{File.basename(item["link"])}/"
      { :content => content,
        :attributes => attributes,
        :identifier => identifier,
        :mtime => mtime }
    end
    last_mod = data.headers["last-modified"].first rescue nil
    { :entries => entries,
      :last_modified => last_mod }
  end

  def fix_date(rss_item)
    key = nil
    for cand in ["dc_date", "date"]
      key = cand if rss_item.has_key? cand
    end

    Time.parse(rss_item[key] + " GMT").utc
  rescue
    Time.new.utc
  end

  def expired?(cache_data)
    if cache_data && cache_data[:cache_time]
      ttl = self.config[:ttl].to_i * 60 # to seconds
      cache_data[:cache_time] + ttl < Time.now.utc
    else
      true
    end
  end

  def read_cache
    if File.exist? cache_file
      YAML.load_file(cache_file)
    else
      { :cache_time => nil,
        :last_modified => nil,
        :entries => []
      }
    end
  end
  
  def write_cache(cache_data)
    cache_data[:cache_time] = Time.now.utc
    open(cache_file, "w") do |f|
      f.write(cache_data.to_yaml)
    end
  end

  def if_modified_since(t)
    t ? { "if-modified-since" => t } : { }
  end

end
