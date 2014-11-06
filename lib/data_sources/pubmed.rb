# encoding: utf-8
require 'nokogiri'
require 'httparty'
require 'nanoc3'

class PubmedPapers < Nanoc3::DataSource
  identifier :pubmed_papers    

  def up
    # set default options here
    @opts = {
      :cache_file => "tmp/pubmed_cache_file.yaml",
      :retmax => 20,
      :ttl => 24,
      :baseurl => "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
      :db => "pubmed",
      :term => "bioconductor"
    }
    # read configuration
    @opts.each do |key, val|
      @opts[key] = self.config[key] unless self.config[key].nil?
    end
    FileUtils.mkdir_p(File.dirname(@opts[:cache_file]))
  end
    
  def items
    fetch()
  end
    
    
  private 
    
  def fetch
    cache_data = read_cache()      
    entries = cache_data[:entries]
      
    if expired?(cache_data)
      puts "PubMed cache expired, querying NCBI"
      res = query_pubmed
      if res
	write_cache(res)
	entries = res[:entries]
	puts "PubMed 200 UPDATE saved"
      end
    else
      puts "Using cached PubMed data"
    end
    
    entries
  end
  
  def query_pubmed    
    baseurl = @opts[:baseurl]
    db =  @opts[:db]
    retmax = @opts[:retmax]
    term = @opts[:term]

    ## search
    search = "#{baseurl}esearch.fcgi?db=#{db}&term=#{term}&retmax=#{retmax}"
    doc = getXML(search)
    return nil if doc.nil?
    
    id_list = doc.xpath("/eSearchResult/IdList/Id")
    return nil if id_list.empty?

    ## query for results 
    query = "#{baseurl}esummary.fcgi?db=#{db}&id=#{join(id_list)}"
    doc = getXML(query)
    return nil if doc.nil?

    items = doc.xpath("/eSummaryResult/DocSum")
    return nil if items.empty?
    
    # XML attribute name mapping  
    mapping = {
      :date => "PubDate",
      :epub => "EPubDate",
      :pubmed => ["History", "pubmed"],
      :title => "Title",
      :author => ["AuthorList", "Author"],
      :journal => "Source",
      :volume => "Volume",
      :issue => "Issue",
      :pages => "Pages",
      :doi => "DOI"
    }
    
    ## iterate over results
    entries = []
    items.each do |item|
      id = "/#{item.xpath("./Id").text}/"
      attributes = {}
      mapping.each {
	|key, val| 
        content = extract(item, val)
        attributes[key] = content unless content.length == 0 # ommit missing entries
      }
      
      begin
	date = Time.strptime(attributes[:epub], "%Y %b %d")
      rescue
	begin
	  date = Time.strptime(attributes[:date], "%Y %b %d")
	rescue
	  begin
	    date = Time.strptime(attributes[:pubmed],  "%Y/%m/%d %H:%M")
	  rescue
	    date = Time.new(Time.now.year-3, Time.now.month, Time.now.day) # if no clue about the actual date, just put it 3 years in the past
	  end
	end
      end
      
      attributes[:date] = date

      entries.push Nanoc3::Item.new("unused", attributes, id, nil)
    end
    
    ## sort by date
    entries.sort! { |x, y| y.attributes[:date] <=> x.attributes[:date] }
    
    return {
      :timestamp => Time.now.utc, # mark the retrieval time
      :entries => entries
    }     
  end
      
  def getXML(url)
    begin
      data = HTTParty.get(url)
      rescue Timeout::Error
	puts "PubMed timeout error"
	return nil
      rescue TCPSocket::SocketError
	puts "Socket error"
	return nil
    end
    if data.code != 200
      puts "Pubmed error: HTTP #{data.code}"
      return nil
    else
      return Nokogiri::XML(data.body)
    end
  end
    
  ## join the content of NodeSet elements
  def join(x, sep = ",")
    y = []
    x.each do |element|
      y.push element.text
    end
    y.join(sep)
  end

  ## function used to extract XML content
  def extract(item, names, sep = ", ")
    names = Array(names)
    paths = names.map { |name| "/Item[@Name='#{name}']" }
    join(item.xpath(".#{paths.join}"), sep = sep)
  end
  
  ## cache management
  
  def read_cache
    if File.exist? @opts[:cache_file]
      YAML.load_file(@opts[:cache_file])
    else
      { 
	:timestamp => nil,
        :entries => []
      }
    end
  end
  
  def write_cache(cache_data)
    open(@opts[:cache_file], "w") do |f|
      f.write(cache_data.to_yaml)
    end
  end
  
  def expired?(cache_data)
    if cache_data[:timestamp]
      ttl = @opts[:ttl].to_i * 60**2 # expiration time in hours
      cache_data[:timestamp] + ttl < Time.now.utc
    else
      true
    end
  end
  
end
