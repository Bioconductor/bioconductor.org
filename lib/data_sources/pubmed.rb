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
      :term => "bioconductor",
      :db => ["pubmed", "pmc"],
      :sort => ["", "electronic+pub+date"]
    }
    # read configuration
    @opts.each do |key, val|
      @opts[key] = self.config[key] unless self.config[key].nil?
    end
    # convert to array if necessary
    [:db, :sort].each do |key|
      @opts[key] = [@opts[key]] unless @opts[key].is_a?(Array)
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
      puts "Publication cache expired, querying NCBI databases"
      res = []
      # query individual databases 
      @opts[:db].each_with_index do |db, i|
	res = res + query_ncbi(db,  @opts[:sort][i])
      end
      # process results
      if !res.empty?
	# remove dups
	res.uniq! { |x| x[:doi] }
	# sort by date
	res.sort! { |x, y| y[:date] <=> x[:date] }
	# take the top ones
	entries = res[0, @opts[:retmax]]
	write_cache({
	  :timestamp => Time.now.utc, # mark the retrieval time
	  :entries => entries
	})
	puts "Publication data saved"
      end
    else
      puts "Using cached publication data"
    end
    
    entries
  end
  
  def query_ncbi(db, sort)
    print(db.upcase+"... ")
    baseurl = @opts[:baseurl]

    ## search
    search = "#{baseurl}esearch.fcgi?db=#{db}&term=#{@opts[:term]}&retmax=#{@opts[:retmax]}&sort=#{sort}"
    doc = getXML(search)
    return [] if doc.nil?
    
    id_list = doc.xpath("/eSearchResult/IdList/Id")
    return [] if id_list.empty?

    ## query for results 
    query = "#{baseurl}esummary.fcgi?db=#{db}&id=#{join(id_list)}"
    doc = getXML(query)
    return [] if doc.nil?

    items = doc.xpath("/eSummaryResult/DocSum")
    return [] if items.empty?
    
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
        attributes[key] = content unless content.length == 0 # omit missing entries
      }
      
      begin
	date = Time.strptime(attributes[:epub], "%Y %b %d")
      rescue
	begin
	  date = Time.strptime(attributes[:epub], "%Y/%m/%d")
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
      end
      
      attributes[:date] = date

      entries.push Nanoc3::Item.new("unused", attributes, id, nil)
    end
    puts("done")
    return entries
  end
      
  def getXML(url)
    begin
      data = HTTParty.get(url)
      rescue Timeout::Error
	puts "Timeout error"
	return nil
      rescue TCPSocket::SocketError
	puts "Socket error"
	return nil
    end
    if data.code != 200
      puts "Error: HTTP #{data.code}"
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
