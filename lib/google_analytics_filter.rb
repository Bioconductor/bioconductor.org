class GoogleAnalyticsFilter < Nanoc3::Filter
  identifier :google_analytics

  # Adds a javascript snippet to links to non-html files so that they can be tracked in google analytics
  
  
  
  
  def fixme(link)
    
    link =~ /href=['|"]([^'|"]*)/i
    href = $1
    return link if href.nil?
    return link if href =~ /^http:|^https|^mailto:/i

    return link if href =~ /\#/
    segs = href.split(".")
    return link if segs.last =~ /\/|%|\?/
    

    begin
      num = Integer(segs.last)
      return link
    rescue ArgumentError => ex
    rescue TypeError => te # ruby 1.9.2
    end
      
    if segs.length > 1 and segs.last.downcase != "html"
        link.gsub!(/<a/, %Q(<a onClick="javascript: pageTracker._trackPageview('#{href}'); "))
    end

    link
    
  end
  
  def run(content, params={})
    #<a href="http://www.example.com/files/map.pdf" onClick="javascript: pageTracker._trackPageview('/downloads/map'); ">
    regex = Regexp.compile(/(<a\s[^>]*>)/im)
    return content unless content =~ regex
    content.gsub(regex) {|i| fixme($1)}
  end
end