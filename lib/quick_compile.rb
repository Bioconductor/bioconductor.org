

if ENV.has_key?("QUICK_NANOC_COMPILE") && ENV["QUICK_NANOC_COMPILE"] == "true"
  puts "suppressing GMANE and biocViews compilation"
  base_command = nil 
  ObjectSpace.each_object() do |obj| 
    next unless obj.respond_to? :site
    base_command = obj 
    begin
      sources = base_command.site.config[:data_sources]
      base_command.site.config[:data_sources] = [sources.first]
    rescue
    end
  end
end


