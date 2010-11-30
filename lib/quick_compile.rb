

if ENV.has_key?("QUICK_NANOC_COMPILE") && ENV["QUICK_NANOC_COMPILE"] == "true"
  puts "suppressing GMANE and biocViews compilation"
  base_command = nil 
  ObjectSpace.each_object(Nanoc3::CLI::Base) do |obj| 
    base_command = obj 
    sources = base_command.site.config[:data_sources]#.collect{|i|i[:type] == "filesystem_unified"}
    base_command.site.config[:data_sources] = [sources.first]
  end
end