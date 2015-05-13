require 'httparty'
require 'fileutils'

def generate_build_shields(outdir, build_db)
  FileUtils.rm_rf outdir
  FileUtils.mkdir_p  outdir
  data = File.readlines(build_db)
  packages = data.map{|i| i.split('#').first}.uniq

  colors = {"OK" => "green", "WARNINGS" => "yellow",
    "ERROR" => "red", "TIMEOUT" => "AA0088"}

  srcdir = File.join("assets", "images", "shields", "builds")

  for package in packages
    relevant = data.find_all{|i| i =~ /^#{package}#/}
    statuses = relevant.map {|i| i.split(' ').last.strip}
    statuses = statuses.reject{|i| i == "NotNeeded"}
    statuses = statuses.uniq
    if statuses.length == 1 and statuses.first == "OK"
      final_status = "OK"
    elsif statuses.include? "ERROR"
      final_status = "ERROR"
    elsif statuses.include? "TIMEOUT"
      final_status = "TIMEOUT"
    elsif statuses.include? "WARNINGS"
      final_status = "WARNINGS"
    else
      raise "Logic fail!"
    end
    color = colors[final_status]
    puts "creating shield for #{package} in #{outdir}..." # remove me?
    FileUtils.cp(File.join(srcdir, "#{final_status}.svg"), File.join(outdir, "#{package}.svg"))
  end
end