#!/usr/bin/env ruby

require 'yaml'
require 'httparty'
require 'date'
require 'fileutils'
require 'descriptive_statistics'
require 'sequel'
require 'net/http'
require 'json'

################################################################
#
# This script contains main functions and helper functions
# to generate the badges on the package landing pages
#
################################################################


############################################
#
# 1. Platforms:
#    Also known as availability
#
#  See Rake task: get_availability_shields
#
############################################

def get_availability(item, numeric_version, reldev)

  img = platform_availability(item)
  availabilityBadge(item['Package'], img, numeric_version, reldev)

end

def availabilityBadge(pkg, img, numeric_version, reldev)

  puts "Creating badge for #{pkg} :  #{img}"
  srcdir = File.join('assets', 'images', 'shields', 'availability')
  destdir = File.join('assets', 'shields', 'availability', reldev)
  FileUtils.mkdir_p destdir
  src = File.join(srcdir, "#{img}.svg")
  dest = File.join(destdir, "#{pkg}.svg")
  res = FileUtils.copy(src, dest)
  puts("    copied #{src} to #{dest}")

end

def platform_availability(item)

  pkg = item['Package']
  unsupported = item['UnsupportedPlatforms']
  status = item['PackageStatus']
  img = "unknown-build"
  if status == "Deprecated"
    img = "none"
  else
    if unsupported == "None" || unsupported.nil?
      img = "all"
    else
      img = "some"
    end
  end
  return(img)

end

#########################################
#
# 2. Rank:
#    lower number more downloads
#
#  See rake task:  process_downloads_data
#
#########################################


def downloadBadge(repo, destdir, release=false)

  filtered_data = getRanking(repo, release)
  len = filtered_data.length.to_s

  filtered_data.each_with_index { |(key, value), index|
    pkg =  key
    shield = File.join(destdir, "#{pkg}.svg")
    rank = "#{value} / #{len}"
    puts pkg
    puts rank

    template = File.read(File.join('assets', 'images', 'shields', 'downloads', 'download-template.svg'))
    newbadge = template.gsub(/99999\/99999/,rank)
    newbadge = newbadge.gsub(/x=\"(765)\"/, 'x="700"')
    newbadge = newbadge.gsub(/width=\"(120)\"/, 'width="110"')
    newbadge = newbadge.gsub(/textLength=\"(750)\"/, '')
    File.open(shield, "w") { |file| file.write(newbadge) }

    puts "done"
  }

end

def getRanking(repo, release=false)

  site_config = YAML.load_file("./config.yaml")
  if release
    ver = site_config["release_version"]
  else
    ver = site_config["devel_version"]
  end

  if ["bioc", "workflows"].include? repo
     url = File.join("https://bioconductor.org/packages/stats/", repo, (repo+"_pkg_scores.tab"))
     json_file = File.join("assets/packages/json/", ver, repo, "packages.json")
  else
     url = File.join("https://bioconductor.org/packages/stats/",("data-"+repo), (repo+"_pkg_scores.tab"))
     json_file = File.join("assets/packages/json/", ver, "data", repo, "packages.json")
  end
  urls = [url]

  raw_data = Hash.new(0)

  urls.each do |url|
    url2 = URI.parse(url)
    req = Net::HTTP.new(url2.host, url2.port)
    req.use_ssl = true
    res = req.request_head(url2.path)
    if res.code == "200"
      lines = HTTParty.get(url).split("\n")
      for line in lines
        next if line =~ /^Package\tDownload_score/ # skip header
        package, distinct_ips = line.strip.split(/\t/)
        raw_data[package] = Integer(distinct_ips)
      end
    else
      if File.exists? json_file
        json = JSON.parse(File.read(json_file))
        json.keys.each do |pkg|
          raw_data[pkg] = json[pkg]["Rank"]
        end
      end
    end
  end

  sorted_data = Hash[raw_data.sort_by(&:last).to_a.reverse]

  # filter on above helpers for active packages
  case repo
  when "workflows"
      pkgs = get_list_of_workflows(release=release)
  when "annotation"
      pkgs = get_annotation_package_list(release=release)
  when "experiment"
      pkgs = get_list_of_packages(bioc=false, release=release)
  when "bioc"
      pkgs = get_list_of_packages(bioc=true, release=release)
  end

  filtered_data = sorted_data.select{ |k,v| pkgs.include?(k) }
  # add packages with no download stats yet
  nostats = pkgs.reject{|x| filtered_data.keys.include? x}
  nostats.each do |pkg|
      filtered_data[pkg] = 0
  end

  # add sorting ranking for ties
  # ties will have highest rank, i.e  if 1:3 are all the same
  # 1:3 get ranked 1/4 then 4 would rank 4/4 etc ...
  preVal = 0
  rankHash = Hash.new(0)
  filtered_data.each_with_index { |(key, value), index|
    if preVal != value
      rankHash["#{key}"] = index+1
      preVal = value
    else
      rankHash["#{key}"] = rankHash[rankHash.keys[index-1]]
    end
  }
  rankHash

end

######################################
#
# 3. Posts:
#    support site tags
#
#  See rake task: get_support_tag_info
#
######################################

# tried to move code to this file
# however then requires knowledge of
# Bioconductor sql password for the support
# site.
# moved back to separate file which requires
# the script only when trying to generate the
# badge


# sql password requirement is legacy
# consider moving code here 

# See get_support_tag_info.rb

############################################
#
# 4. In Bioc:
#    in bioconductor since
#
#  See rake task: get_years_in_bioc_shields
#
############################################

# See lib/helpers.rb
# get_year_shield / years_in_bioc / since
# That code remains there so it is accessible when used for
# package landing page



######################################
#
# 5. build:
#    build status
#
#   See rake task: get_build_dbs
#
######################################

def generate_build_shields(outdir, build_db, version)
  FileUtils.rm_rf outdir
  FileUtils.mkdir_p  outdir
  data = File.readlines(build_db)
  packages = data.map{|i| i.split('#').first}.uniq
  site_config = YAML.load_file("./config.yaml")
  if (version == "devel")
    activebuilders = site_config["active_devel_builders"].values
  else
    activebuilders = site_config["active_release_builders"].values
  end

  colors = {"OK" => "green", "WARNINGS" => "yellow",
    "ERROR" => "red", "TIMEOUT" => "AA0088"}

  srcdir = File.join("assets", "images", "shields", "builds")

  for package in packages
    relevantAll = data.find_all{|i| i =~ /^#{package}#/}
    builderName = relevantAll.map {|i| i.split('#')[1]}
    relevant = Array.new
    indexKeepLog = builderName.map {|i| activebuilders.include?(i) }
    indexKeepLog.each_index.select{|i|
      if (indexKeepLog[i] == true)
        relevant.push(relevantAll[i])
      end
    }
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
    elsif statuses.include? "NA"
      final_status = "OK"
    elsif statuses.include? "skipped"
      final_status = "OK"
    else
      raise "Logic fail!"
    end
    color = colors[final_status]
    puts("creating shield for #{package} in #{outdir}...") # remove me?
    FileUtils.cp(File.join(srcdir, "#{final_status}.svg"), File.join(outdir, "#{package}.svg"))
    sleep(1)
  end
end


###############################################
#
# 6. updates:
#    last commit date
#
#  See rake task: get_last_commit_date_shields
#
###############################################

#
# All badge generation task executed
# by Rake task
#


###########################################
#
# 7. dependencies:
#    dependency count
#
#  See rake task: process_dependency_badge
#
###########################################

def dependencyBadge(repo, destdir, release=false)
  site_config = YAML.load_file("./config.yaml")
  numeric_version = (release) ? site_config['release_version'] : site_config['devel_version']
  json_file = (repo == 'experiment') ? File.join("assets", "packages", "json", numeric_version, "data", "experiment","packages.json") : File.join("assets", "packages", "json", numeric_version, repo, "packages.json")
  json_obj = JSON.parse(File.read(json_file))
  counts = []
  json_obj.each_key { |key|
    dep = json_obj[key]["dependencyCount"]
    if not dep.nil?
      counts.push(dep.to_i)
    end
  }
  for pkg in json_obj.keys
    puts "#{pkg}"
    info = json_obj[pkg]
    shield = File.join(destdir, "#{pkg}.svg")
    if (info.key?("dependencyCount"))
      cnt = info["dependencyCount"]
      puts "#{pkg} : #{cnt}"
      template = File.read(File.join('assets', 'images', 'shields', 'dependencies', 'dependency-template.svg'))
      newbadge = template.gsub(/9999/, cnt)
      if counts.percentile(80).floor <= cnt.to_i
        newbadge = newbadge.gsub(/#007ec6/, '#dd8822')
      end
      if counts.percentile(95).floor <= cnt.to_i
        newbadge = newbadge.gsub(/#dd8822/, '#bb3333')
      end
      newbadge = newbadge.gsub(/textLength=\"(270)\"/, '')
      File.open(shield, "w") { |file| file.write(newbadge) }
    else
      puts "#{pkg} : dependencies not found"
      FileUtils.cp(File.join('assets', 'images', 'shields', 'dependencies', 'unknown-dependencies.svg'), shield)
    end
  end
  puts "done"

end




######################################
#
#
# Other helper functions
#
#
######################################

def get_list_of_packages(bioc=true, release=false)
    config = YAML.load_file("./config.yaml")
    path = "../manifest"

    if release
      version = config['release_version']
      branch = "RELEASE_#{version.sub('.', '_')}"
    else
      version = config['devel_version']
      branch = 'devel'
    end

    if bioc
      manifest_file = "software.txt"
    else
      manifest_file = "data-experiment.txt"
    end

    filepath = "#{path}/#{manifest_file}"
    # check about git pull and forcing overwrite
    system("git -C #{path} checkout #{branch} && git -C #{path} pull")

    file = File.open(filepath, "rb")
    contents = file.read
    pkgs = contents.to_s.split("\n").find_all{|i| i =~ /^Package:/}.map{|i| i.sub("Package:", "").strip}.sort_by{|i|i.downcase}
    file.close
    system("git -C #{path} checkout devel")
    pkgs
end

def get_annotation_package_list(release=false)
  pkgs = []
  if release
    version = 'release'
  else
    version = 'devel'
  end
  url = "http://bioconductor.org/packages/#{version}/data/annotation/src/contrib/PACKAGES"
  res = HTTParty.get(url).to_s
  res.split("\n").each do |line|
    if line =~ /^Package: /
      pkgs << line.strip.sub(/^Package: /, "")
    end
  end
  pkgs
end

def get_list_of_workflows(release=false)
    config = YAML.load_file("./config.yaml")
    path = "../manifest"

    if release
      version = config['release_version']
      branch = "RELEASE_#{version.sub('.', '_')}"
    else
      version = config['devel_version']
      branch = 'devel'
    end

    manifest_file = "workflows.txt"

    filepath = "#{path}/#{manifest_file}"
    # check about git pull and forcing overwrite
    system("git -C #{path} checkout #{branch} && git -C #{path} pull")

    file = File.open(filepath, "rb")
    contents = file.read
    pkgs = contents.to_s.split("\n").find_all{|i| i =~ /^Package:/}.map{|i| i.sub("Package:", "").strip}.sort_by{|i|i.downcase}
    file.close
    system("git -C #{path} checkout devel")
    pkgs
end
