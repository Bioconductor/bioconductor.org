require 'yaml'
require 'httparty'
require 'date'
require 'fileutils'
require 'descriptive_statistics'

def get_list_of_packages(bioc=true, release=false)
    config = YAML.load_file("./config.yaml")
    path = "../manifest"

    if release
      version = config['release_version']
      branch = "RELEASE_#{version.sub('.', '_')}"
    else
      version = config['devel_version']
      branch = 'master'
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
    system("git -C #{path} checkout master")
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
      branch = 'master'
    end

    manifest_file = "workflows.txt"

    filepath = "#{path}/#{manifest_file}"
    # check about git pull and forcing overwrite
    system("git -C #{path} checkout #{branch} && git -C #{path} pull")

    file = File.open(filepath, "rb")
    contents = file.read
    pkgs = contents.to_s.split("\n").find_all{|i| i =~ /^Package:/}.map{|i| i.sub("Package:", "").strip}.sort_by{|i|i.downcase}
    file.close
    system("git -C #{path} checkout master")
    pkgs
end

def downloadBadge(repo, destdir, release=false)

  filtered_data = getRanking(repo, release)
  len = filtered_data.length.to_s

  filtered_data.each_with_index { |(key, value), index|
    pkg =  key
    shield = File.join(destdir, "#{pkg}.svg")
    rank = "#{value} / #{len}"
    puts pkg
    puts rank
    resp = HTTParty.get("https://img.shields.io/badge/rank-#{URI::encode(rank)}-blue.svg")
    if (resp.code == 200)
      fh = File.open(shield, 'w')
      fh.write(resp.to_s)
      fh.close
    else
      puts "ERROR: "+resp.code.to_s
      FileUtils.cp(File.join('assets', 'images', 'shields',
       'downloads', "unknown-downloads.svg"), shield)
    end
    puts "done"
  }

end

def getRanking(repo, release=false)
  if ["bioc", "workflows"].include? repo
     url = File.join("https://bioconductor.org/packages/stats/", repo, (repo+"_pkg_scores.tab"))
  else
     url = File.join("https://bioconductor.org/packages/stats/",("data-"+repo), (repo+"_pkg_scores.tab"))
  end
  urls = [url]

  raw_data = Hash.new(0)

  urls.each do |url|
    lines = HTTParty.get(url).split("\n")
    for line in lines
      next if line =~ /^Package\tDownload_score/ # skip header
      package, distinct_ips = line.strip.split(/\t/)
      raw_data[package] = Integer(distinct_ips)
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
      clr = "blue"
      if counts.percentile(80).floor <= cnt.to_i
        clr = "orange"
      end
      if counts.percentile(95).floor <= cnt.to_i
         clr =  "red"
      end
      resp = HTTParty.get("https://img.shields.io/badge/dependencies-#{cnt}-#{clr}.svg")
      if (resp.code == 200)
        fh = File.open(shield, 'w')
        fh.write(resp.to_s)
        fh.close
      else
        puts "ERROR: "+resp.code.to_s
        FileUtils.cp(File.join('assets', 'images', 'shields', 'unknown-dependencies.svg'), shield)
      end
    else
      puts "#{pkg} : dependencies not found"
      FileUtils.cp(File.join('assets', 'images', 'shields', 'unknown-dependencies.svg'), shield)
    end
  end
  puts "done"

end


def platform_availability(item)

  pkg = item['Package']
  unsupported = item['UnsupportedPlatforms']
  status = item['PackageStatus']
  img = "unknown-build"
  if status == "Deprecated"
    img = "none"
  else
    if unsupported == "None"
      img = "all"
    else
      img = "some"
    end
  end
  return(img)

end

def get_availability(item, numeric_version)

  img = platform_availability(item)
  availabilityBadge(item['Package'], img, numeric_version)

end

def availabilityBadge(pkg, img, numeric_version)

  puts "Creating badge for #{pkg} :  #{img}"
  srcdir = File.join('assets', 'images', 'shields', 'availability')
  destdir = File.join('assets', 'shields', 'availability', numeric_version)
  FileUtils.mkdir_p destdir
  src = File.join(srcdir, "#{img}.svg")
  dest = File.join(destdir, "#{pkg}.svg")
  res = FileUtils.copy(src, dest)
  puts("    copied #{src} to #{dest}")

end
