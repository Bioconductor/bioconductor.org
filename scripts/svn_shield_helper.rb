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

def downloadBadge(repo, destdir)

  if ["bioc", "workflows"].include? repo
     url = File.join("https://bioconductor.org/packages/stats/", repo, (repo+"_pkg_stats.tab"))
  else
     url = File.join("https://bioconductor.org/packages/stats/",("data-"+repo), (repo+"_pkg_stats.tab"))
  end
  urls = [url]

  d = Date.parse(Time.now.to_s)
  last6 = []
  for i in 1..6 do
    x = d << i
    last6 << [x.year.to_s, Date::ABBR_MONTHNAMES[x.month]]
  end

  raw_data = Hash.new(0)
  percentiles = {}

  urls.each do |url|
    lines = HTTParty.get(url).split("\n")
    for line in lines
      next if line =~ /^Package\tYear/ # skip header
      package, year, month, distinct_ips, downloads = line.strip.split(/\t/)
      if last6.find{|i| i == [year, month]} # was it in the last 6 full months?
        raw_data[package] = (raw_data[package] + Integer(downloads))
      end
    end
  end

  sorted_data = Hash[raw_data.sort_by(&:last).to_a.reverse]
  len = sorted_data.length.to_s

  sorted_data.each_with_index { |(key, value), index|
    dx = (index + 1).to_s
    pkg =  key
    shield = File.join(destdir, "#{pkg}.svg")
    rank = "#{dx}/#{len}"
    puts pkg
    puts rank
    resp = HTTParty.get("https://img.shields.io/badge/downloads-#{URI::encode(rank)}-blue.svg")
    if (resp.code == 200)
      fh = File.open(shield, 'w')
      fh.write(resp.to_s)
      fh.close
    else
      FileUtils.cp(File.join('assets', 'images', 'shields',
       'downloads', "unknown-downloads.svg"), shield)
    end
  }

end
