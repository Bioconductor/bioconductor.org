require 'yaml'
require 'httparty'

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
