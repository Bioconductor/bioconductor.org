require 'yaml'
require 'httparty'

def get_list_of_packages(bioc=true, release=false)
    config = YAML.load_file("./config.yaml")

    if release
      version = config['release_version']
      branch = "branches/RELEASE_#{version.sub('.', '_')}"
    else
      version = config['devel_version']
      branch = 'trunk'
    end

    if bioc
      repos = 'bioconductor'
      extra = 'madman/Rpacks'
      manifest_file = "bioc_#{version}.manifest"
    else
      repos = 'bioc-data'
      extra = 'experiment/pkgs'
      manifest_file = "bioc-data-experiment.#{version}.manifest"
    end


    url = "https://hedgehog.fhcrc.org/#{repos}/#{branch}/#{extra}/#{manifest_file}"
    auth = {:username => "readonly", :password => "readonly"}
    resp = HTTParty.get(url, :basic_auth => auth)
    resp.to_s.split("\n").find_all{|i| i =~ /^Package:/}.map{|i| i.sub("Package:", "").strip}.sort_by{|i|i.downcase}
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