require 'rest_client'
require 'json'
require 'yaml'

require_relative './badge_generation.rb'

def create_package_list()

  site_config = YAML.load_file("./config.yaml")
  rel_ver = site_config["release_version"]
  dev_ver = site_config["devel_version"]

  pkgs = []

  ## master get from manifest so any new packages yet to build are included

  File.foreach("../manifest/software.txt"){ |line|
    if (line.start_with? "Package:")
      line.sub! "Package:", ""
      line.strip!
      pkgs.push line
    end
  }
  File.foreach("../manifest/data-experiment.txt"){ |line|
    if (line.start_with? "Package:")
      line.sub! "Package:", ""
      line.strip!
      pkgs.push line
    end
  }
  File.foreach("../manifest/workflows.txt"){ |line|
    if (line.start_with? "Package:")
      line.sub! "Package:", ""
      line.strip!
      pkgs.push line
    end
  }
  File.foreach("../manifest/books.txt"){ |line|
    if (line.start_with? "Package:")
      line.sub! "Package:", ""
      line.strip!
      pkgs.push line
    end
  }

  ## Release -- could either switch to git release branch
  ## use VIEWS
  ## or simplest option to use the existing json file

  rel_file = File.join("assets/packages/json/", rel_ver, "data/experiment/packages.json")
  if File.exists? rel_file
    json = JSON.parse(File.read(rel_file))
    pkgs = pkgs | json.keys
  end
  rel_file = File.join("assets/packages/json/", rel_ver, "bioc/packages.json")
  if File.exists? rel_file
    json = JSON.parse(File.read(rel_file))
    pkgs = pkgs | json.keys
  end
  rel_file = File.join("assets/packages/json/", rel_ver, "workflows/packages.json")
  if File.exists? rel_file
    json = JSON.parse(File.read(rel_file))
    pkgs = pkgs | json.keys
  end

  ## Annotation does not have complete manifest

  rel_file = File.join("assets/packages/json/", rel_ver, "data/annotation/packages.json")
  dev_file = File.join("assets/packages/json/", dev_ver, "data/annotation/packages.json")

  if File.exists? rel_file
    json = JSON.parse(File.read(rel_file))
    pkgs = pkgs | json.keys
  end
  if File.exists? dev_file
    json = JSON.parse(File.read(dev_file))
    pkgs = pkgs | json.keys
  end

  ## ensure unique package list
  pkgs.uniq!

  file = File.open("tmp/packageNameList.txt", "w")
  pkgs.each{ |line|
    file.puts line
  }
  file.close()

end


def get_support_tag_info()

## temp = RestClient.post("https://support.bioconductor.org/api/tags/list/", :tags => File.new("/home/shepherd/Justpkgs.txt"))
## my_hash = JSON.parse(temp.body)
##


end
