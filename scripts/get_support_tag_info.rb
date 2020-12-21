require 'rest_client'
require 'json'
require 'yaml'

require_relative './badge_generation.rb'

def create_package_list()

  site_config = YAML.load_file("./config.yaml")
  rel_ver = site_config["release_version"]
  dev_ver = site_config["devel_version"]

  ## Need both devel and release branches!!!
  ## Create an empty array
  ## add elements from master and release branch
  ## array = []
  ## array.push element
  ## array.uniq!
  
  file = File.open("tmp/packageNameList.txt", "w")
  File.foreach("../manifest/software.txt"){ |line|
    if (line.start_with? "Package:")
      line.sub! "Package:", ""
      line.strip!
      file.puts line
    end
  }
  File.foreach("../manifest/data-experiment.txt"){ |line|
    if (line.start_with? "Package:")
      line.sub! "Package:", ""
      line.strip!
      file.puts line
    end
  }
  File.foreach("../manifest/workflows.txt"){ |line|
    if (line.start_with? "Package:")
      line.sub! "Package:", ""
      line.strip!
      file.puts line
    end
  }
  File.foreach("../manifest/books.txt"){ |line|
    if (line.start_with? "Package:")
      line.sub! "Package:", ""
      line.strip!
      file.puts line
    end
  }
  ## Annotation does not have complete manifest
  ## Create from 
  ## http://bioconductor.org/packages/devel/data/annotation/VIEWS
  ## or from preexisting json
  rel_file = File.join("assets/packages/json/", rel_ver, "data/annotation/packages.json")
  dev_file = File.join("assets/packages/json/", dev_ver, "data/annotation/packages.json")

  pkgs = []
  if File.exists? rel_file
    json = JSON.parse(File.read(rel_file))
    pkgs = pkgs + json.keys
  end
  if File.exists? dev_file
    json = JSON.parse(File.read(dev_file))
    pkgs = pkgs | json.keys
  end
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
