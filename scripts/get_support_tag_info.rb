## temp = RestClient.post("https://support.bioconductor.org/api/tags/list/", :tags => File.new("/home/shepherd/Justpkgs.txt"))
## my_hash = JSON.parse(temp.body)
##


require 'rest_client'
require 'json'

require_relative './badge_generation.rb'

def create_package_list()
  file = File.open("tmp/packageNameList.txt", "w")
  File.foreach("../manifest/software.txt"){ |line|
    if (line.start_with? "Package:")
      line.sub! "Package: ", ""
      line.strip!
      file.puts line
    end
  }
  File.foreach("../manifest/data-experiment.txt"){ |line|
    if (line.start_with? "Package:")
      line.sub! "Package: ", ""
      line.strip!
      file.puts line
    end
  }
  File.foreach("../manifest/workflows.txt"){ |line|
    if (line.start_with? "Package:")
      line.sub! "Package: ", ""
      line.strip!
      file.puts line
    end
  }
  File.foreach("../manifest/books.txt"){ |line|
    if (line.start_with? "Package:")
      line.sub! "Package: ", ""
      line.strip!
      file.puts line
    end
  }
  ## Annotation does not have complete manifest
  ## Create from view  
  ## http://bioconductor.org/packages/devel/data/annotation/VIEWS

  
  file.close()
  
end


def get_support_tag_info()

  
  
end
