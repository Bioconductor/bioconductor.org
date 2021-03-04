require 'rest_client'
require 'json'
require 'yaml'

require_relative './badge_generation.rb'

def create_package_list()

  site_config = YAML.load_file("./config.yaml")
  rel_ver = site_config["release_version"]
  dev_ver = site_config["devel_version"]

  puts "getting package list"
  pkgs = []
  [true, false].each do |state|
    pkgs += get_list_of_packages(state)
    pkgs += get_list_of_packages(bioc=false, state)
    pkgs += get_annotation_package_list(state)
    pkgs += get_list_of_workflows(state)
  end
  pkgs = pkgs.uniq

  puts "writing package file"
  file = File.open("tmp/packageNameList.txt", "w")
  pkgs.each{ |line|
    file.puts line
  }
  file.close()

end


def get_support_tag_info()

  if File.exists? "tmp/packageNameList.txt"
    res = RestClient.post("https://support.bioconductor.org/api/tags/list/", :tags => File.new("tmp/packageNameList.txt"))
    my_hash = JSON.parse(res.body)


    zero_shield = File.join("assets", "images", "shields", "posts",
                            "zero.svg")
    dest_dir = File.join("assets", "shields", "posts")
    ## remove dir first?
    FileUtils.mkdir_p dest_dir

    pkgs = []
    [true, false].each do |state|
      pkgs += get_list_of_packages(state)
      pkgs += get_list_of_packages(bioc=false, state)
      pkgs += get_annotation_package_list(state)
      pkgs += get_list_of_workflows(state)
    end
    pkgs = pkgs.uniq

    for pkg in pkgs
      puts "getting shield for #{pkg}"
      if my_hash.has_key? pkg.downcase
        value = my_hash[pkg.downcase]
        answer = value["answer_count"]
        total = value["total"]
        shield_text = "#{answer} / #{total}"
        shield = File.join(dest_dir, "#{pkg}.svg")
        template = File.read(File.join('assets', 'images', 'shields', 'posts', 'posts-template.svg'))
        newbadge = template.gsub(/0000 \/ 0000/, shield_text)
        File.open(shield, "w") { |file| file.write(newbadge) }
      else
        puts "zero_shield"
        FileUtils.cp zero_shield, File.join(dest_dir, "#{pkg}.svg")
      end
    end

  else
    puts "ERROR: tmp/packageNameList.txt was not created"
  end

end
