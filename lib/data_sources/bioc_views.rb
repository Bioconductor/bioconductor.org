class BiocViews < Nanoc::DataSource

  require 'pp'
  require 'rubygems'
  require 'json'

  identifier :bioc_views

  def hard_coded_repos()
    {
      "bioc/" => "Software", "data/annotation/" => "AnnotationData",
      "data/experiment/" => "ExperimentData", "workflows/" => "Workflow"
    }
  end

  def get_repos(version, repos)
    return hard_coded_repos() unless version == @site_config["devel_version"]
    h = {}

    for repo in @site_config["devel_repos"]
      key = "#{repo}/"
      h[key] = repos[key]
    end
    h
  end

  # TODO - Find out if there is a way to skip items() altogether if
  # things are not found in up()
  def up
    @bad_packages = ["snpMatrix2"] # don't process these
    @repos = hard_coded_repos()
    @good_to_go = true

    @site_config = YAML.load_file("./config.yaml")

    for version in @site_config["versions"]
      @repos = get_repos(version, @repos)
      @repos.each_pair do |k,v|

        dir = "#{config[:json_dir]}/#{version}/#{k}"
        unless test(?f, "#{dir}/packages.json")
#  if one repo fails they all fail
#          @good_to_go = false
          puts "'packages.json' missing in #{dir}"
        end

      end
    end

    @all_packages = {}

    if @good_to_go
      for version in @site_config["versions"]
        hsh = {}

        @repos = get_repos(version, @repos)

        @repos.each_pair do |k,v|
          key = k.gsub(/\/$/, "")
          dir = "#{config[:json_dir]}/#{version}/#{k}"
          if File.exist?("#{dir}/packages.json")
              json_file = File.open("#{dir}/packages.json")
              obj = JSON.parse(json_file.readlines.join("\n"))
              #@all_packages[key] = obj
              hsh[key] = obj
              for bad in @bad_packages
                hsh.delete(bad)
              end
          else
            puts "ERROR: BiocViews DataSource: no JSON files; no initialization"
          end
        end
        @all_packages[version] = hsh
      end
    else
      puts "ERROR: BiocViews DataSource: no JSON files; no initialization"
    end
  end

  def get_index_page(packages, repo, version)

    info = []

    packages.each_pair do |k,v|
      hsh = {}
      hsh[:name] = k
      hsh[:Maintainer] = v["Maintainer"]
      hsh[:Title] = v["Title"]
      info.push hsh
    end

    newinfo = info.sort{|a,b| a[:name].downcase <=> b[:name].downcase}

    bioc_version_num = version
    if (version == @site_config["release_version"])
      bioc_version_str = "Release"
    elsif (version == @site_config["devel_version"])
      bioc_version_str = "Development"
    else
      bioc_version_str = nil
    end

    title = "#{version} #{repo} Packages"

    subnav = []
    subnav.push({:include => "/_bioc_release_packages/"})
    subnav.push({:include => "/_bioc_devel_packages/"})
    subnav.push({:include => "/_bioc_older_packages/"})

    attributes = {
        :rebase => true,
        :package_index_page => true,
        :info => newinfo,
        :bioc_version_num => bioc_version_num,
        :bioc_version_str => bioc_version_str,
        :repo => repo,
        :title => title,
        :subnav => subnav
    }

    item = new_item("", attributes, Nanoc::Identifier.new("all-#{repo}-#{version}", type: :legacy))
    rep = Nanoc::Int::ItemRep.new(item, :package_index_page)


    item
  end

  def items

    unless @good_to_go
      puts "ERROR: BiocViews DataSource: no JSON files; detail pages not built"
      return []
    end

    items = []

    link_list = [:Depends, :Imports, :Suggests, :Enhances,
      :LinkingTo, :dependsOnMe, :importsMe, :suggestsMe, :linksToMe]


    for version in @site_config["versions"]

      @repos = get_repos(version, @repos)
      @repos.each_pair do |k,v|
        dir = "#{config[:json_dir]}/#{version}/#{k}"

        if File.exist?("#{dir}/packages.json")
            json_file = File.open("#{dir}/packages.json")

            packages = JSON.parse(
              File.read(json_file,
                :external_encoding => 'utf-8',
                :internal_encoding => 'utf-8'
              )
            )
            for bad in @bad_packages
              packages.delete(bad)
            end

            items.push(get_index_page(packages, v, version))

            for package in packages.keys
              repo = k
              id = "/#{version}/#{repo}#{package}/"
              pkgs = packages[package]
              subnav = []
 
              title = pkgs["Package"]
              if  version == @site_config["devel_version"]
                title += " (development version)"
              end
              if (version == @site_config["release_version"])
                bioc_version_str = "Release"
              elsif (version == @site_config["devel_version"])
                bioc_version_str = "Development"
              else
                bioc_version_str = nil
              end

              add_sym = {}
              for sym in link_list
                new_sym_name = "#{sym.to_s}_repo".to_sym
                new_sym = []
                for x in to_array(pkgs["#{sym}"])
                  new_sym.push(is_bioc_package?(x, version))
                end
                add_sym[new_sym_name] = new_sym
              end

              temp = pkgs.merge({
                  :rebase => true,
#                  :subnav => subnav,
                  :title => title,
                  :repo => repo,
                  :bioc_version_num => version,
                  :bioc_version_str => bioc_version_str
              })
              attributes = temp.merge(add_sym)

              identifier = Nanoc::Identifier.new(id, type: :legacy)

              item = new_item(" ", attributes, identifier)
              rep = Nanoc::Int::ItemRep.new(item, :unique_name)

              items.push item
            end
        else
          puts "ERROR: BiocViews DataSource: no JSON files; detail pages not built"
        end
      end
    end

    items
  end

  # if "packageName" is a bioC package, return its repo directory,
  # otherwise return false
  def is_bioc_package?(packageName, version)
    return false if packageName.nil?
    packageName = packageName.split(/[ |(]/).first
    @all_packages[version].each_pair do |k,v|
      return k if v.has_key?(packageName)
    end
    false
  end

end
