require 'rubygems'
require 'nanoc3/tasks'
require 'yaml'
require 'fileutils'
require './lib/data_sources/gmane_list.rb'
require './scripts/search_indexer.rb'
require './scripts/parse_bioc_views.rb'
require './scripts/get_json.rb'
require 'open3'
require 'find'
require 'pathname'
require 'json'
require 'pp'
require 'httparty'

include Open3


@clear_search_index_commands = [
  "curl -s http://localhost:8983/solr/update --data-binary '<delete><query>*:*</query></delete>' -H 'Content-type:text/xml; charset=utf-8'",
  "curl -s http://localhost:8983/solr/update --data-binary '<optimize/>' -H 'Content-type:text/xml; charset=utf-8'",
  "curl -s  http://localhost:8983/solr/update --data-binary '<commit/>' -H 'Content-type:text/xml; charset=utf-8'"
  ]



desc "write version info to doc root for javascript to find"
task :write_version_info do
  site_config = YAML.load_file("./config.yaml")
  js = %Q(var develVersion = "#{site_config["devel_version"]}";\nvar releaseVersion="#{site_config["release_version"]}";\n)
  js += %Q(var versions = [")
  js += site_config["versions"].join(%Q(",")) + %Q("];\n)
  f = File.open("#{site_config["output_dir"]}/js/versions.js", "w")
  f.puts js
  f.close
end

desc "copy assets to output directory"
task :copy_assets do
  site_config = YAML.load_file("./config.yaml")
  output_dir = site_config["output_dir"]
  system "rsync -gprt --partial --exclude='.svn' assets/ #{output_dir}"
end

desc "Run nanoc compile"
task :compile => [:pre_compile, :real_compile, :post_compile]

desc "Pre-compilation tasks"
task :pre_compile do
  FileUtils.mkdir_p "content/packages"
  site_config = YAML.load_file("./config.yaml")
  for version in site_config["versions"]
    destdir = "content/packages/#{version}"
    FileUtils.mkdir_p destdir
    puts "copying bioc_views.html to #{destdir}/BiocViews.html"
    FileUtils.rm_f "#{destdir}/BiocViews.html"
    FileUtils.rm_f "output/packages/#{version}/BiocViews.html"
    FileUtils.mkdir_p "output/packages/#{version}/BiocViews"
    unless(ENV.has_key?("QUICK_NANOC_COMPILE") && ENV["QUICK_NANOC_COMPILE"] == "true")
      FileUtils.cp "assets/help/bioc-views.html", "#{destdir}/BiocViews.html", {:preserve => false}
      FileUtils.cp "assets/help/bioc-views.yaml", "#{destdir}/BiocViews.yaml", {:preserve => false}
    end
  end
  
end

task :real_compile do
  fail "compilation failed" unless system "nanoc co"
end

task :post_compile do
  puts "running post-compilation tasks..."
  site_config = YAML.load_file("./config.yaml")
  for entry in Dir.entries("#{site_config["output_dir"]}/packages")
    next if entry =~ /^\./
    next unless entry =~ /^[0-9]/
    entry = "#{site_config["output_dir"]}/packages/#{entry}"
    next unless Kernel.test(?d, entry)
    puts "copying index.html to #{entry}..."
    FileUtils.cp("assets/extra/index.html", 
      "#{entry}/index.html")
  end
  cwd = FileUtils.pwd
  FileUtils.cd "#{site_config["output_dir"]}/packages"

  FileUtils.rm_f "devel"
  FileUtils.rm_f"release"
  begin
    FileUtils.ln_s "#{site_config["release_version"]}", "release"
    FileUtils.ln_s "#{site_config["devel_version"]}", "devel"
  rescue NotImplementedError => ex
    puts "skipping symlink as ln_s is not supported on this platform"
  end
  puts "Generated symlinks for release and devel"
  FileUtils.cd cwd
end

desc "Nuke output directory !! uses rm -rf !!"
task :real_clean do
  site_config = YAML.load_file("./config.yaml")
  output_dir = site_config["output_dir"]
  FileUtils.rm_rf(output_dir)
  FileUtils.mkdir_p(output_dir)
end

desc "copy config.yaml to assets"
task :copy_config do
    FileUtils.copy("config.yaml", "assets/")
end

desc "Build the bioconductor.org site (default)"
task :build => [ :compile, :copy_config, :copy_assets, 
    :write_version_info, :write_version_number ]

task :default => :build

desc "deploy (sync) to staging (run on merlot2)"
task :deploy_staging do
  dst = '/loc/www/bioconductor-test.fhcrc.org'
  site_config = YAML.load_file("./config.yaml")
  output_dir = site_config["output_dir"]
  system "rsync -av --links --partial --partial-dir=.rsync-partial --exclude='.svn' #{output_dir}/ #{dst}"
  chmod_cmd = "chmod -R a+r /loc/www/bioconductor-test.fhcrc.org/packages/json"
  system chmod_cmd
end

desc "deploy (sync) to production"
task :deploy_production do
  site_config = YAML.load_file("./config.yaml")
  src = '/loc/www/bioconductor-test.fhcrc.org'
  dst = site_config["production_deploy_root"]
  system "rsync -av --links --partial --partial-dir=.rsync-partial --exclude='.svn' #{src}/ #{dst}/"
end


desc "Clear search index (and local cache)"
task :clear_search_index do
  for command in @clear_search_index_commands
    puts command
    system command
  end
  dir = (`hostname`.chomp == 'merlot2') ? "/home/biocadmin" : "."
  
  FileUtils.rm_f "#{dir}/search_indexer_cache.yaml"
end


#todo fix this
desc "Clear search index (and local cache) on production. This will cause searches to fail until indexing is re-done!"
task :clear_search_index_production do
  for command in @clear_search_index_commands
    puts %Q(ssh webadmin@bioconductor.org "#{command}")
    system %Q(ssh webadmin@bioconductor.org "#{command}")
  end
  puts "ssh webadmin@bioconductor.org rm -f /home/webadmin/search_indexer_cache.yaml"
  system "ssh webadmin@bioconductor.org rm -f /home/webadmin/search_indexer_cache.yaml"
end

desc "Re-index the site for the search engine"
task :search_index do
  if (SearchIndexer.is_solr_running?)
    hostname = `hostname`.chomp
    args = [] # directory to index, location of cache file, location of output shell script, url of site
    if (hostname =~ /^dhcp/i) # todo - don't test this, instead see if nanoc is installed
      args = ['./output', './', 'scripts', 'http://localhost:3000'] 
    elsif (hostname == 'merlot2')
      args = ['/loc/www/bioconductor-test.fhcrc.org', '/home/biocadmin', '/home/biocadmin', 'http://bioconductor-test.fhcrc.org']
    end
    pwd = FileUtils.pwd
    si = SearchIndexer.new(args)
    FileUtils.cd pwd # just in case
    cmd = "#{args[2]}/index.sh"

    chmod_cmd = "chmod u+x #{cmd}"
    system chmod_cmd
    # todo fork the indexer to a background task because it could take a while in some cases

    stdin, stdout, stderr = Open3.popen3("sh ./scripts/index.sh")
    puts "running indexer, stdout = "
    puts stdout.readlines
    puts "stderr = "
    puts stderr.readlines
  else
    puts "solr is not running, not re-indexing site."
  end
end

desc "Re-run search indexing on production"
task :index_production do
  system("scp config.yaml webadmin@bioconductor.org:~")
  system("scp scripts/search_indexer.rb webadmin@bioconductor.org:~")
  system("scp scripts/get_links.rb webadmin@bioconductor.org:~")
  system(%Q(ssh webadmin@bioconductor.org "cd /home/webadmin && ./get_links.rb /extra/www/bioc > links.txt 2>not_found.txt"))
  system(%Q(ssh webadmin@bioconductor.org "cd /home/webadmin && /home/webadmin/do_index.rb"))
  #system("ssh webadmin@bioconductor.org chmod +x /home/webadmin/index.sh")
  system(%Q(ssh webadmin@bioconductor.org "/bin/sh /home/webadmin/index.sh"))
end

desc "Re-run search indexing cran package home pages on production"
task :index_cran_production do
  system("scp scripts/cran_search_indexer.rb webadmin@bioconductor.org:~")
  system("ssh webadmin@bioconductor.org /home/webadmin/do_index_cran.rb")
  system("ssh webadmin@bioconductor.org /bin/sh /home/webadmin/index_cran.sh")
end

desc "Runs nanoc's dev server on localhost:3000"
task :devserver => [:build] do
  system "nanoc aco"
end


## If this doesn't work, do:
##   rake prepare_json && rake json2js
task :get_json => [:prepare_json, :json2js]

task :json2js do
  # convert json to javascript expression suitable for inclusion via a script tag
  Find.find "assets/packages/json" do |file|
    next unless file =~ /\.json$/
    dir = Pathname.new(file).basename.to_s
    fn = file.split("/").last
    jsfile = file.sub(/\.json$/, ".js")
    if fn == "tree.json"
      f = File.open(file)
      lines = f.readlines
      f.close
      json = lines.join
      obj = JSON.parse json
      hsh = {"data" => obj}
      json = hsh.to_json
      s = "var dataTree = #{json};"
      outf = File.open(jsfile, "w")
      outf.print s
      outf.close
    elsif fn == "packages.json"
      var = nil
      if file =~ /\/bioc\//
        var = "bioc_packages"
      elsif file =~ /\/data\/annotation/
        var = "data_annotation_packages"
      elsif file =~ /\/data\/experiment/
        var = "data_experiment_packages"
      end
      unless var.nil?
        f = File.open(file)
        lines = f.readlines
        json = lines.join.strip
        obj = JSON.parse json
        
        ret = []
        obj.values.each do |v|
          line = []
          line.push v['Package']
          line.push v['Maintainer']
          line.push v['Title']
          ret.push line
        end
        
        holder = {"content" => ret}
        outf = File.open(jsfile, "w")
        outf.print("var #{var} = #{holder.to_json};")
        outf.close
      end
    end
  end
end


desc "Get JSON files required for BiocViews pages"
task :prepare_json do
  json_dir = "assets/packages/json"
  FileUtils.mkdir_p json_dir
  site_config = YAML.load_file("./config.yaml")
  versions = site_config["versions"]
  devel_version = site_config["devel_version"]
  devel_repos = site_config["devel_repos"]
  version_str = '"' + versions.join('","') + '"'
  devel_repos_str = '"' + devel_repos.join('","') + '"'


  for version in versions
    if version == devel_version
      repos = devel_repos
    else
      repos = ["data/annotation", "data/experiment", "bioc"]
    end
   
    gj = GetJson.new(version, "assets/packages/json")
  end
end

desc "Create CloudFormation templates"
task :generate_cf_templates do
  FileUtils.mkdir_p "cloud_formation/output"
  config = YAML.load_file("./config.yaml")
  dir = Dir.new("cloud_formation")
  for file in dir
    next unless file =~ /\.json$/
    f = File.open("cloud_formation/#{file}")
    lines = f.readlines
    str = lines.join()
    repl = str.gsub(/<%=[^%]*%>/) do |s|
      s = s.gsub("<%=", "").gsub("%>", "").strip()
      eval(s)
    end
    outfile = File.open("cloud_formation/output/#{file}", "w")
    outfile.write(repl)
    outfile.close
  end
  puts "Don't forget to copy templates (in cloudformation/output)"
  puts "to S3 and mark them as public!"
end

desc "Get Docbuilder Workflows"
task :get_workflows do
  site_config = YAML.load_file("./config.yaml")
  home = Dir.pwd
  #FileUtils.rm_rf "workflows_tmp"
  FileUtils.mkdir_p "workflows_tmp"
  dest_dir = "help/workflows" # eventually this will change to help/workflows
  f = File.open("content/#{dest_dir}.yaml", "w")
  f.puts "---"
  f.puts "title: Workflows"
  f.close
  ##indexfile = File.open("content/#{dest_dir}.md", "w")
  ## You must have the appropriate private key in order for this to work.
  unless ENV["SKIP_WORKFLOW_RSYNC"] == "true"
    system(%Q(rsync --delete -ave "ssh -i #{ENV['HOME']}/.ssh/docbuilder" jenkins@docbuilder.bioconductor.org:~/repository/ workflows_tmp))
    system(%Q(chmod -R a+r workflows_tmp))
  end    

  FileUtils.mkdir_p "assets/packages/#{site_config['release_version']}/workflows"
  system(%Q(rsync -av workflows_tmp/CRANrepo/#{site_config['release_version']}/ assets/packages/#{site_config['release_version']}/workflows/))

  auth = {:username => "readonly", :password => "readonly"}
  json = HTTParty.get("https://hedgehog.fhcrc.org/bioconductor/trunk/madman/workflows/manifest.json",
    :basic_auth => auth, :verify => false).body
  manifest = JSON.parse(json)  

  dir = Dir.new("workflows_tmp")
  for entry in dir.entries
    next if entry =~ /^\./
    next if ["CRANrepo", "manifest.txt"].include? entry

    unless ENV["IGNORE_WORKFLOW_MANIFEST"] == "true"
      unless manifest.keys.include? entry and manifest[entry].include? "web"
        next
      end
    else
      puts "ignoring workflow manifest"
    end

    fullpath = "workflows_tmp/#{entry}"
    if test ?d, fullpath # if it exists and is a directory
      # indexfile.puts
      # indexfile.puts "## Workflows in #{entry.capitalize}:"
      # indexfile.puts
      asset_dir = "assets/#{dest_dir}/#{entry}"
      md_dir = "content/#{dest_dir}/#{entry}"
      FileUtils.rm_rf asset_dir
      FileUtils.rm_rf md_dir
      FileUtils.mkdir_p asset_dir
      dir2 = Dir.new(fullpath)
      vignettes = dir2.entries.find_all {|i| i =~ /\.md$/i}
      multivig = false
      multivig = true if vignettes.length() > 1
      FileUtils.mkdir_p md_dir if multivig
      for vignette in vignettes
        vigname = vignette.gsub(/\.md$/i, "")
        pdf = "#{vigname}.pdf"
        FileUtils.mkdir_p "#{asset_dir}/#{vigname}" if multivig
        if multivig
          dest = "#{asset_dir}/#{vigname}"
          md_dest = md_dir
        else
          dest = "#{asset_dir}/#{entry}.R"
          md_dest = "content/#{dest_dir}"
        end
        yaml = YAML::load(File.open("#{fullpath}/#{vigname}.yaml"))
        FileUtils.mv "#{fullpath}/#{vigname}.R", dest
        ["md", "yaml"].each do |suffix|
          if multivig
            FileUtils.mv "#{fullpath}/#{vigname}.#{suffix}", md_dest
          else
            FileUtils.mv "#{fullpath}/#{vigname}.#{suffix}", \
              "#{md_dest}/#{entry}.#{suffix}"
          end
        end

        # if multivig
        #   indexfile.puts "* [#{yaml['title']}](/#{dest_dir}/#{entry}/#{vigname})"
        # else
        #   indexfile.puts "* [#{yaml['title']}](/#{dest_dir}/#{entry})"
        # end

        ["tar.gz", "tgz", "zip"].each do |suffix|
          file = Dir.glob("#{fullpath}/*.#{suffix}").first
          next if file.nil?
          FileUtils.mv "#{file}", asset_dir
        end
        FileUtils.rm_f "#{fullpath}/#{vigname}.Rmd"
      end
      dir = Dir.new(fullpath)
      for entry in dir.entries
        next if entry =~ /^\./
        for vig in vignettes
          vigname = vig.gsub(/\.md$/i, "")
          if multivig
            dest = "#{asset_dir}/#{vigname}"
          else
            dest = asset_dir
          end
          FileUtils.cp_r "#{fullpath}/#{entry}", dest
        end
      end
    end
  end
  #indexfile.close
  #FileUtils.rm_rf "workflows_tmp"
end

desc "write version number to endpoint"
task :write_version_number do
    config = YAML.load_file("./config.yaml")
    f = File.open("assets/bioc-version", "w")
    f.print(config["release_version"])
end

task :my_task, :arg1 do |t, args|
  hargs = args.to_hash
  puts "hargs were:"
  pp hargs
  puts "is it nil? #{hargs.nil?}"
  puts hargs.class
  puts "is it empty? #{hargs.empty?}"
  puts "keys:"
  pp hargs.keys()
end

desc "Get build result summaries to build RSS feeds (arg: bioc or data-experiment)" 
# requires connection to internal hutch network 
task :get_build_result_dcfs, :repo do |t, args|
    hargs = args.to_hash
    if hargs.empty? or !hargs.has_key? :repo
      repo = "bioc"
    else
      repo = hargs[:repo]
    end
    unless ["bioc", "data-experiment"].include? repo
      puts "Argument must be either 'bioc' or 'data-experiment'."
      next
    end
    config = YAML.load_file("./config.yaml")
    tmpdir = repo=="bioc" ? "tmp/build_dcfs" : "tmp/data_build_dcfs"
    FileUtils.mkdir_p tmpdir
    ary = []
    for version in ["release", "devel"]
        if version == "release"
            machine = config["active_release_builders"]["linux"]
            biocversion = config["release_version"]
        else
            machine = config["active_devel_builders"]["linux"]
            biocversion = config["devel_version"]
        end
        unless (config["devel_repos"].include? repo.gsub("/", "-"))
          next
        end
        cmd = (%Q(rsync --delete --include="*/" --include="**/*.dcf" --exclude="*" -ave "ssh -o StrictHostKeyChecking=no -i #{ENV['HOME']}/.ssh/bioconductor.org.rsa" biocbuild@#{machine}:~/public_html/BBS/#{biocversion}/#{repo}/nodes #{tmpdir}/#{version}))
        system(cmd)
    end
end