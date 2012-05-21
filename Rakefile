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

desc "Run nanoc3 compile"
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
  system "nanoc3 co"
end

task :post_compile do
  puts "running post-compilation tasks..."
  site_config = YAML.load_file("./config.yaml")
  src = "#{site_config["output_dir"]}/packages/#{site_config["release_version"]}/BiocViews.html"
  other_versions = site_config["versions"] - [site_config["release_version"]]
  cwd = FileUtils.pwd
  FileUtils.cd "#{site_config["output_dir"]}/packages"

  FileUtils.rm_f "devel"
  FileUtils.rm_f"release"
  FileUtils.ln_s "#{site_config["release_version"]}", "release"
  FileUtils.ln_s "#{site_config["devel_version"]}", "devel"
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

desc "Build the bioconductor.org site (default)"
task :build => [ :compile, :copy_assets, :write_version_info ]

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
  system(%Q(ssh webadmin@bioconductor.org "cd /home/webadmin && ./get_links.rb > links.txt 2>not_found.txt"))
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
  system "nanoc3 aco"
end


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
        outf = File.open(jsfile, "w")
        outf.print("var #{var} = #{json};")
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
  config = YAML.load_file("./config.yaml")
  dir = Dir.new("cloud_formation")
  for file in dir
    next unless file =~ /_with_symbols.json$/
    outname = file.gsub(/with_symbols/, "")
    f = File.open("cloud_formation/#{file}")
    lines = f.readlines
    str = lines.join()
    repl = str.gsub(/<%=[^%]*%>/) do |s|
      s = s.gsub("<%=", "").gsub("%>", "").strip()
      eval(s)
    end
  end
  outfile = File.open("cloud_formation/#{outname}", "w")
  outfile.write(repl)
  outfile.close
  puts "Don't forget to copy templates to S3 and mark them as public!"
end