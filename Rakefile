require 'rubygems'
#require 'nanoc/tasks'
require 'yaml'
require 'fileutils'
require './lib/data_sources/gmane_list.rb'
require './scripts/search_indexer.rb'
require './scripts/parse_bioc_views.rb'
require './scripts/get_json.rb'
require './scripts/generate_build_shields.rb'
require './scripts/svn_shield_helper.rb'
require './scripts/get_post_tag_info.rb'
require 'open3'
require 'find'
require 'pathname'
require 'json'
require 'pp'
require 'httparty'
require 'nokogiri'
require 'descriptive_statistics'
require './lib/helpers.rb'
require 'csv'

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
task :compile => [:pre_compile,
  :real_compile, :post_compile]

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
  system("scp config.yaml webadmin@master.bioconductor.org:~")
  system("scp scripts/search_indexer.rb webadmin@master.bioconductor.org:~")
  system("scp scripts/get_links.rb webadmin@master.bioconductor.org:~")
  system(%Q(ssh webadmin@master.bioconductor.org "cd /home/webadmin && ./get_links.rb /extra/www/bioc > links.txt 2>not_found.txt"))
  system(%Q(ssh webadmin@master.bioconductor.org "cd /home/webadmin && /home/webadmin/do_index.rb"))
  #system("ssh webadmin@master.bioconductor.org chmod +x /home/webadmin/index.sh")
  system(%Q(ssh webadmin@master.bioconductor.org "/bin/sh /home/webadmin/index.sh"))
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
        obj = JSON.parse(
          File.read(file)#,
          #   :external_encoding => 'iso-8859-1',
          # )
        )

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
  puts "This might work:"
  puts "aws s3 cp --acl=public-read --recursive cloud_formation/output  s3://bioc-cloudformation-templates"
end

desc "Get Docbuilder Workflows"
task :get_workflows do
  site_config = YAML.load_file("./config.yaml")
  home = Dir.pwd
  #FileUtils.rm_rf "workflows_tmp"
  FileUtils.mkdir_p "workflows_tmp"
  dest_dir = "help/workflows"
  # f = File.open("content/#{dest_dir}.yaml", "w")
  # f.puts "---"
  # f.puts "title: Workflows"
  # f.close
  ##indexfile = File.open("content/#{dest_dir}.md", "w")
  ## You must have the appropriate private key in order for this to work.
  unless ENV["SKIP_WORKFLOW_RSYNC"] == "true"
    system(%Q(rsync --delete -ave "ssh -i #{ENV['HOME']}/.ssh/docbuilder" jenkins@docbuilder.bioconductor.org:~/repository/ workflows_tmp))
    system(%Q(chmod -R a+r workflows_tmp))
    Find.find('workflows_tmp') do |path|
      if (!File.directory?(path)) and (File.basename(path).start_with? '.')
        FileUtils.rm path
      end
    end
  end

  FileUtils.mkdir_p "assets/packages/#{site_config['release_version']}/workflows"
  system(%Q(rsync -av workflows_tmp/CRANrepo/#{site_config['release_version']}/ assets/packages/#{site_config['release_version']}/workflows/))
  Find.find("assets/packages/#{site_config['release_version']}/workflows/") do |path|
    if (!File.directory?(path)) and (File.basename(path).start_with? '.')
      FileUtils.rm path
    end
  end

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
# requires internet connection
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
        FileUtils.mkdir_p(File.join(tmpdir, version))
        if version == "release"
            machine = config["active_release_builders"]["linux"]
            biocversion = config["release_version"]
        else
            machine = config["active_devel_builders"]["linux"]
            biocversion = config["devel_version"]
        end
        unless (config["devel_repos"].include? repo.gsub("-", "/"))
          next
        end

        res = HTTParty.get("http://bioconductor.org/checkResults/#{version}/#{repo}-LATEST/STATUS_DB.txt")
        f = File.open(File.join(tmpdir, version, "STATUS_DB.txt"), "w")
        f.write(res)
        f.close

        #cmd = (%Q(rsync --delete --include="*/" --include="**/*.dcf" --exclude="*" -ave "ssh -o StrictHostKeyChecking=no -i #{ENV['HOME']}/.ssh/bioconductor.org.rsa" biocbuild@#{machine}:~/public_html/BBS/#{biocversion}/#{repo}/nodes #{tmpdir}/#{version}))
        #system(cmd)
    end
end

# should be run via cron every 15 minutes
desc "download build system databases and get build shields"
task :get_build_dbs do
  build_dbs_dir = File.join(%w(tmp build_dbs))
  FileUtils.mkdir_p(build_dbs_dir)
  %w(release devel).each do |version|
    %w(bioc data-experiment).each do |repo|
      url = "http://master.bioconductor.org/checkResults/#{version}/#{repo}-LATEST/STATUS_DB.txt"
      dest_file_name = File.join build_dbs_dir, "#{version}-#{repo}.dcf"
      dest_etag_name = dest_file_name.sub("dcf", "etag")
      etag = HTTParty.head(url).headers["etag"]
      if (!File.exists? dest_etag_name) or File.readlines(dest_etag_name).first != etag
        shield_dir = File.join("assets", "shields", "build", version, repo)
        FileUtils.mkdir_p shield_dir
        efh = File.open(dest_etag_name, 'w')
        efh.write etag
        efh.close
        body = HTTParty.get(url).to_s
        fh = File.open(dest_file_name, "w")
        fh.write(body)
        fh.close
        url2 = url.sub "STATUS_DB", 'meat-index'
        body2 = HTTParty.get(url2).to_s
        fh2 = File.open(dest_file_name.sub(/dcf$/, "meat-index.txt"), 'w')
        fh2.write(body2)
        fh2.close
        generate_build_shields(shield_dir, dest_file_name)

      end
    end
  end
end


desc "get and preprocess svn logs" # run this via cron every day (or more often?)
task :get_svn_logs do
  destdir = File.join('tmp', 'svnlogs')
  FileUtils.mkdir_p destdir
  shield_dir = File.join('assets', 'shields', 'commits')
  FileUtils.mkdir_p shield_dir
  today = Date.today
  now = DateTime.new(today.year, today.month, today.day)
  sixmonthsago = now
  months = [now]

  6.times do
    tmp = sixmonthsago.prev_month
    months << tmp
    sixmonthsago = tmp
  end
  months.reverse!
  ranges = []
  for i in 0..(months.length()-2)
    ranges.push months[i]..months[i+1]
  end
  h = {'bioc.log' => 'https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/',
    'data-experiment.log' => 'https://hedgehog.fhcrc.org/bioc-data/trunk/experiment/pkgs/'}
  h.each_pair do |destfile,url|
      full_shield_dir = File.join(shield_dir, destfile.split('.').first)
      FileUtils.mkdir_p full_shield_dir
      command = "svn log --username=readonly --password=readonly --non-interactive --no-auth-cache -v --xml -r {#{sixmonthsago.strftime("%Y-%m-%d")}}:{#{now.strftime("%Y-%m-%d")}} #{url}"

      stdin, stdout, stderr, wait_thr = Open3.popen3(command)

      fh = File.open(File.join(destdir, destfile), "w")
      while !stdout.eof
        line = stdout.readline
        fh.write line
      end
      fh.close
      stdout.close
      stderr.close
      stdin.close
      unless wait_thr.value.exitstatus == 0
        raise "svn log command failed!"
      end
      fh = File.open(File.join(destdir, destfile), 'r')
      doc = Nokogiri::XML(fh)
      fh.close
      logentries = doc.xpath("//log/logentry")
      dates = []
      for logentry in logentries
        slop = Nokogiri::Slop(logentry.to_xml)
        datestr = slop.logentry.date.children.first.text
        date = DateTime.iso8601(datestr)
        dates << date
      end

      mindate = dates.min.prev_day
      maxdate = dates.max.next_day
      ranges[0] =  mindate..ranges[0].last
      ranges[ranges.length-1] = ranges.last.first..maxdate
      res = Hash.new { |hash, key| hash[key] = {} }

      for logentry in logentries
        thisone = {}
        slop = Nokogiri::Slop(logentry.to_xml)
        datestr = slop.logentry.date.children.first.text
        date = DateTime.iso8601(datestr)
        paths = slop.logentry.paths.children.text.split("\n")
        for path in paths
          next if path.empty?
            segs = path.split '/'
            next if segs.length < 6
            thisone[segs[4]] = 1
        end
        # === operator does not work as expected, so....
        daterange = ranges.find{|i| i === date}
        daterange = ranges.find{|i| date > i.first and date < i.last} if daterange.nil?
        for item in thisone.keys
          unless ranges.include? daterange
            raise "there's a problem!"
          end
          if res[daterange].has_key? item
            res[daterange][item] += 1
          else
            res[daterange][item]= 1
          end
        end
      end
      packages = get_list_of_packages(bioc=(destfile=='bioc.log'))
      for package in packages
        hits = []
        for range in ranges
          if res[range].has_key? package
            hits << res[range][package]
          else
            hits << 0
          end
        end
        avg = hits.inject(0.0) { |sum, el| sum + el } / hits.size
        avg_s = sprintf("%0.2f", avg)
        puts "#{package}: #{avg_s}" # comment me out
        # parallelize this to make it faster?
        response = HTTParty.get("https://img.shields.io/badge/commits-#{avg_s}-1881c2.svg")
        if response.code == 200
          fh = File.open(File.join(full_shield_dir, "#{package}.svg"), 'w')
          fh.write(response.to_s)
          fh.close
        end
      end
  end


end

# run me with cron every day or so...
desc "process downloads data"
task :process_downloads_data do
  lines = HTTParty.get("http://s3.amazonaws.com/bioc-download-summaries/download_summary.csv")
  raw_data = []
  for line in lines
    raw_data << {num: line.first.to_i, pkg: line.last}
  end
  pkgs = []
  [true, false].each do |state|
    pkgs += get_list_of_packages(state)
    pkgs += get_annotation_package_list(state)
  end
  pkgs = pkgs.uniq
  avgs = {}
  for pkg in pkgs
    relevant = raw_data.find_all{|i| i[:pkg] == pkg}
    hits = relevant.map{|i| i[:num]}
    if hits.length < 6
      diff = 6 - hits.length
      diff.times {hits << 0}
    end
    avg = hits.inject(0.0) { |sum, el| sum + el } / hits.size
    avg = 0 if avg.nan?
    avgs[pkg] = avg
  end
  percentiles = {}
  data = avgs.values
  for pkg in pkgs
    percentiles[pkg] = data.percentile_rank(avgs[pkg])
  end
  srcdir = File.join('assets', 'images', 'shields', 'downloads')
  destdir = File.join('assets', 'shields', 'downloads')
  FileUtils.rm_rf destdir
  FileUtils.mkdir_p destdir
  percentiles.each_pair do |k, v|
    img = nil
    case v
    when 95..100
      img = 'top5.svg'
    when 80...95
      img = 'top20.svg'
    when 50...80
      img = 'top50.svg'
    else
      img = 'available.svg'
    end
    # puts "#{k}: #{img} (#{v})"
    FileUtils.cp(File.join(srcdir, img), File.join(destdir, "#{k}.svg"))
  end
end

# set this to run in crontab
desc "get info about post tags"
task :get_post_tag_info do
  get_post_tag_info()
end

desc "do push tasks"
task :push => [:copy_assets, :deploy_staging, :deploy_production]

# run me in crontab daily
desc "get years-in-bioc shields"
task :get_years_in_bioc_shields do
  sconfig = YAML.load_file("./config.yaml")
  sconfig[:release_dates] = sconfig["release_dates"]
  sconfig[:release_dates] = sconfig[:release_dates].inject({}){|memo,(k,v)| memo[k.to_sym] = v; memo}
  puts "Getting manifest information..."
  url = "https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/"
  branch_url = "https://hedgehog.fhcrc.org/bioconductor/branches/RELEASE_"
  auth = {:username => "readonly", :password => "readonly"}
  page = HTTParty.get("https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/",
    :basic_auth => auth, :verify => false).body
  manifests = {}
  for line in page.split("\n")
    if line =~ /bioc_([0123456789.]+)\.manifest/
      ver = $1
      manifests[ver] = []
      if ver == sconfig["devel_version"]
        manifest_url = "#{url}bioc_#{ver}.manifest"
      else
        manifest_url = "#{branch_url}#{ver.gsub(".", "_")}/madman/Rpacks/bioc_#{ver}.manifest"
      end
      manifest = HTTParty.get(manifest_url, :basic_auth => auth, :verify=>false).body
      for l in manifest.split("\n")
        next unless l =~ /^Package: /
        manifests[ver].push l.chomp.sub(/^Package: /, "").strip
      end
    end
  end
  sconfig[:manifests] = manifests

  sconfig[:manifest_keys] = manifests.keys.sort do |a,b|
    amaj, amin = a.split(".")
    bmaj, bmin = b.split(".")
    amaj = Integer(amaj)
    amin = Integer(amin)
    bmaj = Integer(bmaj)
    bmin = Integer(bmin)
    if amaj == bmaj
      amin <=> bmin
    else
      amaj <=> bmaj
    end
  end
  pkgs = get_list_of_packages()
  for pkg in pkgs
    get_year_shield(pkg, true, sconfig)
  end
end

# run me in crontab
desc "get test coverage shields"
task :get_coverage_shields do
  config = YAML.load_file("./config.yaml")

  branches = ["release-#{config['release_version']}", "devel"]

  branches.each do |branch|
    dirname = branch
    dirname = "release" if branch =~ /^release/
    codecov_branch = branch
    codecov_branch = "master" if branch == "devel"
    shield_dir = "assets/shields/coverage/#{dirname}"

    packages = get_list_of_packages(true, branch=="release")
    unless File.exists? shield_dir
      FileUtils.mkdir_p shield_dir
    end
    for package in packages
      url = "https://codecov.io/github/Bioconductor-mirror/#{package}/coverage.svg?branch=#{codecov_branch}"
      cov = HTTParty.head(url).headers["x-coverage"]
      cov_color = coverage_color(cov)
      cov += "%" unless cov == "unknown"
      puts "Downloading test-coverage shield for #{package} in #{dirname}..."
      resp = HTTParty.get("https://img.shields.io/badge/test_coverage-#{URI::encode(cov)}-#{cov_color}.svg")
      if resp.code == 200
        shield = File.join(shield_dir, "#{package}.svg")
        fh = File.open(shield, "w")
        fh.write(resp.to_s)
        fh.close
      end
    end

  end
end

desc "get all shields"
task :get_all_shields => [:get_build_dbs, :get_svn_logs,
  :process_downloads_data, :get_post_tag_info,
  :get_years_in_bioc_shields, :get_coverage_shields, :copy_assets]

# should be run every time mirror info in config.yaml changes
# that's hard to remember do to, so
# make sure this is run via crontab every hour
desc "extract mirror information to csv file"
task :mirror_csv do
    config = YAML.load_file("./config.yaml")
    CSV.open(File.join("assets", "BioC_mirrors.csv"), "w") do |csv|
      csv << ["Name","Country","City","URL","Host","Maintainer","OK","CountryCode"]
      for mirror_outer in config['mirrors']
        country = mirror_outer.keys.first
        country_mirrors = mirror_outer.values
        for mirrors in country_mirrors
          for mirror in mirrors
            if mirror['contact'].is_a? Array
              mirror['contact'] = mirror['contact'].first
              mirror['contact_email'] = mirror['contact_email'].first
            end
            maintainer = "#{mirror['contact']} <#{mirror['contact_email']}>".sub("@", " # ")
            data = ["#{country} (#{mirror['city']})", country, mirror['city'], mirror['mirror_url'],
              mirror['institution'], maintainer,
              check_mirror_url(mirror['mirror_url']), mirror['country_code']]
            csv << data
            if mirror.has_key? "https_mirror_url"
              data[3] = mirror['https_mirror_url']
              data[0] = data[0] + " [https]"
              data[6] = check_mirror_url(mirror['https_mirror_url'])
              csv << data
            end
          end
        end
      end
    end
end
