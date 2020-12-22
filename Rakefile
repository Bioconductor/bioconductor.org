require 'rubygems'
require 'nanoc'
require 'yaml'
require 'fileutils'
require './lib/data_sources/gmane_list.rb'
require './scripts/search_indexer.rb'
require './scripts/parse_bioc_views.rb'
require './scripts/get_json.rb'
require './scripts/badge_generation.rb'
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
require 'date'
require 'dcf'

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
  system "rsync -gprt --partial --exclude='.git' assets/ #{output_dir}"
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
  fail "compilation failed" unless system "bundle exec nanoc co"
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
  FileUtils.rm_f "release"
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

desc "deploy (sync) to staging server root on staging"
task :deploy_staging do
  dst = '/loc/www/bioconductor-test.fhcrc.org'
  site_config = YAML.load_file("./config.yaml")
  output_dir = site_config["output_dir"]
  system "rsync -av --links --partial --partial-dir=.rsync-partial --exclude='.git' #{output_dir}/ #{dst}"
  chmod_cmd = "chmod -R a+r /loc/www/bioconductor-test.fhcrc.org/packages/json"
  system chmod_cmd
end

desc "deploy (sync) to production"
task :deploy_production do
  site_config = YAML.load_file("./config.yaml")
  src = '/loc/www/bioconductor-test.fhcrc.org'
  dst = site_config["production_deploy_root"]
  system "rsync -av --links --partial --partial-dir=.rsync-partial --exclude='.git' #{src}/ #{dst}/"
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
  system "bundle exec nanoc aco"
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
      elsif file =~ /\/workflows\//
	var = "workflow_packages"
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
	  line.push v['Rank']
	  line.push v['dependencyCount']
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

  for version in versions
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

desc "write version number to endpoint"
task :write_version_number do
    config = YAML.load_file("./config.yaml")
    f = File.open("assets/bioc-version", "w")
    f.print(config["release_version"])
    f = File.open("assets/bioc-devel-version", "w")
    f.print(config["devel_version"])
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

desc "do push tasks"
task :push => [:copy_assets, :deploy_staging, :deploy_production]

desc "pull manifests"
task :pull_manifests do
  system("git -C ../manifest/ pull --all")
end

# originally run every 15 minutes - but it will only change with build report so
# now once a day
desc "download build system databases and get build shields"
task :get_build_dbs do
  build_dbs_dir = File.join(%w(tmp build_dbs))
  FileUtils.mkdir_p(build_dbs_dir)
  %w(release devel).each do |version|
    %w(bioc data-experiment workflows).each do |repo|
      puts "Working On: #{version} #{repo}"
      url = "http://master.bioconductor.org/checkResults/#{version}/#{repo}-LATEST/STATUS_DB.txt"
      dest_file_name = File.join build_dbs_dir, "#{version}-#{repo}.dcf"
      dest_etag_name = dest_file_name.sub("dcf", "etag")
      etag = HTTParty.head(url).headers["etag"]
      urlcode = HTTParty.head(url).response.code
      if ((!File.exists? dest_etag_name) or File.readlines(dest_etag_name).first != etag) and urlcode == "200"
	shield_dir = File.join("assets", "shields", "build", version, repo)
	FileUtils.mkdir_p shield_dir
	efh = File.open(dest_etag_name, 'w')
	efh.write etag
	efh.close
	body = HTTParty.get(url).to_s
	fh = File.open(dest_file_name, "w")
	fh.write(body)
	fh.close
	url2 = url.sub "STATUS_DB.txt", 'meat-index.dcf'
	body2 = HTTParty.get(url2).to_s
	fh2 = File.open(dest_file_name.sub(/dcf$/, "meat-index.txt"), 'w')
	fh2.write(body2)
	fh2.close
	puts shield_dir
	puts dest_file_name
	generate_build_shields(shield_dir, dest_file_name)

      end
    end
  end
end

# off download stats which are generated tue/fri so run wed/sat
desc "process downloads data"
task :process_downloads_data do

  destdir = File.join('assets', 'shields', 'downloads', 'devel')
  FileUtils.rm_rf destdir
  FileUtils.mkdir_p destdir

  repos = ["bioc", "annotation", "experiment", "workflows"]
  repos.each do |repo|
    puts "GENERATING BADGES: " + repo + ", devel"
    downloadBadge(repo, destdir, false)
  end

  destdir = File.join('assets', 'shields', 'downloads', 'release')
  FileUtils.rm_rf destdir
  FileUtils.mkdir_p destdir

  repos.each do |repo|
    puts "GENERATING BADGES: " + repo + ", release"
    downloadBadge(repo, destdir, true)
  end

end


# set this to run in crontab
desc "get pkg availability info"
task :get_availability_shields  do
  site_config = YAML.load_file("./config.yaml")
  for reldev in ['release', 'devel']
    puts "Working on #{reldev}"
    numeric_version = (reldev == 'release') ? site_config['release_version'] : site_config['devel_version']
    meat_index_file = File.join("tmp", "build_dbs","#{reldev}-bioc.meat-index.txt")
    json_file = File.join("assets", "packages", "json", numeric_version, "bioc", "packages.json")

    indexList=[]
    if File.exists? meat_index_file

      mitxt = File.readlines(meat_index_file).join
      meat_index = Dcf.parse(mitxt)
      json_obj = JSON.parse(File.read(json_file))

      if (not meat_index.nil?)
	for item in meat_index
	  get_availability(item, numeric_version)
	  indexList.push(item['Package'])
	end
	unknown = json_obj.keys.sort - indexList.sort
	if unknown.length != 0
	  for item in unknown
	    availabilityBadge(item, "unknown-build", numeric_version)
	  end
	end
      else
	puts "ERROR meat_index nil"
      end
    else
      puts "ERROR meat_index doesn't exist"
    end
  end
end


desc "get info about support site activity"
task :get_supportsite_info_shield do
  require './scripts/get_support_tag_info.rb'
  create_package_list()
  get_support_tag_info()
end

# shouldn't be run daily - will update minimally
desc "get years-in-bioc shields"
task :get_years_in_bioc_shields do
  sconfig = YAML.load_file("./config.yaml")
  sconfig[:release_dates] = sconfig["release_dates"]
  sconfig[:release_dates] = sconfig[:release_dates].inject({}){|memo,(k,v)| memo[k.to_sym] = v; memo}
  rel_ver = sconfig["release_version"]
  dev_ver = sconfig["devel_version"]
  all_ver = sconfig["release_dates"].keys
  all_ver.push dev_ver
  man_path = "../manifest/"

  manifests = {}

  all_ver.each { |v|
    ver = v.gsub(/\./,"_")
    manifests[v] = []
    if v == dev_ver
      system("git -C #{man_path} checkout master")
    else
      system("git -C #{man_path} checkout RELEASE_#{ver}")
    end
    if $? == 0
      soft="#{man_path}/software.txt"
      File.open(soft, "r") do |f|
	f.each_line do |line|
	  next unless line =~ /^Package: /
	  manifests[v].push line.chomp.sub(/^Package: /, "").strip
	end
      end
    end
  }
  sconfig[:manifests] = manifests
  system("git -C #{man_path} checkout master")

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

# run crontab daily - only off VIEWS file creation so will not update
# more frequent than build report
desc "get last commit date shields"
task :get_last_commit_date_shields do

  site_config = YAML.load_file("./config.yaml")
  for reldev in ['release', 'devel']
     for repo in ['bioc', 'data-experiment', 'workflows']
	puts "writing badges for #{reldev} #{repo}"
	numeric_version = (reldev == 'release') ? site_config['release_version'] : site_config['devel_version']
	unsupported_platforms = {}
	json_file = (repo == 'data-experiment') ? File.join("assets", "packages", "json", numeric_version, "data", "experiment","packages.json") : File.join("assets", "packages", "json", numeric_version, repo, "packages.json")
	json_obj = JSON.parse(File.read(json_file))
	release_date = Date.strptime(site_config['release_dates'][site_config['release_version'].to_s], "%m/%d/%Y")

	day = []
	week = []
	month = []
	three = []
	release = []
	none = []
	unknown = []

	for k in json_obj.keys
	    date = json_obj[k]["git_last_commit_date"]
	    if (date.nil?)
		if (reldev == 'release')
		   none.push k
		else
		   unknown.push k
		end
	    else
		date = Date.parse(date)
		if (date < release_date)
		    none.push k
		elsif (date < (Date.today - 90) && date >= release_date)
		    release.push k
		elsif (date < (Date.today - 30) && date >= (Date.today - 90))
		    three.push k
		elsif (date < (Date.today - 7) && date >= (Date.today - 30))
		    month.push k
		elsif (date < (Date.today - 1) && date >= (Date.today - 7))
		    week.push k
		elsif (date >= (Date.today - 1))
		    day.push k
		else
		    unknown.push k
		end
	     end
	end

	# make badges
	shield_dir = File.join("assets", "shields", "lastcommit", reldev, repo)
	FileUtils.mkdir_p shield_dir
	if (not day.empty?)
	    puts "DAY"
	    img = File.join("assets","images","shields","lastcommit", "Day.svg")
	    day.each{|pkg|
		FileUtils.cp img, File.join(shield_dir, (pkg+".svg"))
		puts "    #{pkg}"
	    }
	end
	if (not week.empty?)
	    puts "WEEK"
	    img = File.join("assets","images","shields","lastcommit", "Week.svg")
	    week.each{|pkg|
		FileUtils.cp img, File.join(shield_dir, (pkg+".svg"))
		puts "    #{pkg}"
	    }
	end
	if (not month.empty?)
	    puts "MONTH"
	    img = File.join("assets","images","shields","lastcommit", "Month.svg")
	    month.each{|pkg|
		FileUtils.cp img, File.join(shield_dir, (pkg+".svg"))
		puts "    #{pkg}"
	    }
	end
	if (not three.empty?)
	    "THREE MONTHS"
	    img = File.join("assets","images","shields","lastcommit", "ThreeMonths.svg")
	    three.each{|pkg|
		FileUtils.cp img, File.join(shield_dir, (pkg+".svg"))
		puts "    #{pkg}"
	    }
	end
	if (not release.empty?)
	    "SINCE RELEASE"
	    img = File.join("assets","images","shields","lastcommit", "LastRelease.svg")
	    release.each{|pkg|
		FileUtils.cp img, File.join(shield_dir, (pkg+".svg"))
		puts "    #{pkg}"
	    }
	end
	if (not none.empty?)
	    "NONE"
	    img = File.join("assets","images","shields","lastcommit", "None.svg")
	    none.each{|pkg|
		FileUtils.cp img, File.join(shield_dir, (pkg+".svg"))
		puts "    #{pkg}"
	    }
	end
	if (not unknown.empty?)
	    "UNKNOWN"
	    img = File.join("assets","images","shields","lastcommit", "Unknown.svg")
	    unknown.each{|pkg|
		FileUtils.cp img, File.join(shield_dir, (pkg+".svg"))
		puts "    #{pkg}"
	    }
	end

     end # repo
  end  # reldev
end # task


# this shouldn't update very often, run crontab once a month
desc "process dependency badge"
task :process_dependency_badge do

  destdir = File.join('assets', 'shields', 'dependencies', 'devel')
  FileUtils.rm_rf destdir
  FileUtils.mkdir_p destdir

  repos = ["bioc", "experiment", "workflows"]
  repos.each do |repo|
    puts "GENERATING BADGES: " + repo + ", devel"
    dependencyBadge(repo, destdir, false)
  end

  destdir = File.join('assets', 'shields', 'dependencies', 'release')
  FileUtils.rm_rf destdir
  FileUtils.mkdir_p destdir

  repos.each do |repo|
    puts "GENERATING BADGES: " + repo + ", release"
    dependencyBadge(repo, destdir, true)
  end

end

desc "get all shields"
task :get_all_shields => [:get_build_dbs,
  :process_downloads_data, :get_post_tag_info,
  :get_years_in_bioc_shields, :copy_assets,
  :get_availability_shields, :get_last_commit_date_shields,
  :process_dependency_badge]

# should be run every time mirror info in config.yaml changes
# that's hard to remember do to, so
# make sure this is run via crontab every hour
desc "extract mirror information to csv file"
task :mirror_csv do
    config = YAML.load_file("./config.yaml")
    CSV.open(File.join("assets", "BioC_mirrors.csv"), "w") do |csv|
      csv << ["Name","Country","City","URL","Host","Maintainer","OK",
	      "CountryCode","Comment"]
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
	    # As of BioC 3.5 https is required for all sites
	    # First row is https (CRAN request)
	    data = ["#{country} (#{mirror['city']}) [https]",
	      country,
	      mirror['city'],
	      mirror['https_mirror_url'],
	      mirror['institution'],
	      maintainer,
	      check_mirror_url(mirror['https_mirror_url']),
	      mirror['country_code'],
	      mirror['rsync']]
	    csv << data
	    # Second row is http
	    data[3] = mirror['mirror_url']
	    data[0] = "#{country} (#{mirror['city']})"
	    data[6] = check_mirror_url(mirror['mirror_url'])
	    csv << data
	  end
	end
      end
    end
end
