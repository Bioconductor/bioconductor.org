require 'nanoc3/tasks'
require 'yaml'
require 'fileutils'
require 'lib/data_sources/gmane_list.rb'
#require 'lib/data_sources/bioc_views.rb'
require 'scripts/search_indexer.rb'
require 'open3'

include Open3


desc "copy assets to output directory"
task :copy_assets do
  site_config = YAML.load_file("./config.yaml")
  output_dir = site_config["output_dir"]
  system "rsync -gprt --partial --exclude='.svn' assets/ #{output_dir}"
end

desc "Run nanoc3 compile"
task :compile do
  system "nanoc3 co"
end

desc "Nuke output directory !! uses rm -rf !!"
task :real_clean do
  site_config = YAML.load_file("./config.yaml")
  output_dir = site_config["output_dir"]
  FileUtils.rm_rf(output_dir)
  FileUtils.mkdir_p(output_dir)
end

desc "Build the bioconductor.org site (default)"
task :build => [ :compile, :copy_assets ]

task :default => :build

task :deploy_staging do
  dst = '/loc/www/bioconductor-test.fhcrc.org'
  site_config = YAML.load_file("./config.yaml")
  output_dir = site_config["output_dir"]
  system "rsync -av --partial --partial-dir=.rsync-partial --exclude='.svn' #{output_dir}/ #{dst}"
end

task :deploy_production do
  site_config = YAML.load_file("./config.yaml")
  src = '/loc/www/bioconductor-test.fhcrc.org'
  dst = site_config["production_deploy_root"]
  system "rsync -av --partial --partial-dir=.rsync-partial --exclude='.svn' #{src}/ #{dst}/"
end

desc "Re-index the site for the search engine"
task :search_index do
  if (SearchIndexer.is_solr_running?)
    hostname = `hostname`.chomp
    args = [] # directory to index, location of cache file, location of output shell script
    if (hostname =~ /^dhcp/i)
      args = ['./output', './', 'scripts'] 
    elsif (hostname == 'merlot2')
      args = ['/loc/www/bioconductor-test.fhcrc.org', '/home/biocadmin', '/home/biocadmin']
    end
    pwd = FileUtils.pwd
    si = SearchIndexer.new(args)
    FileUtils.cd pwd # just in case
    cmd = "#{args.last}/index.sh"

    chmod_cmd = "chmod u+x #{cmd}"
    system chmod_cmd
    # todo fork the indexer to a background task because it could take a while in some cases
    #stdin, stdout, stderr = Open3.popen3("#{cmd}")
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
  system("scp scripts/search_indexer.rb webadmin@krait:~")
  system("ssh webadmin@krait /home/webadmin/do_index.rb")
  system("ssh webadmin@krait chmod +x /home/webadmin/index.sh")
  system("ssh webadmin@krait /home/webadmin/index.sh")
end

desc "Runs nanoc's dev server on localhost:3000"
task :devserver => [:build] do
  system "nanoc3 aco"
end
