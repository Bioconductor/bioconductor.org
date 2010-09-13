require 'nanoc3/tasks'
require 'yaml'
require 'fileutils'
require 'lib/data_sources/gmane_list.rb'
require 'lib/data_sources/bioc_views.rb'
require 'scripts/search_indexer.rb'


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
    args = []
    if (hostname =~ /^dhcp/i)
      args = ['./output', './', './scripts']
    elsif (hostname == 'merlot2')
      args = ['/loc/www/bioconductor-test.fhcrc.org', '/home/biocadmin', '/home/biocadmin']
    end
    pwd = FileUtils.pwd
    si = SearchIndexer.new(args)
    FileUtils.cd pwd # just in case

    chmod_cmd = "chmod o+x #{args.last}/index.sh"
    system chmod_cmd
    # we could run the indexer here but it could take a while, maybe better to fork it to a background task
    result = `ls&`
    puts result
  else
    puts "solr is not running, not re-indexing site."
  end
end

desc "Runs nanoc's dev server on localhost:3000"
task :devserver => [:build] do
  system "nanoc3 aco"
end
