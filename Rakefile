require 'nanoc3/tasks'
require 'yaml'
require 'fileutils'
require 'lib/data_sources/gmane_list.rb'

task :copy_assets do
  site_config = YAML.load_file("./config.yaml")
  output_dir = site_config["output_dir"]
  system "rsync -gprt --partial --exclude='.svn' assets/ #{output_dir}"
end

task :compile do
  system "nanoc3 co"
end

task :build => [ :compile, :copy_assets ]

task :real_clean do
  site_config = YAML.load_file("./config.yaml")
  output_dir = site_config["output_dir"]
  FileUtils.rm_rf(output_dir)
  FileUtils.mkdir_p(output_dir)
end

task :default => :build

task :deploy_merlot2_local do
  dst = '/usr/local/nginx/sites/bioconductor.org'
  site_config = YAML.load_file("./config.yaml")
  output_dir = site_config["output_dir"]
  system "rsync -gvprt --partial --exclude='.svn' #{output_dir}/ #{dst}"
end
