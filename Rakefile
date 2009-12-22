require 'nanoc3/tasks'
require 'yaml'
require 'fileutils'

task :copy_assets do
  site_config = YAML.load_file("./config.yaml")
  output_dir = site_config["output_dir"]
  system "rsync -gprt --partial --exclude='.svn' assets/ #{output_dir}"
end

task :compile do
  system "nanoc3 co"
end

task :build => [ :compile, :copy_assets ]

task :default => :build
