require 'json'

module GenerateRedirects

  def pkg_is_only_in_devel(pkg, db, config)
    release_version = config["release_version"]
    devel_version = config["devel_version"]
    result = db[pkg]
    result[:versions].length == 1 and 
      result[:versions].first == devel_version
  end

  def pkg_is_in_release(pkg, db, config)
    release_version = config["release_version"]
    devel_version = config["devel_version"]
    result = db[pkg]
    result[:versions].include? release_version
  end

  def generate_redirects()
    config = YAML.load_file("./config.yaml")
    release_version = config["release_version"]
    devel_version = config["devel_version"]
    basedir = File.join("assets", "packages", "json")
    d = Dir.new(basedir)
    db = {}
    for dir in d.entries
      next if [".", ".."].include? dir
      fulldir = File.join(basedir, dir)
      next unless File.directory? fulldir
      for subdir in ['bioc', File.join('data', 'annotation'), File.join('data', 'experiment')]
        repos = subdir.sub(File::SEPARATOR, "/")
        jsonfile = File.join(fulldir, subdir, "packages.json")
        next unless File.exists? jsonfile
        fh = File.read(jsonfile)
        obj = JSON.parse(fh)
        for key in obj.keys
          if db.has_key? key
            db[key][:versions] << dir
          else
            db[key] = {:repos => repos, :versions => [dir]}
          end
        end
      end        
    end
    # better to write to assets or directly to output???
    outfile = File.join('output', 'packages', '.htaccess')
    out_fh = File.open(outfile, "w")
    out_fh.puts '<IfModule mod_rewrite.c>'
    out_fh.puts 'RewriteEngine On'
    out_fh.puts 'RewriteBase /packages/'
    out_fh.puts
      db.each_pair do |key, value|
        repos = value[:repos]
        versions = value[:versions]

        #RewriteRule ^Biobase/$ /packages/release/bioc/html/Biobase.html [PT]
        if pkg_is_only_in_devel(key, db, config)
          out_fh.puts "RewriteRule ^#{key}/$ /packages/devel/#{repos}/html/#{key}.html [PT]"
        elsif pkg_is_in_release(key, db, config)
          out_fh.puts "RewriteRule ^#{key}/$ /packages/release/#{repos}/html/#{key}.html [PT]"
        end


        for version in versions
          out_fh.puts "RewriteRule ^#{version}/#{key}/$ /packages/#{version}/#{repos}/html/#{key}.html [PT]"
          if version == release_version
            out_fh.puts "RewriteRule ^release/#{key}/$ /packages/release/#{repos}/html/#{key}.html [PT]"
          end
          if version == devel_version
            out_fh.puts "RewriteRule ^devel/#{key}/$ /packages/devel/#{repos}/html/#{key}.html [PT]"
          end
        end
        out_fh.puts
      end
    out_fh.puts
    out_fh.puts '</IfModule>'
    out_fh.close
  end
end