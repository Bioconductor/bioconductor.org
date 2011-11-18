#!/usr/bin/env ruby

if ARGV.size != 2
  puts "supply hostname and bioC version number"
  exit
end

host, biocVer = ARGV

url = "https://hedgehog.fhcrc.org/bioconductor/"

branch_url = url + "RELEASE_" + biocVer.gsub(".", "_") + "_branch"

manifest_url = branch_url + "/madman/Rpacks/bioc_#{biocVer}.manifest"

result = `curl -s -u readonly:readonly #{manifest_url}`

if (result =~ /403 Forbidden/)
  branch_url = url + "trunk"

  manifest_url = branch_url + "/madman/Rpacks/bioc_#{biocVer}.manifest"
  result = `curl -s -u readonly:readonly #{manifest_url}`
  
end

lines = result.split("\n")
manifest_pkgs = []


for line in lines
  next unless line =~ /^Package: /
  manifest_pkgs.push line.chomp.strip().gsub(/^Package: /, "")
end

if manifest_pkgs.empty?
  puts "nothing in manifest? Bailing..."
  exit
end

find_output = `ssh biocadmin@#{host} find "~/PACKAGES/#{biocVer}/bioc"`

lines = find_output.split("\n")

bad_pkgs = []

for line in lines
  line.chomp!
  if line =~ /\.tar\.gz$|\.tgz$|\.zip$/
    segs = line.split("/")
    fullPkg = segs.last
    pkgName = fullPkg.split("_").first
    if (!manifest_pkgs.include? pkgName)
      puts "BAD PACKAGE! removing #{fullPkg}"
      puts "to delete it, edit this script."
      ##result = `ssh biocadmin@#{host} rm #{line}` # uncomment this line to activate deletion
    end
    
  end
end
