#!/usr/bin/env ruby


require 'rubygems'
require 'sqlite3'
require 'time'
require 'pp'

f = File.open("manifest")
lines = f.readlines
manifest = []
for line in lines
  next if line.chomp.empty?
  next if line =~ /^#/
  manifest.push line.sub("Package: ", "").chomp
end

db = SQLite3::Database.new("pkgdownloads_db.dtenenba_backup.sqlite")


sql = "select package, count(*) as count, biocrepo from access_log where biocrepo = 'bioc' group by package order by count asc"

rows = db.execute(sql)

filtered  = []
for row in rows
  next unless row.last == 'bioc'
  filtered.push row
end

pp filtered
__END__

counts = {}

for package in manifest
#  puts "#{package}..."
  #sql = "select day_month_year, month_year from access_log where package = ?"
  sql = "select count(*) from access_log where package = ?"
  rows = db.execute(sql, package)
  #counts[package] = rows.size
  counts[package] = rows.first
end



counts.sort {|a,b|b.first<=>a.first}.each do |elem|
  puts "#{elem.first}\t#{elem[1]}"
end
