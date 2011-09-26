#!/usr/bin/env ruby

require 'rubygems'
require 'sqlite3'
require 'time'
require 'pp'

us = ["d.tenenbaum", "p.aboyoun", "n.gopalakrishnan", "hpages@fhcrc.org", "mtmorgan@fhcrc.org",
  "v.obenchain", "c.wong", "m.carlson"]

db = SQLite3::Database.new( "data.db" )

packages = db.execute(%Q(select distinct package from entries where package != "" order by package))

h = {}

rr = []

for package in packages
  package = package.first
  rows = db.execute(%Q(select date, author from entries where package = ? order by date desc limit 20), package)
  r = []
  ok = true
  for row in rows
    t = Time.at(row.first)
    datestr = t.strftime("%m/%e/%Y").gsub(" ", "0")
    next if us.include? row.last
    rr.push [package, datestr, row.last, row.first] if ok
    ok = false
    break unless ok
  end
  h[package] = r
end


good = []
for item in rr
  next if item.empty?
  good.push item
end

sorted = good.sort_by{|e|e.last}

#pp sorted

f = File.open("manifest")
lines = f.readlines
manifest = []
for line in lines
  next if line.chomp.empty?
  next if line =~ /^#/
  manifest.push line.sub("Package: ", "").chomp
end


filtered = []

for item in sorted
  filtered.push item if manifest.include? item.first
end

for item in filtered
  puts "#{item.first}\t#{item[1]}\t#{item[2]}"
end


