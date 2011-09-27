#!/usr/bin/env ruby


require 'rubygems'
require 'sqlite3'
require 'time'
require 'date'
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

# todo add index on 'method' and 'errorcode' columns
sql = "select package, day_month_year, month_year, method, errorcode from access_log where package = ? /* and method = 'GET'*/  /*and errorcode = '200'*/"



counts = {}

cutoff = 10
i = 0
for package in manifest.reverse
  #break if i == cutoff
  i += 1
  $stderr.puts "#{package}..."
  rows = db.execute(sql, package)
  
  filtered = []
  for row in rows
    next unless row.last == "200"
    next unless row[3] == 'GET'
    filtered.push row
  end
  
  rows = []
  
  for row in filtered
    dmy = row[1]
    dmydate = DateTime.strptime(dmy, "%d/%b/%y")
    dmytime = Time.parse(dmydate.to_s)
    dmyi = dmytime.to_i
    if counts.has_key? package
      if counts[package]["dmyi"] > dmyi
        counts[package] = {"size" => filtered.size, "dmyi" => dmyi, "dmy" => dmy}
      end
    else
      counts[package] = {"size" => filtered.size, "dmyi" => dmyi, "dmy" => dmy}
    end
  end
end


$stderr.puts "------"
#pp counts

counts.each_pair do |k, v|
  puts "#{k}\t#{v["size"]}\t#{v["dmy"]}\t#{v["dmyi"]}"
end

