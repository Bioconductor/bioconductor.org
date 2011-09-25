#!/usr/bin/env ruby

require 'rubygems'
require 'sqlite3'
require 'rexml/document'
require 'time'

include REXML


FileUtils.rm_rf "data.db"
system(%Q(echo "CREATE TABLE entries (id int, package text, revision int, author text, date int, logmsg text);"|sqlite3 data.db))

db = SQLite3::Database.new( "data.db" )

# 2011-09-25T15:11:01.592582Z
# 2011-09-25T15:11:01.Z


def remove_millis(str)
    segs = str.split(".")
    segs.pop
    "#{segs.first}.Z"
end

#date = DateTime.strptime("2011-09-25T15:11:01.592582Z", "%Y-%m-%dT%H:%M:%S.%Z")


# create log.xml as follows:
# svn log -v --xml https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/ > log.xml

file = File.new("log.xml")
doc = Document.new(file)


id = 1
ok = true
doc.elements.each("log/logentry") do |logentry|
  puts "done" if ok
  ok = false
  revision = logentry.attributes["revision"]
  tmp = []
  logentry.elements.each("author") {|i| tmp.push i}
  author = tmp.last.text
  logentry.elements.each("msg") {|i| tmp.push i}
  logmsg = tmp.last.text
  logentry.elements.each("date") {|i| tmp.push i}
  datestr = remove_millis(tmp.last.text)
  date = Time.parse(datestr, "%Y-%m-%dT%H:%M:%S.%Z").to_i
  packages = {}
  logentry.elements.each("paths/path") do |path|
    pathstr = path.text
    pathstr.sub! "/trunk/madman/Rpacks/", ""
    segs = pathstr.split "/"
    next if segs.length < 2
    package = segs.first
    packages[package] = 1
  end
  packages.each_key do |package|
    sql = "insert into entries values (?, ?, ?, ?, ?, ?)"
    db.execute sql, id, package, revision, author, date, logmsg
    id += 1
    #puts "inserting row #{id}"
  end
end