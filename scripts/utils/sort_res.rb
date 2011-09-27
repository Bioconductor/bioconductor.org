#!/usr/bin/env ruby

f = File.open("res")
lines = f.readlines

##lines.each{|line| line.chomp!}

lines.sort! do |a,b|
  sa = a.split("\t")
  sb = b.split("\t")
  ia = Integer(sa[1])
  ib = Integer(sb[1])
  ia <=> ib
end

for line in lines
  puts line
end