#!/usr/bin/env ruby

# This script should eventually be available at
# https://hedgehog.fhcrc.org/bioconductor/trunk/bioconductor.org/scripts

require 'date'

# assuming we are in the same dir as the DB...
dir = Dir.pwd

now = Date.today

months = []

t = now

for i in 0..5
  t = t<<1
  months << t
end

months.reverse!

month_names = []

months.each {|i| month_names << i.strftime("%b/%Y")}

sql = ".mode csv\nselect  count(distinct ips), package from access_log where month_year in ('#{month_names.join "','"}') group by month_year, package;"

puts sql
