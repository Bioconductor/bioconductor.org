#!/usr/bin/env ruby

require 'sequel'
require 'fileutils'
require 'httparty'

require_relative './svn_shield_helper.rb'

DB = Sequel.connect("postgres://biostar:#{ENV['POSTGRESQL_PASSWORD']}@habu:5432/biostar")



def get_post_tag_info()

  pkgs = []

  [true, false].each do |state|
    pkgs += get_list_of_packages(state)
    pkgs += get_annotation_package_list(state)
  end

  pkgs = pkgs.uniq

  posts_post = DB[:posts_post]


  today = Date.today
  now = DateTime.new(today.year, today.month, today.day)
  sixmonthsago = now
  months = [now]

  6.times do
    tmp = sixmonthsago.prev_month
    months << tmp
    sixmonthsago = tmp
  end
  months.reverse!
  ranges = []
  for i in 0..(months.length()-2)
    ranges.push months[i]..months[i+1]  
  end
  new_range = ranges.last.first..ranges.last.last.next_day
  ranges.pop
  ranges.push new_range


  res = posts_post.where("lastedit_date > ?", sixmonthsago).select(:id, :tag_val, :status, :type, :has_accepted, :root_id, :parent_id, :reply_count).all


  hsh = Hash.new { |h, k| h[k] = [] }

  for item in res
    id = item[:id].to_i
    tags = item[:tag_val].split(',')
    for tag in tags
      tag.strip!
      tag.downcase!
      hsh[tag] << id
    end
  end

  hsh.each_pair {|k,v| hsh[k] = v.sort.uniq}

  # Support activity: tagged questions, answers / comments per question;
  # % closed, 6 month rolling average

  zero_shield = File.join("assets", "images", "shields", "posts",
    "zero.svg")
  dest_dir = File.join("assets", "shields", "posts")
  # remove dir first?
  FileUtils.mkdir_p dest_dir

  for pkg in pkgs
    if hsh.has_key? pkg.downcase
      num = hsh[pkg.downcase].length
      relevant = res.find_all{|i| hsh[pkg.downcase].include? i[:id]}
      questions = relevant.find_all{|i| i[:id] == i[:parent_id]}

      q = questions.length
      closed = questions.find_all{|i| i[:has_accepted] == true}.length
      answers = []
      comments = []

      for question in questions
        answers << res.find_all{|i| question[:id] == i[:root_id] and i[:type] == 1}.length
        comments << res.find_all{|i| question[:id] == i[:root_id] and i[:type] == 6}.length
      end


      a_avg =  sprintf("%0.1g", answers.inject(0.0) { |sum, el| sum + el } / answers.size)
      c_avg =  sprintf("%0.1g", comments.inject(0.0) { |sum, el| sum + el } / comments.size)
      shield_text = "#{q} / #{a_avg} / #{c_avg} / #{closed}".gsub(' ', '%20')
      puts "getting shield for #{pkg}"
      response = HTTParty.get("https://img.shields.io/badge/posts-#{shield_text}-87b13f.svg")
      sf = File.open(File.join(dest_dir, "#{pkg}.svg"), "w")
      sf.write(response.to_s)
      sf.close

    else
      FileUtils.cp zero_shield, File.join(dest_dir, "#{pkg}.svg")
    end
  end

end

if __FILE__ == $0
  do_it()
end
