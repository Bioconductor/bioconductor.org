#!/usr/bin/env ruby

# This script should loop through all active software packages
# (found at http://bioconductor.org/packages/devel/bioc/VIEWS) and for each one,
# see if 1) the maintainer has an account on the support site and 2) that
# they are 'watching' the tag for the package.

# If they are not on the support site, we send them an email.
# If they are on the support site, we add the package to their
# watched tags if it is not there already.

# Don't do anything if the maintainer is maintainer@bioconductor.org.

# In order to minimize annoying emails from us, we should probably
# only send one email to each maintainer even if that maintainer
# maintains more than one package.


require 'sequel'
require 'fileutils'
require 'httparty'

DB = Sequel.connect("postgres://biostar:#{ENV['POSTGRESQL_PASSWORD']}@habu:5432/biostar")
users = DB[:users_user]
profiles = DB[:users_profile]

body = HTTParty.get('http://bioconductor.org/packages/devel/bioc/VIEWS').to_s

lines = body.split("\n")

pm = {}

curpkg = nil
cont_mode = false
mtr_mode = false
mtr_lines = []

def get_emails(input)
  emails = []

  emails
end


for line in lines
  if line.start_with? "Package: "
    curpkg = line.sub("Package:", "").strip
  end
  if line.start_with? "Maintainer: "
    mtr_mode = true
    mtr_lines << line
    next
  end
  if mtr_mode
    if line =~ /^\s/
      mtr_lines << line
    else
      mtr_mode = false
      raw = mtr_lines.join "\n"
      emails = raw.scan(/<([^>]*)>/).flatten
      pm[curpkg] = emails
      #puts curpkg if pm[curpkg].count("<") > 1
      mtr_lines = []
    end
  end
end


need_to_register = []

pm.each_pair do |k, v|
  for email in v
    next if email.downcase == "maintainer@bioconductor.org"
    # select(:id, :watched_tags). # doesn't work
    record = profiles.join(:users_user,
      :id => :user_id).where(Sequel.ilike(:email,
        "%#{email}%")).first

    if record.nil?
      # puts "#{email} isn't signed up, so #{k} is SOL!"
      need_to_register.push({email: email, package: k})
    else
      # make sure all watched tags are lowercase
      watched_tags = record[:watched_tags]
      if watched_tags =~/[A-Z]/
        puts "converting #{watched_tags} to lowercase"
        watched_tags.downcase!
        profiles.where(:id =>
          record[:id]).update(watched_tags: watched_tags)
      end
      if watched_tags != "" and watched_tags !~ /,/ and watched_tags =~ /\s/
        puts "whoa, #{watched_tags} is badly formatted!"
      end
      # puts "watched tags for #{email} are #{watched_tags}."
      pkgs = watched_tags.downcase.split(",")
      pkgs = pkgs.map{|i| i.strip}
      if pkgs.include? k.downcase
        # puts "hells yeah, good lookin out, #{email} has #{k} in watched tags!"
      else
        # add it...
        pkgs << k.downcase
        puts "adding #{k.downcase} to #{email}'s watched tags"
        puts "to end up with"
        puts pkgs.join(", ")
        result = pkgs.join(", ")
        if result.length >= 255

          puts "ERROR! field is too long to be inserted (#{result.length}) !!!!!!"
        else
          profiles.where(id:
            record[:id]).update(watched_tags: pkgs.join(", "))
        end
      end
    end
  end
end

require 'pry';binding.pry


puts 'bye'
