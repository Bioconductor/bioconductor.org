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
require 'net/smtp'
require 'yaml'


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
      mtr_lines = []
    end
  end
end


need_to_register = Hash.new { |h, k| h[k] = [] }

pm.each_pair do |k, v|
  for email in v
    next if email.downcase == "maintainer@bioconductor.org"
    # select(:id, :watched_tags). # doesn't work
    record = profiles.join(:users_user,
      :id => :user_id).where(Sequel.ilike(:email,
        "%#{email}%")).first

    if record.nil?
      need_to_register[email] << k
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
      pkgs = watched_tags.downcase.split(",")
      pkgs = pkgs.map{|i| i.strip}
      if pkgs.include? k.downcase
        # great, all is right with the world.
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


if ARGV.length > 0 and ARGV.first == 'nag'
  blacklist = []
  blacklist_file = File.join(File.expand_path(File.dirname(__FILE__)), "maintainer_blacklist.txt")

  if File.exists?(blacklist_file)
    blacklist = File.readlines(blacklist_file).map{|i|i.strip.downcase}
  end
  configfile = File.join(File.expand_path(File.dirname(__FILE__)),
    "mailconfig.yml")
  mailconfig = YAML::load_file(configfile)
  need_to_register.each_pair do |k,v|
    if blacklist.include? k.downcase
      puts "not sending email to #{k} (#{v.join ', '}) because email is in the blacklist"
      next
    end
    message = <<"MESSAGE_END"
From: Dan Tenenbaum <dtenenba@fredhutch.org>
To: #{k}
Subject: Please register for the Bioconductor Support site

Hi, this is an automated message (sent on behalf of Dan Tenenbaum)
asking you to please register for the Bioconductor support site
so that you can more effectively maintain your package(s)
(#{v.join(", ")}). 

If you have registered on the support site, you have not done so
under the same email address that is found in the Maintainer
field of the DESCRIPTION file in your package (#{k}). You should either
register under that email address, or change the email address in
the Maintainer field to match the address that is already subscribed
to the support site (and increment your version number so that
the change is propagated).

To register, visit the link:

https://support.bioconductor.org/accounts/signup/

Once you have registered, visit the link
https://support.bioconductor.org/profile/#div_id_watched_tags
in order to add the packages you maintain to the "Watched tags" list.
For each (comma-separated, lowercase) entry in this list, you'll
receive an email whenever someone posts a question tagged with your
package name. (If you forget this step, it's OK, your packages
will automatically be added to this field once a day).

Thank you for your understanding of this process.
We are just trying to bring people with questions closer to
the people with answers, and make sure that questions
asked don't fall into a hole. 

This automated email will be sent out once a week. 
In order to avoid further nagging, please register today.

Thanks for contributing to Bioconductor and for helping the
community by supporting your package. 

MESSAGE_END
    puts "emailing #{k} about #{v.join(", ")}"
    smtp = Net::SMTP.new(mailconfig['server'], mailconfig['port'])
    smtp.enable_starttls
    smtp.start('localhost', mailconfig['username'],
      mailconfig['password'])
    smtp.send_message(message, "dtenenba@fredhutch.org", k)
  end
end



puts 'Done.'
