#!/usr/bin/env ruby

# to be run by a crontab

# checks to see if there is a message saying
# that a docbuilder build has completed. If so, runs
# rake get_workflows 

require 'rubygems'
require 'aws-sdk'

credentials = YAML::load(File.open('/home/biocadmin/bioc-test-web/etc/credentials.yaml'))

AWS.config({:access_key_id => credentials['access_key_id'],
    :secret_access_key => credentials['secret_key'],
    :region => credentials['region']})
sqs = AWS::SQS.new()
queues = sqs.queues.to_a

q = queues.find {|i| i.url =~ /buildcomplete$/}

msg = q.receive_message

unless msg.nil?
    if msg.body == "build complete"
        Dir.chdir("/home/biocadmin/bioc-test-web/bioconductor.org") do
            system("rake get_workflows >> /home/biocadmin/bioc-test-web/get_workflows.log 2>&1")
        end
    end
end
