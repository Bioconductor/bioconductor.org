#!/usr/bin/env python


import sys
import re
import time
import calendar
import ConfigParser
from boto import ec2
import subprocess
import urllib2

config_file = "%s/credentials.cfg" % sys.path[0]

config = ConfigParser.RawConfigParser()
config.read(config_file)

access_key = config.get("credentials", "access_key")
secret_key = config.get("credentials", "secret_key")



pat = re.compile("[0-9]{3}Z$")

ec2conn = ec2.connection.EC2Connection(access_key, secret_key)
reservations = ec2conn.get_all_instances()
instances = [i for r in reservations for i in r.instances]
for i in instances:
    d = i.__dict__
    if 'Name' in d['tags'] and d['tags']['Name'] == 'tryitnow':
        # figure out how long it has been running
        lts = d['launch_time'] #format: 2009-10-27T17:10:22.000Z
        lts = re.sub(pat, "UTC", lts)
        lt = time.strptime(lts, "%Y-%m-%dT%H:%M:%S.%Z")
        now = time.gmtime()
        nowsecs = calendar.timegm(now)
        ltsecs = calendar.timegm(lt)
        diff = nowsecs - ltsecs
        if (diff > (60 * 60 * 2)): # 2 hour timeout
            print "Terminating %s launched at %s." % (d['id'], d['launch_time'])
            i.terminate()
