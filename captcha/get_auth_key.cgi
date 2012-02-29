#!/usr/bin/env python

# workaround until http requests to some_host:8787 will work on mamba

import cgi
import urllib2


print "Content-type: text/plain"
print

form = cgi.FieldStorage()

host = form['host'].value

url = "http://%s:8787/auth-public-key" % host

#print url

while True:
    try:
        f = urllib2.urlopen(url, timeout=1)
        print f.read()
        f.close()
        break
    except urllib2.URLError:
        pass
        #print "exception..."

