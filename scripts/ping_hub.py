#!/usr/bin/env python

from pubsubhubbub_publish import *


f = open("tmp/rss_urls.txt", "r")

urls = []
for line in f:
    urls.append(line.rstrip())

f.close()

try:
    publish('http://pubsubhubbub.appspot.com', urls)
    print("Done.")
except PublishError, e:
    print("Caught Exception!!")

