#!/usr/bin/env python

from pubsubhubbub_publish import *
import sys

if (len(sys.argv) == 1):
    repo = "bioc"
else:
    if sys.argv[1] == "bioc" or sys.argv[1] == "data-experiment":
        repo = sys.argv[1]
    else:
        print("Argument must be 'bioc' or 'data-experiment'.")
        exit(1)

if repo == "bioc":
    filename = "tmp/rss_urls.txt"
else:
    filename = "tmp/data_rss_urls.txt"

f = open(filename, "r")

urls = []
for line in f:
    urls.append(line.rstrip())

f.close()

try:
    publish('http://pubsubhubbub.appspot.com', urls)
    print("Done.")
except PublishError, e:
    print("Caught Exception!!")

