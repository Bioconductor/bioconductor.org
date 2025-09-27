#!/usr/bin/env python

from pubsubhubbub_publish import *
import sys

if (len(sys.argv) == 1):
    buildtype = "bioc"
else:
    buildtype = sys.argv[1]

if buildtype == "bioc":
    rssfile = "tmp/rss_urls.txt"
elif buildtype == "data-annotation":
    rssfile = "tmp/data_annnotation_rss_urls.txt"
elif buildtype == "data-experiment":
    rssfile = "tmp/data_experiment_rss_urls.txt"
elif buildtype == "workflows":
    rssfile = "tmp/workflows_rss_urls.txt"
elif buildtype == "books":
    rssfile = "tmp/books_rss_urls.txt"
elif buildtype == "bioc-longtests":
    rssfile = "tmp/longtests_rss_urls.txt"
else:
    print("Argument (%s) is not a valid buildtype" % buildtype)
    exit(1)

f = open(rssfile, "r")

urls = []
for line in f:
    urls.append(line.rstrip())

f.close()

try:
    publish('http://pubsubhubbub.appspot.com', urls)
    print("Done.")
except PublishError, e:
    print("Caught Exception!!")

