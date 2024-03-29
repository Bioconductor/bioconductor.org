# Build System RSS Feeds

The Bioconductor Build System updates RSS feeds for each package
if there were any issues (warnings, errors, timeouts) with the package
build, in both release and devel.

Package maintainers are encouraged to subscribe to 
these feeds to be notified immediately if there are any build problems.

There is a separate RSS URL for each package. The feed URLs take this form:

*Software builds*:

    http://bioconductor.org/rss/build/packages/PACKAGENAME.rss

*Experiment Data builds*:

    http://bioconductor.org/rss/build/data-experiment/packages/PACKAGENAME.rss

*Workflow builds*:

    http://bioconductor.org/rss/build/workflows/packages/PACKAGENAME.rss

*Book builds*:

    http://bioconductor.org/rss/build/books/packages/PACKAGENAME.rss

*Long Tests builds*:

    http://bioconductor.org/rss/build/longtests/packages/PACKAGENAME.rss

To get the feed URL for your package, replace "PACKAGENAME" with the name of
your package. Feed URLs can be pasted into any RSS reader (they are
not meant to be viewed by web browsers). 

(Note that these are actually Atom feeds, but will work in any modern
RSS reader.)

At present, feeds are updated:
- daily (except on Sundays) at 2:45pm EST for the Software builds
- on Tuesdays, Thursdays, and Saturdays, at 4:45pm EST for the Experiment Data builds
- on Tuesdays and Fridays at 3:45pm EST for the Workflow builds
- on Mondays, Wednesdays, and Fridays, at 3:45pm EST for the Book builds
- on Sundays at 12:45am EST for the Long Tests builds

The feeds support PubSubHubbub (also known as PuSH) which means updates
will be pushed immediately to PuSH-enabled readers, eliminating the need to poll
the feeds for new content.

<!--
If a package had no build issues, its feed is not updated.
-->

## RSS To Email

If you are not in the habit of using an RSS reader, you can receive the feed notifications in your email client.

In Thunderbird, RSS feeds can be
[displayed in the email client](http://kb.mozillazine.org/Thunderbird_:_FAQs_:_RSS_Basics). This is possible
in other desktop email clients as well.

There are 
[a number of services](http://blog.themeforest.net/resources/7-rss-to-emailsms-services-you-can-use-for-your-item-feed/)
which will email you the contents of RSS feeds.

You can also set up the [rss2email](http://www.allthingsrss.com/rss2email/)
script to do this for you.

## About the Feeds

The RSS feeds contain a link to your package's daily build report which contains
clickable icons labeled with the build status (OK, ERROR, WARNINGS, TIMEOUT),
that look something like this:

<img src="buildreport.jpg" width="267" height="110"/>

Click on the icon to see the full report. Note that the builds
run daily, and if the issue has cleared up, you may not
see the condition reported by the RSS feed.

Build issues may be transient (for example, a package tries to
download a file from a site which may be down temporarily),
or they may be the result of a problem with the build system itself.
If you have questions about a build issue, contact the
[bioc-devel](http://bioconductor.org/help/mailing-list/)
mailing list.

