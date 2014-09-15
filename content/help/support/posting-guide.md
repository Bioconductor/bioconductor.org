# ![](/images/icons/magnifier.gif)Posting guide #

Revised: 18 August, 2014

## Why? ##

* Good questions attract good responses.

* Formulating a good question can help point to the answer.

* Good questions help avoid embarrassment.

## Key points ##

* Use the latest Bioconductor release (or 'devel') version. Ensure
  that your packages are [up-to-date](/install/).

* Provide a minimal, self-contained example (data and code) that
  someone else can cut-and-paste into a new R session to reproduce
  your problem.

* Include the output of `sessionInfo()`. If the example produces an
  error, provide the output of `traceback()` after the error occurs.

* Include the maintainer of the package in the email, e.g.,
  `packageDescription('GenomicRanges')$Maintainer`.

## Preparing ##

* Read the relevant R documentation.  Use `help.start()` to start the
  HTML search engine. If you are having trouble with function
  `somefunc`, try `?somefunc`. If you are searching for a function,
  try `help.search("somefunc")`. Read the vignette(s) for the
  package(s) using `vignette(package="somepkg")`

* Search the [FAQ](/docs/faq/) and
  [support site](https://support.bioconductor.org/local/search/page/)
  and [R-help](http://tolstoy.newcastle.edu.au/R/) archives for
  similar posts. Try a [Google](http://www.google.com/) search.

## Composing ##

* Compose a new message with a new subject line; only reply to
  an existing post if you are elaborating on or answering a previous
  question.

* Use an informative subject line. This will help attract responders
  and helps others when they search the archives.  For example, a
  subject line "time course experiment using limma" is better than
  just "limma".

* Tag your post appropriately. Choose one or more tags so that people
  following those tags will see your posts. If you're posting about a
  specific package, tag your post with the name of that package.

* Identify yourself.  If you are logged in with a non-professional email
  account, like gmail or hotmail, then
  give your full name and professional affiliation.

* If you are asking for advice on how to use a particular function or
  package, then explain fully what documentation you have already read
  and why this hasn't yet fully answered your question.  This allows a
  responder to answer your specific question instead of simply
  referring you to the existing documentation.

* If your question is about package development, send email to the
  [bioc-devel list](/support/#bioc-devel); otherwise post to the
  [Bioconductor support site](https://support.bioconductor.org).

* Package developers always appreciate being alerted to possible bugs,
  but be very sure that you have used the package correctly.  In most
  cases, best practice is (i) to double-check the documentation and
  then (ii) report the behavior that was unexpected to you.

* Read over your post. Is it polite and easy to understand? Remember
  that responses on the Bioconductor support site are from volunteers.

## Replying ##

* If possible, write an answer that can be understood by readers with
  different scientific backgrounds and skill levels / skill
  profiles.

* When answering a question, consider including some explanation of
  how you arrived at your solution. This way, you help people not only
  to solve their problem at hand, but also to help themselves in the
  future.

* Rudeness and ad hominem comments are not acceptable.

* Brevity is OK. Consider, though, that information that is obvious to
  you may be very helpful to others.

* Provide URLs to other relevant threads or web sites.

## Acknowledgments ##

This posting guide has been heavily adapted from the
[R-help posting guide](http://www.r-project.org/posting-guide.html),
various comments and suggestions by Bioconductor posters and inspired
by Eric Steven Raymond's essay on
[How To Ask Questions The Smart Way](http://www.catb.org/%7Eesr/faqs/smart-questions.html).

