# ![](/images/icons/magnifier.gif)Posting guide #

## Why? ##

* Good questions attract good responses.

* Formulating a good question can help point to the answer.

* Good questions help avoid embarrassment.

## Preparing ##

* Read the relevant R documentation.  Use `help.start()` to start the
  HTML search engine. If you are having trouble with function
  `somefunc`, try `?somefunc`. If you are searching for a function,
  try `help.search("somefunc")`. Read the vignette(s) for the
  package(s) using `browseVignettes(package="somepkg")`

* Search the [FAQ](http://bioconductor.org/docs/faq/) and
  [Bioconductor](http://dir.gmane.org/gmane.science.biology.informatics.conductor)
  and [R-help](http://tolstoy.newcastle.edu.au/R/) archives for
  similar posts. Try a [Google](http://www.google.com/) search.
  
* Please include the R code that is causing the problem and enough data
  for someone else to run the code and get the same problem.
  Simplify code to a minimal, self-contained example. If reporting an
  error, be sure to reproduce the error in a new R session, started
  with the `--vanilla` option to avoid loading .Rprofile or .RData
  files. 
  

  
* Ensure that you are using the latest Bioconductor release and that
  your installed packages are [up-to-date](/install/).

## Composing ##

* You may reply to an existing email if you are elaborating on or
  answering a previous question.  However, if you are asking a new
  question, then start a new email message with a new subject line.
  Don't simply add a new question on a different topic as part of a
  reply to an existing thread.

* Use an informative subject line that is as specific as
  possible. This will help attract responders and also helps others in
  the future when they search the archives.  For example, a subject
  line "time course experiment using limma" might be better than just
  "limma".

* Identify yourself.  If you are using a non-professional email
  account, like gmail or hotmail, or the guest posting facility, then
  give your full name and professional affiliation.  Anonymous
  postings are much less likely to get responses.

* If you are asking for advice as to how to use a particular function
  or package, then explain fully what documentation you have already
  read and why this hasn't yet fully answered your question.  This
  allows a responder to answer your specific question instead of
  simply referring you to the existing documentation.

* If you are reporting a code problem, then include a simple and
  reproducible example along with the output of `sessionInfo()`. The
  example should be reproducible by others, meaning that someone else
  reading your post should be able to run the code themselves and get
  the same output as you did.  If the example produces an error,
  provide the error message and the output of `traceback()`.  Provide
  giving actual output to demonstrate what you mean when you wish to
  indicate that the code does not work as expected.

* Send email to the appropriate list.  Use the Bioconductor mailing
  list for questions about specific package, or conceptual
  questions. The R-help mailing list is for questions about the
  underlying R program. The Bioc-devel and R-devel lists are for
  discussing code development and other technical issues.

* Ensure that your code is readable, and can be cut and pasted into an
  R session.

* Use plain text instead of HTML; it is smaller in size and easier to
  read.

* The following attachment types are accepted: png, pdf, rda/Rdata. Total 
  message size cannot exceed 1MB. If larger attachments are essential, post
  them to a publicly accessible location and include the link in
  your email.

* Package developers always appreciate being alerted to possible bugs,
  but be very sure that you have used the package correctly.  In most
  cases, best practice is (i) to double-check the documentation and
  then (ii) report the behavior that was unexpected to you.

* Read over your mail. Is it polite and easy to understand? Remember
  that responses on the Bioconductor mailing list are from volunteers.

## Replying ##

* Respond to everyone in the list, which ensures that your response is
  archived.

* If possible, write an answer that can be understood by readers with
  different scientific backgrounds and skill levels / skill
  profiles.

* When answering a question, consider including some explanation of
  how you arrived at your solution (rather than just writing a
  solution). This way, you help people not only to solve their problem
  at hand, but also to help themselves in the future.

* Rudeness and ad hominem comments are not acceptable.

* Brevity is OK. Consider, though, that information that is obvious
  and not worth mentioning to you may be very helpful to others.

* If you believe the issue has been discussed before, please give the
  URL of the relevant thread or web site.

## Acknowledgments ##

This posting guide has been heavily adapted from the
[R-help posting guide](http://www.r-project.org/posting-guide.html),
various comments and suggestions by Bioconductor posters and inspired
by Eric Steven Raymond's essay on
[How To Ask Questions The Smart Way](http://www.catb.org/%7Eesr/faqs/smart-questions.html).
