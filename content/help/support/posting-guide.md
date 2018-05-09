# ![](/images/icons/magnifier.gif)Posting guide #

Revised: 8 May, 2018

## Why? ##

* Well-formed questions attract better responses and help save
  everyone's time. 
* Please remember that the helpers on the support site are volunteers
  who are often quite knowledgeable about Bioconductor, but won't know
  what you have done already, nor what you are trying to do. The key
  points below are intended to help you give enough information that
  an experienced person could provide a useful answer without
  requesting more information.

## Key points ##

* Use the latest Bioconductor [release version](/packages/release/BiocViews.html#___Software). 
  Ensure that your packages are [up-to-date](/install#update-bioconductor-packages).
* Post all of your R code.
* Include a copy of any error or warning messages that appeared in R.
* If your question involves experimental data, include an example of
  the [sample-level covariate data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4509590/figure/F2/)
  (one row per sample, one column per covariate). If it would help
  answer your technical question, and can be shared, explain the
  motivation behind your experiment.
* Always paste the output of `sessionInfo()` at the end of your
  post. Alternatively, the output of `session_info()` from the
  *devtools* or *sessioninfo* packages can be posted. 
* If possible, provide a minimal, self-contained example that someone
  else can cut-and-paste into a new R session to reproduce your
  problem. 
* If the example produces an error, provide the output of
  `traceback()` after the error occurs.

## Communication with package maintainers ##

* Tag your post with the names of the appropriate packages. Putting
  the correct tag is important, as this triggers an automatic email
  from the system to the package maintainer(s). Capitalization doesn't
  matter, but spelling does matter for the package tags. 
* Do not email the package maintainer directly. As much as possible,
  keep all conversations on the support site, and do not revert to
  email unless indicated by the other party. Some package maintainers
  have an alternate URL for bug reports, found under *BugReports* on
  the package landing page. 
* Why? Because the support site is public and searchable, it acts as
  an archive providing benefit to other users who may face similar
  errors or issues. 

## Before posting ##


* Read the vignettes for the package in question using
  `browseVignettes("somepkg")`. Search the vignette text for sections
  relevant to your question.
* Read the relevant R documentation. If you are having trouble with
  function `somefunc`, read the relevant sections of `?somefunc`,
  including the sections about function arguments and function
  output. 
* Search the support site for similar posts. The best way to search
  the support site is using [Google](https://google.com) by starting
  with `site:support.bioconductor.org` followed by a space and a set 
  of keywords and the package name. 
* If searching the support site with Google doesn't help, then try searching
  the [FAQ](/docs/faq/)
  and [R-help](https://tolstoy.newcastle.edu.au/R/) archives for
  similar posts.
* Ensure that your question is on topic. The Bioconductor support site
  is intended primarily for helping people with questions about using
  Bioconductor software. Questions about non-Bioconductor R packages
  or general statistics should be asked on 
  [R-help](https://www.r-project.org/help.html) or
  [Stackoverflow](https://stackoverflow.com/). 
  Questions about bioinformatics in general, but not
  about Bioconductor packages can be posted to 
  [Biostars](https://www.biostars.org/).

## Composing ##

* Use an informative post title. This will help attract responders and
  helps others when they search the support site. For example, a subject
  line "time course experiment using limma" is better than just
  "limma". 
* If you are asking for advice on how to use a particular function or
  package, then explain fully what documentation you have already read
  and why this hasn’t yet fully answered your question. This allows a
  responder to answer your specific question instead of simply referring
  you to the existing documentation. 
* If your question is about package development, send email to the
  [bioc-devel list](/support/#bioc-devel); otherwise post to the
  [Bioconductor support site](https://support.bioconductor.org).
* As mentioned above, tag your post with the names of the appropriate
  packages.
* Read over your post. Is it polite and easy to understand? Remember
  that responses on the Bioconductor support site are from volunteers
  trying to help.

## Replying ##

* Only use the *Add Answer* button if you are answering the original
  Question.
* Use the *Add Comment* button to create a threaded conversation off
  of an existing Answer. 
* If possible, write an answer that can be understood by readers with
  different scientific backgrounds, skill levels / skill profiles. Use
  simple English whenever possible to make your answer accessible to
  the widest audience. 
* When answering a question, consider including some explanation of
  how you arrived at your solution. This way, you help people not only
  to solve their problem at hand, but also to help themselves in the
  future. 
* Rudeness and personal attacks are not acceptable. 
* Brevity is OK. Consider, though, that information that is obvious to
  you may be very helpful to others. 
* Provide URLs to other relevant threads or web sites. 

## Acknowledgments ##

This posting guide has been heavily adapted from the
[R-help posting guide](http://www.r-project.org/posting-guide.html),
various comments and suggestions by Bioconductor posters and inspired
by Eric Steven Raymond's essay on
[How To Ask Questions The Smart Way](http://www.catb.org/%7Eesr/faqs/smart-questions.html).

