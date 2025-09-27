# ![](/images/icons/magnifier.gif) Posting guide

<p class="text-large">
General advise to posting on the Bioconductor help sites. 
Please remember when posting a question or response to abide by the
<a class="text-large" href="/about/code-of-conduct/">Bioconductor Code of
Conduct</a>.
</p>

Revised: 1 Dec, 2023

<a name="top"></a>

- [Why](#why)
- [Key points](#keypoints)
- [Communication with Package Maintainers](#maintainers)
- [Before Posting](#before)
- [Composing](#composing)
- [Replying](#replying)



<a name="why"></a>

## Why?

* Well-formed questions attract better responses and help save
  everyone's time.

* Please remember that the helpers on the support site are volunteers
  who are often quite knowledgeable about _Bioconductor_, but won't
  know what you have done already, nor what you are trying to do. The
  key points below are intended to help you give enough information
  that an experienced person could provide a useful answer without
  requesting more information.

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>
<a name="keypoints"></a>

## Key points

* Use the latest _Bioconductor_ [release version][].  Ensure that your
  packages are [up-to-date][].

* Post all of your _R_ code.

* Include a copy of all output including any error or warning messages that
  appeared in _R_.
  
* If your question involves experimental data, include an example of
  the [sample-level covariate data][] (one row per sample, one column
  per covariate). If it would help answer your technical question, and
  can be shared, explain the motivation behind your experiment.

* Always paste the output of `sessionInfo()` at the end of your post.

* If possible, provide a minimal, self-contained example that someone
  else can cut-and-paste into a new _R_ session to reproduce your
  problem.

* If the example produces an error, provide the output of
  `traceback()` after the error occurs.

[release version]: /packages/release/BiocViews.html#___Software
[up-to-date]: /install#update-bioconductor-packages
[sample-level covariate data]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4509590/figure/F2/

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>
<a name="maintainers"></a>

## Communication with package maintainers

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

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>
<a name="before"></a>

## Before posting


* Read the vignettes for the package in question using
  `browseVignettes("somepkg")`. Search the vignette text for sections
  relevant to your question.

* Read the relevant _R_ documentation. If you are having trouble with
  function `somefunc()`, read the relevant sections of `?somefunc`,
  including the sections about function arguments and function output.

* Search the support site for similar posts. The best way to search
  the support site is using [Google](https://google.com) by starting
  with `site:support.bioconductor.org` followed by a space and a set
  of keywords and the package name.

* Ensure that your question is on topic. The _Bioconductor_ support site
  is intended primarily for helping people with questions about using
  _Bioconductor_ software. Questions about non-_Bioconductor_ _R_ packages
  or general statistics should be asked on
  [_R_-help](https://www.r-project.org/help.html) or
  [Stackoverflow](https://stackoverflow.com/).  Questions about
  bioinformatics in general, but not about _Bioconductor_ packages can
  be posted to [Biostars](https://www.biostars.org/).

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>
<a name="composing"></a>

## Composing

* Use an informative post title. This will help attract responders and
  helps others when they search the support site. For example, a
  subject line "time course experiment using limma" is better than
  just "limma".

* If you are asking for advice on how to use a particular function or
  package, then explain fully what documentation you have already read
  and why this hasn’t yet fully answered your question. This allows a
  responder to answer your specific question instead of simply
  referring you to the existing documentation.

* As mentioned above, tag your post with the names of the appropriate
  packages.

* Include code chunks with all output. Use a leading and trailing triple
  backtick to indicate a code chunk (<code>```</code>). You can directly copy and paste code
  from your console inbetween the backticks. 

* Always include a code chunk with the results of `sessionInfo()`. This will
  indicate what os platform version you are running as well as your versions of
  R and R/_Bioconductor_ package versions. 

* Read over your post. Is it polite and easy to understand? Remember
  that responses on the _Bioconductor_ support site are from volunteers
  trying to help.

* There is also a helpful tutorial to using/formatting with the markdown
  editor. This can be found [here](https://support.bioconductor.org/p/117436/)

* Package development questions are best asked by sending an email to the
  [bioc-devel list](https://stat.ethz.ch/mailman/listinfo/bioc-devel).


<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>
<a name="replying"></a>

## Replying

* Use the *Add Comment* button to ask for further explanation of
  answers, or to respond to requests for more information.

* Only use the *Add Answer* button if you are answering the original
  question.

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

<p class="back_to_top">[ <a href="#top">Back to top</a> ]</p>



<p class="publications-footer">
This posting guide has been heavily adapted from the
<a href="http://www.r-project.org/posting-guide.html">R-help posting guide</a>,
various comments and suggestions by _Bioconductor_ posters and
inspired by Eric Steven Raymond's essay on
<a href="http://www.catb.org/%7Eesr/faqs/smart-questions.html">How To Ask Questions The Smart Way</a>.
</p>