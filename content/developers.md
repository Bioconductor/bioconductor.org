_Bioconductor_ is an open development project, meaning that all
developers from the scientific community are able to contribute
software. Listed below are helpful links which will guide developers
at different stages of their package. 

<h2 id="make">Design your package</h2>

* [Primer for Making Packages](/developers/how-to/buildingPackagesForBioc/)   
  Understanding the structure of a typical R package and helpful pointers to 
  keep in mind when writing a package for _Bioconductor_. 

* [Package Guidelines](/developers/package-guidelines/)    
  In order to foster coherence within the project, developers
  are encouraged to follow package guidelines to make it easier for
  others to use and extend the software.

* [Creating Workflow Vignettes](/developers/how-to/workflows/)   
  Workflows are becoming extremely common in Biology. Learn how you can 
  submit your own to _Bioconductor_.

Some Helpful Coding Resources

<ul class="inline_list">
    <li><a href="/developers/how-to/coding-style/">Coding&nbsp;Style</a></li>
    <li><a href="/developers/how-to/useDevel/">Use&nbsp;Bioc&nbsp;Devel&nbsp;Version</a></li>
    <li><a href="/developers/how-to/version-numbering/">Version&nbsp;Numbering</a></li>
    <li><a href="/developers/how-to/deprecation/">Function&nbsp;Deprecation</a></li>
    <li><a href="/developers/how-to/mavericks-howto/">C++/Mavericks&nbsp;Best&nbsp;Practises</a></li>
    <li><a href="/developers/how-to/c-debugging/">Debug&nbsp;C/C++&nbsp;Code</a></li>
    <li><a href="/developers/how-to/unitTesting-guidelines/">Unit&nbsp;Tests</a></li>
    <li><a href="/developers/how-to/biocViews/">Using&nbsp;biocViews</a></li>
</ul>


<h2 id="submit">Submit to <i>Bioconductor</i></h2>

* [Package Submission](/developers/package-submission)    
  New packages are added to the project through a package submission
  process after a review. Developers can submit software, experiment 
  data and annotation data packages to us. 

* [New Packages](/developers/new_packages/)   
  Last 100 packages added to devel branch of Bioconductor are listed here. 

<h2 id="maintenance">Resources for Package Maintenance</h2>

* [Source Control](/developers/how-to/source-control/)    
  Packages contributed to _Bioconductor_ are added to the project
  subversion repository. Contributing developers are then 
  provided with password protected write access to their 
  portion of the codebase. Anonymous access to
  this repository is also available.

* [Subversion Log](/developers/svnlog/)   
  Recent commits to the development branch of _Bioconductor_ subversion
  repository. 

* [Build Reports](/checkResults/)   
  Changes in R and _Bioconductor_ could result in the malfunction of
  software packages. Therefore, package maintainers should periodically
  check their packages to ensure that they are still working as
  expected. To facilitate this, _Bioconductor_ maintains a daily build
  system to check that every package in the release and development 
  branches can be built via `R CMD build` and checked via `R CMD check`.

* [Release Schedule](/developers/release-schedule/)   
  Bioconductor has two releases in a year and package developers
  are expected to make sure that their package passes build and check without
  any errors and warnings before each release. Packages that fail to
  pass these checks will regretfully be dropped from the next release
  of _Bioconductor_.

* [Developer Mailing List](https://stat.ethz.ch/mailman/listinfo/bioc-devel)    
  The developer mailing list facilitates communication amongst _Bioconductor_ 
  developers. 

* [Using Bioconductor Git Mirrors](/developers/how-to/git-mirrors/)    
  Use Git and Github if you prefer not to use Subversion.

* [RSS Feeds](/developers/rss-feeds/)   
  _Bioconductor_ build system updates RSS feeds for each package 
  if there were any issues (warnings, errors, timeouts) with 
  the package build, in both release and devel. Package maintainers 
  are encouraged to subscribe to these feeds to be 
  notified immediately if there are any build problems.

* [End-of-life for a Bioconductor Package](/developers/package-end-of-life/)   
  You may choose to stop supporting your package at _Bioconductor_. Please
  follow these instructions to deprecate your package. 

<h2 id="contact">Contact us</h2>

For questions / concerns, please contact us at 
email: bioc-devel NEAR r-project POINT org





