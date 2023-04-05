# ![](/images/icons/magnifier.gif)Source Control

The _Bioconductor_ project is maintained in a Git source control
system. Package maintainers update their packages by pushing changes
to their git repositories.

Please see the [Git Version
Control](https://contributions.bioconductor.org/git-version-control.html) for
different scenerios and workflows from setting up remotes, transitioning from
svn, syncing github and Bioconductor git.bioconductor.org repositories, and
troubleshooting. 

## Annotation packages

Traditional Annotation packages are not stored in GIT due to the size of
annotation files. To update an existing Annotation package please send an email
to maintainer@bioconductor.org. A member of the Bioconductor team will be in
contact to receive the updated package.

Newer annotation packages can be stored in GIT as it is a requirement to use the
[AnnotationHub](https://bioconductor.org/packages/AnnotationHub/) or similar
server hosted data. The larger sized files are not included directly in the package.
To contribute a new Annotation package please contact hubs@bioconductor.org
for guidance and read the documentation on [How to Create A Hub
package](https://bioconductor.org/packages/devel/bioc/vignettes/HubPub/inst/doc/CreateAHubPackage.html).

Currently direct updates to annotation packages, even those stored on git, are
not supported. If you wish to updated an annotation package, make required
changes and push to git.bioconductor.org. Then send an email to
hubs@bioconductor.org or maintainer@bioconductor.org requesting the package be propagated.


## More help

Need more help? Ask on the [bioc-devel](/help/mailing-list/) mailing
list.
