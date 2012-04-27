# caBIG&reg;

[caBIG](https://cabig.nci.nih.gov/)&reg; is a National Cancer Institute
initiative to facilitate data sharing between members of the cancer research
community.  [caGrid](http://cagrid.org/) is infrastructure to support
data and analysis 'services'.

## Work flows

Several Bioconductor work flows are available as caGrid analytic services;
instances of these services are discoverable on the
[caGrid portal](http://cagrid-portal.nci.nih.gov) as those services that are
hosted at the Fred Hutchinson Cancer Research Center here in Seattle.

Service documentation is
[available](https://gforge.nci.nih.gov/docman/?group_id=175) under '02
Bioconductor caGrid Services'. The functionality provided by the services is
documented in the Service Documentation archive; the Installation Guide
provides instruction on establishing your own instances.

The code for the services is in
[gforge](https://gforge.nci.nih.gov/scm/?group_id=175).

The
[RWebServices](/packages/release/bioc/html/RWebServices.html)
package has a vignette on enabling web services on additional packages,
including how these are combined with the caGrid `introduce` infrastructure,
to create custom caGrid services from Bioconductor scripts.

## Using services

caGrid services are SOAP-based, and are most easily accessed through a
java-based client built using caGrid `introduce` infrastructure. Each service
in gforge contains a demonstration client. Our local instance is no longer
available.
