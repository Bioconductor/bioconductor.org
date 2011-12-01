#!/usr/bin/env sh

cd /srv/www/htdocs

rsync -Have ssh --delete --checksum \
 --exclude='.svn/' \
 --exclude=responses/ \
 --exclude=schema/ \
 --exclude=scripts/ \
 --exclude=docs/ \
 --exclude=include/admin.inc \
 --exclude=include/attendees.inc \
 --exclude=admin.php  \
 --exclude=report.php \
 --exclude=batchsheet.php \
 --exclude=attendees.php \
/srv/www/htdocs/mailform/ webadmin@mamba:/extra/www/event_reg/mailform

