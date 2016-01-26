# Solr Backup and Restore Procedures

For the purposes of recovery & development, we need a backup of our Solr data. 
One way to backup the Solr data is described below [^1].

#### Initial Steps
Since some of the operations (e.g. using the web / admin UI) are simpler via http://localhost:8983/solr/ , I recommend using an ssh tunnel from your local
machine to the master node.
```
# Connect and bind port 8983 
# The port on the left represents our local machine
# The host:port combination on the right is in the context of master
ssh webadmin@master.bioconductor.org -L8983:localhost:8983
```
Inside this ssh connection, everything behaves as you'd expect.  One thing
you'll need is a directory for Solr to store / retrieve it's backup data:
```
# If it doesn't already exist, create a backup directory
mkdir ~/SOLR-BACKUP-DIR
```

#### Backup
This will backup the `default` core in `/home/webadmin/SOLR-BACKUP-DIR`. 
It's important to provide the `location` and `name` arguments.

```
# Invoke the backup endpoint
curl 'http://localhost:8983/solr/default/replication?command=backup&location=/home/webadmin/SOLR-BACKUP-DIR&name=2016-JAN-26'
```
You can view the status of the backup by executing the following command.  The 
relevant information will be at the bottom of the output :
```
curl 'http://localhost:8983/solr/default/replication?command=details&indent=on'
```
In my test, this creates a directory `/home/webadmin/SOLR-BACKUP-DIR/snapshot.2016-JAN-26` on master.

#### Restoration
We'll restore the `default` core to a new core called `restore-test`.  To do so, we need to create it:
```
# Inside the master node, use the Solr executable 
$HOME/solr-5.2.1/bin/solr create -c restore-test
```
If that succeeds, you should see output like this :
```
Setup new core instance directory:
/home/webadmin/solr-5.2.1/server/solr/restore-test

Creating new core 'restore-test' using command:
http://localhost:8983/solr/admin/cores?action=CREATE&name=restore-test&instanceDir=restore-test

{
  "responseHeader":{
    "status":0,
    "QTime":470},
  "core":"restore-test"}

```
Before we perform the restore, we can do a sanity check on our core : 

```
curl 'http://localhost:8983/solr/restore-test/replication?command=restorestatus&indent=on'
<?xml version="1.0" encoding="UTF-8"?>
<response>

<lst name="responseHeader">
  <int name="status">0</int>
  <int name="QTime">0</int>
</lst>
<lst name="restorestatus">
  <str name="status">No restore actions in progress</str>
</lst>
</response>
```


Next, you'll restore the previous backup (`2016-JAN-26`) to our new
core (`restore-test`).  Remember to provide the `location` and `name` arguments.  You'll 
receive feedback immediately :
```
curl 'http://localhost:8983/solr/restore-test/replication?command=restore&location=/home/webadmin/SOLR-BACKUP-DIR&name=2016-JAN-26&indent=on'
<?xml version="1.0" encoding="UTF-8"?>
<response>

<lst name="responseHeader">
  <int name="status">0</int>
  <int name="QTime">0</int>
</lst>
<str name="status">OK</str>
</response>
```

The status of the restoration can be viewed by navigating to http://localhost:8983/solr/restore-test/replication?command=restorestatus&indent=on .  Or, 
using curl as we did before : 
```
curl 'http://localhost:8983/solr/restore-test/replication?command=restorestatus&indent=on'
<?xml version="1.0" encoding="UTF-8"?>
<response>

<lst name="responseHeader">
  <int name="status">0</int>
  <int name="QTime">0</int>
</lst>
<lst name="restorestatus">
  <str name="snapshotName">snapshot.2016-JAN-26</str>
  <str name="status">In Progress</str>
</lst>
</response>
```
Once the restoration is complete, the status will indicate success : 
```
curl 'http://localhost:8983/solr/restore-test/replication?command=restorestatus&indent=on'
<?xml version="1.0" encoding="UTF-8"?>
<response>

<lst name="responseHeader">
  <int name="status">0</int>
  <int name="QTime">0</int>
</lst>
<lst name="restorestatus">
  <str name="snapshotName">snapshot.2016-JAN-26</str>
  <str name="status">success</str>
</lst>
</response>
```


[^1]: Based on Solr backup documentation: https://cwiki.apache.org/confluence/display/solr/Making+and+Restoring+Backups+of+SolrCores
