# Deleting a Solr Core

## This is very dangerous!
Be sure you're deleting the correct core before proceeding.

#### Initial Steps
Since some of the operations (e.g. using the web / admin UI) are simpler via http://localhost:8983/solr/ , I recommend using an ssh tunnel from your local
machine to the master node.
```
# Connect and bind port 8983 
# The port on the left represents our local machine
# The host:port combination on the right is in the context of master
ssh webadmin@master.bioconductor.org -L8983:localhost:8983
```
Inside this ssh connection, everything behaves as you'd expect.  

#### Delete a named core
Inside the ssh session, assuming you want to delete a core named `restore-test` , do the following : 
```
$HOME/solr-5.2.1/bin/solr delete -c restore-test
```
