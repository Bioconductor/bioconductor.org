# Remove Large Data Files and Clean Git Tree

__Goal:__ Git remembers. Sometimes large data files are added to git
repository (intentionally or unintentionally) causing the size of the repository
to become large. It is necessary to remove the files and clean the git tree from
tracking in order to reduce the size.

These steps should be run on your local copy and (if necessary) pushed to your
own github repository. The steps below assume origin is a user maintained github
repository.

**NOTE:** Anyone that is maintaining the package repository (with a local copy)
  should run steps 1-3.


1. Download [BFG Repo-Cleaner](https://rtyley.github.io/bfg-repo-cleaner/)

2. Run BFG Repo-Cleaner on your package directory

   In the location of your package, run the following command
   ```
   java -jar <path to download>/bfg-1.13.0.jar --strip-blobs-bigger-than 100M <your
   package>
   ```
   **Note:** The above command would remove any file that is 100Mb or
   larger. Adjust this argument based on the size of the files you are cleaning up
   after. It should be lower than the offending file size.

3. Run clean up
   ```
   cd <your package>
   git reflog expire --expire=now --all && git gc --prune=now --aggressive
   ```

4. Push Changes
   ```
   git push -f origin
   ```

5. Request updates on the git.bioconductor.org repository location.

   The Bioconductor git server does not allow `-f` or to force push to the
   git.bioconductor.org location. Please email <bioc-devel@r-project.org>
   explaining the package has been cleaned for large data files and needs to be
   reset.
