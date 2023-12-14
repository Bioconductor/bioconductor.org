<h1 class="format-bold">📖 Bioconductor Publications</h1>
<div class="page-container">
  <div class="upper-container">
    <section class="page-intro">
      <h2>Publications</h2>
      <p class="text-large">
        This is a select list of featured and recent publications that cite
	Bioconductor, along with select books, book chapters, and original
	Bioconductor technical reports. 
      </p>  
      <div class="sidebar-box">
      <div class="sidebar-box-header">
      <p class="format-bold">Literature Search</p>
      </div>
      <nav>
        <div>
          <a
            href="http://www.ncbi.nlm.nih.gov/pubmed/?term=bioconductor&sort=date"
            >PubMed</a
          >
        </div>
        <div>
          <a href="http://www.ncbi.nlm.nih.gov/pmc/?term=bioconductor&sort=ePubDate"
            >PubMed Central</a
          >
        </div>
        <div>
          <a href="http://scholar.google.com/scholar?q=bioconductor&btnG=Search"
            >Google Scholar</a
          >
        </div>
      </nav>
      </div>

    </section>
    <div class="upper-content">
      <section>
         <h4>Featured</h4>
       
       <div class="white_box">
       Alexander TA, Irizarry RA, Bravo HC. <a
       href="https://doi.org/10.1093/biostatistics/kxab030">Capturing discrete latent
       structures: choose LDs over PCs.</a> Biostatistics. 2021 Sep 1:kxab030. doi:
       10.1093/biostatistics/kxab030. Epub ahead of print. PMID: 34467372.
       </div>

       <div class="grey_box">
       Oh S, Geistlinger L, Ramos M, Blankenberg D, van den Beek M, Taroni JN, Carey
       VJ, Greene CS, Waldron L, Davis S. <a
       href="https://www.nature.com/articles/s41467-022-31411-3">GenomicSuperSignature facilitates
       interpretation of RNA-seq experiments through robust, efficient comparison to
       public databases.</a> Nat Commun. 2022 Jun 27;13(1):3695. doi:
       10.1038/s41467-022-31411-3. PMID: 35760813; PMCID: PMC9237024
       </div>

       <div class="white_box">
       Schatz MC, Philippakis AA, Afgan E, Banks E, Carey VJ, Carroll RJ, Culotti A,
       Ellrott K, Goecks J, Grossman RL, Hall IM, Hansen KD, Lawson J, Leek JT, Luria
       AO, Mosher S, Morgan M, Nekrutenko A, O'Connor BD, Osborn K, Paten B, Patterson
       C, Tan FJ, Taylor CO, Vessio J, Waldron L, Wang T, Wuichet K. <a
       href="https://www.sciencedirect.com/science/article/pii/S2666979X21001063?via%3Dihub">Inverting the
       model of genomics data sharing with the NHGRI Genomic Data Science Analysis,
       Visualization, and Informatics Lab-space.</a> Cell Genom. 2022 Jan
       12;2(1):100085. doi: 10.1016/j.xgen.2021.100085. Epub 2022 Jan 13. PMID:
       35199087; PMCID: PMC8863334
       </div>

       <div class="grey_box">
       He D, Zakeri M, Sarkar H, Soneson C, Srivastava A, Patro R. <a
       href="https://www.nature.com/articles/s41592-022-01408-3">Alevin-fry unlocks
       rapid, accurate and memory-frugal quantification of single-cell RNA-seq
       data.</a> Nat Methods. 2022 Mar;19(3):316-322. doi: 10.1038/s41592-022-01408-3. Epub
       2022 Mar 11. PMID: 35277707; PMCID: PMC8933848
       </div>

      </section>
    </div>
  </div>
  <hr class="page-divider" />
  <div class="lower-container">
    <div class="lower-content">
      <section>
          <h4>Recent</h4>

	  10 most recent PubMed and PubMed Central citations mentioning &quot;*Bioconductor*&quot;.
	  For a complete list, go to 
	  <a href="http://www.ncbi.nlm.nih.gov/pubmed/?term=bioconductor" target="_blank">PubMed</a> and
	  <a href="http://www.ncbi.nlm.nih.gov/pmc/?term=bioconductor&sort=ePubDate" target="_blank">PubMed Central</a>.
	  Last updated <abbr class="timeago" title="<%=get_pubmed_cache_date()%>"><%=get_pubmed_cache_date()%></abbr>.

	  <p></p>
	  <%= render('/_pubmed/', :pb_item => find_item(@items, '/help/publications/papers/pubmed/')) %>

      </section>

      <section>
           <h4>Books</h4>
	   A selection of Bioconductor releated books, book chapters, and online
           books.

	   <p><br></p>
	   
           <h6>Books</h6>

	   <a
	   href="books/bioinformatics-and-computational-biology-solutions/"
	   title="Landmark Bioconductor Book">
	   <img src="books/bioinformatics-and-computational-biology-solutions/bcbs_small.jpg" height="110"
	   alt="Landmark Bioconductor Book" width="85" /></a>&nbsp;&nbsp;&nbsp;&nbsp;

	   <a
	   href="books/r-programming-for-bioinformatics/" title="R Programming for Bioinformatics">
	   <img
	   src="books/r-programming-for-bioinformatics/rbioinf-small.jpg"
	   height="110" alt="R Programming for Bioinformatics" width="85" /></a>&nbsp;&nbsp;&nbsp;&nbsp;

	   <a href="books/bioconductor-case-studies/" title="Bioconductor Case Studies">
	   <img
	   src="books/bioconductor-case-studies/case-studies-small.jpg"
	   height="110" alt="Bioconductor Case Studies" width="85" /></a>

	   <p><br></p>

           <h6>Book Chapters</h6>

	   <a href="book-chapters/MiMB">Analyzing biological data using R: methods for graphs and networks</a>

	   <p><br></p>
	   
	   <h6>Online Books</h6>

	   <div class="white_box">
	   <a href="https://www.huber.embl.de/msmb/">Modern Statistics for Modern Biology</a> 
	   </div>

	   <div class="grey_box">
	   <a href="/books/release/OSCA">Orchestrating Single-Cell Analysis with Bioconductor</a>
	   </div>

	   <div class="white_box">
	   <a href="https://jserizay.com/OHCA/">Orchestrating Hi-C Analysis with Bioconductor</a>
	   </div>

	   <div class="grey_box">
	   <a href="https://rformassspectrometry.github.io/book/">R for Mass Spectrometry</a>
	   </div>

	   <div class="white_box">
	   <a href="/books/release/SingleRBook">SingleRBook</a>
	   </div>

	   <div class="grey_box">
	   <a href="/books/release/csawBook">csawBook</a>
	   </div>

	   <p><br></p>
	   
           <h4>Other</h4>
	   F1000 Channel,Technical Reports, and Compendia
	   
	   <p><br></p>
	   
           <h6>F1000 Channel</h6>
	   
 	   <a href="http://f1000research.com/channels/bioconductor">F1000 Research Channel</a>

  	   <p><br></p>

	   <h6>Initial Bioconductor Technical Reports</h6>

	   Bioconductor Core (November 2002) <a href="tech-reports/relProjTR.pdf"><em>An Overview of Projects in Computing for Genomic Analysis</em></a>
	    <p><br></p>
	   Bioconductor Core (November 2002) <a href="tech-reports/currProgTR.pdf"><em>Bioconductor: Assessment of Current Progress</em></a>

      </section>

    </div>
  </div>

</div>