# CSAMA 2017: Statistical Data Analysis for Genome-Scale Biology

June 11-16, 2017<br />
Bressanone-Brixen, Italy<br />
URL: [http://www.huber.embl.de/csama2017/](2017/CSAMA/http://www.huber.embl.de/csama2017/)

Lecturers: Jennifer Bryan, RStudio and UBC; Vincent J. Carey, Harvard
Medical School; Laurent Gatto, University of Cambridge; Wolfgang
Huber, European Molecular Biology Laboratory (EMBL), Heidelberg;
Martin Morgan, Roswell Park Cancer Institute, Buffalo; Johannes
Rainer, European Academy of Bozen (EURAC); Charlotte Soneson,
University of Zurich; Levi Waldron, CUNY School of Public Health at
Hunter College, New York.

Teaching Assistants: Simone Bell, EMBL, Heidelberg; Vladislav Kim,
EMBL, Heidelberg; Lori Shepherd, RPCI, Buffalo; Mike L. Smith, EMBL,
Heidelberg.

## Resources

Source: [Github](https://github.com/Bioconductor/CSAMA/blob/2017/README.md)

**Monday, June 12**

Lectures

- Introduction to R and Bioconductor
  ([html](lectures/1-monday/lecture-01-r-intro/lecture-01-r-intro.html)).
- Computing with Sequences and Ranges
  ([html](lectures/1-monday/lecture-02-bioc-intro/lecture-02-bioc-intro.html)).
- Tabular data management
  ([pdf](lectures/1-monday/lecture-03-rectangular-data/lect03_bryan_tabular-data-management.pdf)).
- Annotation resources
  ([html](lectures/1-monday/lecture-04-a-annotation-intro/lecture-04a-annotation-intro.html));
  EnsemblDb ([pdf](lectures/1-monday/lecture-04b-ensembldb/04-AnnotationResources-ensembldb-beamer.pdf)).

Labs

- Introduction to Bioc. Introduction
  ([html](labs/1-monday/lab-01-intro-to-r-bioc/L1.1-r-intro-morgan.html)),
  _R_ [html](labs/1-monday/lab-01-intro-to-r-bioc/lab-1-intro-to-r-bioc.html)),
  _Bioconductor_
  ([html](labs/1-monday/lab-01-intro-to-r-bioc/L1.2-bioc-intro-morgan.html)),
  Data representations
  ([html](labs/1-monday/lab-01-intro-to-r-bioc/L1.3-bioc-data-representation-morgan.html)),
  Annotation
  ([html](labs/1-monday/lab-01-intro-to-r-bioc/L1.4-bioc-annotation-morgan.html))
  (data:
  [ALLphenoData.tsv](labs/1-monday/lab-01-intro-to-r-bioc/ALLphenoData.tsv),
  [BRFSS-subset.csv](labs/1-monday/lab-01-intro-to-r-bioc/BRFSS-subset.csv),
  [symgo.csv](labs/1-monday/lab-01-intro-to-r-bioc/symgo.csv)).
- Use of Git and GitHub with R, RStudio, and R Markdown
  ([html](labs/1-monday/lab-02-git-github-r-rmd-rstudio/README.html))

**Tuesday, June 13**

Lectures

- Basics of sequence alignment and aligners
  ([pdf](lectures/2-tuesday/lec05-alignmentbased-rnaseq.pdf)).
- RNA-Seq data analysis and differential expression
  ([pdf](lectures/2-tuesday/lec06-deseq2-huber-compressed.pdf)).
- New workflows for RNA-seq
  ([pdf](lectures/2-tuesday/lec07-alignmentfree-rnaseq.pdf)).
- Hypothesis testing
  ([pdf](lectures/2-tuesday/lec08-testing-huber-compressed.pdf)).

Labs

- End-to-end RNA-Seq workflow
  ([html](labs/2-tuesday/lab-03-rnaseq/rnaseqGene_CSAMA2017.html))
- Independent hypothesis weighting
  ([html](labs/2-tuesday/lab-03b-ihw/IHW.html))

**Wednesday, June 14**

Lectures

- Multiple testing
  ([txt](lectures/3-wednesday/lecture09-multiple-testing/lecture09-readme.txt)).
- Linear models (basic intro))
  ([html](lectures/3-wednesday/lecture10-linear-models-basic-intro/Waldron_linearmodels.html)).
- Experimental design, batch effects and confounding
  ([pdf](lectures/3-wednesday/lecture11-experimental-design/lec11-experimental-design.pdf)).
- Robust statistics: median, MAD, rank test, Spearman, robust linear model
  ([html](lectures/3-wednesday/lecture12-robust-statistics/robust-statistics.html)).

**Thursday, June 15**

Lectures

- Visualization, the grammar of graphics and ggplot2
  ([pdf](lectures/4-thursday/lecture-13-graphics-and-visualisation/lecture-13-viz-huber-compressed.pdf)).
- Mass spec proteomics
  ([html](lectures/4-thursday/lecture-14a-proteomics/lab-14a-proteomics.html))
  and metabolomics
  ([pdf](lectures/4-thursday/lecture-14b-Metabolomics/Metabolomics-beamer.pdf)).
- Clustering and classification
  ([html](lectures/4-thursday/lecture-15-clust/lect-15-clust-class.html)).
- Resampling: cross-validation, bootstrap, and permutation tests
  ([html](lectures/4-thursday/lecture-16-resampling/Waldron_CSAMA2017_resampling.html)).
- Analysis of microbiome marker gene data
  ([pdf](lectures/4-thursday/lecture-17-microbial-genomics/lecture-17-microbiome.pdf)).

Labs

- Mass spec proteomics & metabolomics
  Proteomics ([html](labs/4-thursday/lab-04-Mass_spec_proteomics_and_metabolomics/01-proteomics/lab.html)),
  Metabolomics ([html](labs/4-thursday/lab-04-Mass_spec_proteomics_and_metabolomics/02-metabolomics-preprocessing/metabolomics-preprocessing.html)).
- MultiAssayExperiment
  ([html](labs/4-thursday/lab-05-MultiAssayExperiment/MultiAssayExperiment-lab.html)),
  cheatsheet
  ([pdf](labs/4-thursday/lab-05-MultiAssayExperiment/MultiAssayExperiment_cheatsheet.pdf)).

**Friday, June 16**

Lectures

- Gene set enrichment analysis
  ([pdf](lectures/5-friday/lecture-18-gene-set-enrichment/lecture-18-gene-set-enrichment.pdf)).
- Working with large-scale
  ([pdf](lectures/5-friday/lecture-19-big-data/lecture-19-big-data.pdf))
  and remote
  ([html](lectures/5-friday/lecture-19a-remote-data/lecture-19a-oom.html))
  data.
- Developer Practices, Writing functions
  ([html](lectures/5-friday/lecture-20-writing-functions/README.html)).
- Developer Practices, Writing packages
  ([html](lectures/5-friday/lecture-21-writing-packages/README.html)).

Labs

- Graphics ([pdf](labs/5-friday/lab-06-graphics/lab6-graphics.pdf))
  (data: [diabetes.csv](labs/5-friday/lab-06-graphics/diabetes.csv)).
- Machine learning (supervised)
  ([pdf](labs/5-friday/lab-07-machine-learning/lab7-machine-learning.pdf)).
- Large data, performance, and parallelization; large-scale efficient
  computation with genomic intervals
  ([html](labs/5-friday/lab-08-efficient-and-parallel-r/L8.1-efficient-and-parallel-r.html)).
