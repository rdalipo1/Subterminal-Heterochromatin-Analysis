# Subterminal-Heterochromatin-Analysis


### Currently a work in progress...


## Before Starting Set-Up and Required Programs

### Required Programs

- Perl
- Jellyfish
- R
  - [karyoploteR](https://www.bioconductor.org/packages/release/bioc/html/karyoploteR.html)
  - [HMM](https://cran.r-project.org/web/packages/HMM/index.html)
  - [dunn.test](https://cran.r-project.org/web/packages/dunn.test/index.html)
- Bedtools


### Genome Set Up

- Genomes used in the paper are shown in the supplemental table:
  
Species	      Name	    Sex	  BioSample ID	BioProject ID
Chimpanzee	  AG18354	  M	    SAMN30216104	PRJNA916736
Gorilla	      Jim	      M	    SAMN04003007	PRJNA942267
X. laevis	    J_2021	  F	    SAMN18035302	PRJNA313213
X. tropicalis	Nigerian	F	    SAMN13041969	PRJNA577946

- These can be obtained from either NCBI or Xenbase (add links). NCBI downloads can be done through the [datasets](https://github.com/ncbi/datasets) tool.

- Once the fa file has been downloaded, I like to rename my chromosomes to this format: `ChrXX` and remove any scaffolds since those will not be used for this analysis. Here's an example with chimp:

`>CM054434.2 Pan troglodytes isolate AG18354 chromosome 1, whole genome shotgun sequence`

`>Chr1`



