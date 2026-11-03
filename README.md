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
  
| Species       | Name     | Sex | BioSample ID | BioProject ID |
|---------------|----------|-----|--------------|---------------|
| Chimpanzee    | AG18354  | M   | SAMN30216104 | PRJNA916736   |
| Gorilla       | Jim      | M   | SAMN04003007 | PRJNA942267   |
| X. laevis     | J_2021   | F   | SAMN18035302 | PRJNA313213   |
| X. tropicalis | Nigerian | F   | SAMN13041969 | PRJNA577946   |


- These can be obtained from either NCBI or Xenbase (add links). NCBI downloads can be done through the [datasets](https://github.com/ncbi/datasets) tool.

- Once the fa file has been downloaded, I like to rename my chromosomes to this format: `ChrXX` and remove any scaffolds since those will not be used for this analysis. Here's an for chromosome renaming example with chimp:

  *>CM054434.2 Pan troglodytes isolate AG18354 chromosome 1, whole genome shotgun sequence*

  *>Chr1*

  - This can be done with awk, sed or any other regex method of choice.

- You will then need to generate a chromosome size file. You can use the 'genome.size.from.fa.sh' script to generate this for you. The two arguments it needs is the fa file and the name of the resulting chr.size file:

    `./genome.size.from.fa.sh chimp.t2t.chronly.renamed.headers.fa chimp.t2t.chronly.chr.size.txt`

  - Head of an example file shown below:

| Chr1    | 231299448 |
| Chr2    | 201579334 |
| Chr3    | 191670129 |
| Chr4    | 178964043 |
| Chr5    | 183320697 |
| Chr6    | 174648847 |




