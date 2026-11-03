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
- BLAST+


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
 
- Next is to create a database of the fa file with blast. This can be done with one line of code as long as blast is installed in the system:

  `makeblastdb -in chimp.t2t.chronly.renamed.headers.fa -dbtype nucl -parse_seqids -title chimp.t2t -out chimp.t2t.db`

  - The `-parse_seqids` will apply the renamed headers and you can use the `chimp.t2t.db` name to refer to the database in the future

- Lastly, generate a chromosome a size file. The 'genome.size.from.fa.sh' script can be used to generate this. The two arguments it needs is the fa file and the name of the resulting chr.size file:

    `./genome.size.from.fa.sh chimp.t2t.chronly.renamed.headers.fa chimp.t2t.chronly.chr.size.txt`

  - Head of an example file shown below:
    
|      |           |
|------|-----------|
| Chr1 | 231299448 |
| Chr2 | 201579334 |
| Chr3 | 191670129 |
| Chr4 | 178964043 |
| Chr5 | 183320697 |
| Chr6 | 174648847 |


### Step 1. K-mer counting and filtering

- To start, we are using the `SubTeloJelly.variable.pl` perl script which calls Jellyfish and blastdbcmd. The arguments are as follows

  - **Argument 1**: Genome size file (instructions on how to make it can be found in **Genome Set Up** section above
  - **Argument 2**: 

