# Subterminal-Heterochromatin-Analysis


### Currently a work in progress...


## Before Starting Set-Up and Required Programs

### Required Programs

- Perl
- Jellyfish
- python3
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

    **Note:** For our analysis, we did not use acrocentric chromosomes for Steps 1 and 2 below. For chimp and gorilla, you can manually remove acrocentrics as         they are described in the supplemental figure II.8 of [Yoo D., et al (2025)](https://www.nature.com/articles/s41586-025-08816-3). We still attempt to call         subterminal regions in these chromosomes, but we do not use k-mer counts in these regions to differentiate subterminal and interstitial k-mers statistically. 
    


### Step 1. K-mer counting and filtering

- To start, we are using the `SubTeloJelly.variable.pl` perl script which calls Jellyfish and blastdbcmd. The arguments are as follows

  - **Argument 1**: Genome(chromosome) size file (instructions on how to make it can be found in **Genome Set Up** section above
  - **Argument 2**: Blast generated database
  - **Argument 3**: Total genome size. This doesn't need to be exact. I used 3G for primates and *X. laevis* and 2G for *X. tropicalis*
  - **Argument 4**: K-mer size. For *Xenopus* we used k=13 and for primates it was k=7
  - **Argument 5**: Subterminal Heterochromatin Search Space Size. This argument is key as the rest of the analysis heavily depends on the estimated size of the                       search space. For primates we used 5MB. For *Xenopus* we used 20MB (Arguments need to be full numbers). (Details for why are in the paper)
  - **Argument 6**: Label for search space. 5M and 20M in our case
    
- Example for chimp:

    `perl ~/scripts/SubTeloJelly.variable.pl ~/chimp.genome/t2t/chimp.t2t.chr.length.noacros.txt ~/chimp.genome/t2t/chimp.t2t.db 3000000000 13 5000000 5MB &`

  - This job may take a long time and you should use `srun` to not overload your clusters headnode
 
  - This script produces 3 types of files. The sorted. files will have k-mer counts by chromosome and region. We want to merge each file that starts with sorted.      including the file headers. Using the script `multi_join.sh` in this way creates the `STcount.merged.txt` file:
 
      `./multi_join.sh sorted.Chr* > STcount.merged.txt`

  - The resulting file should have every unique k-mer (row) count for each chromosomal region (column). The regions are defined by the 5th argument in the script             `SubTeloJelly.variable.pl`, so for each chromosome, there should be the 5p, 3p and interstitial counts as columns.
    
    #### K-mer filtering
  - Your choice of k (Argument 4) will determine how many unique k-mers are in the file. Higher k means more unique k-mers. This means that there will be some         counts in this table that will be very low or 0, making it difficult to run statisticall tests on them.
 
  - The script `filter_JellyRows.pl` takes the `STcount.merged.txt` file as the first argument and an integer that defines the minimum required count for a k-mer      in each region. Example usage:

      `perl ./filter_JellyRows.pl STcount.merged.txt 100 > STcount.merged.filtered.txt`

      **Note:** This script removes the header from `STcount.merged.txt` You will need to paste it from the original in an extra step.
    
  - In the paper, we used k=7 for the primates and k=13 for *Xenopus*. K=7 produces 8192 unique k-mers none of which require filtering (at least for the genomes       used in this paper). K=13 produces much more unique k-mers, depending on the genome, thus we filtered the count by a minimum of 100 as shown in the example        above.
 
 ### Step 2. Running ANOVA or Dunn's tests | Creating norm.tbl and tukey.tbl or dunn.tbl

  - In this step, we are running the applicable statistical test on each k-mer to ask for its ability to differentiate between the designated 5p and 3p terminal       ends and the interstitial region of chromosomes. Given the normal distributions 13-mers produce in *Xenopus*, we use ANOVA + TukeyHSD. For the non-normal          distributions of 7-mers in primates, we use Kruscal-Wallis + Dunn's

  - The scripts `Kmap.diploid.ST.dunns.R`, `Kmap.diploid.ST.R` and `Kmap.tetraploid.ST.R` are used for this analysis depending on the species and k. They all          generally do the same thing: normalizing counts based on the size of each chromosome and creating the `norm.tbl file`, and applying either ANOVA or Kruskal-       Wallis tests followed by TukeyHSD or Dunn's tests respecitivelly, creating the `tukey.tbl` or `dunn.tbl` files.


    #### Parallelizing Step 2. for sbatch job submissions

  ### Step 3. Finding k-mer locations along each chromosome

  ### Step 4. Determining emission and applying HMM

  
