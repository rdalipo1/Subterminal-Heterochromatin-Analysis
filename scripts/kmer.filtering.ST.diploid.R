

args<- commandArgs(trailingOnly = TRUE)

tukey_tbl <- read.table(args[1], row.names = 1)

norm_tbl <- read.table(args[2], row.names = 1, header = T)

out_suffix <- args[3]


ST <- c()
IT <- c()


for(kmer in rownames(tukey_tbl)){
    if(tukey_tbl[kmer,13] < 0.05){ # ST - IT Bonf pval
        if(tukey_tbl[kmer,9] > 0){ # ST - IT difference if > then ST
            ST <- append(ST, kmer)
        }
        if(tukey_tbl[kmer,9] < 0){ # ST - IT difference if < then IT
            IT <- append(IT, kmer)
        }
    }
}


write.table(ST, paste("ST.", out_suffix,".txt", sep = ""), col.names = F, row.names = F, quote = F)
write.table(IT, paste("IT.", out_suffix,".txt", sep = ""), col.names = F, row.names = F, quote = F)

