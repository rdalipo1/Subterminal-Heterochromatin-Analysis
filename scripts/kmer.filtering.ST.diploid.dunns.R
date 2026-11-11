

args<- commandArgs(trailingOnly = TRUE)

dunns_tbl <- read.table(args[1], row.names = 1)

norm_tbl <- read.table(args[2], row.names = 1, header = T)

out_suffix <- args[3]


ST <- c()
IT <- c()


for(kmer in rownames(dunns_tbl)){
    if(dunns_tbl[kmer,6] < 0.05){ # IT - ST Bonf pval
        if(dunns_tbl[kmer,4] < 0){ # IT - ST z-score if < 0 then ST
            ST <- append(ST, kmer)
        }
        if(dunns_tbl[kmer,4] > 0){ # IT - ST z-score if > 0 then IT
            IT <- append(IT, kmer)
        }
    }
}


write.table(ST, paste("ST.", out_suffix,".txt", sep = ""), col.names = F, row.names = F, quote = F)
write.table(IT, paste("IT.", out_suffix,".txt", sep = ""), col.names = F, row.names = F, quote = F)

