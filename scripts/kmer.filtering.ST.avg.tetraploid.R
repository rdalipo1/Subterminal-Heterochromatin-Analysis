
args<- commandArgs(trailingOnly = TRUE)

tukey_tbl <- read.table(args[1], row.names = 1)

norm_tbl <- as.matrix(read.table(args[2], row.names = 1, header = T))

out_suffix <- args[3]

ST_colnames <- colnames(norm_tbl)[!grepl("IT$", colnames(norm_tbl))]
L_colnames <- ST_colnames[grepl("L", ST_colnames)]
S_colnames <- ST_colnames[!grepl("L", ST_colnames)]

L.ST <- c()
S.ST <- c()
ST.Shared <- c()
L.IT <- c()
S.IT <- c()
IT.Shared <- c()

L.ST.contaminant <- c()
S.ST.contaminant <- c()

for(kmer in rownames(tukey_tbl)){
    ####        LST Specific 
    if( tukey_tbl[kmer,13] < 0.05 & tukey_tbl[kmer,9] > 0 &  # L.ST - L.IT
        tukey_tbl[kmer,18] > 0.05 &                          # S.ST - L.IT
        tukey_tbl[kmer,23] > 0.05 &                          # S.ST - S.IT
        tukey_tbl[kmer,38] < 0.05 & tukey_tbl[kmer,34] < 0 & # S.ST - L.ST
        tukey_tbl[kmer,33] < 0.05 & tukey_tbl[kmer,29] < 0   # S.IT - L.ST
    ){
        if (mean(norm_tbl[kmer,S_colnames]) > 1e-6){
            ### KMER is contaminating L specific
            L.ST.contaminant <- append(L.ST.contaminant, kmer)
            }
        else{
            L.ST <- append(L.ST, kmer)
            }
        }
    ###         SST Specific
    if( tukey_tbl[kmer,13] > 0.05 &                          # L.ST - L.IT
        tukey_tbl[kmer,18] < 0.05 & tukey_tbl[kmer,14] > 0 & # S.ST - L.IT
        tukey_tbl[kmer,23] < 0.05 & tukey_tbl[kmer,19] > 0 & # S.ST - S.IT
        tukey_tbl[kmer,38] < 0.05 & tukey_tbl[kmer,34] > 0 & # S.ST - L.ST
        tukey_tbl[kmer,33] > 0.05                            # S.IT - L.ST
    ){
        if (mean(norm_tbl[kmer,L_colnames]) > 1e-6){
            ### KMER is contaminating S specific
            S.ST.contaminant <- append(S.ST.contaminant, kmer)
            }
        else{ 
            S.ST <- append(S.ST, kmer)
            }
        }

    ###         ST Shared
    if( tukey_tbl[kmer,13] < 0.05 & tukey_tbl[kmer,9]  > 0 & # L.ST - L.IT
        tukey_tbl[kmer,18] < 0.05 & tukey_tbl[kmer,14] > 0 & # S.ST - L.IT
        tukey_tbl[kmer,23] < 0.05 & tukey_tbl[kmer,19] > 0 & # S.ST - S.IT
        tukey_tbl[kmer,38] > 0.05 &                          # S.ST - L.ST
        tukey_tbl[kmer,33] < 0.05 & tukey_tbl[kmer,29] < 0   # S.IT - L.ST
    ){
        ST.Shared <- append(ST.Shared, kmer)
        }

    ###         L IT
    if( tukey_tbl[kmer,13] < 0.05 & tukey_tbl[kmer,9]  < 0 & # L.ST - L.IT
        tukey_tbl[kmer,18] < 0.05 & tukey_tbl[kmer,14] < 0 & # S.ST - L.IT
        tukey_tbl[kmer,23] > 0.05 &                          # S.ST - S.IT
        tukey_tbl[kmer,38] > 0.05 &                          # S.ST - L.ST
        tukey_tbl[kmer,33] > 0.05 &                          # S.IT - L.ST
        tukey_tbl[kmer,28] < 0.05 & tukey_tbl[kmer,24] < 0   # S.IT - L.IT
    ){
        L.IT <- append(L.IT, kmer)

        }

    ###         S IT
    if( tukey_tbl[kmer,13] > 0.05 &                          # L.ST - L.IT
        tukey_tbl[kmer,18] > 0.05 &                          # S.ST - L.IT
        tukey_tbl[kmer,23] < 0.05 & tukey_tbl[kmer,19] < 0 & # S.ST - S.IT
        tukey_tbl[kmer,38] > 0.05 &                          # S.ST - L.ST
        tukey_tbl[kmer,33] < 0.05 & tukey_tbl[kmer,29] > 0 & # S.IT - L.ST
        tukey_tbl[kmer,28] < 0.05 & tukey_tbl[kmer,24] > 0   # S.IT - L.IT
    ){
        S.IT <- append(S.IT, kmer)

        }
    ###         Shared IT
    if( tukey_tbl[kmer,13] < 0.05 & tukey_tbl[kmer,9] < 0  & # L.ST - L.IT
        tukey_tbl[kmer,18] < 0.05 & tukey_tbl[kmer,14] < 0 & # S.ST - L.IT
        tukey_tbl[kmer,23] < 0.05 & tukey_tbl[kmer,19] < 0 & # S.ST - S.IT
        tukey_tbl[kmer,38] > 0.05 &                          # S.ST - L.ST
        tukey_tbl[kmer,33] < 0.05 & tukey_tbl[kmer,29] > 0 & # S.IT - L.ST
        tukey_tbl[kmer,28] > 0.05                            # S.IT - L.IT
    ){
        IT.Shared <- append(IT.Shared, kmer)

        }
}


write.table(L.ST, paste("L.ST.", out_suffix,".txt", sep = ""), col.names = F, row.names = F, quote = F)
write.table(S.ST, paste("S.ST.", out_suffix,".txt", sep = ""), col.names = F, row.names = F, quote = F)
write.table(L.IT, paste("L.IT.", out_suffix,".txt", sep = ""), col.names = F, row.names = F, quote = F)
write.table(S.IT, paste("S.IT.", out_suffix,".txt", sep = ""), col.names = F, row.names = F, quote = F)
write.table(ST.Shared, paste("ST.Shared.", out_suffix,".txt", sep = ""), col.names = F, row.names = F, quote = F)
write.table(IT.Shared, paste("IT.Shared.", out_suffix,".txt", sep = ""), col.names = F, row.names = F, quote = F)
write.table(L.ST.contaminant, paste("LST.contaminant.", out_suffix,".txt", sep = ""), col.names = F, row.names = F, quote = F)
write.table(S.ST.contaminant, paste("SST.contaminant.", out_suffix,".txt", sep = ""), col.names = F, row.names = F, quote = F)

