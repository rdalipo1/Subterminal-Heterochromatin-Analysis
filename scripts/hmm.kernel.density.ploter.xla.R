
chr_lenghts <- read.table("/data/home/rdalipo1/v10.XlaXtr.homology/Xla/Xenla10.chr.size.txt", header = FALSE, col.names = c("Chr", "Length"))

regions <- read.table("/data/home/rdalipo1/subtelo_proj/final.window.region.tsv", header = T)



pdf("kernal.density.by.region.plots.pdf", useDingbats = F)

#all L subtelo
lmer_ST_density <- density(log10(regions[regions$Category != 'IT' & grepl('L', regions$Chr), 'Lmer']+1))
smer_ST_density <- density(log10(regions[regions$Category != 'IT' & grepl('L', regions$Chr), 'Smer']+1))
LST_ST_density <- density(log10(regions[regions$Category != 'IT' & grepl('L', regions$Chr), 'L.ST']+1))
SST_ST_density <- density(log10(regions[regions$Category != 'IT' & grepl('L', regions$Chr), 'S.ST']+1))
LIT_ST_density <- density(log10(regions[regions$Category != 'IT' & grepl('L', regions$Chr), 'L.IT']+1))
SIT_ST_density <- density(log10(regions[regions$Category != 'IT' & grepl('L', regions$Chr), 'S.IT']+1))
ST.shared_ST_density <- density(log10(regions[regions$Category != 'IT' & grepl('L', regions$Chr), 'ST.shared']+1))
IT.shared_ST_density <- density(log10(regions[regions$Category != 'IT' & grepl('L', regions$Chr), 'IT.shared']+1))

plot(smer_ST_density$x, length(regions[regions$Category != 'IT' & grepl('L', regions$Chr), 'Smer']) * smer_ST_density$y, type = 'l', lty = 1, col = '#d62728', ylab = 'Number of 100kb windows', xlab = 'log10(Number of bp overlap for each kmer category)', xlim = c(-0.5, 4.5),ylim = c(0, 6000),  main = "L Subtelomere") # RED
lines(lmer_ST_density$x, length(regions[regions$Category != 'IT' & grepl('L', regions$Chr), 'Lmer']) * lmer_ST_density$y, col = '#1f77b4', lty = 1) #blue
lines(LST_ST_density$x, length(regions[regions$Category != 'IT' & grepl('L', regions$Chr), 'L.ST']) * LST_ST_density$y, col = '#1f77b4', lty = 6) #BLUE
lines(SST_ST_density$x, length(regions[regions$Category != 'IT' & grepl('L', regions$Chr), 'S.ST']) * SST_ST_density$y, col = '#d62728', lty = 6) #RED
lines(LIT_ST_density$x, length(regions[regions$Category != 'IT' & grepl('L', regions$Chr), 'L.IT']) * LIT_ST_density$y, col = '#1f77b4', lty = 3) #teal '#17becf'
lines(SIT_ST_density$x, length(regions[regions$Category != 'IT' & grepl('L', regions$Chr), 'S.IT']) * SIT_ST_density$y, col = '#d62728', lty = 3) #orange
lines(ST.shared_ST_density$x, length(regions[regions$Category != 'IT' & grepl('L', regions$Chr), 'ST.shared']) * ST.shared_ST_density$y, col = 'green', lty = 1) #green
lines(IT.shared_ST_density$x, length(regions[regions$Category != 'IT' & grepl('L', regions$Chr), 'IT.shared']) * IT.shared_ST_density$y, col = 'black', lty = 1) #black


legend("topright", legend = c("lmer", "smer", "LST", "SST", "LIT", "SIT", "ST.shared", "IT.shared"), col = c('#1f77b4', '#d62728', '#1f77b4', '#d62728', '#1f77b4', "#d62728", 'green', 'black'), lty = c(1,1,6,6,3,3,1,1), cex = 0.8)

par(mfrow = c(4,2))

plot(smer_ST_density$x, length(regions[regions$Category != 'IT' & grepl('L', regions$Chr), 'Smer']) * smer_ST_density$y, type = 'l', lty = 1, col = '#d62728', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4.5),ylim = c(0, 6000),  main = "L Subtelo | Smer")

plot(lmer_ST_density$x, length(regions[regions$Category != 'IT' & grepl('L', regions$Chr), 'Lmer']) * lmer_ST_density$y, type = 'l', lty = 1, col = '#1f77b4', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4.5),ylim = c(0, 6000),  main = "L Subtelo | Lmer")

plot(LST_ST_density$x, length(regions[regions$Category != 'IT' & grepl('L', regions$Chr), 'L.ST']) * LST_ST_density$y, type = 'l', lty = 1, col = '#1f77b4', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4.5),ylim = c(0, 6000),  main = "L Subtelo | LST")

plot(SST_ST_density$x, length(regions[regions$Category != 'IT' & grepl('L', regions$Chr), 'S.ST']) * SST_ST_density$y, type = 'l', lty = 1, col = '#d62728', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4.5),ylim = c(0, 6000),  main = "L Subtelo | SST")

plot(LIT_ST_density$x, length(regions[regions$Category != 'IT' & grepl('L', regions$Chr), 'L.ST']) * LIT_ST_density$y, type = 'l', lty = 1, col = '#1f77b4', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4.5),ylim = c(0, 6000),  main = "L Subtelo | LIT")

plot(SIT_ST_density$x, length(regions[regions$Category != 'IT' & grepl('L', regions$Chr), 'S.ST']) * SIT_ST_density$y, type = 'l', lty = 1, col = '#d62728', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4.5),ylim = c(0, 6000),  main = "L Subtelo | SIT")

plot(ST.shared_ST_density$x, length(regions[regions$Category != 'IT' & grepl('L', regions$Chr), 'ST.shared']) * ST.shared_ST_density$y, type = 'l', lty = 1, col = 'green', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4.5),ylim = c(0, 6000),  main = "L Subtelo | ST.shared")

plot(IT.shared_ST_density$x, length(regions[regions$Category != 'IT' & grepl('L', regions$Chr), 'IT.shared']) * IT.shared_ST_density$y, type = 'l', lty = 1, col = 'black', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4.5),ylim = c(0, 6000),  main = "L Subtelo | IT.shared")



#all S subtelo

par(mfrow = c(1,1))

lmer_ST_density <- density(log10(regions[regions$Category != 'IT' & grepl('S', regions$Chr), 'Lmer']+1))
smer_ST_density <- density(log10(regions[regions$Category != 'IT' & grepl('S', regions$Chr), 'Smer']+1))
LST_ST_density <- density(log10(regions[regions$Category != 'IT' & grepl('S', regions$Chr), 'L.ST']+1))
SST_ST_density <- density(log10(regions[regions$Category != 'IT' & grepl('S', regions$Chr), 'S.ST']+1))
LIT_ST_density <- density(log10(regions[regions$Category != 'IT' & grepl('S', regions$Chr), 'L.IT']+1))
SIT_ST_density <- density(log10(regions[regions$Category != 'IT' & grepl('S', regions$Chr), 'S.IT']+1))
ST.shared_ST_density <- density(log10(regions[regions$Category != 'IT' & grepl('S', regions$Chr), 'ST.shared']+1))
IT.shared_ST_density <- density(log10(regions[regions$Category != 'IT' & grepl('S', regions$Chr), 'IT.shared']+1))

plot(smer_ST_density$x, length(regions[regions$Category != 'IT' & grepl('S', regions$Chr), 'Smer']) * smer_ST_density$y, type = 'l', lty = 1, col = '#d62728', ylab = 'Number of 100kb windows', xlab = 'log10(Number of bp overlap for each kmer category)', xlim = c(-0.5, 4.5),ylim = c(0, 5000),  main = "S Subtelomere") # RED
lines(lmer_ST_density$x, length(regions[regions$Category != 'IT' & grepl('S', regions$Chr), 'Lmer']) * lmer_ST_density$y, col = '#1f77b4', lty = 1) #blue
lines(LST_ST_density$x, length(regions[regions$Category != 'IT' & grepl('S', regions$Chr), 'L.ST']) * LST_ST_density$y, col = '#1f77b4', lty = 6) #BLUE
lines(SST_ST_density$x, length(regions[regions$Category != 'IT' & grepl('S', regions$Chr), 'S.ST']) * SST_ST_density$y, col = '#d62728', lty = 6) #RED
lines(LIT_ST_density$x, length(regions[regions$Category != 'IT' & grepl('S', regions$Chr), 'L.IT']) * LIT_ST_density$y, col = '#1f77b4', lty = 3) #teal '#17becf'
lines(SIT_ST_density$x, length(regions[regions$Category != 'IT' & grepl('S', regions$Chr), 'S.IT']) * SIT_ST_density$y, col = '#d62728', lty = 3) #orange
lines(ST.shared_ST_density$x, length(regions[regions$Category != 'IT' & grepl('S', regions$Chr), 'ST.shared']) * ST.shared_ST_density$y, col = 'green', lty = 1) #green
lines(IT.shared_ST_density$x, length(regions[regions$Category != 'IT' & grepl('S', regions$Chr), 'IT.shared']) * IT.shared_ST_density$y, col = 'black', lty = 1) #black



legend("topright", legend = c("lmer", "smer", "LST", "SST", "LIT", "SIT", "ST.shared", "IT.shared"), col = c('#1f77b4', '#d62728', '#1f77b4', '#d62728', '#1f77b4', "#d62728", 'green', 'black'), lty = c(1,1,6,6,3,3,1,1), cex = 0.8)

par(mfrow = c(4,2))

plot(lmer_ST_density$x, length(regions[regions$Category != 'IT' & grepl('S', regions$Chr), 'Lmer']) * lmer_ST_density$y, type = 'l', lty = 1, col = '#1f77b4', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4.5),ylim = c(0, 5000),  main = "S Subtelo | Lmer")

plot(smer_ST_density$x, length(regions[regions$Category != 'IT' & grepl('S', regions$Chr), 'Smer']) * smer_ST_density$y, type = 'l', lty = 1, col = '#d62728', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4.5),ylim = c(0, 5000),  main = "S Subtelo | Smer")

plot(LST_ST_density$x, length(regions[regions$Category != 'IT' & grepl('S', regions$Chr), 'L.ST']) * LST_ST_density$y, type = 'l', lty = 1, col = '#1f77b4', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4.5),ylim = c(0, 5000),  main = "S Subtelo | LST")

plot(SST_ST_density$x, length(regions[regions$Category != 'IT' & grepl('S', regions$Chr), 'S.ST']) * SST_ST_density$y, type = 'l', lty = 1, col = '#d62728', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4.5),ylim = c(0, 5000),  main = "S Subtelo | SST")

plot(LIT_ST_density$x, length(regions[regions$Category != 'IT' & grepl('S', regions$Chr), 'L.ST']) * LIT_ST_density$y, type = 'l', lty = 1, col = '#1f77b4', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4.5),ylim = c(0, 5000),  main = "S Subtelo | LIT")

plot(SIT_ST_density$x, length(regions[regions$Category != 'IT' & grepl('S', regions$Chr), 'S.ST']) * SIT_ST_density$y, type = 'l', lty = 1, col = '#d62728', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4.5),ylim = c(0, 5000),  main = "S Subtelo | SIT")

plot(ST.shared_ST_density$x, length(regions[regions$Category != 'IT' & grepl('S', regions$Chr), 'ST.shared']) * ST.shared_ST_density$y, type = 'l', lty = 1, col = 'green', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4.5),ylim = c(0, 5000),  main = "S Subtelo | ST.shared")

plot(IT.shared_ST_density$x, length(regions[regions$Category != 'IT' & grepl('S', regions$Chr), 'IT.shared']) * IT.shared_ST_density$y, type = 'l', lty = 1, col = 'black', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4.5),ylim = c(0, 5000),  main = "S Subtelo | IT.shared")


#all L inner


par(mfrow = c(1,1))

lmer_IT_density <- density(log10(regions[regions$Category == 'IT' & grepl('L', regions$Chr), 'Lmer']+1))
smer_IT_density <- density(log10(regions[regions$Category == 'IT' & grepl('L', regions$Chr), 'Smer']+1))
LST_IT_density <- density(log10(regions[regions$Category == 'IT' & grepl('L', regions$Chr), 'L.ST']+1))
SST_IT_density <- density(log10(regions[regions$Category == 'IT' & grepl('L', regions$Chr), 'S.ST']+1))
LIT_IT_density <- density(log10(regions[regions$Category == 'IT' & grepl('L', regions$Chr), 'L.IT']+1))
SIT_IT_density <- density(log10(regions[regions$Category == 'IT' & grepl('L', regions$Chr), 'S.IT']+1))
ST.shared_IT_density <- density(log10(regions[regions$Category == 'IT' & grepl('L', regions$Chr), 'ST.shared']+1))
IT.shared_IT_density <- density(log10(regions[regions$Category == 'IT' & grepl('L', regions$Chr), 'IT.shared']+1))

plot(smer_IT_density$x, length(regions[regions$Category == 'IT' & grepl('L', regions$Chr), 'Smer']) * smer_IT_density$y, type = 'l', lty = 1, col = '#d62728', ylab = 'Number of 100kb windows', xlab = 'log10(Number of bp overlap for each kmer category)', xlim = c(-0.5, 4.5),ylim = c(0, 30000),  main = "L Intertelomere") # RED
lines(lmer_IT_density$x, length(regions[regions$Category == 'IT' & grepl('L', regions$Chr), 'Lmer']) * lmer_IT_density$y, col = '#1f77b4', lty = 1) #blue
lines(LST_IT_density$x, length(regions[regions$Category == 'IT' & grepl('L', regions$Chr), 'L.ST']) * LST_IT_density$y, col = '#1f77b4', lty = 6) #BLUE
lines(SST_IT_density$x, length(regions[regions$Category == 'IT' & grepl('L', regions$Chr), 'S.ST']) * SST_IT_density$y, col = '#d62728', lty = 6) #RED
lines(LIT_IT_density$x, length(regions[regions$Category == 'IT' & grepl('L', regions$Chr), 'L.IT']) * LIT_IT_density$y, col = '#1f77b4', lty = 3) #teal '#17becf'
lines(SIT_IT_density$x, length(regions[regions$Category == 'IT' & grepl('L', regions$Chr), 'S.IT']) * SIT_IT_density$y, col = '#d62728', lty = 3) #orange
lines(ST.shared_IT_density$x, length(regions[regions$Category == 'IT' & grepl('L', regions$Chr), 'ST.shared']) * ST.shared_IT_density$y, col = 'green', lty = 1) #green
lines(IT.shared_IT_density$x, length(regions[regions$Category == 'IT' & grepl('L', regions$Chr), 'IT.shared']) * IT.shared_IT_density$y, col = 'black', lty = 1) #black

legend("topright", legend = c("lmer", "smer", "LST", "SST", "LIT", "SIT", "ST.shared", "IT.shared"), col = c('#1f77b4', '#d62728', '#1f77b4', '#d62728', '#1f77b4', "#d62728", 'green', 'black'), lty = c(1,1,6,6,3,3,1,1), cex = 0.8)



par(mfrow = c(4,2))

plot(lmer_IT_density$x, length(regions[regions$Category == 'IT' & grepl('L', regions$Chr), 'Lmer']) * lmer_IT_density$y, type = 'l', lty = 1, col = '#1f77b4', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4.5),ylim = c(0, 30000),  main = "L Intertelo | Lmer")

plot(smer_IT_density$x, length(regions[regions$Category == 'IT' & grepl('L', regions$Chr), 'Smer']) * smer_IT_density$y, type = 'l', lty = 1, col = '#d62728', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4.5),ylim = c(0, 30000),  main = "L Intertelo | Smer")

plot(LST_IT_density$x, length(regions[regions$Category == 'IT' & grepl('L', regions$Chr), 'L.ST']) * LST_IT_density$y, type = 'l', lty = 1, col = '#1f77b4', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4.5),ylim = c(0, 30000),  main = "L Intertelo | LST")

plot(SST_IT_density$x, length(regions[regions$Category == 'IT' & grepl('L', regions$Chr), 'S.ST']) * SST_IT_density$y, type = 'l', lty = 1, col = '#d62728', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4.5),ylim = c(0, 30000),  main = "L Intertelo | SST")

plot(LIT_IT_density$x, length(regions[regions$Category == 'IT' & grepl('L', regions$Chr), 'L.ST']) * LIT_IT_density$y, type = 'l', lty = 1, col = '#1f77b4', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4.5),ylim = c(0, 30000),  main = "L Intertelo | LIT")

plot(SIT_IT_density$x, length(regions[regions$Category == 'IT' & grepl('L', regions$Chr), 'S.ST']) * SIT_IT_density$y, type = 'l', lty = 1, col = '#d62728', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4.5),ylim = c(0, 30000),  main = "L Intertelo | SIT")

plot(ST.shared_IT_density$x, length(regions[regions$Category == 'IT' & grepl('L', regions$Chr), 'ST.shared']) * ST.shared_IT_density$y, type = 'l', lty = 1, col = 'green', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4.5),ylim = c(0, 30000),  main = "L Intertelo | ST.shared")

plot(IT.shared_IT_density$x, length(regions[regions$Category == 'IT' & grepl('L', regions$Chr), 'IT.shared']) * IT.shared_IT_density$y, type = 'l', lty = 1, col = 'black', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4.5),ylim = c(0, 30000),  main = "L Intertelo | IT.shared")

#all S inner

par(mfrow = c(1,1))

lmer_IT_density <- density(log10(regions[regions$Category == 'IT' & grepl('S', regions$Chr), 'Lmer']+1))
smer_IT_density <- density(log10(regions[regions$Category == 'IT' & grepl('S', regions$Chr), 'Smer']+1))
LST_IT_density <- density(log10(regions[regions$Category == 'IT' & grepl('S', regions$Chr), 'L.ST']+1))
SST_IT_density <- density(log10(regions[regions$Category == 'IT' & grepl('S', regions$Chr), 'S.ST']+1))
LIT_IT_density <- density(log10(regions[regions$Category == 'IT' & grepl('S', regions$Chr), 'L.IT']+1))
SIT_IT_density <- density(log10(regions[regions$Category == 'IT' & grepl('S', regions$Chr), 'S.IT']+1))
ST.shared_IT_density <- density(log10(regions[regions$Category == 'IT' & grepl('S', regions$Chr), 'ST.shared']+1))
IT.shared_IT_density <- density(log10(regions[regions$Category == 'IT' & grepl('S', regions$Chr), 'IT.shared']+1))

plot(smer_IT_density$x, length(regions[regions$Category == 'IT' & grepl('S', regions$Chr), 'Smer']) * smer_IT_density$y, type = 'l', lty = 1, col = '#d62728', ylab = 'Number of 100kb windows', xlab = 'log10(Number of bp overlap for each kmer category)', xlim = c(-0.5, 4.5),ylim = c(0, 25000),  main = "S Intertelomere") # RED
lines(lmer_IT_density$x, length(regions[regions$Category == 'IT' & grepl('S', regions$Chr), 'Lmer']) * lmer_IT_density$y, col = '#1f77b4', lty = 1) #blue
lines(LST_IT_density$x, length(regions[regions$Category == 'IT' & grepl('S', regions$Chr), 'L.ST']) * LST_IT_density$y, col = '#1f77b4', lty = 6) #BLUE
lines(SST_IT_density$x, length(regions[regions$Category == 'IT' & grepl('S', regions$Chr), 'S.ST']) * SST_IT_density$y, col = '#d62728', lty = 6) #RED
lines(LIT_IT_density$x, length(regions[regions$Category == 'IT' & grepl('S', regions$Chr), 'L.IT']) * LIT_IT_density$y, col = '#1f77b4', lty = 3) #teal '#17becf'
lines(SIT_IT_density$x, length(regions[regions$Category == 'IT' & grepl('S', regions$Chr), 'S.IT']) * SIT_IT_density$y, col = '#d62728', lty = 3) #orange
lines(ST.shared_IT_density$x, length(regions[regions$Category == 'IT' & grepl('S', regions$Chr), 'ST.shared']) * ST.shared_IT_density$y, col = 'green', lty = 1) #green
lines(IT.shared_IT_density$x, length(regions[regions$Category == 'IT' & grepl('S', regions$Chr), 'IT.shared']) * IT.shared_IT_density$y, col = 'black', lty = 1) #black

legend("topright", legend = c("lmer", "smer", "LST", "SST", "LIT", "SIT", "ST.shared", "IT.shared"), col = c('#1f77b4', '#d62728', '#1f77b4', '#d62728', '#1f77b4', "#d62728", 'green', 'black'), lty = c(1,1,6,6,3,3,1,1), cex = 0.8)



par(mfrow = c(4,2))

plot(lmer_IT_density$x, length(regions[regions$Category == 'IT' & grepl('S', regions$Chr), 'Lmer']) * lmer_IT_density$y, type = 'l', lty = 1, col = '#1f77b4', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4),ylim = c(0, 25000),  main = "S Intertelo | Lmer")

plot(smer_IT_density$x, length(regions[regions$Category == 'IT' & grepl('S', regions$Chr), 'Smer']) * smer_IT_density$y, type = 'l', lty = 1, col = '#d62728', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4),ylim = c(0, 25000),  main = "S Intertelo | Smer")

plot(LST_IT_density$x, length(regions[regions$Category == 'IT' & grepl('S', regions$Chr), 'L.ST']) * LST_IT_density$y, type = 'l', lty = 1, col = '#1f77b4', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4),ylim = c(0, 25000),  main = "S Intertelo | LST")

plot(SST_IT_density$x, length(regions[regions$Category == 'IT' & grepl('S', regions$Chr), 'S.ST']) * SST_IT_density$y, type = 'l', lty = 1, col = '#d62728', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4),ylim = c(0, 25000),  main = "S Intertelo | SST")

plot(LIT_IT_density$x, length(regions[regions$Category == 'IT' & grepl('S', regions$Chr), 'L.ST']) * LIT_IT_density$y, type = 'l', lty = 1, col = '#1f77b4', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4),ylim = c(0, 25000),  main = "S Intertelo | LIT")

plot(SIT_IT_density$x, length(regions[regions$Category == 'IT' & grepl('S', regions$Chr), 'S.ST']) * SIT_IT_density$y, type = 'l', lty = 1, col = '#d62728', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4),ylim = c(0, 25000),  main = "S Intertelo | SIT")

plot(ST.shared_IT_density$x, length(regions[regions$Category == 'IT' & grepl('S', regions$Chr), 'ST.shared']) * ST.shared_IT_density$y, type = 'l', lty = 1, col = 'green', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4),ylim = c(0, 25000),  main = "S Intertelo | ST.shared")

plot(IT.shared_IT_density$x, length(regions[regions$Category == 'IT' & grepl('S', regions$Chr), 'IT.shared']) * IT.shared_IT_density$y, type = 'l', lty = 1, col = 'black', ylab = 'Frequency', xlab = 'log10(bp overlap)', xlim = c(-0.5, 4),ylim = c(0, 25000),  main = "S Intertelo | IT.shared")

dev.off()

# pdf("kernal.density.by.chr.plots.pdf", useDingbats = F)
# for (Chr in unique(regions$Chr)){
#     lmer_ST_density <- density(log10(regions[regions$Category != 'IT' & regions$Chr == Chr, 'Lmer']+1))
#     smer_ST_density <- density(log10(regions[regions$Category != 'IT' & regions$Chr == Chr, 'Smer']+1))
#     LST_ST_density <- density(log10(regions[regions$Category != 'IT' & regions$Chr == Chr, 'L.ST']+1))
#     SST_ST_density <- density(log10(regions[regions$Category != 'IT' & regions$Chr == Chr, 'S.ST']+1))
#     LIT_ST_density <- density(log10(regions[regions$Category != 'IT' & regions$Chr == Chr, 'L.IT']+1))
#     SIT_ST_density <- density(log10(regions[regions$Category != 'IT' & regions$Chr == Chr, 'S.IT']+1))
#     ST.shared_ST_density <- density(log10(regions[regions$Category != 'IT' & regions$Chr == Chr, 'ST.shared']+1))
#     IT.shared_ST_density <- density(log10(regions[regions$Category != 'IT' & regions$Chr == Chr, 'IT.shared']+1))
    
#     lmer_IT_density <- density(log10(regions[regions$Category == 'IT' & regions$Chr == Chr, 'Lmer']+1))
#     smer_IT_density <- density(log10(regions[regions$Category == 'IT' & regions$Chr == Chr, 'Smer']+1))
#     LST_IT_density <- density(log10(regions[regions$Category == 'IT' & regions$Chr == Chr, 'L.ST']+1))
#     SST_IT_density <- density(log10(regions[regions$Category == 'IT' & regions$Chr == Chr, 'S.ST']+1))
#     LIT_IT_density <- density(log10(regions[regions$Category == 'IT' & regions$Chr == Chr, 'L.IT']+1))
#     SIT_IT_density <- density(log10(regions[regions$Category == 'IT' & regions$Chr == Chr, 'S.IT']+1))
#     ST.shared_IT_density <- density(log10(regions[regions$Category == 'IT' & regions$Chr == Chr, 'ST.shared']+1))
#     IT.shared_IT_density <- density(log10(regions[regions$Category == 'IT' & regions$Chr == Chr, 'IT.shared']+1))

#     #Subtelo

#     plot(lmer_ST_density$x, length(regions[regions$Category != 'IT' & regions$Chr == Chr, 'Lmer']) * lmer_ST_density$y, col = '#1f77b4', type = 'l', lty = 1 , ylim = c(0,700),xlim = c(-1,4.5), main = paste(Chr, "Subtelo", sep = ' '), ylab = 'Number of 100kb windows', xlab = 'log10(bp overlap per 100kb window)') #BLUE
#     lines(smer_ST_density$x, length(regions[regions$Category != 'IT' & regions$Chr == Chr, 'Smer']) * smer_ST_density$y, col = '#d62728', lty = 1) #RED
#     lines(LST_ST_density$x, length(regions[regions$Category != 'IT' & regions$Chr == Chr, 'L.ST']) * LST_ST_density$y, col = '#1f77b4', lty = 6) #BLUE
#     lines(SST_ST_density$x, length(regions[regions$Category != 'IT' & regions$Chr == Chr, 'S.ST']) * SST_ST_density$y, col = '#d62728', lty = 6) #RED
#     lines(LIT_ST_density$x, length(regions[regions$Category != 'IT' & regions$Chr == Chr, 'L.IT']) * LIT_ST_density$y, col = '#1f77b4', lty = 3) #BLUE
#     lines(SIT_ST_density$x, length(regions[regions$Category != 'IT' & regions$Chr == Chr, 'S.IT']) * SIT_ST_density$y, col = "#d62728", lty = 3) #RED
#     lines(ST.shared_ST_density$x, length(regions[regions$Category != 'IT' & regions$Chr == Chr, 'ST.shared']) * ST.shared_ST_density$y, col = 'green', lty = 1) #green
#     lines(IT.shared_ST_density$x, length(regions[regions$Category != 'IT' & regions$Chr == Chr, 'IT.shared']) * IT.shared_ST_density$y, col = 'black', lty = 1) #black

#     legend("topright", legend = c("lmer", "smer", "LST", "SST", "LIT", "SIT", "ST.shared", "IT.shared"), col = c('#1f77b4', '#d62728', '#1f77b4', '#d62728', '#17becf', "orange", 'green', 'black'), lty = c(1,1,6,6,3,3,1,1), cex = 0.8)

#     #Intertelo

#     plot(lmer_IT_density$x, length(regions[regions$Category == 'IT' & regions$Chr == Chr, 'Lmer']) * lmer_IT_density$y, col = '#1f77b4', type = 'l', lty = 1 , ylim = c(0,4800), xlim = c(-1,4.5), main = paste(Chr, "Intertelo", sep = ' '), ylab = 'Number of 100kb windows', xlab = 'log10(bp overlap per 100kb window)') #BLUE
#     lines(smer_IT_density$x, length(regions[regions$Category == 'IT' & regions$Chr == Chr, 'Smer']) * smer_IT_density$y, col = '#d62728', lty = 1) #RED
#     lines(LST_IT_density$x, length(regions[regions$Category == 'IT' & regions$Chr == Chr, 'L.ST']) * LST_IT_density$y, col = '#1f77b4', lty = 6) #BLUE
#     lines(SST_IT_density$x, length(regions[regions$Category == 'IT' & regions$Chr == Chr, 'S.ST']) * SST_IT_density$y, col = '#d62728', lty = 6) #RED
#     lines(LIT_IT_density$x, length(regions[regions$Category == 'IT' & regions$Chr == Chr, 'L.IT']) * LIT_IT_density$y, col = '#1f77b4', lty = 3) #teal
#     lines(SIT_IT_density$x, length(regions[regions$Category == 'IT' & regions$Chr == Chr, 'S.IT']) * SIT_IT_density$y, col = "#d62728", lty = 3) #orange
#     lines(ST.shared_IT_density$x, length(regions[regions$Category == 'IT' & regions$Chr == Chr, 'ST.shared']) * ST.shared_IT_density$y, col = 'green', lty = 1) #green
#     lines(IT.shared_IT_density$x, length(regions[regions$Category == 'IT' & regions$Chr == Chr, 'IT.shared']) * IT.shared_IT_density$y, col = 'black', lty = 1) #black

#     legend("topright", legend = c("lmer", "smer", "LST", "SST", "LIT", "SIT", "ST.shared", "IT.shared"), col = c('#1f77b4', '#d62728', '#1f77b4', '#d62728', '#1f77b4', "#d62728", 'green', 'black'), lty = c(1,1,6,6,3,3,1,1), cex = 0.8)

# }
# dev.off()
