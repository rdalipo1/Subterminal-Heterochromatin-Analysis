
args <- commandArgs(trailingOnly = TRUE)
# 1 - bpoverlap.by.chrwindow.with.regions.emit.tsv
# 2 - chr_lengths file
# 3 - pdf file name
# 4 - low cut off
# 5 - high cut off
# 6 - species
# 7 - terminal region size, e.g. 5

regions <- read.table(args[1], header = TRUE)


chr_lengths <- read.table(args[2], header = FALSE, col.names = c("Chr","Length"))

pdf(paste0(args[3],'.kernel.density.pdf'), useDingbats = F)
ST_density <- density(log10(regions[regions$Category != 'IT', 'ST.kmer'] + 1))
IT_density <- density(log10(regions[regions$Category == 'IT', 'ST.kmer'] + 1))

plot(ST_density$x, length(regions[regions$Category != 'IT','ST.kmer']) * ST_density$y, col = '#1f77b4', type = 'l', lty = 1 , main = 'Terminal 20MB', ylab = 'Number of 100kb windows', xlab = 'log10(bp overlap per 100kb window)')

plot(IT_density$x, length(regions[regions$Category == 'IT','ST.kmer']) * IT_density$y, col = '#1f77b4', type = 'l', lty = 1 , main = 'Interstitial', ylab = 'Number of 100kb windows', xlab = 'log10(bp overlap per 100kb window)')
dev.off()


pdf(paste0(args[3],".density.historgrams.pdf"), useDingbats = F)
par(mfrow = c(2,1))
hist(log10(regions[regions$Category == 'IT',"ST.kmer"]), breaks = 200, xlim=c(1,5), main = paste0(args[6]," Intersitial"), xlab = paste0("log10(bp overlap per 100kb window) | Green Lines: ", args[4], ", ", args[5]))
abline(v = log10(as.numeric(args[4])), col = 'green')
abline(v = log10(as.numeric(args[5])), col = 'green')
hist(log10(regions[regions$Category != 'IT',"ST.kmer"]), breaks = 200, xlim=c(1,5),main = paste0(args[6]," Terminal ", args[7],"MB"), xlab = paste0("log10(bp overlap per 100kb window) | Green Lines: ", args[4], ", ", args[5]))
abline(v = log10(as.numeric(args[4])), col = 'green')
abline(v = log10(as.numeric(args[5])), col = 'green')
dev.off()


