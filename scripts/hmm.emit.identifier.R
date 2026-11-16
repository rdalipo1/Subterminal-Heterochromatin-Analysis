library(dplyr)


args <- commandArgs(trailingOnly = TRUE)
# 1 - bpoverlap.by.chrwindow.with.regions.tsv
# 2 - lower cut off
# 3 - higher cut off
# 4 - name 

data <- read.table(args[1], header = T)

last_rows <- tapply(1:nrow(data), data$Chr, max)
data <- data[-last_rows,]

for( i in 1:nrow(data)) {
if(data$ST.kmer[i] > as.numeric(args[2])) {		
	data$emit[i] <- "ST"
	}
if(data$ST.kmer[i] < as.numeric(args[3])){		
	data$emit[i] <- "IT"
}
if(data$ST.kmer[i] > as.numeric(args[2]) & data$ST.kmer[i] < as.numeric(args[3])){
	data$emit[i] <- "Unknown"}
}

table(data$emit)

write.table(data, paste0(args[4],".bpoverlap.by.chrwindow.with.regions.emit.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

