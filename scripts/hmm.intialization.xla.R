
require(HMM)

all_emits <- read.table("/data/home/rdalipo1/subtelo_proj/final.window.with.emit.tsv", header = T)


states <- c('lit','lst', 'sst','sit')
symbols <- c('LIT','LST', 'SIT','SST', 'Unknown')

start <- c(0.05,0.45,0.45, 0.05)

transprob <- matrix(nrow = 4, ncol = 4)
rownames(transprob) <- states
colnames(transprob) <- states


transprob[1,] <- c(.95,.04,.005,.005)
transprob[2,] <- c(.03,.95,.015,.005)
transprob[3,] <- c(.005,.015,.95,.03)
transprob[4,] <- c(.005,.005,.04,.95)


transprob

### Choosing training chromosome pairs and setting up emission probability table

emitprob <- matrix(nrow =4 , ncol = 5)
rownames(emitprob) <- states
colnames(emitprob) <- symbols


chr2l.st <- subset(all_emits, all_emits$Chr == "Chr2L" & all_emits$Category != "IT")
chr2s.st <- subset(all_emits, all_emits$Chr == "Chr2S" & all_emits$Category != "IT")
chr2l.it <- subset(all_emits, all_emits$Chr == "Chr2L" & all_emits$Category == "IT")
chr2s.it <- subset(all_emits, all_emits$Chr == "Chr2S" & all_emits$Category == "IT")


emitprob[1,]<- table(chr2l.it$emit)/sum(table(chr2l.it$emit))
emitprob[2,]<- table(chr2l.st$emit)/sum(table(chr2l.st$emit))
emitprob[3,]<- table(chr2s.st$emit)/sum(table(chr2s.st$emit))
emitprob[4,]<- table(chr2s.it$emit)/sum(table(chr2s.it$emit))

emitprob

chr2l <- subset(all_emits, all_emits$Chr == "Chr2L")
chr2s <- subset(all_emits, all_emits$Chr == "Chr2S")
# table(chr2l$emit)
# table(chr2s$emit)
# table(all_emits$emit)

inithmm <- initHMM(states,symbols,startProbs = start, transProbs = transprob, emissionProbs = emitprob)


all_emits$path <- NA


for(i in unique(all_emits$Chr)){
    current_chr <- subset(all_emits, all_emits$Chr == i)
    current_chr$path <- viterbi(inithmm, current_chr$emit)
    all_emits$path[rownames(all_emits) %in% rownames(current_chr)] <- current_chr$path
}

all_paths <- c("lit", "lst", "sst", "sit") 

write.table(all_emits, "final.window.with.emit.path.tsv", quote = F, col.names = T)


### Visualizing hmm assignments 

library(karyoploteR)


genome_size <- read.table("/data/home/rdalipo1/v10.XlaXtr.homology/Xla/Xenla10.chr.size.txt")
genome_size <- genome_size[1:18, ]
rownames(genome_size) <- genome_size$V1
genome_xla <- read.table("/data/home/rdalipo1/v10.XlaXtr.homology/Xla/Xenla10.chr.size.txt")
genome_xla$Start <- 1
genome_xla <- genome_xla[c('V1', 'Start', 'V2')]
genome_xla <- genome_xla[grepl("Chr", genome_xla$V1),]
genome_xla <- genome_xla[c(1,3,5,7,9,11,13,15,17,2,4,6,8,10,12,14,16,18),]
genome_xla <- toGRanges(genome_xla)


all_emits$val <- NA
all_emits[,21] <- 1
lst <- subset(all_emits[,c(1,2,3)], all_emits[,20] == "lst")
lst <- toGRanges(lst, genome = genome_xla)

sst <- subset(all_emits[,c(1,2,3)], all_emits[,20] == "sst")
sst <- toGRanges(sst, genome = genome_xla)

lit <- subset(all_emits[,c(1,2,3)], all_emits[,20] == "lit")
lit <- toGRanges(lit, genome = genome_xla)

sit <- subset(all_emits[,c(1,2,3)], all_emits[,20] == "sit")
sit <- toGRanges(sit, genome = genome_xla)


chr_list <- c('Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7','Chr8','Chr9_10')


pdf("hmm.path.karyoplots.xla.pdf", useDingbats = F)
kp <- plotKaryotype(genome = genome_xla, chromosomes = 'all', cex = 0.8, plot.type = 1)

for(Chr in chr_list){
    L <- paste(Chr,"L", sep = "")
    S <- paste(Chr,"S", sep = "")

    kpRect(kp, data = lst, col = "cyan", y0 = 0, y1 = 1, border = NA)
    kpRect(kp, data = sst, col = "yellow", y0 = 0, y1 = 1, border = NA)
    kpRect(kp, data = lit, col = "blue", y0 = 0, y1 = 1, border = NA)
    kpRect(kp, data = sit, col = "red", y0 = 0, y1 = 1, border = NA)
}
dev.off()


