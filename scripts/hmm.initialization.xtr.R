
require(HMM)
args <- commandArgs(trailingOnly = TRUE)
# 1 - bpoverlap.by.chrwindow.with.regions.emit.tsv
# 2 - file identifier (e.x: chimp)
# 3 - chr.size file (2 columns)

all_emits <- read.table(args[1], header = T)



states <- c('it', 'st')
symbols <- c('IT', 'ST', 'Unknown')

start <- c(.05,.95)

transprob <- matrix(nrow = 2, ncol = 2)
rownames(transprob) <- states
colnames(transprob) <- states

transprob[1,] <- c(.99,.01)
transprob[2,] <- c(.05,.95)

transprob

### Choosing training chromosome and setting up emission table

emitprob <- matrix(nrow =2 , ncol = 3)
rownames(emitprob) <- states
colnames(emitprob) <- symbols


chr6.st <- subset(all_emits, all_emits$Chr == "Chr6" & all_emits$Category != "IT")
chr6.it <- subset(all_emits, all_emits$Chr == "Chr6" & all_emits$Category == "IT")



emitprob[1,]<- table(chr6.it$emit)/sum(table(chr6.it$emit))
emitprob[2,]<- table(chr6.st$emit)/sum(table(chr6.st$emit))

emitprob


inithmm <- initHMM(states,symbols,startProbs = start, transProbs = transprob, emissionProbs = emitprob)


all_emits$path <- NA


for(i in unique(all_emits$Chr)){
    current_chr <- subset(all_emits, all_emits$Chr == i)
    current_chr$path <- viterbi(inithmm, current_chr$emit)
    all_emits$path[rownames(all_emits) %in% rownames(current_chr)] <- current_chr$path
}

write.table(all_emits, paste0(args[2],'.bpoverlap.by.chrwindow.with.regions.emit.path.tsv'), quote = F, col.names = T, row.names = F)

### Visualizing hmm assignments using karyoploteR

library(karyoploteR)
library(dplyr)

genome_size <- read.table(args[3])
genome_size <- genome_size[1:10, ]
rownames(genome_size) <- genome_size$V1


### Hard coded Bredeson et al 2024 subterminal calls for comparison

all_emits <- all_emits %>%
  mutate(
    Bredeson = case_when(
    Chr == "Chr1" & (Start < 34000000 | End > (genome_size[1,2] - 28470000)) ~ "ST",
    Chr == "Chr2" & (Start < 27000000 | End > (genome_size[2,2] - 28030000)) ~ "ST",
    Chr == "Chr3" & (Start < 14000000 | End > (genome_size[3,2] - 21870000)) ~ "ST",
    Chr == "Chr4" & (Start < 21000000 | End > (genome_size[4,2] - 23960000)) ~ "ST",
    Chr == "Chr5" & (Start < 24000000 | End > (genome_size[5,2] - 30030000)) ~ "ST",
    Chr == "Chr6" & (Start < 26000000 | End > (genome_size[6,2] - 23480000)) ~ "ST",
    Chr == "Chr7" & (Start < 25000000 | End > (genome_size[7,2] - 21560000)) ~ "ST",
    Chr == "Chr8" & (Start < 17000000 | End > (genome_size[8,2] - 30240000)) ~ "ST",
    Chr == "Chr9" & (Start < 20000000 | End > (genome_size[9,2] - 21210000)) ~ "ST",
    Chr == "Chr10" & (Start < 13000000 | End > (genome_size[10,2] - 10430000)) ~ "ST",
    TRUE ~ "IT"
  ))


genome_xtr <- read.table(args[3])
genome_xtr$Start <- 0
genome_xtr <- genome_xtr[c('V1', 'Start', 'V2')]
genome_xtr <- genome_xtr[grepl("Chr", genome_xtr$V1),]
genome_xtr <- toGRanges(genome_xtr)


## Collapsing hmm calls so it is easier to load the image
collapsed_path <- all_emits %>%
  group_by(Chr) %>%
  mutate(grp = cumsum(path != lag(path, default = first(path)))) %>%
  group_by(Chr, grp, path) %>%
  summarise(Start = min(Start), End = max(End), .groups = "drop") %>%
  select(Chr, Start, End, path)

write.table(collapsed_path, paste0(args[2],'.hmm.path.bed'), quote = F, col.names = T, row.names = F)
st <- toGRanges(as.data.frame(subset(collapsed_path[,1:3], collapsed_path$path == "st")), genome = genome)
it <- toGRanges(as.data.frame(subset(collapsed_path[,1:3], collapsed_path$path == "it")), genome = genome)


all_emits$val <- NA
all_emits[,9] <- 1

stj <- subset(all_emits[,c(1,2,3)], all_emits[,8] == "ST")
stj <- toGRanges(stj, genome = genome_xtr)

itj <- subset(all_emits[,c(1,2,3)], all_emits[,8] == "IT")
itj <- toGRanges(itj, genome = genome_xtr)


hmm <- toGRanges(all_emits[,c(1,2,3,9)], genome = genome_xtr)

head(all_emits)

chr_list <- c('Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7','Chr8','Chr9','Chr10')


pdf(paste0(args[2],".hmm.path.karyoplots.withbredeson.pdf"), useDingbats = F)
kp <- plotKaryotype(genome = genome_xtr, chromosomes = 'all', cex = 0.8, plot.type = 2)

for(Chr in chr_list){
    kpRect(kp, data = it, col = "blue", y0 = 0, y1 = 1, border = NA, data.panel = 1)
    kpRect(kp, data = st, col = "cyan", y0 = 0, y1 = 1, border = NA, data.panel = 1)

    kpRect(kp, data = itj, col = "purple", y0 = 0, y1 = 1, border = NA, data.panel = 2)
    kpRect(kp, data = stj, col = "plum", y0 = 0, y1 = 1, border = NA, data.panel = 2)
}
# kpAddLabels(kp, labels="HMM", data.panel = 1, cex = 0.6, col = 'blue', side = 'right')
# kpAddLabels(kp, labels="Bredeson (2024)", data.panel = 2, cex = 0.6, col = 'blue', side = 'right')


dev.off()


