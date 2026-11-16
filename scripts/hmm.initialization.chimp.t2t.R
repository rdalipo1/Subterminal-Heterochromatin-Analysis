
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

emitprob <- matrix(nrow =2 , ncol = 3)
rownames(emitprob) <- states
colnames(emitprob) <- symbols

### Choosing training chromosome and setting up emmission table

chr7.st <- subset(all_emits, all_emits$Chr == "Chr7" & all_emits$Category != "IT")
chr7.it <- subset(all_emits, all_emits$Chr == "Chr7" & all_emits$Category == "IT")

table(chr7.it$emit)
table(chr7.st$emit)

emitprob[1,]<- table(chr7.it$emit)/sum(table(chr7.it$emit))
emitprob[2,]<- table(chr7.st$emit)/sum(table(chr7.st$emit)

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
                                       
head(all_emits)

library(karyoploteR)
library(dplyr)

genome <- read.table(args[3])
genome_size <- read.table(args[3])
genome$Start <- 0
genome <- genome[c('V1', 'Start', 'V2')]
genome <- toGRanges(genome)


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

chr_list <- c(
  "Chr1", "Chr2", "Chr3", "Chr4", "Chr5",
  "Chr6", "Chr7", "Chr8", "Chr9", "Chr10", "Chr11",
  "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17",
  "Chr18", "Chr19", "Chr20", "Chr21", "Chr22","Chr23", "ChrX"
)

# If only showing the hmm results:
                                       
# pdf(paste0(args[2],".hmm.path.karyoplots.pdf"), useDingbats = F)
# kp <- plotKaryotype(genome = genome, chromosomes = 'all', cex = 0.8, plot.type = 1)

# for(Chr in chr_list){
#     kpRect(kp, data = it, col = "blue", y0 = 0, y1 = 1, border = NA)
#     kpRect(kp, data = st, col = "cyan", y0 = 0, y1 = 1, border = NA)
#     kpAbline(kp, v = 5000000, col = 'red')
#     kpAbline(kp, v = genome_size[,'V2'] - 5000000,col = 'red')
# }
# dev.off()

# Hard coded comparison values from Yoo et al. 2025
                                       
all_emits <- all_emits %>%
  mutate(
    yoo2025 = case_when(
    Chr == "Chr2" & (Start >= 9260 & End <= 3701568 ) ~ "ST",
    Chr == "Chr5" & (Start >= 5259 & End <= 5851120 ) ~ "ST",
    Chr == "Chr6" & (Start >= 112185101 & End <= 118079859) ~ "ST",
    Chr == "Chr6" & (Start >= 7115 & End <= 5762040) ~ "ST",
    Chr == "Chr7" & (Start >= 5226 & End <= 7031729) ~ "ST",
    Chr == "Chr8" & ( Start >= 3508 & End <= 12635833) ~ "ST",
    Chr == "Chr9" & (Start >= 12881 & End <= 4261115) ~ "ST",
    Chr == "Chr10" & (Start >= 9395 & End <= 7093694) ~ "ST",
    Chr == "Chr11" & (Start >= 115835047 & End <= 118858361) ~ "ST",
    Chr == "Chr12" & (Start >= 4835 & End <= 9729033) ~ "ST",
    Chr == "Chr13" & (Start >= 8913 & End <= 14168301) ~ "ST",
    Chr == "Chr14" & (Start >= 45617284 & End <= 52003170) ~ "ST",
    Chr == "Chr14" & (Start >= 3652408 & End <= 5093942) ~ "ST",
    Chr == "Chr14" & (Start >= 9368928 & End <= 10490943) ~ "ST",
    Chr == "Chr15" & (Start >=  8098516 & End <= 9428111) ~ "ST",
    Chr == "Chr15" & (Start >=  12748574 & End <= 14015352) ~ "ST",
    Chr == "Chr16" & (Start >= 4985 & End <= 5082537) ~ "ST",
    Chr == "Chr17" & (Start >= 8036009 & End <= 9192116) ~ "ST",
    Chr == "Chr18" & (Start >= 7369 & End <= 2318617) ~ "ST",
    Chr == "Chr18" & (Start >= 88416064 & End <= 94519508) ~ "ST",
    Chr == "Chr19" & (Start >= 9278 & End <= 10557213) ~ "ST",
    Chr == "Chr19" & (Start >= 96653435 & End <= 99588077) ~ "ST",
    Chr == "Chr20" & (Start >= 4923 & End <= 2444115) ~ "ST",
    Chr == "Chr20" & (Start >= 61577641 & End <= 68073140) ~ "ST",
    Chr == "Chr21" & (Start >= 17389 & End <= 6785509) ~ "ST",
    Chr == "Chr21" & (Start >= 70766427 & End <= 78110418) ~ "ST",
    Chr == "Chr22" & (Start >= 3338094 & End <= 4527108) ~ "ST",
    Chr == "Chr22" & (Start >= 6339443 & End <= 9810504) ~ "ST",
    Chr == "Chr22" & (Start >= 47454653 & End <= 51345201) ~ "ST",
    Chr == "Chr23" & (Start >= 8711135 & End <= 10707668) ~ "ST",
    Chr == "Chr23" & (Start >= 53045785 & End <= 59516181) ~ "ST",
    Chr == "ChrX" & (Start >= 27462 & End <= 685536) ~ "ST",
    TRUE ~ "IT"
  ))


sty <- subset(all_emits[,c(1,2,3)], all_emits[,8] == "ST")
sty <- toGRanges(sty, genome = genome)

ity <- subset(all_emits[,c(1,2,3)], all_emits[,8] == "IT")
ity <- toGRanges(ity, genome = genome)




pdf(paste0(args[2],".hmm.path.karyoplots.withyoo2025.pdf"), useDingbats = F)
kp <- plotKaryotype(genome = genome, chromosomes = 'all', cex = 0.8, plot.type = 2)

for(Chr in chr_list){
    kpRect(kp, data = it, col = "blue", y0 = 0, y1 = 1, border = NA, data.panel = 1)
    kpRect(kp, data = st, col = "cyan", y0 = 0, y1 = 1, border = NA, data.panel = 1)

    kpRect(kp, data = ity, col = "red", y0 = 0, y1 = 1, border = NA, data.panel = 2)
    kpRect(kp, data = sty, col = "yellow", y0 = 0, y1 = 1, border = NA, data.panel = 2)
}
# kpAddLabels(kp, labels="HMM chimp.t2t", data.panel = 1, cex = 0.6, col = 'blue', side = 'right')
# kpAddLabels(kp, labels="Yoo (2025)", data.panel = 2, cex = 0.6, col = 'blue', side = 'right')
dev.off()


