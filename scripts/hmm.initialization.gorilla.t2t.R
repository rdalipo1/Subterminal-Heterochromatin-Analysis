
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

### Choosing training chromosome and setting up emmission table 

emitprob <- matrix(nrow =2 , ncol = 3)
rownames(emitprob) <- states
colnames(emitprob) <- symbols

chr7.st <- subset(all_emits, all_emits$Chr == "Chr7" & all_emits$Category != "IT")
chr7.it <- subset(all_emits, all_emits$Chr == "Chr7" & all_emits$Category == "IT")

table(chr7.it$emit)
table(chr7.st$emit)

emitprob[1,]<- table(chr7.it$emit)/sum(table(chr7.it$emit))
emitprob[2,]<- table(chr7.st$emit)/sum(table(chr7.st$emit))


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
    Chr == "Chr1" & (Start >= 13397  & End <= 5535160) ~ "ST",
    Chr == "Chr1" & (Start >= 241579359 & End <= 243835808 ) ~ "ST",

    Chr == "Chr2" & (Start >= 9938  & End <=9084796) ~ "ST",
    Chr == "Chr2" & (Start >= 211876084 & End <= 215245360) ~ "ST",

    Chr == "Chr3" & (Start >= 6637 & End <=7644626) ~ "ST",
    Chr == "Chr3" & (Start >= 212524079  & End <=215364096 ) ~ "ST",

    Chr == "Chr4" & (Start >= 11322 & End <=8694292 ) ~ "ST",
    Chr == "Chr4" & (Start >= 185478912 & End <=195324532 ) ~ "ST",

    Chr == "Chr5" & (Start >= 10798   & End <= 17587290   ) ~ "ST",
    Chr == "Chr5" & (Start >= 195574868  & End <= 199439332  ) ~ "ST",

    Chr == "Chr6" & (Start >= 10267  & End <= 6601210 ) ~ "ST",
    Chr == "Chr6" & (Start >= 169065139  & End <= 170793360 ) ~ "ST",

    Chr == "Chr7" & (Start >= 7487  & End <= 4470846 ) ~ "ST",
    Chr == "Chr7" & (Start >= 151474829 & End <= 158527604  ) ~ "ST",

    Chr == "Chr8" & ( Start >= 13439 & End <= 8319882 ) ~ "ST",
    Chr == "Chr8" & ( Start >= 146489047  & End <= 150795270  ) ~ "ST",

    Chr == "Chr9" & (Start >= 10736  & End <= 5373931 ) ~ "ST",
    Chr == "Chr9" & (Start >= 144214193 & End <= 147863740 ) ~ "ST",

    Chr == "Chr10" & (Start >= 9359  & End <= 15507589 ) ~ "ST",
    Chr == "Chr10" & (Start >= 149520535   & End <= 154217949  ) ~ "ST",

    Chr == "Chr11" & (Start >= 11076  & End <= 15059675 ) ~ "ST",
    Chr == "Chr11" & (Start >= 145498381  & End <= 150235715 ) ~ "ST",

    Chr == "Chr12" & (Start >= 140882027  & End <= 145891758 ) ~ "ST",

    Chr == "Chr14" & (Start >= 129303957  & End <= 132005803 ) ~ "ST",

    Chr == "Chr15" & (Start >=  7764  & End <= 11481686 ) ~ "ST",
    Chr == "Chr15" & (Start >=  124919076  & End <= 131679779 ) ~ "ST",

    Chr == "Chr16" & (Start >= 104033428  & End <= 107906169 ) ~ "ST",

    Chr == "Chr17" & (Start >= 110461100  & End <= 116047618 ) ~ "ST",

    Chr == "Chr18" & (Start >= 9650  & End <= 2841911 ) ~ "ST",
    Chr == "Chr18" & (Start >= 121510407  & End <= 131342505) ~ "ST",

    Chr == "Chr19" & (Start >= 10532  & End <= 15249603 ) ~ "ST",
    Chr == "Chr19" & (Start >= 112924778  & End <= 120447406 ) ~ "ST",

    Chr == "Chr20" & (Start >= 6626  & End <= 5402075 ) ~ "ST",
    Chr == "Chr20" & (Start >= 71036979  & End <= 80303887 ) ~ "ST",

    Chr == "Chr21" & (Start >= 8424  & End <= 9952369 ) ~ "ST",
    Chr == "Chr21" & (Start >= 76662328  & End <= 84434812 ) ~ "ST",

    Chr == "Chr22" & (Start >= 48678645  & End <= 59498265 ) ~ "ST",

    Chr == "Chr23" & (Start >= 48894677  & End <= 54762103 ) ~ "ST",
    
    Chr == "ChrX" & (Start >= 8875  & End <= 9544179 ) ~ "ST",
    Chr == "ChrX" & (Start >= 168016846  & End <= 177550564 ) ~ "ST",
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
# kpAddLabels(kp, labels="HMM gorilla.t2t", data.panel = 1, cex = 0.6, col = 'blue', side = 'right')
# kpAddLabels(kp, labels="Yoo (2025)", data.panel = 2, cex = 0.6, col = 'blue', side = 'right')


dev.off()
