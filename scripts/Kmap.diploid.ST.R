
args<- commandArgs(trailingOnly = TRUE);
if (length(args)!=6) {
  stop("Correct usage is Rscript Kmap.diploid.ST.R ChrSizeFile STCountTable STSizeNumber KmerCountNumber TukeyOutFile STNormCountOutFile.n", call.=FALSE)
} 

chr.size<-as.matrix(read.table(args[1],row=1,header=F));
ST.count<-as.matrix(read.table(args[2],row=1,header=T));
ST.size<-as.numeric(args[3]);
chr.st.size<-chr.size;
chr.st.size<-chr.st.size-(2*as.numeric(ST.size));
my.IT.list<-c();
my.ST.list<-c();
ST.norm<-ST.count;

for(i in colnames(ST.norm)){if(endsWith(i,'IT')=="FALSE"){ST.norm[,i]<-ST.norm[,i]/ST.size;my.ST.list<-c(i,my.ST.list);}else{my.chr<-as.character(strsplit(i,split="_")[[1]][1]);ST.norm[,i]<-ST.norm[,i]/chr.st.size[my.chr,1];my.IT.list<-c(i,my.IT.list);}}

ST.tukey<-matrix(nrow=length(rownames(ST.norm)),ncol=13);
colnames(ST.tukey)<-c("Group.df","Group.Sumsq","Group.MeanSq","Residual.df","Residual.Sumsq","Residual.MeanSq","Fvalue","Pr","ST-IT.diff","ST-IT.p","ST-IT.95l","ST-IT.95u","ST-IT.bonf");
rownames(ST.tukey)<-rownames(ST.norm);



for(i in rownames(ST.tukey)){

my.iter1.dat<-data.frame(y=character(),Group=character());my.iter1.dat<-rbind(my.iter1.dat,data.frame(y=ST.norm[i,my.IT.list],Group="IT"));my.iter1.dat<-rbind(my.iter1.dat,data.frame(y=ST.norm[i,my.ST.list],Group="ST"));

fit1<-lm(y ~ Group, my.iter1.dat);

my.anova<-anova(fit1)

ST.tukey[i,1]<-my.anova["Group","Df"];ST.tukey[i,2]<-my.anova["Group","Sum Sq"];ST.tukey[i,3]<-my.anova["Group","Mean Sq"];ST.tukey[i,4]<-my.anova["Residuals","Df"];ST.tukey[i,5]<-my.anova["Residuals","Sum Sq"];ST.tukey[i,6]<-my.anova["Residuals","Mean Sq"];ST.tukey[i,7]<-my.anova["Group","F value"];ST.tukey[i,8]<-my.anova["Group","Pr(>F)"];

my.aov<-aov(fit1);

my.Tukey<-TukeyHSD(my.aov);

ST.tukey[i,9]<-my.Tukey$Group["ST-IT","diff"];ST.tukey[i,10]<-my.Tukey$Group["ST-IT","p adj"];ST.tukey[i,11]<-my.Tukey$Group["ST-IT","lwr"];ST.tukey[i,12]<-my.Tukey$Group["ST-IT","upr"];


}
ST.tukey[,13]<-p.adjust(p=as.numeric(ST.tukey[,10]),method="bonferroni",n=as.numeric(args[4]));
write.table(ST.tukey,file=args[5],row=T,col=T,quote=F,sep="\t");
write.table(ST.norm,file=args[6],row=T,col=T,quote=F,sep="\t");
#ST.sig<-subset(ST.tukey,ST.tukey[,13]<=0.05);

