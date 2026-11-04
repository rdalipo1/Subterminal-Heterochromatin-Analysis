#!/data/home/rdalipo1/../asession/anaconda3/envs/newspiedie/bin/Rscript

args<- commandArgs(trailingOnly = TRUE);
if (length(args)!=7) {
  stop("Correct usage is Rscript Kmap.tetraploid.ST.R ChrSizeFile STCountTable STSizeNumber KmerCountNumber SubgenomeTable TukeyOutFile STNormCountOutFile.n", call.=FALSE)
} 

chr.size<-as.matrix(read.table(args[1],row=1,header=F));
ST.count<-as.matrix(read.table(args[2],row=1,header=T));
ST.size<-as.numeric(args[3]);
subgenome1<-as.character(read.table(args[5])$V1);
subgenome2<-as.character(read.table(args[5])$V2);
chr.st.size<-chr.size;
chr.st.size<-chr.st.size-(2*as.numeric(ST.size));
#my.IT.list<-c();
#my.ST.list<-c();
my.SG1.ST.list<-c();my.SG2.ST.list<-c();my.SG1.IT.list<-c();my.SG2.IT.list<-c();
ST.norm<-ST.count;

for(i in colnames(ST.norm)){if(endsWith(i,'IT')=="FALSE"){ST.norm[,i]<-ST.norm[,i]/ST.size;if(grepl("TRUE",paste(startsWith(i,prefix=subgenome1),collapse=""))=="TRUE"){my.SG1.ST.list<-c(i,my.SG1.ST.list);}else if(grepl("TRUE",paste(startsWith(i,prefix=subgenome2),collapse=""))=="TRUE"){my.SG2.ST.list<-c(i,my.SG2.ST.list);}}else{my.chr<-paste(head(as.character(strsplit(i,split="_")[[1]]),-1),collapse="_");ST.norm[,i]<-ST.norm[,i]/chr.st.size[my.chr,1];if(grepl("TRUE",paste(startsWith(i,prefix=subgenome1),collapse=""))=="TRUE"){my.SG1.IT.list<-c(i,my.SG1.IT.list);}else if(grepl("TRUE",paste(startsWith(i,prefix=subgenome2),collapse=""))=="TRUE"){my.SG2.IT.list<-c(i,my.SG2.IT.list);};}}


ST.tukey<-matrix(nrow=length(rownames(ST.norm)),ncol=38);
colnames(ST.tukey)<-c("Group.df","Group.Sumsq","Group.MeanSq","Residual.df","Residual.Sumsq","Residual.MeanSq","Fvalue","Pr","SG1.ST-SG1.IT.diff","SG1.ST-SG1.IT.p","SG1.ST-SG1.IT.95l","SG1.ST-SG1.IT.95u","SG1.ST-SG1.IT.bonf","SG2.ST-SG1.IT.diff","SG2.ST-SG1.IT.p","SG2.ST-SG1.IT.95l","SG2.ST-SG1.IT.95u","SG2.ST-SG1.IT.bonf","SG2.ST-SG2.IT.diff","SG2.ST-SG2.IT.p","SG2.ST-SG2.IT.95l","SG2.ST-SG2.IT.95u","SG2.ST-SG2.IT.bonf","SG2.IT-SG1.IT.diff","SG2.IT-SG1.IT.p","SG2.IT-SG1.IT.95l","SG2.IT-SG1.IT.95u","SG2.IT-SG1.IT.bonf","SG2.IT-SG1.ST.diff","SG2.IT-SG1.ST.p","SG2.IT-SG1.ST.95l","SG2.IT-SG1.ST.95u","SG2.IT-SG1.ST.bonf","SG2.ST-SG1.ST.diff","SG2.ST-SG1.ST.p","SG2.ST-SG1.ST.95l","SG2.ST-SG1.ST.95u","SG2.ST-SG1.ST.bonf");
rownames(ST.tukey)<-rownames(ST.norm);



for(i in rownames(ST.tukey)){



my.iter1.dat<-data.frame(y=character(),Group=character());my.iter1.dat<-rbind(my.iter1.dat,data.frame(y=ST.norm[i,my.SG1.IT.list],Group="SG1.IT"));my.iter1.dat<-rbind(my.iter1.dat,data.frame(y=ST.norm[i,my.SG1.ST.list],Group="SG1.ST"));my.iter1.dat<-rbind(my.iter1.dat,data.frame(y=ST.norm[i,my.SG2.IT.list],Group="SG2.IT"));my.iter1.dat<-rbind(my.iter1.dat,data.frame(y=ST.norm[i,my.SG2.ST.list],Group="SG2.ST"));

fit1<-lm(y ~ Group, my.iter1.dat);

my.anova<-anova(fit1)

ST.tukey[i,1]<-my.anova["Group","Df"];ST.tukey[i,2]<-my.anova["Group","Sum Sq"];ST.tukey[i,3]<-my.anova["Group","Mean Sq"];ST.tukey[i,4]<-my.anova["Residuals","Df"];ST.tukey[i,5]<-my.anova["Residuals","Sum Sq"];ST.tukey[i,6]<-my.anova["Residuals","Mean Sq"];ST.tukey[i,7]<-my.anova["Group","F value"];ST.tukey[i,8]<-my.anova["Group","Pr(>F)"];

my.aov<-aov(fit1);

my.tukey<-TukeyHSD(my.aov);

ST.tukey[i,9]<-my.tukey$Group["SG1.ST-SG1.IT","diff"];ST.tukey[i,10]<-my.tukey$Group["SG1.ST-SG1.IT","p adj"];ST.tukey[i,11]<-my.tukey$Group["SG1.ST-SG1.IT","lwr"];ST.tukey[i,12]<-my.tukey$Group["SG1.ST-SG1.IT","upr"];ST.tukey[i,13]<-p.adjust(p=as.numeric(ST.tukey[i,10]),method="bonferroni",n=as.numeric(args[4]));
ST.tukey[i,14]<-my.tukey$Group["SG2.ST-SG1.IT","diff"];ST.tukey[i,15]<-my.tukey$Group["SG2.ST-SG1.IT","p adj"];ST.tukey[i,16]<-my.tukey$Group["SG2.ST-SG1.IT","lwr"];ST.tukey[i,17]<-my.tukey$Group["SG2.ST-SG1.IT","upr"];ST.tukey[i,18]<-p.adjust(p=as.numeric(ST.tukey[i,15]),method="bonferroni",n=as.numeric(args[4]));
ST.tukey[i,19]<-my.tukey$Group["SG2.ST-SG2.IT","diff"];ST.tukey[i,20]<-my.tukey$Group["SG2.ST-SG2.IT","p adj"];ST.tukey[i,21]<-my.tukey$Group["SG2.ST-SG2.IT","lwr"];ST.tukey[i,22]<-my.tukey$Group["SG2.ST-SG2.IT","upr"];ST.tukey[i,23]<-p.adjust(p=as.numeric(ST.tukey[i,20]),method="bonferroni",n=as.numeric(args[4]))
ST.tukey[i,24]<-my.tukey$Group["SG2.IT-SG1.IT","diff"];ST.tukey[i,25]<-my.tukey$Group["SG2.IT-SG1.IT","p adj"];ST.tukey[i,26]<-my.tukey$Group["SG2.IT-SG1.IT","lwr"];ST.tukey[i,27]<-my.tukey$Group["SG2.IT-SG1.IT","upr"];ST.tukey[i,28]<-p.adjust(p=as.numeric(ST.tukey[i,25]),method="bonferroni",n=as.numeric(args[4]))
ST.tukey[i,29]<-my.tukey$Group["SG2.IT-SG1.ST","diff"];ST.tukey[i,30]<-my.tukey$Group["SG2.IT-SG1.ST","p adj"];ST.tukey[i,31]<-my.tukey$Group["SG2.IT-SG1.ST","lwr"];ST.tukey[i,32]<-my.tukey$Group["SG2.IT-SG1.ST","upr"];ST.tukey[i,33]<-p.adjust(p=as.numeric(ST.tukey[i,30]),method="bonferroni",n=as.numeric(args[4]))
ST.tukey[i,34]<-my.tukey$Group["SG2.ST-SG1.ST","diff"];ST.tukey[i,35]<-my.tukey$Group["SG2.ST-SG1.ST","p adj"];ST.tukey[i,36]<-my.tukey$Group["SG2.ST-SG1.ST","lwr"];ST.tukey[i,37]<-my.tukey$Group["SG2.ST-SG1.ST","upr"];ST.tukey[i,38]<-p.adjust(p=as.numeric(ST.tukey[i,35]),method="bonferroni",n=as.numeric(args[4]))


}

write.table(ST.tukey,file=args[6],row=T,col=T,quote=F,sep="\t");
write.table(ST.norm,file=args[7],row=T,col=T,quote=F,sep="\t");
#ST.sig<-subset(ST.tukey,ST.tukey[,13]<=0.05);

