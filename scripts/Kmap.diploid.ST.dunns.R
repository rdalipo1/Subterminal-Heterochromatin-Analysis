
library(dunn.test)

args<- commandArgs(trailingOnly = TRUE);
if (length(args)!=6) {
  stop("Correct usage is Rscript Kmap.diploid.ST.R ChrSizeFile STCountTable STSizeNumber KmerCountNumber DunnsOutFile STNormCountOutFile.n", call.=FALSE)
} 

chr.size<-as.matrix(read.table(args[1],row=1,header=F));
ST.count<-as.matrix(read.table(args[2],row=1,header=T));
ST.size<-as.numeric(args[3]);
chr.st.size<-chr.size;
chr.st.size<-chr.st.size-(2*as.numeric(ST.size));
my.IT.list<-c();
my.ST.list<-c();
ST.norm<-ST.count;

for(i in colnames(ST.norm)){
    if(endsWith(i,'IT')=="FALSE"){
        ST.norm[,i]<-ST.norm[,i]/ST.size;
        my.ST.list<-c(i,my.ST.list);
        }
    else{
        my.chr<-as.character(strsplit(i,split="_")[[1]][1]);
        ST.norm[,i]<-ST.norm[,i]/chr.st.size[my.chr,1];
        my.IT.list<-c(i,my.IT.list);
        }
}


ST.dunns<-matrix(nrow=length(rownames(ST.norm)),ncol=6);
colnames(ST.dunns)<-c("Group.df","Chi.squared","P","IT-ST.z","IT-ST.p","IT-ST.bonf");
rownames(ST.dunns)<-rownames(ST.norm);


for(i in rownames(ST.dunns)){

    my.iter1.dat<-data.frame(y=character(),Group=character());
    my.iter1.dat<-rbind(my.iter1.dat,data.frame(y=ST.norm[i,my.IT.list],Group="IT"));
    my.iter1.dat<-rbind(my.iter1.dat,data.frame(y=ST.norm[i,my.ST.list],Group="ST"));
    
    my.kruskal <- kruskal.test(y ~ Group, my.iter1.dat);

    ST.dunns[i,1]<-my.kruskal$parameter
    ST.dunns[i,2]<-my.kruskal$statistic
    ST.dunns[i,3]<-my.kruskal$p.value

    capture.output(my.dunns <- dunn.test(my.iter1.dat$y, my.iter1.dat$Group, kw = FALSE), file = NULL)
    my.dunns <- dunn.test(my.iter1.dat$y, my.iter1.dat$Group)

    ST.dunns[i,4] <- my.dunns$Z
    ST.dunns[i,5] <- my.dunns$P
    

}

ST.dunns[,6]<-p.adjust(p=as.numeric(ST.dunns[,5]),method="bonferroni",n=as.numeric(args[4]));
write.table(ST.dunns,file=args[5],row=T,col=T,quote=F,sep="\t");
write.table(ST.norm,file=args[6],row=T,col=T,quote=F,sep="\t");

