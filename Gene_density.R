## This script generate the graphs of genes densities accross chromosomes

library(ggplot2)
test = read.table("GeneDensityZamb",header=TRUE,sep=",")
Chrom2L <-subset(test,Chromosome=="X_Chromosome")
PLotDensity = ggplot(data=Chrom2L, aes(x=Chrom2L$Segment, y=Chrom2L$MaxDs, fill=Chrom2L$NbGenes)) + theme_classic() +
  geom_bar(stat="identity")+geom_line(aes(y=Chrom2L$Ds), size=1, color="black")
PLotDensity
mid<-mean(Chrom2L$NbGenes)
PLotDensity+scale_fill_gradient2(midpoint=mid, low="yellow", mid="pink", high="red") 
mid<-mean(Chrom2L$NbGenes)
PLotDensity+scale_fill_gradient2(midpoint=mid, low="yellow", high="red") 
PLotDensity+scale_fill_gradientn(colours=rainbow(5))
mid<-mean(Chrom2L$NbGenes)
PLotDensity+scale_fill_gradient2(midpoint=mid, low="aquamarine", mid="chartreuse", high="red") 

