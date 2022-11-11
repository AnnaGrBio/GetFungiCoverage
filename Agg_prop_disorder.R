library (ggplot2)

## This script handle the output of aggregation propensity in proteins and calculate the lm model
Data1 = read.table("FinalDataAggregation.txt", sep = ",", header = TRUE)
NbPops = Data1$NbPops
Aggregation = Data1$Aggregation
PlotAggreg <-ggplot(Data1, aes(NbPops,Aggregation))+
  theme_classic()+
  labs(title = "Aggregation propensity")+theme(plot.title = element_text(hjust=0.5))+geom_boxplot(color="blue4", outlier.shape=NA)
PlotAggreg
summary(lm(Aggregation~NbPops))

Pop1 = subset(Data1,NbPops=="ORF_1_Pop")
AggPop1 <- Pop1$Aggregation
mean(AggPop1)
Pop2 = subset(Data1,NbPops=="ORF_2_Pop")
AggPop2 <- Pop2$Aggregation
mean(AggPop2)
Pop3 = subset(Data1,NbPops=="ORF_3_Pop")
AggPop3 <- Pop3$Aggregation
mean(AggPop3)
Pop4 = subset(Data1,NbPops=="ORF_4_Pop")
AggPop4 <- Pop4$Aggregation
mean(AggPop4)
Pop5 = subset(Data1,NbPops=="ORF_5_Pop")
AggPop5 <- Pop5$Aggregation
mean(AggPop5)
Pop6 = subset(Data1,NbPops=="ORF_6_Pop")
AggPop6 <- Pop6$Aggregation
mean(AggPop6)
Pop7 = subset(Data1,NbPops=="ORF_7_Pop")
AggPop7 <- Pop7$Aggregation
mean(AggPop7)

##############################################################

## This script handle the output of disorder in proteins and calculate the lm model
Data2 = read.table("FinalDataDisorder.txt", sep = ",", header = TRUE)
NbPops = Data2$NbPops
Disorder = Data2$Disorder
PlotDisorder <-ggplot(Data2, aes(NbPops,Disorder))+
  theme_classic()+
  labs(title = "Intrinsic disorder")+theme(plot.title = element_text(hjust=0.5))+geom_boxplot(color="blue4", outlier.shape=NA)
summary(lm(Disorder~NbPops))
PlotDisorder

Pop1 = subset(Data2,NbPops=="ORF_1_Pop")
DisorderPop1 <- Pop1$Disorder
mean(DisorderPop1)
Pop2 = subset(Data2,NbPops=="ORF_2_Pop")
DisorderPop2 <- Pop2$Disorder
mean(DisorderPop2)
Pop3 = subset(Data2,NbPops=="ORF_3_Pop")
DisorderPop3 <- Pop3$Disorder
mean(DisorderPop3)
Pop4 = subset(Data2,NbPops=="ORF_4_Pop")
DisorderPop4 <- Pop4$Disorder
mean(DisorderPop4)
Pop5 = subset(Data2,NbPops=="ORF_5_Pop")
DisorderPop5 <- Pop5$Disorder
mean(DisorderPop5)
Pop6 = subset(Data2,NbPops=="ORF_6_Pop")
DisorderPop6 <- Pop6$Disorder
mean(DisorderPop6)
Pop7 = subset(Data2,NbPops=="ORF_7_Pop")
DisorderPop7 <- Pop7$Disorder
mean(DisorderPop7)

#########################################"
## This script handle the output of HCA new domains in proteins and calculate the lm model

Data3 = read.table("path/zFinalDataHCA.txt", sep = ",", header = TRUE)
NbPops = Data3$NbPops
HCA = Data3$HCA
PlotHCA <-ggplot(Data3, aes(NbPops,HCA))+
  theme_classic()+
  labs(title = "Hydrophobic clusters")+theme(plot.title = element_text(hjust=0.5))+geom_boxplot(color="blue4", outlier.shape=NA)
PlotHCA


