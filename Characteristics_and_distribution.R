######## Data1. Nb ORFs in orthogroups vs nb orthogroups
# This script produces one graph of the figure 2

Data1 = read.table("Data1",header=TRUE,sep=" ")
NbORF = Data1$Nb_ORF_in_orthogroup
NbOrthogroups = Data1$Nb_Orthogroups
library(ggplot2)
p1 <- ggplot(Data1, aes(NbORF,NbOrthogroups)) +
  geom_bar(stat="identity") +theme_classic()+
  labs(title = "Nb ORF in orthogroups")+theme(plot.title = element_text(hjust = 0.5))+scale_color_manual(values=c("blue"))
p1

######## Data2. Data 2 Nb lines sharing an orthogroup
# This script produces one graph of the figure 2
Data2 = read.table("Data2",header=TRUE,sep=" ")
NbPops = Data2$Nb_pops_in_orthogroup
NbOrthogroups = Data2$Nb_Orthogroups
p2 <- ggplot(Data2, aes(NbPops, NbOrthogroups)) +
  geom_point() +theme_classic()+geom_line()+
  labs(title = "Nb orthogroups shared by populations")+
  theme(plot.title = element_text(hjust = 0.5))+scale_x_continuous(breaks=NbPops)
p2

######## Data3. Nb proto-genes specific and non-specific per lines
# This script produces one graph of the figure 2
DataPlot = read.table("Data3",header=TRUE,sep=",")
Pop=DataPlot$Population
Status = DataPlot$Status
Nb = DataPlot$Number
pORF <- ggplot(DataPlot, aes(x=Pop, y=Nb, fill=Status)) +
  geom_bar(stat="identity") +theme_classic()+
  scale_fill_manual(values=c("lemonchiffon2", "red", "chocolate4","darkorange",  "yellow", "chartreuse", "cyan", "blue4"))+
  labs(title = "Number of de novo ORF per population")+
  theme(plot.title = element_text(hjust = 0.5))
pORF

######## Data4. Percentage of proto-gene at each genomic position
# This script produces one graph of the figure 2
DataPlot2 = read.table("Data4",header=TRUE,sep=",")
Population = DataPlot2$Population
GenomicPosition = DataPlot2$GenomicPosition
Number = DataPlot2$Number
pORF2 <- ggplot(DataPlot2, aes(x=Population, y=Number, fill=GenomicPosition)) +
  geom_bar(stat="identity") +theme_classic()+
  scale_fill_manual(values=c( "ivory4","khaki","khaki4", "honeydew3","khaki3"))+
  labs(title = "Genomic position of the de novo ORFs")+
  theme(plot.title = element_text(hjust = 0.5))
pORF2

########
# This script access proto-genes properties
Data3 = read.table("proto_genes_properties",header=TRUE,sep=",")
NbPops = Data3$NbPopsSharingOrthogroup
ORFsize = Data3$SplicedORFsize
TPM = Data3$TPM
NbExon = Data3$NbExon
GenomicPosition = Data3$GenomicPosition
library(ggplot2)

## Sizes
G1 <- ggplot(Data3, aes(NbPops, ORFsize)) +
  theme_classic() +
  labs(title = "Sizes of de novo ORFs")+theme(plot.title = element_text(hjust = 0.5))+geom_boxplot(color="blue4", outlier.shape = NA)
G1 + ylim(0,500)

Size = Data3$SplicedORFsize
mean(Size)

Pop1 = subset(Data3,NbPopsSharingOrthogroup=="ORF_1_Pop")
mean(Pop1$SplicedORFsize)
SizePop1 =Pop1$SplicedORFsize
Pop2 = subset(Data3,NbPopsSharingOrthogroup=="ORF_2_Pop")
mean(Pop2$SplicedORFsize)
SizePop2 =Pop2$SplicedORFsize
Pop3 = subset(Data3,NbPopsSharingOrthogroup=="ORF_3_Pop")
mean(Pop3$SplicedORFsize)
SizePop3 =Pop3$SplicedORFsize
Pop4 = subset(Data3,NbPopsSharingOrthogroup=="ORF_4_Pop")
mean(Pop4$SplicedORFsize)
SizePop4 =Pop4$SplicedORFsize
Pop5 = subset(Data3,NbPopsSharingOrthogroup=="ORF_5_Pop")
mean(Pop5$SplicedORFsize)
SizePop5 =Pop5$SplicedORFsize
Pop6 = subset(Data3,NbPopsSharingOrthogroup=="ORF_6_Pop")
mean(Pop6$SplicedORFsize)
SizePop6 =Pop6$SplicedORFsize
Pop7 = subset(Data3,NbPopsSharingOrthogroup=="ORF_7_Pop")
mean(Pop7$SplicedORFsize)
SizePop7 =Pop7$SplicedORFsize
SizePop1_6 = c(SizePop1,SizePop2,SizePop3,SizePop4,SizePop5,SizePop6)
mean(SizePop1_6)
t.test(SizePop7,SizePop1_6)


## TPM
G2 <- ggplot(Data3, aes(NbPops, TPM)) +
  theme_classic() +
  stat_boxplot(fill = NA) +
  labs(title = "TPM of de novo ORFs")+theme(plot.title = element_text(hjust = 0.5))+geom_boxplot(color="black")
G2+ylim(0,100)

mean(Pop1$TPM)
mean(Pop2$TPM)
mean(Pop3$TPM)
mean(Pop4$TPM)
mean(Pop5$TPM)
mean(Pop6$TPM)
mean(Pop7$TPM)
Tpm7 = Pop7$TPM
TpmOther = c(Pop1$TPM, Pop2$TPM, Pop3$TPM, Pop4$TPM, Pop6$TPM)
mean(TpmOther)
t.test(Tpm7, TpmOther)
library(ggplot2)
G3 <- ggplot(Data3, aes(GenomicPosition, NbExon)) +
  theme_classic() +
  labs(title = "Number of exons in the transcript")+theme(plot.title = element_text(hjust = 0.5))+geom_boxplot(color="blue4", outlier.shape= NA)
G3

ExonInside = subset(Data3,GenomicPosition=="ExonInside")
ExonLonger = subset(Data3,GenomicPosition=="ExonLonger")
Intergenic = subset(Data3,GenomicPosition=="Intergenic")
Intronic = subset(Data3,GenomicPosition=="Intronic")
ReverseGenic = subset(Data3,GenomicPosition=="ReverseGenic")
mean(ExonInside$NbExon)
mean(ExonLonger$NbExon)
mean(Intergenic$NbExon)
mean(Intronic$NbExon)
mean(ReverseGenic$NbExon)
AvAllButInside = c(ExonLonger$NbExon, Intergenic$NbExon, Intronic$NbExon, ReverseGenic$NbExon)
t.test(AvAllButInside, ExonInside$NbExon)
Longer_Reverse = c(ExonLonger$NbExon, ReverseGenic$NbExon)
Intergenic_Intronic = c(Intergenic$NbExon, Intronic$NbExon)
t.test(Longer_Reverse, Intergenic_Intronic)


ExonInExonInside <- ExonInside$NbExon
ExonInExonLonger <- ExonLonger$NbExon
ExonInIntergenic <- Intergenic$NbExon
ExonInIntronic <- Intronic$NbExon
ExonInReverseGenic <- ReverseGenic$NbExon
t.test(ExonInReverseGenic,ExonInExonInside)

#introns
p4 <- ggplot(proto_genes_properties, aes(Position, NbIntrons, color = Position)) +
  geom_jitter() + theme_classic()+
  stat_boxplot(fill = NA) +
  labs(title = "Number of intron according to the position")+theme(plot.title = element_text(hjust = 0.5))+geom_boxplot(color="black")+scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
p4

p6 <- ggplot(Data4, aes(Position, IntronSize, color = Position)) +
  geom_jitter() +
  stat_boxplot(fill = NA) +
  labs(title = "Size of the intron according to the position")+theme(plot.title = element_text(hjust = 0.5))+geom_boxplot(color="black")+scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
p6

######## Data5 TE and proto-genes
# TEs overlap to proto-genes
TE = read.table("TEinORF.txt",header=TRUE,sep=",")
PoplationName = TE$Population
Status = TE$Status
Number = TE$Number
library(ggplot2)
BarplotTE <- ggplot(TE, aes(x=PoplationName, y=Number, fill=Status)) +
  geom_bar(stat="identity") +theme_classic()+
  scale_fill_manual(values=c( "tan4","tan2","darkolivegreen1"))+
  labs(title = "De novo ORF overlap to transposable elements")+
  theme(plot.title = element_text(hjust = 0.5))
BarplotTE

