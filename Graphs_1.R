Data1 = read.table("/global/students/homes/agrandch/DeNovoGenes/NewVersion/MakeOrthogroups/Step4MakeSomeStatistics/Data1", header = TRUE, sep=" ")
Nb_ORF_in_orthogroup = Data1$Nb_ORF_in_orthogroup
Nb_Orthogroups = Data1$Nb_Orthogroups
plot(Nb_ORF_in_orthogroup, Nb_Orthogroups, main = "Size of orthogroups", xlab = "Nb ORF in orthogroup", ylab = "Nb orthogroups")


Data2 = read.table("/global/students/homes/agrandch/DeNovoGenes/NewVersion/MakeOrthogroups/Step4MakeSomeStatistics/Data2", header = TRUE, sep=" ")
Nb_pops_in_orthogroup = Data2$Nb_pops_in_orthogroup
Nb_Orthogroups = Data2$Nb_Orthogroups
plot(Nb_pops_in_orthogroup, Nb_Orthogroups, main = "Size of orthogroups", xlab = "Nb pops in orthogroup", ylab = "Nb orthogroups")

Data3 = read.table("/global/students/homes/agrandch/DeNovoGenes/NewVersion/MakeOrthogroups/Step4MakeSomeStatistics/Data3", header = TRUE, sep=",")
NbPop = Data3$NbPopsSharingOrthogroup
Shared_1_Pop <- subset(Data3, NbPop == "ORF_1_Pop")
Shared_2_Pop <- subset(Data3, NbPop == "ORF_2_Pop")
Shared_3_Pop <- subset(Data3, NbPop == "ORF_3_Pop")
Shared_4_Pop <- subset(Data3, NbPop == "ORF_4_Pop")
Shared_5_Pop <- subset(Data3, NbPop == "ORF_5_Pop")
Shared_6_Pop <- subset(Data3, NbPop == "ORF_6_Pop")
Shared_7_Pop <- subset(Data3, NbPop == "ORF_7_Pop")

TPM_1 = Shared_1_Pop$TPM
TPM_2 = Shared_2_Pop$TPM
TPM_3 = Shared_3_Pop$TPM
TPM_4 = Shared_4_Pop$TPM
TPM_5 = Shared_5_Pop$TPM
TPM_6 = Shared_6_Pop$TPM
TPM_7 = Shared_7_Pop$TPM

boxplot(TPM_1, TPM_2, TPM_3, TPM_4, TPM_5, TPM_6, TPM_7)
mean(TPM_1)
mean(TPM_2)
mean(TPM_3)
mean(TPM_4)
mean(TPM_5)
mean(TPM_6)
mean(TPM_7)

Size_1 = Shared_1_Pop$SplicedORFsize
Size_2 = Shared_2_Pop$SplicedORFsize
Size_3 = Shared_3_Pop$SplicedORFsize
Size_4 = Shared_4_Pop$SplicedORFsize
Size_5 = Shared_5_Pop$SplicedORFsize
Size_6 = Shared_6_Pop$SplicedORFsize
Size_7 = Shared_7_Pop$SplicedORFsize

boxplot(Size_1, Size_2, Size_3, Size_4, Size_5, Size_6, Size_7)
mean(Size_1)
mean(Size_2)
mean(Size_7)
















