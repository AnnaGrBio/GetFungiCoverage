import cv2
import os
os.chdir("path to folder")
import random
from collections import Counter


""" docstring

    This program creates formated file for accessing the number of lines and proto-genes per orthogroup, as well as the genomic location

    input
    ------------
    - orthogroups file
    
    output
    ------------
    - three formated file to be processed by R

"""


def open_file(NameFile):
    F = open(NameFile, "r")
    L = F.readlines()
    return L  


### Step1 Get number of genes by categories
def f1_nb_orthogroups_vs_nb_gene(File):
    # This function access the size oforthogroups and the number of genes inside
    Dico = {}
    c = 0
    for i in File:
        ligne = i.split("\n")[0]
        Group = ligne.split(" ")
        SizeOrthogroup = len(Group)-1
        if SizeOrthogroup in Dico.keys():
            Dico[SizeOrthogroup] += 1
        else:
            Dico[SizeOrthogroup] = 1
    return Dico


def build_f1_nbgene_nborthogroup(Dico):
    # This function write a file that contain information of genes number in orthogroups
    F = open("Data1", "w")
    F.write("Nb_ORF_in_orthogroup"+" "+"Nb_Orthogroups"+"\n")
    for i in range(1,1000):
        if i in Dico.keys():
            F.write(str(i)+" "+str(Dico[i])+"\n")
    F.close()


# Step2 Nb orthogroup shared by Nb species
def f2_shared_gene(File):
    # This function analyse how many genes are shared from a same line in one orthogroup
    Dico = {}
    for i in File:
        ligne = i.split("\n")[0]
        Group = ligne.split(" ")
        ListePops = []
        for i in Group[1:]:
            PopName = i.split("_")[0]
            if PopName not in ListePops:
                ListePops.append(PopName)
        NbSpecies = len(ListePops)
        if NbSpecies in Dico.keys():
            Dico[NbSpecies] += 1
        else:
            Dico[NbSpecies] = 1
    return Dico


def build_f2_nb_line_nb_orthogroup(Dico):
    # This function write the final file of the number of populations found in one orthogroup
    F = open("Data2", "w")
    F.write("Nb_pops_in_orthogroup"+" "+"Nb_Orthogroups"+"\n")
    for i in range(1,8):
        if i in Dico.keys():
            F.write(str(i) + " " + str(Dico[i]) + "\n")
    F.close()


# Step 3 DataByGroups
def get_orf_info(InfoFile,ORF):
    # This function retrieves the informations about proto-genes in an info file
    L = []
    for i in InfoFile:
        ligne = i.split("\n")[0]
        if ORF in ligne:
            Data = ligne.split(",")
            NbExon = Data[3]
            Chrom = Data[4]
            TPM = Data[10]
            Start = Data[13]
            Stop = Data[14]
            GenomicPosition = Data[15]
            L = [ORF, NbExon, Chrom, TPM, GenomicPosition, Start, Stop]
            break
    return L
            
def retrive_all_lines_info(Orthogroups, InfoFile):
    # This function calculate the info of proto-genes in orthogroups
    Liste1 = []
    Liste2 = []
    Liste3 = []
    Liste4 = []
    Liste5 = []
    Liste6 = []
    Liste7 = []
    for i in Orthogroups:
        ligne = i.split("\n")[0]
        Group = ligne.split(" ")
        ListePops = []
        NbORFinOrthogroup = 0
        NameOrthogroup = Group[0]
        SubLists = []
        for j in Group[1:]:
            NbORFinOrthogroup+=1
            Info_ORF = get_orf_info(InfoFile,j)
            PopName = j.split("_")[0]
            SplicedName = j.split("_")
            ORFsize = int(SplicedName[4]) - int((SplicedName[3])) + 4
            Info_ORF.append(ORFsize)
            if PopName not in ListePops:
                ListePops.append(PopName)
            SubLists.append(Info_ORF)
        NbSpecies = len(ListePops)
        for j in SubLists:
            j.append(NbORFinOrthogroup)
        if NbSpecies == 1:
            for j in SubLists:
                Liste1.append(j)
        if NbSpecies == 2:
            for j in SubLists:
                Liste2.append(j)
        if NbSpecies == 3:
            for j in SubLists:
                Liste3.append(j)
        if NbSpecies == 4:
            for j in SubLists:
                Liste4.append(j)
        if NbSpecies == 5:
            for j in SubLists:
                Liste5.append(j)
        if NbSpecies == 6:
            for j in SubLists:
                Liste6.append(j)
        if NbSpecies == 7:
            for j in SubLists:
                Liste7.append(j)
    return Liste1, Liste2, Liste3, Liste4, Liste5, Liste6, Liste7


def build_f3_all_info(Liste1, Liste2, Liste3, Liste4, Liste5, Liste6, Liste7):
    # This function writes the final file of info of proto-genes
    F = open("Data3", "w")
    F.write("NbPopsSharingOrthogroup" + "," + "ORFname,NbExon,Chrom,TPM,GenomicPosition,Start,Stop,SplicedORFsize,NbORFinOrthogroup" + "\n")
    for i in Liste1:
        F.write("ORF_1_Pop")
        for j in i:
            F.write(","+str(j))
        F.write("\n")
    for i in Liste2:
        F.write("ORF_2_Pop")
        for j in i:
            F.write(","+str(j))
        F.write("\n")
    for i in Liste3:
        F.write("ORF_3_Pop")
        for j in i:
            F.write(","+str(j))
        F.write("\n")
    for i in Liste4:
        F.write("ORF_4_Pop")
        for j in i:
            F.write(","+str(j))
        F.write("\n")
    for i in Liste5:
        F.write("ORF_5_Pop")
        for j in i:
            F.write(","+str(j))
        F.write("\n")
    for i in Liste6:
        F.write("ORF_6_Pop")
        for j in i:
            F.write(","+str(j))
        F.write("\n")
    for i in Liste7:
        F.write("ORF_7_Pop")
        for j in i:
            F.write(","+str(j))
        F.write("\n")
    F.close()


def main_function():
    Orthogroups = open_file("OrthogroupsSorted.txt")
    # Make Data1
    DicoData1 = f1_nb_orthogroups_vs_nb_gene(Orthogroups)
    MakeFile1(DicoData1)
    # Make Data2
    DicoData2 = f2_shared_gene(Orthogroups)
    build_f2_nb_line_nb_orthogroup(DicoData2)
    # Make Data3
    InfoFile = open_file("DeNovo_InfoFile")
    Liste1, Liste2, Liste3, Liste4, Liste5, Liste6, Liste7 = retrive_all_lines_info(Orthogroups, InfoFile)
    build_f3_all_info(Liste1, Liste2, Liste3, Liste4, Liste5, Liste6, Liste7)







