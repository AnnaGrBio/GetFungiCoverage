#import cv2
import os
os.chdir("path to folder")
import random


""" docstring

    This program search all established genes in 7 genomes, retrieve the one in commun, and store all of them in orthology file
    between the 7 lines

    input
    ------------
    - gtf file of de novo assemblies of 7 genomes
    - cds files of all proteins
    
    outpur
    ------------  
    - FASTA files with orthogroups of established genes
"""


# Build the list with corresponding proteins
def open_file(NameFile):
    F = open(NameFile, "r")
    L = F.readlines()
    return L  


def access_protein_name(File):
    # THis function find name of established proteins in a gtf file
    ListeLongName = []
    ListeShortName = []
    for i in File: 
        LongName = i.split(" ")[0]
        if LongName[0] == ">":
            LongName = LongName[1:]
            ShortName = LongName.split(".")[0]
            ListeLongName.append(LongName)
            ListeShortName.append(ShortName)
    return ListeLongName,ListeShortName


def find_genes_in_commun(L1,L2,L3,L4,L5,L6,L7,L8):
    # This function search for established genes presents in the 7 lines (normaly most of them)
    ListeNamesCommun = []
    for i in L1:
        if i in L2 and i in L3 and i in L4 and i in L5 and i in L6 and i in L7 and i in L8:
            ListeNamesCommun.append(i)
    print (len(ListeNamesCommun))
    return ListeNamesCommun


# Make All files for alignments
def retrieve_proteins_and_rename_them(File):
    # This function retrieve the protein sequence of all stablished genes in each line and rename them
    Dico = {}
    GeneName = ""
    Prot = ""
    for i in File:
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            if GeneName!="":
                Dico[GeneName] = Prot[0:len(Prot)-1]
                Prot = ""
            LongName = ligne.split(" ")[0]
            LongName = LongName[1:]
            GeneName = LongName
        else:
            Prot+=ligne
    Dico[GeneName] = Prot[0:len(Prot)-1]
    return Dico


def build_final_files(DicoProtRef, DicoProtAK5, DicoProtDK5, DicoProtGI5, DicoProtSW5, DicoProtUM, DicoProtYE, DicoProtZamb, ListeCommunLongNames):
    # This function build the final groups of orthology
    ListeNameFiles = []
    for i in ListeCommunLongNames:
        NameFile = "TheCDS_"+i+".fa"
        F = open(NameFile, "w")
        F.write(">"+i+"__"+"REF"+"\n")
        F.write(DicoProtRef[i]+"\n")
        
        F.write(">"+i+"__"+"AK5"+"\n")
        F.write(DicoProtAK5[i]+"\n")
        
        F.write(">"+i+"__"+"DK5"+"\n")
        F.write(DicoProtDK5[i]+"\n")

        F.write(">"+i+"__"+"GI5"+"\n")
        F.write(DicoProtGI5[i]+"\n")
        
        F.write(">"+i+"__"+"SW5"+"\n")
        F.write(DicoProtSW5[i]+"\n")
        
        F.write(">"+i+"__"+"UM"+"\n")
        F.write(DicoProtUM[i]+"\n")
        
        F.write(">"+i+"__"+"YE"+"\n")
        F.write(DicoProtYE[i]+"\n")
        
        F.write(">"+i+"__"+"Zamb"+"\n")
        F.write(DicoProtZamb[i]+"\n")
        F.close()
        ListeNameFiles.append(NameFile)
    F = open("AllFileNames.txt", "w")
    for i in ListeNameFiles:
        F.write(i+"\n")
    F.close()
            

def modify_reference(DataRef):
    NewFile = []
    for i in DataRef:
        if i[0] == ">":
            ligne = i.split(" ")
            Debut = ligne[0]
            End = ligne[1]
            NewDebut = Debut+"_R0"
            NewLine = NewDebut+" "+End
            NewFile.append(NewLine)
            print (NewLine)
        else:
            NewFile.append(i)
    return NewFile
        

def main_function():
    DataRef = open_file("Ref_CDS.fa")
    DataRefGood = modify_reference(DataRef)
    ListeLongNameRef,ListeShortNameRef = access_protein_name(DataRefGood)
    DataAK5 = open_file("FI_CDS.fa")
    ListeLongNameAK5,ListeShortNameAK5 = access_protein_name(DataAK5)
    DataDK5 = open_file("DK_CDS.fa")
    ListeLongNameDK5,ListeShortNameDK5 = access_protein_name(DataDK5)
    DataGI5 = open_file("ES_CDS.fa")
    ListeLongNameGI5,ListeShortNameGI5 = access_protein_name(DataGI5)
    DataSW5 = open_file("SE_CDS.fa")
    ListeLongNameSW5,ListeShortNameSW5 = access_protein_name(DataSW5)
    DataUM = open_file("UA_CDS.fa")
    ListeLongNameUM,ListeShortNameUM = access_protein_name(DataUM)
    DataYE = open_file("TR_CDS.fa")
    ListeLongNameYE,ListeShortNameYE = access_protein_name(DataYE)
    DataZamb = open_file("ZI_CDS.fa")
    ListeLongNameZamb,ListeShortNameZamb = access_protein_name(DataZamb)
    ListeCommunLongNames = find_genes_in_commun(ListeLongNameRef,ListeLongNameAK5,ListeLongNameDK5,ListeLongNameGI5,ListeLongNameSW5,ListeLongNameUM,ListeLongNameYE,ListeLongNameZamb)
    #find_genes_in_commun(ListeShortNameAK5,ListeShortNameDK5,ListeShortNameGI5,ListeShortNameSW5,ListeShortNameUM,ListeShortNameYE,ListeShortNameZamb)
    DicoProtRef = retrieve_proteins_and_rename_them(DataRefGood)
    DicoProtAK5 = retrieve_proteins_and_rename_them(DataAK5)
    DicoProtDK5 = retrieve_proteins_and_rename_them(DataDK5)
    DicoProtGI5 = retrieve_proteins_and_rename_them(DataGI5)
    DicoProtSW5 = retrieve_proteins_and_rename_them(DataSW5)
    DicoProtUM = retrieve_proteins_and_rename_them(DataUM)
    DicoProtYE = retrieve_proteins_and_rename_them(DataYE)
    DicoProtZamb = retrieve_proteins_and_rename_them(DataZamb)
    build_final_files(DicoProtRef, DicoProtAK5, DicoProtDK5, DicoProtGI5, DicoProtSW5, DicoProtUM, DicoProtYE, DicoProtZamb, ListeCommunLongNames)





