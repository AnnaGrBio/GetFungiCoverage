import cv2
import os
os.chdir("path to folder")
import random
from collections import Counter


""" docstring

    THis program calculates intrinsic disorder in proteins

    input
    ------------
    - sequences of de novo genes

    output
    ------------
    - intrinsic disorder of de novo genes
"""

def open_file(NameFile):
    F = open(NameFile, "r")
    L = F.readlines()
    return L  


def build_dic_gene_sequences(File):
    # This function build a dictionary that contain the name of proto-genes associated to their sequence
    D = {}
    for i in File:
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            Name = ligne[1:]
        else:
            Seq = ligne
            D[Name] = ligne
    return D


def calculate_disorder_from_result_file():
    # This function calculate the disorder of a sequence according to the result provided by IUPRED
    Liste = []
    File = open_file("Rst.txt")
    for i in File:
        ligne = i.split()
        if ligne[0] != "#":
            Value = float(ligne[2])
            Liste.append(Value)
    Average = float(sum(Liste))/float(len(Liste))
    os.system('rm Rst.txt')
    return Average


def run_IUPRED(Dico):
    # This function runs IUPRED on each protein
    FinalDico = {}
    for i in Dico.keys():
        F = open("Seq", "w")
        Name = i
        Seq = Dico[i]
        F.write(">"+Name+"\n")
        F.write(Seq)
        F.close()
        os.system('python iupred2a/iupred2a.py Seq short > Rst.txt')
        AverageDisorder = calculate_disorder_from_result_file()
        FinalDico[Name] = AverageDisorder
    return FinalDico

def build_final_file1(FinalDico):
    # This function writes the first file with all disorder results
    F = open("IupredResults.txt", "w")
    for i in FinalDico.keys():
        F.write(i+","+str(FinalDico[i])+"\n")
    F.close()


def build_new_dico(Data3):
    # This function builds a new dictionarw that store intrinsic disotder and nb pop that the proto-gene share
    Translation = {"ORF_1_Pop":1, "ORF_2_Pop":2,"ORF_3_Pop":3,"ORF_4_Pop":4,"ORF_5_Pop":5,"ORF_6_Pop":6,"ORF_7_Pop":7}
    D = {}
    for i in Data3[1:]:
        ligne = i.split(",")
        NbPop = ligne[0]#str(Translation[ligne[0]])
        NameProtoGene = ligne[1]
        D[NameProtoGene] = NbPop
    return D


def build_final_file2(DataDisorder, Dico):
    # Thsi function build the second output file that gone be used for statistical analyses
    F = open("FinalDataDisorder.txt", "w")
    F.write("NbPops,Disorder"+"\n")
    for i in DataDisorder:
        ligne = i.split("\n")[0]
        Name = ligne.split(",")[0]
        Disorder = ligne.split(",")[1]
        if Name in Dico.keys():
            F.write(Dico[Name]+","+Disorder+"\n")
    F.close()
        

def main_function():
    File = open_file("AllPopsProtoGenesProt.fa")
    DicoSeqs = build_dic_gene_sequences(File)
    FinalDico = run_IUPRED(DicoSeqs)
    build_final_file1(FinalDico)
    DataDisorder = open_file("IupredResults.txt")
    Data3 = open_file("Data3")
    DicoProtoGene_NbPop = build_new_dico(Data3)
    build_final_file2(DataDisorder, DicoProtoGene_NbPop)
