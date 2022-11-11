import cv2
import os
os.chdir("path to folder")
import random
from collections import Counter


""" docstring

    This program claculates aggregation propensity of each proto gene and store this information also with nb of line that share the gene

    input
    ------------
    - protein sequences of proto-genes

    output
    ------------
    - file with aggregation propensity

"""


def open_file(NameFile):
    F = open(NameFile, "r")
    L = F.readlines()
    return L  


def build_sec_dic(File):
    D = {}
    for i in File:
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            Name = ligne[1:]
        else:
            Seq = ligne
            D[Name] = ligne
    return D


def read_file_aggregation_propensity():
    Liste = []
    File = open_file("VTS1.txt")
    for i in File[1:]:
        ligne = i.split()
        Data = float(ligne[5])
        Liste.append(Data)
    AggProp = float(sum(Liste))/float(len(Liste))
    os.system('rm VTS1.txt')
    return AggProp


def run_TANGO(Dico):
    FinalDico = {}
    for i in Dico.keys():
        Name = i
        Seq = Dico[i]
        Command = "./tango_x86_64_release VTS1 ct=" + "N" + " nt="+"N" +" ph="+"7.2" +" te="+"303"+" io="+"0.02"+"seq="+Seq+">> Lala.txt"
        os.system(Command)
        AggProp = read_file_aggregation_propensity()
        FinalDico[Name] = AggProp
    return FinalDico


def build_final_file(FinalDico):
    F = open("TANGOResults.txt", "w")
    for i in FinalDico.keys():
        F.write(i+","+str(FinalDico[i])+"\n")
    F.close()


def build_new_dico(Data3):
    Translation = {"ORF_1_Pop":1, "ORF_2_Pop":2,"ORF_3_Pop":3,"ORF_4_Pop":4,"ORF_5_Pop":5,"ORF_6_Pop":6,"ORF_7_Pop":7}
    D = {}
    for i in Data3[1:]:
        ligne = i.split(",")
        NbPop = ligne[0]#str(Translation[ligne[0]])
        NameProtoGene = ligne[1]
        D[NameProtoGene] = NbPop
    return D


def make_formated_file(DataTango, Dico):
    F = open("FinalDataAggregation.txt", "w")
    F.write("NbPops,Aggregation"+"\n")
    for i in DataTango:
        ligne = i.split("\n")[0]
        Name = ligne.split(",")[0]
        Aggreg = ligne.split(",")[1]
        if Name in Dico.keys():
            F.write(Dico[Name]+","+Aggreg+"\n")
    F.close()
        
def main_function():
    File = open_file("AllPopsProtoGenesProt.fa") 
    DicoSeqs = build_sec_dic(File)
    FinalDico = run_TANGO(DicoSeqs)
    build_final_file(FinalDico)
    DataTango = openFile("TANGOResults.txt")
    Data3 = openFile("Data3")
    DicoProtoGene_NbPop = build_new_dico(Data3)
    make_formated_file(DataTango, DicoProtoGene_NbPop)



