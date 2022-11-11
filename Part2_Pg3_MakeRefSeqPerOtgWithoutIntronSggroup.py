import cv2
import os
os.chdir("path to folder")
import random
from collections import Counter


""" docstring

    This program select a referent sequences for orthogroups in which the ORF is not present in all lines and contain NO intron
    It onlu search orthogroups where one sequence per line is present

    input
    ------------
    - orthogroup file from previous step

    output
    ------------
    - file with referent seq per orthogropus + some info about the orthogroup/sequence

"""

def openFile(NameFile):
    F = open(NameFile, "r")
    L = F.readlines()
    return L  


def search_missing_lines(ListePresentes):
    # this function determines which lines are present and absents from the orthogroup
    list_all_lines = ["FI", "DK", "ES", "SE", "UA", "TR", "ZI"]
    list_missing = []
    for pop in list_all_lines:
        if pop not in ListePresentes:
            list_missing.append(pop)
    return list_missing


def build_final_file(File):
    # This function creates the final file
    FinalFile = open("InfoFile_DeNovo_Enabling_WithoutIntrons.txt", "w")
    FinalFile.write("Orthogroup,NameDeNovo,NamePop,Exons,Chrom,Start,End,PopsToSearchIn"+"\n")
    PopsPresentes = []
    RefSeq = ""
    for i in File:
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            if len(PopsPresentes)>0:
                if len(RefSeq)>0 and len(PopsPresentes)<7:
                    list_missingPops = search_missing_lines(PopsPresentes)
                    FinalFile.write(RefSeq)
                    for j in list_missingPops:
                        FinalFile.write(","+j)
                    FinalFile.write("\n")
            PopsPresentes = []
            RefSeq = ""
            Orthogroup = ligne[1:len(ligne)-1]
        else:
            Data = ligne.split(",")
            NameDeNovo = Data[0]
            NamePop = NameDeNovo.split("_")[0]
            if NamePop not in PopsPresentes:
                PopsPresentes.append(NamePop)
            Chrom = Data[2]
            Exons = Data[1]
            Start = Data[5]
            End = Data[6]
            if RefSeq == "" and Exons == "1":
                RefSeq = Orthogroup+","+NameDeNovo+","+NamePop+","+Exons+","+Chrom+","+Start+","+End
    if len(RefSeq)>0 and len(PopsPresentes)<7:
        list_missingPops = search_missing_lines(PopsPresentes)
        FinalFile.write(RefSeq)
        for j in list_missingPops:
            FinalFile.write(","+j)
    FinalFile.close()
        

def main_function():
    Orthogroups = openFile("Orthogroups_info.txt")
    build_final_file(Orthogroups)






