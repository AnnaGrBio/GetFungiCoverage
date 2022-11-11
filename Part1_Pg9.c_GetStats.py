import cv2
import os
os.chdir("path to folder")
import random


""" docstring

    This program calculates the percentage of each position of proto-genes
    
    input
    ------------
    - previous output file of the position in genomes of proto-genes
    
    output
    - percentage of proto-genes in each position
"""


def open_file(NameFile):
    F = open(NameFile, "r")
    L = F.readlines()
    return L  


def build_dic_positions(FilePositions):
    # this function associate each proto-gene to its genomic position
    Dico = {}
    for i in FilePositions:
        ligne = i.split("\n")[0]
        ligne = ligne.split(",")
        Dico[ligne[0]] = ligne[1]
    return Dico


def retrieve_position(Dico):
    # This function sum all found position and calculate the final percentage
    UndirectionalExonLonger = 0
    UndirectionalExonInside = 0
    Intergenic = 0
    ReverseGenic = 0
    Pseudogene = 0
    NcRNA = 0
    Intronic = 0
    Weird = 0
    ExonLonger = 0
    ExonInside = 0
    for i in Dico.keys():
        PosORF = Dico[i]
        if PosORF == "UndirectionalExonLonger":
            UndirectionalExonLonger += 1
        elif PosORF == "UndirectionalExonInside":
            UndirectionalExonInside += 1
        elif PosORF == "Intergenic":
            Intergenic += 1
        elif PosORF == "ReverseGenic":
            ReverseGenic += 1
        elif PosORF == "Pseudogene":
            Pseudogene += 1
        elif PosORF == "NcRNA":
            NcRNA += 1
        elif PosORF == "ExonLonger":
            ExonLonger += 1
        elif PosORF == "ExonInside":
            ExonInside += 1
        elif PosORF == "Intronic":
            Intronic += 1
        else:
            Weird += 1
    print ("UndirectionalExonLonger : "+str(UndirectionalExonLonger))
    print ("UndirectionalExonInside : "+str(UndirectionalExonInside))
    print ("ExonInside : "+str(ExonInside))
    print ("ExonLonger : "+str(ExonLonger))
    print ("Intergenic : "+str(Intergenic))
    print ("ReverseGenic : "+str(ReverseGenic))
    print ("Pseudogene : "+str(Pseudogene))
    print ("NcRNA : "+str(NcRNA))
    print ("Intronic : "+str(Intronic))
    print ("Weird : "+str(Weird))
    

def main_function():



    FilePositions = open_file("AK5_PositionsInGenome")
    Dico = build_dic_positions(FilePositions)
    retrieve_position(Dico)
