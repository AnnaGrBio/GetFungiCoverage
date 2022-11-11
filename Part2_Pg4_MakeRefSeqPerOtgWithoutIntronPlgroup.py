import cv2
import os
os.chdir("path to folder")
import random
from collections import Counter


""" docstring

    This program select a referent sequences for orthogroups in which the ORF is not present in all lines and contain NO intron
    It onlu search orthogroups where several sequences per line are present

    input
    ------------
    - orthogroup file from previous step

    output
    ------------
    - file with referent seq per orthogropus + some info about the orthogroup/sequence

"""

def openFile(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  


def build_dic_sequences(DataAK5, DataDK5, DataGI5, DataSW5, DataUM, DataYE, DataZamb):
    # This function retrieve one referent sequence per orthogroup
    Dico = {}
    for i in DataAK5:
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            Name = ligne[1:]
        else:
            Seq = ligne
            Dico[Name] = Seq
    for i in DataDK5:
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            Name = ligne[1:]
        else:
            Seq = ligne
            Dico[Name] = Seq
    for i in DataGI5:
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            Name = ligne[1:]
        else:
            Seq = ligne
            Dico[Name] = Seq
    for i in DataSW5:
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            Name = ligne[1:]
        else:
            Seq = ligne
            Dico[Name] = Seq
    for i in DataUM:
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            Name = ligne[1:]
        else:
            Seq = ligne
            Dico[Name] = Seq
    for i in DataYE:
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            Name = ligne[1:]
        else:
            Seq = ligne
            Dico[Name] = Seq
    for i in DataZamb:
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            Name = ligne[1:]
        else:
            Seq = ligne
            Dico[Name] = Seq
    return Dico


def build_final_file(DataInfo, Dico):
    F = open("OrthogroupsRefSeqWithoutIntron.fa", "w")
    for i in DataInfo[1:]:
        Name = i.split(",")[1]
        Seq = Dico[Name]
        F.write(">"+Name+"\n")
        F.write(Seq+"\n")
    F.close()


def main_function():
    DataAK5 = openFile("FI_DeNovoORFwithStop_DNA.fa")
    DataDK5 = openFile("DK_DeNovoORFwithStop_DNA.fa")
    DataGI5 = openFile("ES_DeNovoORFwithStop_DNA.fa")
    DataSW5 = openFile("SE_DeNovoORFwithStop_DNA.fa")
    DataUM = openFile("UA_DeNovoORFwithStop_DNA.fa")
    DataYE = openFile("TR_DeNovoORFwithStop_DNA.fa")
    DataZamb = openFile("ZI_DeNovoORFwithStop_DNA.fa")
    DicoAllSeqs = build_dic_sequences(DataAK5, DataDK5, DataGI5, DataSW5, DataUM, DataYE, DataZamb)
    DataInfo = openFile("InfoFile_DeNovo_Enabling_WithoutIntrons.txt")
    build_final_file(DataInfo, DicoAllSeqs)











