import cv2
import os
os.chdir("/global/students/homes/agrandch/DeNovoGenes/NewVersion/MakeOrthogroups/Step8MakeFileEnablingMutationWITHintrons")
import random
from collections import Counter


""" docstring

    This program select a referent sequences for orthogroups in which the ORF is not present in all lines and contain at least one intron
    It onlu search orthogroups where several sequences per line are present

    input
    ------------
    - orthogroup file from previous step

    output
    ------------
    - file with referent seq per orthogropus + some info about the orthogroup/sequence

"""

def open_file(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  


def build_dic_seq(DataAK5, DataDK5, DataGI5, DataSW5, DataUM, DataYE, DataZamb):
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
    F = open("OrthogroupsRefSeqWithIntron.fa", "w")
    for i in DataInfo[1:]:
        Name = i.split(",")[1]
        Seq = Dico[Name]
        F.write(">"+Name+"\n")
        F.write(Seq+"\n")
    F.close()


def main_function():
    DataAK5 = open_file("FI_DeNovoGene_Intron_Exon.fa")
    DataDK5 = open_file("DK_DeNovoGene_Intron_Exon.fa")
    DataGI5 = open_file("ES_DeNovoGene_Intron_Exon.fa")
    DataSW5 = open_file("SE_DeNovoGene_Intron_Exon.fa")
    DataUM = open_file("UA_DeNovoGene_Intron_Exon.fa")
    DataYE = open_file("TR_DeNovoGene_Intron_Exon.fa")
    DataZamb = open_file("ZI_DeNovoGene_Intron_Exon.fa")
    DicoAllSeqs = build_dic_seq(DataAK5, DataDK5, DataGI5, DataSW5, DataUM, DataYE, DataZamb)
    DataInfo = open_file("InfoFile_DeNovo_Enabling_WithIntrons.txt")
    build_final_file(DataInfo, DicoAllSeqs)











