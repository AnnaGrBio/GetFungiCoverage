import cv2
import os
os.chdir("path to folder")
import random


""" docstring

    This program builds a new file containing all de novo ORFs, 
    but add their respective stop codon to the end of the sequence, 
    as it is not implemented in the file resulting from getORF
    
    input
    ------------
    - position of unspliced ORFs
    - de novo ORF file
    
    output
    ------------
    - new file with de novo orf plus their stop codon
"""


def open_file(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  


def build_dic_stop_codon(File):
    # This function builds the dictionary will all found stop codons
    D = {}
    for i in File[1:]:
        ligne = i.split("\n")[0]
        Content = ligne.split(",")
        Name = Content[0]
        Stop = Content[3]
        D[Name] = Stop
    return D        
        
        
def build_dic_sequences(File):
    # This function builds a dictionary with all sequences
    Dico = {}
    Seq = ""
    Name = ""
    for i in File:
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            if Name != "":
                Dico[Name] = Seq
                Seq = ""
            Name = ligne[1:]
        else:
            Seq += ligne
    Dico[Name] = Seq
    return Dico


def build_final_file(Name, DicoStop, DicoSeqs):
    # This function builds the final dictionary with sequences and stop codons
    F = open(Name, "w")
    for i in DicoSeqs.keys():
        Seq = DicoSeqs[i].upper()
        if i in DicoStop.keys():
            Stop = DicoStop[i]
            Seq = Seq + Stop
            F.write(">" + i + "\n")
            F.write(Seq + "\n")
    F.close()
                

def main_function():
    Stops = open_file("FI_UnsplicedORFpositions")
    DicoStop = build_dic_stop_codon(Stops)
    SequencesNuc = open_file("FI_DeNovoORF_DNA.fa")
    DicoSeqsNuc = build_dic_sequences(SequencesNuc)
    SequencesProt = open_file("FI_DeNovoORF_Protein.fa")
    DicoSeqsProt = build_dic_sequences(SequencesProt)
    build_final_file("FI_DeNovoORFwithStop_DNA.fa", DicoStop, DicoSeqsNuc)
    build_final_file("FI_DeNovoORFwithStop_Protein.fa", DicoStop, DicoSeqsProt)
