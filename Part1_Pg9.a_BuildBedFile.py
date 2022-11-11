import cv2
import os
os.chdir("path to folder")
import random


""" docstring

    This program builds a bed file for the de novo ORFs
    
    input
    ------------
    - info file
    
    output
    ------------
    - bed file of de novo ORFs
"""

def open_file(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  


def get_transcript_list(File):
    # This function retreive the name of each transcripts from the previous info file
    L = []
    for i in File:
        if i[0] == ">":
            Name = i[1:len(i)-1]
            L.append(Name)
    return L


def build_bed_file(Info, NameFinal):
    # This function builds the final bed file for the transcripts and ORFs from my info file
    F = open(NameFinal, "w")
    for i in Info[1:]:
        ligne = i.split("\n")[0]
        ligne = ligne.split(",")
        transcript_name = ligne[12]
        Chrom = ligne[4]
        Beg = ligne[13]
        End = ligne[14]
        TPM = float(ligne[10])
        Direction = ligne[7]
        if Direction == ".":
            Direction = "-"
        if Beg!= "-" and End!= "-":
            if int(Beg)>int(End):
                M = End
                End = Beg
                Beg = M
            if TPM>=0.5:
                F.write(Chrom+"	"+Beg+"	"+End+"	"+Direction+"	"+transcript_name+"\n")
    F.close()
    

def main_function():
    Info = open_file("FI_InformationFile")
    build_bed_file(Info, "FIdenovo_bed")

