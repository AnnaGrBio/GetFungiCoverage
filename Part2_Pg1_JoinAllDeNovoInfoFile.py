import cv2
import os
os.chdir("path to folder")
import random
from collections import Counter


""" docstring

    This program builds a new info file containing only information about de novo genes for each line.

    input
    ------------
    - information file of each line

    output
    ------------
    - one single info file with information about de novo ORFs of each lines

"""


def open_file(NameFile):
    F = open(NameFile, "r")
    L = F.readlines()
    return L 
 

def extract_info_infofile(File):
    # This function extracts information of de novo genes from info file and store these lines in a list that it returns
    liste = []
    c = 0
    for i in File[1:]:
        ligne = i.split("\n")[0]
        Data = ligne.split(",")
        if Data[15] != "-":
            liste.append(i)
    return liste


def build_final_output(L1, L2, L3, L4, L5, L6, L7):
    # This function builds the final info file
    F = open("DeNovo_InfoFile", "w")
  F.write("GeneName,NbTranscripts,TranscriptName,NbExon,Chrom,StartTranscript,EndTranscript,Direction,Cov,FPKM,TPM,DeNovoORF,ORFname,StartUnsplicedORF,StopUnsplicedORF,GenomicPosition"+"\n")
    for i in L1:
        F.write(i)
    for i in L2:
        F.write(i)
    for i in L3:
        F.write(i)
    for i in L4:
        F.write(i)
    for i in L5:
        F.write(i)
    for i in L6:
        F.write(i)
    for i in L7:
        F.write(i)
    F.close()


def main_function():
    FileAK5 = open_file("FI_Final_InformationFile")
    listeDeNovoAK5 = extract_info_infofile(FileAK5)
    FileDK5 = open_file("DK_Final_InformationFile")
    listeDeNovoDK5 = extract_info_infofile(FileDK5)
    FileGI5 = open_file("ES_Final_InformationFile")
    listeDeNovoGI5 = extract_info_infofile(FileGI5)
    FileSW5 = open_file("SE_Final_InformationFile")
    listeDeNovoSW5 = extract_info_infofile(FileSW5)
    FileUM = open_file("UA_Final_InformationFile")
    listeDeNovoUM = extract_info_infofile(FileUM)
    FileYE = open_file("TR_Final_InformationFile")
    listeDeNovoYE = extract_info_infofile(FileYE)
    FileZamb = open_file("ZI_Final_InformationFile")
    listeDeNovoZamb = extract_info_infofile(FileZamb)
    build_final_output(listeDeNovoAK5, listeDeNovoDK5, listeDeNovoGI5, listeDeNovoSW5, listeDeNovoUM, listeDeNovoYE, listeDeNovoZamb)





