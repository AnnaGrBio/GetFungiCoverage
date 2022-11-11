import cv2
import os
os.chdir("path to folder")
import random
from collections import Counter


""" docstring

    This program calculates the presence and overlap to TEs of de novo ORFs.
    To process, it uses the TE annotation and lowering provided done in previous steps

    input
    ------------    
    - de novo ORF sequences

    output
    ------------
    - NUAber of de novo sequences that belong to one of the 3 categories:
        - inside TE
        - overlap to TE
        - out of TEs 

"""


def open_file(NameFile):
    F = open(NameFile, "r")
    L = F.readlines()
    return L  


def store_TE_info_sequences(FileSeq):
    # This function estimate the full presence or overlap of a TE to de novo genes, and store the results in a dic
    TE_list = ["a", "t", "g", "c"]
    not_TE_list = ["A", "T", "G", "C"]
    DicoFinal = {}
    c = 0
    for i in FileSeq:
        ligne = i.split("\n")[0]
        if i[0] == ">":
            Name = ligne[1:]
        else:
            TEpresent = False
            RegularPresent = False
            for i in TE_list:
                if i in ligne:
                    TEpresent = True
                    break
            for i in not_TE_list:
                if i in ligne:
                    RegularPresent = True
            if TEpresent == True and RegularPresent == False:
                DicoFinal[Name] = "InsideTE"
            elif TEpresent == False and RegularPresent == True:
                DicoFinal[Name] = "OutsideTE"
            elif TEpresent == True and RegularPresent == True:
                DicoFinal[Name] = "OverlapTE"
            else:
                print (ligne)
    return DicoFinal
            

def count_distribution(DataORF,DicoFI, DicoDK, DicoES,DicoSE, DicoUA, DicoTR,DicoZI):
    # For each lines, this function counts the distribution of each categories of ORFs to TE overlap
    FI_InsideTE = 0
    FI_OutsideTE = 0
    FI_OverlapTE = 0
    DK_InsideTE = 0
    DK_OutsideTE = 0
    DK_OverlapTE = 0
    ES_InsideTE = 0
    ES_OutsideTE = 0
    ES_OverlapTE = 0
    SE_InsideTE = 0
    SE_OutsideTE = 0
    SE_OverlapTE = 0
    UA_InsideTE = 0
    UA_OutsideTE = 0
    UA_OverlapTE = 0
    TR_InsideTE = 0
    TR_OutsideTE = 0
    TR_OverlapTE = 0
    ZI_InsideTE = 0
    ZI_OutsideTE = 0
    ZI_OverlapTE = 0
    c = 0
    for i in DataORF[1:]:
        ORFName = i.split(",")[1]
        if ORFName in DicoFI.keys():
            Value = DicoFI[ORFName]
            if Value == "InsideTE":
                FI_InsideTE += 1
            elif Value == "OutsideTE":
                FI_OutsideTE += 1
            else:
                FI_OverlapTE += 1
        if ORFName in DicoDK.keys():
            Value = DicoDK[ORFName]
            if Value == "InsideTE":
                DK_InsideTE += 1
            elif Value == "OutsideTE":
                DK_OutsideTE += 1
            else:
                DK_OverlapTE += 1
        if ORFName in DicoES.keys():
            Value = DicoES[ORFName]
            if Value == "InsideTE":
                ES_InsideTE += 1
            elif Value == "OutsideTE":
                ES_OutsideTE += 1
            else:
                ES_OverlapTE += 1
        if ORFName in DicoSE.keys():
            Value = DicoSE[ORFName]
            if Value == "InsideTE":
                SE_InsideTE += 1
            elif Value == "OutsideTE":
                SE_OutsideTE += 1
            else:
                SE_OverlapTE += 1
        if ORFName in DicoUA.keys():
            Value = DicoUA[ORFName]
            if Value == "InsideTE":
                UA_InsideTE += 1
            elif Value == "OutsideTE":
                UA_OutsideTE += 1
            else:
                UA_OverlapTE += 1
        if ORFName in DicoTR.keys():
            Value = DicoTR[ORFName]
            if Value == "InsideTE":
                TR_InsideTE += 1
            elif Value == "OutsideTE":
                TR_OutsideTE += 1
            else:
                TR_OverlapTE += 1
        if ORFName in DicoZI.keys():
            Value = DicoZI[ORFName]
            if Value == "InsideTE":
                ZI_InsideTE += 1
            elif Value == "OutsideTE":
                ZI_OutsideTE += 1
            else:
                ZI_OverlapTE += 1
        else:
            print ORFName
        c += 1   
    F = open("TEinORF.txt", "w")
    F.write("Population,Status,NUAber" + "\n")
    F.write("FI,InsideTE," + str(FI_InsideTE) + "\n")
    F.write("FI,OutsideTE," + str(FI_OutsideTE) + "\n")
    F.write("FI,OverlapTE," + str(FI_OverlapTE) + "\n")
    F.write("DK,InsideTE," + str(DK_InsideTE) + "\n")
    F.write("DK,OutsideTE," + str(DK_OutsideTE) + "\n")
    F.write("DK,OverlapTE," + str(DK_OverlapTE) + "\n")
    F.write("ES,InsideTE," + str(ES_InsideTE) + "\n")
    F.write("ES,OutsideTE," + str(ES_OutsideTE) + "\n")
    F.write("ES,OverlapTE," + str(ES_OverlapTE) + "\n")
    F.write("SE,InsideTE," + str(SE_InsideTE) + "\n")
    F.write("SE,OutsideTE," + str(SE_OutsideTE) + "\n")
    F.write("SE,OverlapTE," + str(SE_OverlapTE) + "\n")
    F.write("UA,InsideTE," + str(UA_InsideTE) + "\n")
    F.write("UA,OutsideTE," + str(UA_OutsideTE) + "\n")
    F.write("UA,OverlapTE," + str(UA_OverlapTE) + "\n")
    F.write("TR,InsideTE," + str(TR_InsideTE) + "\n")
    F.write("TR,OutsideTE," + str(TR_OutsideTE) + "\n")
    F.write("TR,OverlapTE," + str(TR_OverlapTE) + "\n")
    F.write("ZI,InsideTE," + str(ZI_InsideTE) + "\n")
    F.write("ZI,OutsideTE," + str(ZI_OutsideTE) + "\n")
    F.write("ZI,OverlapTE," + str(ZI_OverlapTE) + "\n")
    F.close()
    

def main_function():   
    FileAllSeqencesFI = open_file("FI_NucleotideORF")
    FileAllSeqencesDK = open_file("DK_NucleotideORF")
    FileAllSeqencesES = open_file("ES_NucleotideORF")
    FileAllSeqencesSE = open_file("SE_NucleotideORF")
    FileAllSeqencesUA = open_file("UA_NucleotideORF")
    FileAllSeqencesTR = open_file("TR_NucleotideORF")
    FileAllSeqencesZI = open_file("ZI_NucleotideORF")
    DataORF = open_file("Data3")
    DicoFI = store_TE_info_sequences(FileAllSeqencesFI)
    DicoDK = store_TE_info_sequences(FileAllSeqencesDK)
    DicoES = store_TE_info_sequences(FileAllSeqencesES)
    DicoSE = store_TE_info_sequences(FileAllSeqencesSE)
    DicoUA = store_TE_info_sequences(FileAllSeqencesUA)
    DicoTR = store_TE_info_sequences(FileAllSeqencesTR)
    DicoZI = store_TE_info_sequences(FileAllSeqencesZI)
    print ("Dico Seq Done")
    count_distribution(DataORF,DicoFI, DicoDK, DicoES,DicoSE, DicoUA, DicoTR,DicoZI)  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
