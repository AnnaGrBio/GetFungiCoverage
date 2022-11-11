import cv2
import os
os.chdir("path to folder")
import random
from collections import Counter


""" docstring

This program calculates movement between protoo-genes thathave been found in several positions in on single line

    input
    ------------
    - Orthogroup info file

    output
    ------------
    - matrix of proto-genes moves suitable for circos
    
"""


def open_file(NameFile):
    F = open(NameFile, "r")
    L = F.readlines()
    return L  


def assess_proto_genes_overlap(ListePos, PosStart, PosEnd):
    #
    Overlap = False
    for i in ListePos:
        QueryStart = i[0]
        QueryEnd = i[1]
        if PosStart<=QueryStart and PosEnd>=QueryEnd:
            Overlap = True
            break
        elif PosStart>=QueryStart and PosEnd<=QueryEnd:
            Overlap = True
            break
        elif PosStart<=QueryStart and PosEnd>QueryStart:
            Overlap = True
            break
        elif PosStart<QueryEnd and PosEnd>QueryEnd:
            Overlap = True
            break
    return Overlap


def prepare_matrix_output(Matrix, C_2L, C_2R, C_3L, C_3R, C_4, C_X, C_Y, C_Mito):
    Total = C_2L+C_2R+C_3L+C_3R+C_4+C_X+C_Y+C_Mito
    if Total <=1:
        return Matrix
    else:
        if C_2L>0:
            if C_2L>1:
                Matrix[2][2]+=1
            if C_2R>0:
                Matrix[2][3]+=1
            if C_3L>0:
                Matrix[2][4]+=1
            if C_3R>0:
                Matrix[2][5]+=1
            if C_4>0:
                Matrix[2][6]+=1
            if C_X>0:
                Matrix[2][7]+=1
            if C_Y>0:
                Matrix[2][8]+=1
            if C_Mito>0:
                Matrix[2][9]+=1
        if C_2R>0:
            if C_2L>0:
                Matrix[3][2]+=1
            if C_2R>1:
                Matrix[3][3]+=1
            if C_3L>0:
                Matrix[3][4]+=1
            if C_3R>0:
                Matrix[3][5]+=1
            if C_4>0:
                Matrix[3][6]+=1
            if C_X>0:
                Matrix[3][7]+=1
            if C_Y>0:
                Matrix[3][8]+=1
            if C_Mito>0:
                Matrix[3][9]+=1
        if C_3L>0:
            if C_2L>0:
                Matrix[4][2]+=1
            if C_2R>0:
                Matrix[4][3]+=1
            if C_3L>1:
                Matrix[4][4]+=1
            if C_3R>0:
                Matrix[4][5]+=1
            if C_4>0:
                Matrix[4][6]+=1
            if C_X>0:
                Matrix[4][7]+=1
            if C_Y>0:
                Matrix[4][8]+=1
            if C_Mito>0:
                Matrix[4][9]+=1
        if C_2R>0:
            if C_2L>0:
                Matrix[5][2]+=1
            if C_2R>0:
                Matrix[5][3]+=1
            if C_3L>0:
                Matrix[5][4]+=1
            if C_3R>1:
                Matrix[5][5]+=1
            if C_4>0:
                Matrix[5][6]+=1
            if C_X>0:
                Matrix[5][7]+=1
            if C_Y>0:
                Matrix[5][8]+=1
            if C_Mito>0:
                Matrix[5][9]+=1
        if C_4>0:
            if C_2L>0:
                Matrix[6][2]+=1
            if C_2R>0:
                Matrix[6][3]+=1
            if C_3L>0:
                Matrix[6][4]+=1
            if C_3R>0:
                Matrix[6][5]+=1
            if C_4>1:
                Matrix[6][6]+=1
            if C_X>0:
                Matrix[6][7]+=1
            if C_Y>0:
                Matrix[6][8]+=1
            if C_Mito>0:
                Matrix[6][9]+=1
        if C_X>0:
            if C_2L>0:
                Matrix[7][2]+=1
            if C_2R>0:
                Matrix[7][3]+=1
            if C_3L>0:
                Matrix[7][4]+=1
            if C_3R>0:
                Matrix[7][5]+=1
            if C_4>0:
                Matrix[7][6]+=1
            if C_X>1:
                Matrix[7][7]+=1
            if C_Y>0:
                Matrix[7][8]+=1
            if C_Mito>0:
                Matrix[7][9]+=1
        if C_Y>0:
            if C_2L>0:
                Matrix[8][2]+=1
            if C_2R>0:
                Matrix[8][3]+=1
            if C_3L>0:
                Matrix[8][4]+=1
            if C_3R>0:
                Matrix[8][5]+=1
            if C_4>0:
                Matrix[8][6]+=1
            if C_X>0:
                Matrix[8][7]+=1
            if C_Y>1:
                Matrix[8][8]+=1
            if C_Mito>0:
                Matrix[8][9]+=1
        if C_Mito>0:
            if C_2L>0:
                Matrix[9][2]+=1
            if C_2R>0:
                Matrix[9][3]+=1
            if C_3L>0:
                Matrix[9][4]+=1
            if C_3R>0:
                Matrix[9][5]+=1
            if C_4>0:
                Matrix[9][6]+=1
            if C_X>0:
                Matrix[9][7]+=1
            if C_Y>0:
                Matrix[9][8]+=1
            if C_Mito>1:
                Matrix[9][9]+=1
        return Matrix


def build_final_table(SearchedPopName,File, FileName):
    Matrix = [["Data","Data","199,21,133","255,182,193","210,105,30","205,133,63","244,164,96","222,184,135","188,143,143","255,218,185"],
              ["Data","Data","2L","2R","3L","3R","4","X","Y","M"],
              ["199,21,133","2L",0,0,0,0,0,0,0,0],
              ["255,182,193","2R",0,0,0,0,0,0,0,0],
              ["210,105,30","3L",0,0,0,0,0,0,0,0],
              ["205,133,63","3R",0,0,0,0,0,0,0,0],
              ["244,164,96","4",0,0,0,0,0,0,0,0],
              ["222,184,135","X",0,0,0,0,0,0,0,0],
              ["188,143,143","Y",0,0,0,0,0,0,0,0],
              ["255,218,185","M",0,0,0,0,0,0,0,0]]
    C_2L = 0          
    Pos2L = []
    C_2R = 0
    Pos2R = []
    C_3L = 0
    Pos3L = []
    C_3R = 0
    Pos3R = []
    C_4 = 0
    Pos4 = []
    C_X = 0
    PosX = []
    C_Y = 0
    PosY = []
    C_Mito = 0
    PosMito = []
    Start = False
    for i in File[1:]:
        if i[0] == ">":
            if Start == True:
                Matrix = prepare_matrix_output(Matrix, C_2L, C_2R, C_3L, C_3R, C_4, C_X, C_Y, C_Mito)
                C_2L = 0
                Pos2L = []
                C_2R = 0
                Pos2R = []
                C_3L = 0
                Pos3L = []
                C_3R = 0
                Pos3R = []
                C_4 = 0
                Pos4 = []
                C_X = 0
                PosX = []
                C_Y = 0
                PosY = []
                C_Mito = 0
                PosMito = []
            else:
                Start = True
        else:
            ligne = i.split(",")
            TotalName = ligne[0]
            PopName = TotalName.split("_")[0]
            if PopName == SearchedPopName:
                Chrom = ligne[2]
                PosStart = int(ligne[5])
                PosEnd = int(ligne[6])
                if Chrom == "2L_Chromosome":
                    if len(Pos2L) == 0:
                        Pos2L.append([PosStart, PosEnd])
                        C_2L +=1
                    else:
                        Overlap = assess_proto_genes_overlap(Pos2L, PosStart, PosEnd)
                        if Overlap == False:
                            Pos2L.append([PosStart, PosEnd])
                            C_2L +=1

                if Chrom == "2R_Chromosome":
                    if len(Pos2R) == 0:
                        Pos2R.append([PosStart, PosEnd])
                        C_2R +=1
                    else:
                        Overlap = assess_proto_genes_overlap(Pos2R, PosStart, PosEnd)
                        if Overlap == False:
                            Pos2R.append([PosStart, PosEnd])
                            C_2R +=1
                            
                if Chrom == "3L_Chromosome":
                    if len(Pos3L) == 0:
                        Pos3L.append([PosStart, PosEnd])
                        C_3L +=1
                    else:
                        Overlap = assess_proto_genes_overlap(Pos3L, PosStart, PosEnd)
                        if Overlap == False:
                            Pos3L.append([PosStart, PosEnd])
                            C_3L +=1
                            
                if Chrom == "3R_Chromosome":
                    if len(Pos3R) == 0:
                        Pos3R.append([PosStart, PosEnd])
                        C_3R +=1
                    else:
                        Overlap = assess_proto_genes_overlap(Pos3R, PosStart, PosEnd)
                        if Overlap == False:
                            Pos3R.append([PosStart, PosEnd])
                            C_3R +=1
                            
                if Chrom == "4_Chromosome":
                    if len(Pos4) == 0:
                        Pos4.append([PosStart, PosEnd])
                        C_4 +=1
                    else:
                        Overlap = assess_proto_genes_overlap(Pos4, PosStart, PosEnd)
                        if Overlap == False:
                            Pos4.append([PosStart, PosEnd])
                            C_4 +=1
                            
                if Chrom == "X_Chromosome":
                    if len(PosX) == 0:
                        PosX.append([PosStart, PosEnd])
                        C_X +=1
                    else:
                        Overlap = assess_proto_genes_overlap(PosX, PosStart, PosEnd)
                        if Overlap == False:
                            PosX.append([PosStart, PosEnd])
                            C_X +=1
                            
                if Chrom == "Y_Chromosome":
                    if len(PosY) == 0:
                        PosY.append([PosStart, PosEnd])
                        C_Y +=1
                    else:
                        Overlap = assess_proto_genes_overlap(PosY, PosStart, PosEnd)
                        if Overlap == False:
                            PosY.append([PosStart, PosEnd])
                            C_Y +=1
                            
                if Chrom == "mitochondrion_Chromosome":
                    if len(PosMito) == 0:
                        PosMito.append([PosStart, PosEnd])
                        C_Mito +=1
                    else:
                        Overlap = assess_proto_genes_overlap(PosMito, PosStart, PosEnd)
                        if Overlap == False:
                            PosMito.append([PosStart, PosEnd])
                            C_Mito +=1           
    Matrix = prepare_matrix_output(Matrix, C_2L, C_2R, C_3L, C_3R, C_4, C_X, C_Y, C_Mito)
    F = open(FileName, "w")
    for i in Matrix:
        for j in i:
            F.write(str(j)+" ")
        F.write("\n")
    F.close()
       

def main_function():                     
    Orthogroups = open_file("Orthogroups_info.txt")
    build_final_table("FI",Orthogroups, "AK5_circlePlot")
    build_final_table("DK",Orthogroups, "DK5_circlePlot")
    build_final_table("ES",Orthogroups, "GI5_circlePlot")
    build_final_table("SE",Orthogroups, "SW5_circlePlot")
    build_final_table("UMUA",Orthogroups, "UM_circlePlot")
    build_final_table("TR",Orthogroups, "YE_circlePlot")
    build_final_table("ZI",Orthogroups, "Zamb_circlePlot")










