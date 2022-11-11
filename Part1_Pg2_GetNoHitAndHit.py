import cv2
import os
from Bio.sequence import sequence
os.chdir("path to folder")
import random

""" docstring

    This program assess the file resulting from a multiple BLAST search, 
    according to the p-value that you set, it will retrieve the sequences
    that got an hit and the sequences that did not, and print them in
    output files
    
    input
    ------------
    - File with BLAST results
    
    output
    ------------
    - protein and/or DNA with and without hit
"""

def open_file(nameFile):
    # This function opens a file in a read mode
    F = open(nameFile, "r")
    L = F.readlines()
    return L  


def detect_no_hit(F, Evalue):
    # This function search the E value in a BLAST output file
    nb_total_query = 0
    final_list_no_hit = [] # Store name of ORFs without hits
    final_list_hit = [] # Store name of ORFs with hits
    for i in range(0, len(F)):
        if F[i][0:6] == "Query=":
            nb_total_query += 1
            #print nb_total_query
            if F[i+5] == "***** No hits found *****"+"\n":
                Nom = F[i][7:]
                final_list_no_hit.append(Nom)
            else:
                lala = F[i+6]
                lili = lala.split()
                str_value = lili[len(lili)-1]
                FValue = float(str_value)
                if float(FValue)>float(Evalue):
                    Nom =  F[i][7:]
                    final_list_no_hit.append(Nom)
                else:
                    Nom = F[i][7:]
                    final_list_hit.append(Nom)
    print ("Nb initial ORF : " + str(nb_total_query))
    print ("Nb Total de novo : " + str(len(final_list_no_hit)))
    print ("Nb Total Hits : " + str(len(final_list_hit)))
    return final_list_no_hit, final_list_hit # Returns the final list with hits and no hits


def final_file_no_hit(nameFile, Liste, protein_file):
    # This function retrieves the protein sequences of the ORF which have no hit and store them in a file
    DicoProt = {}
    name = ""
    sequence = ""
    for i in protein_file:
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            if name != "":
                DicoProt[name] = sequence
                sequence = ""
            name = ligne[1:]
        else:
            sequence += ligne
    DicoProt[name] = sequence
    F = open(nameFile, "w")
    for i in Liste:
        Nom = i.split("\n")[0]
        F.write(">" + Nom + "\n")
        F.write(DicoProt[Nom] + "\n")
    F.close()
    

def final_file_hit(nameFile, Liste, protein_file):
    # This function retrieves the protein sequences of the ORF which have at least one hit and store them in a file
    DicoProt = {}
    name = ""
    sequence = ""
    for i in protein_file:
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            if name != "":
                DicoProt[name] = sequence
                sequence = ""
            name = ligne[1:]
        else:
            sequence += ligne
    DicoProt[name] = sequence
    F = open(nameFile, "w")
    for i in Liste:
        Nom = i.split("\n")[0]
        F.write(">" + Nom + "\n")
        F.write(DicoProt[Nom] + "\n")
    F.close()


def main_function():
    F = open_file("FIblastOutput.txt")
    print ("******************")
    print ("FileOpened")
    print ("******************")
    E_VALUE = 0.01
    list_no_hits, list_hits = detect_no_hit(F, E_VALUE)
    protein_file = open_file("FI_PutativeDeNovoORF_Protein.fa")
    nucleotid_file = open_file("FI_PutativeDeNovoORF_DNA.fa")
    final_file_no_hit("FI_DeNovoORF_Protein.fa", list_no_hits,protein_file)
    final_file_no_hit("FI_DeNovoORF_DNA.fa", list_no_hits,nucleotid_file)
    final_file_hit("FI_ProtHits.fa", list_hits, protein_file)
    
    
    
    
