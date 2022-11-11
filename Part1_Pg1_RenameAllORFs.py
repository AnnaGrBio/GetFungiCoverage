import cv2
import os
os.chdir("Folder")
import random


""" Docstring

    This script is used to rename the ORFs that were extracted with the software getORF

    input
    ------------
    - Output file from getORF in nucleotide (eg. AK5_All_ORF_Nucleotide)
    - Output file from getORF in protein (eg. AK5_All_ORF_Protein)
    
    output
    ------------
    - File of nucleotide ORF with new header including transcript name and position in transcript
    - File of protein ORF with new header including transcript name and position in transcript
    
"""

def open_file(nameFile):
    # This function open file in lecture mode
    F = open(nameFile, "r")
    L = F.readlines()
    return L  

## Step 1. Store sequences in Dictionary ##
def sequence_storage(File): 
    # This function stored the ORFs in a dictionary
    Dico = {}
    name = ""
    seq = ""
    for i in File:
        ligne = i.split("\n")[0]
        if ligne[0] == ">": # Search for the header of the sequence
            if name != "":
                Dico[name] = seq
                seq = ""
            name = ligne[1:]
            #print name
        else:
            seq += ligne
    Dico[name] = seq
    return Dico # Returns a dictionary with the name of the ORF in transcript and the ORF sequence


## Step 2. Remove reverse AK5_All_ORF_Nucleotide ##
def remove_reverse_seq(Dico):
    # This function remove all ORFs that were found in reverse direction in the transcripts
    for i in Dico.keys():
        if "REVERSE SENSE" in i:
            Dico.pop(i)
    
    
## Step 3. Rename Transcript ##
def rename_transcripts(dico, namePop):
    # This function renames the transcripts
    # The new name will include the name of the population in which it was found, the transcript name and 
    # Start and Stop in transcript
    new_dico = {}
    for i in dico.keys():
        ligne = i.split(" ")
        start = ligne[1][1:]
        end = ligne[3][0:len(ligne[3])-1]
        new_name = namePop + "_" + ligne[0] + "_" + start + "_" + end
        new_dico[new_name] = dico[i]
    return new_dico # return the dictionary with the new name associated to the sequence


## Step 4. write_output_file output file ##
def write_output_file(DicoNuc, DicoProt, Popname):
    # This function produce the finals output files, one with the ORFs and one with their corresponding proteins
    # with the correct name as a header
    F = open((Popname + "_NucleotideORF"), "w")
    for i in DicoNuc.keys():
        F.write_output_file(">" + i + "\n")
        F.write_output_file(DicoNuc[i] + "\n")
    F.close()
    F = open((Popname+"_ProteinORF"), "w")
    for i in DicoProt.keys():
        F.write_output_file(">" + i + "\n")
        F.write_output_file(DicoProt[i] + "\n")
    F.close()


def main_function():
    DNA = open_file("AK5_All_ORF_Nucleotide")
    Protein = open_file("AK5_All_ORF_Protein")
    ###### Step 1. Store sequences in Dictionary
    dico_DNA = sequence_storage(DNA) 
    dico_protein = sequence_storage(Protein)
    print ("le nombre initial d ORF est "+str(len(dico_protein)))
    ###### Step 2. Remove reverse AK5_All_ORF_Nucleotide
    remove_reverse_seq(dico_DNA)
    remove_reverse_seq(dico_protein)
    print ("le nombre final d ORF est "+str(len(dico_protein)))
    ###### Step 3. Rename ORFs
    new_dico_DNA = rename_transcripts(dico_DNA, "AK5")
    new_dico_protein = rename_transcripts(dico_protein, "AK5")
    ###### Step 4. write_output_file output file
    write_output_file(new_dico_DNA, new_dico_protein, "AK5")
    
    
## Run Code ##
main_function()




















