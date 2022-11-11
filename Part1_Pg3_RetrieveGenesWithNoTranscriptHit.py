import cv2
import os
from Bio.Seq import Seq
os.chdir("path to directory")
import random


""" docstring

    This program aims to remove all ORF with no hit which belong to a transcript were hits were found
    
    input
    ------------
    - a gtf file of the transcriptome assembly
    - the file from pg2 output with no hit ORFs
    
    output
    ------------
    - a corrected file with ORFs without hits
"""

def open_file(NameFile):
    F = open(NameFile, "r")
    L = F.readlines()
    return L  


####### Step1 Retrieve all transcripts
def get_all_transcripts(File):
    # This function retrieves all transcripts name from a gtf file and store them in a list
    liste_transcripts = []
    for i in File[2:]:
        detection = i.split("	")[2] # gtf parsing
        if detection == "transcript":
            ligne = i.split(";")
            first_name = ligne[1].split(" ")[2]
            transcript_name = first_name[1:len(first_name)-1]
            liste_transcripts.append(transcript_name)
    return liste_transcripts # return list of transcripts

        
####### Step2 Remove hits from the list
def remove_hit(liste_transcripts, Hits):
    # this function remove hits from the list of hits when one of them belong to a transcript that shows hits
    for i in Hits:
        if i[0] == ">":
            transcript_name = i.split("_")[1]
            if transcript_name in liste_transcripts:
                liste_transcripts.remove(transcript_name)
    return liste_transcripts 


####### Step3 Make Dico of DNA
def build_dic_prot(no_hit, liste_transcripts_without_hits):
    # This function retrives the longest ORF in a transcript that has several de novo ORFs
    D = {}
    list_de_novo_transcript = []
    for i in no_hit:
        if i[0] == ">":
            transcript_name = i.split("_")[1]
            if transcript_name in liste_transcripts_without_hits:
                Name = i
                if transcript_name not in list_de_novo_transcript:
                    list_de_novo_transcript.append(transcript_name)
            else:
                Name = ""
        else:
            if Name != "":
                Seq = i
                D[Name] = Seq
                Name = ""
                Seq = ""
    print ("The total number of De novo transcript: "+str(len(list_de_novo_transcript)))
    return D


####### Step4 Make Dico Proteins
def build_dic_nuc(File, DicoProt):
    # This function build a dictionary that contain purified list of no hit proteins
    DicoNuc = {}
    for i in File:
        if i[0] == ">":
            if i in DicoProt.keys():
                Name = i
            else:
                Name = ""
        else:
            if Name != "":
                Seq = i
                #print Seq
                DicoNuc[Name] = Seq
                Seq = ""
                Name = ""
    return DicoNuc


## Step5 Make final files
def make_final_output(Name, Dico):
    # this function creates the output file
    F = open(Name, "w")
    for i in Dico.keys():
        F.write(i)
        F.write(Dico[i])
    F.close()


def main_function():
    GTF = open_file("FItranscriptome.gtf")
    ListeAllTranscript = get_all_transcripts(GTF)
    print ("The total number of transcripts is: " + str(len(ListeAllTranscript)))
    Hits = open_file("FI_Hits.fa")
    liste_transcripts_without_hits = remove_hit(ListeAllTranscript, Hits)
    print ("The total number of transcripts without hit: " + str(len(liste_transcripts_without_hits)))
    no_hit = open_file("FI_PutativeDeNovoRound1.fa")
    DicoORFsProt = build_dic_prot(no_hit, liste_transcripts_without_hits)
    AllNucleotides = open_file("FI_NucleotideORF")
    DicoORFsNucleotide = build_dic_nuc(AllNucleotides, DicoORFsProt)
    make_final_output("FI_SortedORF_WithoutHit_DNA", DicoORFsNucleotide)
    make_final_output("FI_SortedORF_WithoutHit_Protein", DicoORFsProt)

