import cv2
import os
os.chdir("path to folder")
import random


""" docstring

    This program retrieve the longest ORF per transcript whose none ORF has a hit.
    It also assess the presence of a stop codon at the end of each ORFs, by retrieving it from the transcript.
    ORFs without STOP codon, or ORFs that end by the only fact that the transcript ends, are not taken into acount.
    
    input
    ------------
    - transcriptome assembly
    - ORF sorted from Pg3
    
    output
    - ORF sorted in such a way only one ORF is associated to a de novo transcript
    
"""


def open_file(nameFile):
    F=open(nameFile, "r")
    L=F.readlines()
    return L  


### Preliminary work. Get Transcript Size
def get_transcrip_size(File):
    # This function calculate transcripts sizes
    D = {}
    name = ""
    seq = ""
    for i in File:
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            if name != "":
                D[name] = len(seq)
                seq = ""
            name = ligne[1:].split(" ")[0]
        else:
            seq+=ligne
    D[name] = len(seq)
    return D


def get_transcrip_sequence(File):
    # This function access the DNA sequence of transcripts, and store it in a dictionary with the name
    D = {}
    name = ""
    seq = ""
    for i in File:
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            if name != "":
                D[name] = seq
                seq = ""
            name = ligne[1:].split(" ")[0]
        else:
            seq+=ligne
    D[name]=seq
    return D


###### Step 1. Store sequences in Dictionary
def sequence_storage(File):
    # This function access the DNA sequence of transcripts, and store it in a dictionary with the name
    Dico = {}
    name = ""
    seq = ""
    for i in File:
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            if name != "":
                Dico[name] = seq
                seq = ""
            name = ligne[1:]
            #print name
        else:
            seq+=ligne
    Dico[name]=seq
    return Dico


###### Step 2. Get longest ORF
def get_maximum(Liste,size_whole_transcript,sequence_whole_transcript):
    # This function search for the longest ORF in a list of ORFs. It discards all of the sequences which do not end by a STOP codon
    list_stop_codon = ["TAG", "TAA", "TGA"]
    Max = 0
    list_max = []
    for i in Liste:
        Value = i[2]-i[1]
        if Value > Max:
            if size_whole_transcript-int(i[2])>3: #3
                if sequence_whole_transcript[i[2]:i[2]+3].upper() in list_stop_codon:
                    Max = Value
                    list_max = i
    return list_max # Returns the longest sequence of the list


def get_longest_ORF(DicoNuc, DicoProt,DicoSizeTranscript,DicoseqTranscript):
    # This functoin as for aim to build a dictionary with the longest ORF per transcript without hit
    # It also make sure that the sequence is consistent with the transcript, and retrieve the stop codon from the transcript directly
    dic_per_transcript = {}
    dic_longest_sequence = {}
    dic_longest_prot = {}
    for i in DicoNuc.keys():
        ligne = i.split("_")
        name = ligne[0]+"_"+ligne[1]
        sublist = [ligne[2], int(ligne[3]), int(ligne[4])]
        if name in dic_per_transcript.keys():
            dic_per_transcript[name].append(sublist)
        else:
            dic_per_transcript[name] = [sublist]
    print ("Le nombre total de transcript est: "+str(len(dic_per_transcript)))
    for i in dic_per_transcript.keys():
        T = i.split("_")[1]
        size_whole_transcript = DicoSizeTranscript[T]
        sequence_whole_transcript = DicoseqTranscript[T]
        max_sublist = get_maximum(dic_per_transcript[i],size_whole_transcript,sequence_whole_transcript)
        if len(max_sublist)>0:
            longest_transcript_renamed = i+"_"+max_sublist[0]+"_"+str(max_sublist[1])+"_"+str(max_sublist[2])
            dic_longest_sequence[longest_transcript_renamed]=DicoNuc[longest_transcript_renamed]
            dic_longest_prot[longest_transcript_renamed]=DicoProt[longest_transcript_renamed]
    print ("Le nombre total de transcript avec une ORF est: "+str(len(dic_longest_sequence)))
    return dic_longest_sequence, dic_longest_prot
                

###### Step 3. write_output_file output file
def write_output_file(DicoNuc, DicoProt, Popname):
    # This function writes the final output file for DNA and protein
    F = open((Popname+"_PutativeDeNovoORF_DNA.fa"), "w")
    for i in DicoNuc.keys():
        F.write_output_file(">"+i+"\n")
        F.write_output_file(DicoNuc[i]+"\n")
    F.close()
    F = open((Popname+"_PutativeDeNovoORF_Protein.fa"), "w")
    for i in DicoProt.keys():
        F.write_output_file(">"+i+"\n")
        F.write_output_file(DicoProt[i]+"\n")
    F.close()


def main_function()
    Transcript = open_file("FItranscriptomeAssembly.fa")
    DicoSizeTranscript = get_transcrip_size(Transcript)
    DicoseqTranscript = get_transcrip_sequence(Transcript)
    DNA = open_file("FI_SortedORF_WithoutHit_DNA")
    Protein = open_file("FI_SortedORF_WithoutHit_Protein")
    DicoDNA = sequence_storage(DNA) 
    DicoProtein = sequence_storage(Protein)
    FinalDicoDNA, FinalDicoProtein = get_longest_ORF(DicoDNA, DicoProtein,DicoSizeTranscript,DicoseqTranscript)
    write_output_file(FinalDicoDNA, FinalDicoProtein, "FI")




















