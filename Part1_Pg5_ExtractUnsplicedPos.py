import cv2
import os
os.chdir("path to folder")
import random

""" docstring

    This program aims to:
    - Rebuild unspliced transcript
    - Rebuild spliced sequence and compare it to the spliced sequence extracted from get ORF to make sure all is working
    - Extract positions of spliced transcript and position in genome of the unspliced ORFs
    
    input
    ------------
    - genome of the line
    - gtf file of the transcriptome
    - file of de novo ORFs
    
    output
    ------------
    - information file
"""


def openFile(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  


############################## Step 1. Store the ORF and their position in the transcript in a dictionary
def build_ORF_dictionary(File):
    # This function builds a dictionary with all ORFs associated to their positions in the transcripts and an
    # other on associated to their sequence
    DicoSeqs = {}
    DicoPos = {}
    Title = ""
    Seq = ""
    for i in File:
        if i[0] == ">":
            if Title!= "":
                DicoSeqs[Title] = Seq
                Seq = ""
            Title = i[1:].split("\n")[0]
            Liste = Title.split("_")
            NameTranscript = Liste[1]
            Start = int(Liste[3])
            End = int(Liste[4])
            ListePos = [Start, End]
            DicoPos[Title] = ListePos
        else:
            ligne = i.split("\n")[0]
            Seq+=ligne
    DicoSeqs[Title] = Seq
    return DicoSeqs, DicoPos


############################## Step 2. Store Unspliced Positions
def build_dictionary_unspliced_positions(File):
    # This function determine the positions in chromosomes of unspliced transcripts and store them in a dictionary
    D = {}
    TranscriptName = ""
    for i in File:
        if i[0] != "#":
            ligne = i.split("	")
            if ligne[2] == "transcript":
                if TranscriptName!="":
                    L = [Chrom,Direction,ListExon]
                    D[TranscriptName] = L
                    #print D
                    #break
                ListExon = []
                Chrom = ligne[0]
                #Start = int(ligne[3])
                #End = int(ligne[4])
                Direction = ligne[6]
                #print Chrom,Start,End,Direction
                SubLigne = ligne[8].split(";")
                SubSubLigne = SubLigne[1].split(" ")
                TranscriptName = SubSubLigne[2][1:len(SubSubLigne[2])-1]
            else:
                StartExon = int(ligne[3])
                EndExon = int(ligne[4])
                SubList = [StartExon,EndExon]
                ListExon.append(SubList)
    L = [Chrom,Direction,ListExon]
    D[TranscriptName] = L        
    return D


############################## Step 3. Rebuild spliced Sequences
def build_dictionary_chromosomes(File):
    DicoChromosomes = {}
    Seq = ""
    Name = ""
    for i in File:
        if i[0] == ">":
            if Name != "":
                DicoChromosomes[Name] = Seq
                Seq = ""
            Name = i[1:].split("\n")[0]
        else:
            Data = i.split("\n")[0]
            Seq += Data
    DicoChromosomes[Name] = Seq
    return DicoChromosomes
    
    
def reverse_transcript(NewTranscript):
    # this function revere the transcripts
    NewTRanscript = ""
    Dico = {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N","a":"t", "t":"a", "g":"c", "c":"g", "n":"n"}
    L = list(NewTranscript.strip())
    for i in reversed(L):
        NewLetter = Dico[i]
        NewTRanscript+=NewLetter
    return NewTRanscript
        
        
        
def build_new_spliced_transcript(Liste,DicoChrom,Dico):
    # this function build the spliced transcript in order to compare it to the extracted spliced transcript
    Chrom = Liste[0]
    Direction = Liste[1]
    ListeExon = Liste[2]
    SeqChrom = DicoChrom[Chrom]
    c = 3
    if Direction == "+":
        NewTranscript = ""
        NewTranscript+=SeqChrom[ListeExon[0][0]-3:ListeExon[0][0]-1]
        DicoPosSplicedTranscript={1:ListeExon[0][0]-3, 2:ListeExon[0][0]-2}
        for i in ListeExon:
            Exon = SeqChrom[i[0]-1:i[1]]
            for j in range(i[0]-1,i[1]):
                DicoPosSplicedTranscript[c]=j
                c+=1
            NewTranscript+=Exon
        if len(SeqChrom)>ListeExon[len(ListeExon)-1][1]:
            NewTranscript+=SeqChrom[ListeExon[len(ListeExon)-1][1]]
        DicoPosSplicedTranscript[c] = ListeExon[len(ListeExon)-1][1]
    elif Direction == "-" or Direction == ".": # Decouverte donc, les "." sont des reverse
        NewTranscript = ""
        NewTranscript+=SeqChrom[ListeExon[0][0]-3:ListeExon[0][0]-1]
        DicoPosSplicedTranscript={1:ListeExon[0][0]-3, 2:ListeExon[0][0]-2}
        for i in ListeExon:
            Exon = SeqChrom[i[0]-1:i[1]]
            for j in range(i[0]-1,i[1]):
                DicoPosSplicedTranscript[c]=j
                c+=1
            NewTranscript+=Exon
        if len(SeqChrom)>ListeExon[len(ListeExon)-1][1]:
            NewTranscript+=SeqChrom[ListeExon[len(ListeExon)-1][1]]
        DicoPosSplicedTranscript[c] = ListeExon[len(ListeExon)-1][1]
        NewTranscript = reverse_transcript(NewTranscript)
    return NewTranscript,DicoPosSplicedTranscript


def find_max_value_dic(D):
    Max =0
    for i in D.keys():
        if i>Max:
            Max = i
    return Max


def rebuild_spliced_seq(File, Dico_Info, DicoSeqs):
    # This function rebuild the ORF and extracts the stop codon
    FinalDicoPos = {}
    DicoChrom = build_dictionary_chromosomes(File)
    ListeOfficialStopCodon = ["TAG","TAA", "TGA"]
    c = 0
    CompteurFalse = 0
    DicoTranslation = {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N","a":"t", "t":"a", "g":"c", "c":"g", "n":"n"}
    for NameTranscript in Dico_Info.keys():
        c+=1
        Liste = Dico_Info[NameTranscript]
        NewSplicedTranscript,DicoStockage = build_new_spliced_transcript(Liste,DicoChrom,Dico_Info)
        Ajouter = True
        for longName in DicoSeqs.keys():
            ligne = longName.split("_")
            NameShortTranscript = ligne[1]
            if NameShortTranscript == NameTranscript:
                #print longName
                StartORF = int(ligne[3])
                EndORF = int(ligne[4])
                SplicedORF = DicoSeqs[longName]
                if  Liste[1] == "+":
                    NewSplicedORF = NewSplicedTranscript[StartORF+1:EndORF+2]
                    ChromosomeSeq = DicoChrom[Liste[0]]
                    StartORFinGenome = DicoStockage[StartORF+2]
                    EndORFinGenome = DicoStockage[EndORF+2]
                    EndORF_StopCodonInGenome = DicoStockage[EndORF+2]+3
                    StopCodon = ChromosomeSeq[EndORFinGenome+1:EndORFinGenome+4].upper()
                    if StopCodon not in ListeOfficialStopCodon:
                        CompteurFalse += 1
                        Ajouter = False
                    if SplicedORF != NewSplicedORF:
                        CompteurFalse += 1
                        Ajouter = False
                elif Liste[1] == "-" or Liste[1] == ".":
                    NewSplicedORF = NewSplicedTranscript[StartORF:EndORF+1]
                    ChromosomeSeq = DicoChrom[Liste[0]]
                    StartORFinGenome = DicoStockage[find_max_value_dic(DicoStockage)-StartORF]
                    EndORFinGenome = DicoStockage[find_max_value_dic(DicoStockage)-EndORF]
                    EndORF_StopCodonInGenome = EndORFinGenome-3
                    StopCodon = reverse_transcript(ChromosomeSeq[EndORF_StopCodonInGenome:EndORFinGenome]).upper()
                    if SplicedORF != NewSplicedORF:
                        CompteurFalse += 1
                        Ajouter = False
                    if StopCodon not in ListeOfficialStopCodon:
                        CompteurFalse += 1
                        Ajouter = False
                if Ajouter == True:
                    FinalDicoPos[longName] = [StartORFinGenome,EndORF_StopCodonInGenome,StopCodon]
                break
    return FinalDicoPos
                
                        
### Make final file
def make_final_file(Dico, Name):
    F = open(Name+"_UnsplicedORFpositions", "w")
    F.write("TranscriptName,StartORFinGenome,EndORFinGenome, StopCodon"+"\n")
    for i in Dico.keys():
        F.write(str(i)+","+str(Dico[i][0])+","+str(Dico[i][1])+","+str(Dico[i][2])+"\n")
    F.close()


def main_function():
    Nucs = openFile("FI_DeNovoORF_DNA.fa")
    DicoSeqs, DicoPosSpliced = build_ORF_dictionary(Nucs)
    GTF = openFile("FItranscriptome.gtf")
    Dico_Info = build_dictionary_unspliced_positions(GTF)
    Genome = openFile("FIfinalGenome.masked.fa")
    print len(DicoSeqs)
    FinalDicoPosUnsplicedORF = rebuild_spliced_seq(Genome, Dico_Info, DicoSeqs)
    print len(FinalDicoPosUnsplicedORF)
    make_final_file(FinalDicoPosUnsplicedORF, "AK5")






