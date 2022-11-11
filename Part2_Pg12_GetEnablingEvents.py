#import pandas as pd
import cv2
import os
os.chdir("path to folder")
import random
from collections import Counter


""" docstring

    This program access further enabling mutations between proto-genes and non-coding homologs

    input
    ------------
    - BLAST resuts
    - proto-genes

    output
    ------------
    - supp enabling mutations
"""


def openFile(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  


#Step 1 Option 1: Get BLAST hits with Best hit
def access_best_hit_rsts(File):
    Liste = []
    print ("le nombre initial de hits est : "+str(len(File)))
    OldQueryName = ""
    OldQueryPop = ""
    for i in File[1:]:
        ligne1 = i.split("\n")[0]
        ligne2 = ligne1.split(",")
        NameQueryDeNovo = ligne2[0]
        NameQueryPop = ligne2[11]
        if NameQueryDeNovo!=OldQueryName:
            OldQueryName = NameQueryDeNovo
            OldQueryPop = NameQueryPop
            Liste.append(ligne1)
        elif NameQueryDeNovo == OldQueryName and NameQueryPop!=OldQueryPop:
            Liste.append(ligne1)
            OldQueryPop = NameQueryPop
    print ("le nombre final de Best hits est : "+str(len(Liste)))
    return Liste


#Step 1 Option 2: Get BLAST hits with All hit
def access_all_hit_rsts(File):
    Liste = []
    for i in File[1:]:
        ligne1 = i.split("\n")[0]
        Liste.append(ligne1)
    return Liste
    

# Step 2 Make File All Enabling Mutations
def determine_enabling_mutation(Liste, SeuilSequenceSize):
    CompteurNoMutation = 0
    ListeFinale = []
    PresenceATG = ""
    PresenceStop = ""
    PresenceCodon = ""
    SameSize = ""
    TranscriptionEvent = ""
    PrematureStopCodon = ""
    Elements = ""
    for i in Liste:
        Data = i.split(",")
        Name = Data[0]
        Hit = int(Data[1])
        Synteny = int(Data[2])
        ATG = int(Data[3])
        Stop = int(Data[4])
        Codon = int(Data[5])
        SizeQuery = Data[6]
        TargetSize = Data[7]
        TranscriptionEvent = Data[22]
        PrematureStopCodon = Data[23]
        Pop = Data[11]
        if Hit!=0:
            Threshold = float(100)*float(TargetSize)/float(SizeQuery)
            if ATG == 0:
                PresenceATG = "NoATG"
                Elements+=PresenceATG
            if Stop == 0:
                PresenceStop = "NoStop"
                if len(Elements)>0:
                    Elements+="_"+PresenceStop
                else:
                    Elements+=PresenceStop
            if Codon ==0:
                PresenceCodon = "NoTripletCodon"
                if len(Elements)>0:
                    Elements+="_"+PresenceCodon
                else:
                    Elements+=PresenceCodon
            if float(Threshold)<float(SeuilSequenceSize):
                SameSize = "DifferentSize"
                if len(Elements)>0:
                    Elements+="_"+SameSize
                else:
                    Elements+=SameSize
            if TranscriptionEvent!="TranscriptFull":
                TranscriptionEvent = "NoTranscript"
                if len(Elements)>0:
                    Elements+="_"+TranscriptionEvent
                else:
                    Elements+=TranscriptionEvent
            if PrematureStopCodon=="PrematureStopBefore75":
                PrematureStopCodon = "PrematureStop"
                if len(Elements)>0:
                    Elements+="_"+PrematureStopCodon
                else:
                    Elements+=PrematureStopCodon
            if Elements == "":
                CompteurNoMutation+=1
            NewListe = [Name, Pop, Elements]
            ListeFinale.append(NewListe)
        PresenceATG = ""
        PresenceStop = ""
        PresenceCodon = ""
        SameSize = ""
        TranscriptionEvent = ""
        PrematureStopCodon = ""
        Elements = ""
    print ("Number of sequences with no mutations and potentialy coding : "+str(CompteurNoMutation))
    return ListeFinale
            
            
# Step 3 Make Final Files
def build_final_file(Liste):
    DicoAllMutations = {}
    Dico6Mutations = {}
    TotalMutations = 0
    for i in Liste:
        Mutations = i[2]
        if Mutations != "":
            if Mutations not in DicoAllMutations.keys():
                DicoAllMutations[Mutations] = 1
            else:
                DicoAllMutations[Mutations]+=1
        if Mutations!= "":
            TotalMutations+=1
            Data = Mutations.split("_")
            for j in Data:
                if j not in Dico6Mutations.keys():
                    Dico6Mutations[j] = 1
                else:
                    Dico6Mutations[j] +=1
    F = open("ResultsAllEnablingMutations.txt", "w")
    F.write("NameMutation,Number,Percentage"+"\n")
    for i in DicoAllMutations.keys():
        Percentage = float(DicoAllMutations[i])*float(100)/float(TotalMutations)
        F.write(i+","+str(DicoAllMutations[i])+","+str(Percentage)+"\n")
    F.close()
    F = open("Results5EnablingMutations.txt", "w")
    F.write("NameMutation,Number,Percentage"+"\n")
    for i in Dico6Mutations.keys():
        Percentage = float(Dico6Mutations[i])*float(100)/float(TotalMutations)
        F.write(i+","+str(Dico6Mutations[i])+","+str(Percentage)+"\n")
    F.close()
        
    
# Supp analyse ...  make graph of percentage of transcirption ListeAllEnablingEvents
def build_supp_file1(Liste):
    Dico = {"UncompleteTranscript":0, "ReverseTranscript":0, "UncompleteReverseTranscript":0, "NoTranscript":0}
    T = 0
    for i in Liste:
        Data = i.split(",")
        TranscriptionEvent = Data[22]
        if TranscriptionEvent == "TranscriptHalf":
            Dico["UncompleteTranscript"]+=1
            T+=1
        elif TranscriptionEvent == "TranscriptReverse":
            Dico["ReverseTranscript"]+=1
            T+=1
        elif TranscriptionEvent == "TranscriptHalfReverse":
            Dico["UncompleteReverseTranscript"]+=1
            T+=1
        elif TranscriptionEvent == "NoTranscript":
            Dico["NoTranscript"]+=1
            T+=1
    F = open("Data_Transcript_Mutation.txt", "w")
    F.write("TranscriptEvent,Percentage"+"\n")
    for i in Dico.keys():
        F.write(i+",")
        F.write(str(float(Dico[i])*float(100)/float(T))+"\n")
    F.close()


def main_function():
    BLAST_AllResults = openFile("FinalBLAST.rst")
    ListeBestHits = access_best_hit_rsts(BLAST_AllResults)
    #ListeAllHits = access_all_hit_rsts(BLAST_AllResults)
    ListeAllEnablingEvents = determine_enabling_mutation(ListeBestHits, 90)
    build_final_file(ListeAllEnablingEvents)
    build_supp_file1(ListeBestHits)




