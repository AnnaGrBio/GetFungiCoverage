#import pandas as pd
import cv2
import os
os.chdir("path to folder")
import random
from collections import Counter


""" docstring

    This program search transcription status and premature stop codons in non-coding othologs

    input
    ------------
    - result file from BLAST
    - genomes
    - gtf transcriptome assemblies

    output
    ------------
    - supp info about enabling mutations
"""


def open_file(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  


# Step 1 Store all Transcript event in Dicos
def fill_transcript_dic(Dico, GTF):
    for i in GTF[2:]:
        ligne = i.split()
        Chrom = ligne[0]
        Type = ligne[2]
        Start = int(ligne[3])
        End = int(ligne[4])
        Direction = ligne[6]
        if Type == "transcript":
            TPM = float(ligne[17][1:len(ligne[17])-2])
            if TPM > 0.5:
                Subliste = [Start, End, Direction]
                if Chrom in Dico.keys():
                    Dico[Chrom].append(Subliste)
                else:
                    Dico[Chrom] = [Subliste]
    

# Step 2 Store all Pops chromosomes in dictionaries
def fill_genomes_dic(Dico, GenomeFile):
    Name = ""
    Seq = ""
    for i in GenomeFile:
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            if Name != "":
                Dico[Name] = Seq
                Seq = ""
            Name = ligne[1:]
        else:
            Seq += ligne.upper()
    Dico[Name] = Seq

        
# Step 3 MakeListe With Direction ATG
def sort_result_file(File):
    NewFile = []
    print ("Total NUAber of hits : "+str(len(File)))
    NewFile.append(File[0])
    for i in File[1:]:
        ligne = i.split("\n")[0]
        Data = ligne.split(",")
        Seq = Data[12]
        if len(Seq)<60:
            continue
        else:
            NewFile.append(i)
    print ("Total NUAber of sorted hits (hit len sup 60) : "+str(len(NewFile)))
    return NewFile


def build_list_best_blast_hit(InnitialResultFile, DicoGenomeFI, DicoGenomeDK, DicoGenomeES, DicoGenomeSE, DicoGenomeUA, DicoGenomeTR, DicoGenomeZI):
    ListeFinale = []
    TotalForward = 0
    TotalReverse = 0
    TotalOther = 0
    Total = 0
    c = 0
    SortedResultFile = sort_result_file(InnitialResultFile)
    for i in SortedResultFile[1:]:
        ligne = i.split("\n")[0]
        Data = ligne.split(",")
        PopTarget = Data[11]
        Start = Data[12][0:3]
        Chromosome = Data[10]
        TargetStart = Data[8]
        if TargetStart != "NaN":
            TargetStart = int(TargetStart)
            if PopTarget == "FI":
                AllSeq = DicoGenomeFI[Chromosome]
                PutativeStart = AllSeq[TargetStart:TargetStart+3]
                if Start == PutativeStart:
                    TotalForward += 1
                    Direction = "+"
                elif Start != PutativeStart and len(Start) == 3:
                    Direction = "-"
                    TotalReverse += 1
                else:
                    Direction = "NaN"
                    TotalOther += 1
            elif PopTarget == "DK":
                AllSeq = DicoGenomeDK[Chromosome]
                PutativeStart = AllSeq[TargetStart:TargetStart+3]
                if Start == PutativeStart:
                    TotalForward += 1
                    Direction = "+"
                elif Start != PutativeStart and len(Start) == 3:
                    Direction = "-"
                    TotalReverse += 1
                else:
                    Direction = "NaN"
                    TotalOther += 1
            elif PopTarget == "ES":
                AllSeq = DicoGenomeES[Chromosome]
                PutativeStart = AllSeq[TargetStart:TargetStart+3]
                if Start == PutativeStart:
                    TotalForward += 1
                    Direction = "+"
                elif Start != PutativeStart and len(Start) == 3:
                    Direction = "-"
                    TotalReverse += 1
                else:
                    Direction = "NaN"
                    TotalOther += 1
            elif PopTarget == "SE":
                AllSeq = DicoGenomeSE[Chromosome]
                PutativeStart = AllSeq[TargetStart:TargetStart+3]
                if Start == PutativeStart:
                    TotalForward += 1
                    Direction = "+"
                elif Start != PutativeStart and len(Start) == 3:
                    Direction = "-"
                    TotalReverse += 1
                else:
                    Direction = "NaN"
                    TotalOther += 1
            elif PopTarget == "UA":
                AllSeq = DicoGenomeUA[Chromosome]
                PutativeStart = AllSeq[TargetStart:TargetStart+3]
                if Start == PutativeStart:
                    TotalForward += 1
                    Direction = "+"
                elif Start != PutativeStart and len(Start) == 3:
                    Direction = "-"
                    TotalReverse += 1
                else:
                    Direction = "NaN"
                    TotalOther += 1
            elif PopTarget == "TR":
                AllSeq = DicoGenomeTR[Chromosome]
                PutativeStart = AllSeq[TargetStart:TargetStart+3]
                if Start == PutativeStart:
                    TotalForward += 1
                    Direction = "+"
                elif Start != PutativeStart and len(Start) == 3:
                    Direction = "-"
                    TotalReverse += 1
                else:
                    Direction = "NaN"
                    TotalOther += 1
            elif PopTarget == "ZI":
                AllSeq = DicoGenomeZI[Chromosome]
                PutativeStart = AllSeq[TargetStart:TargetStart+3]
                if Start == PutativeStart:
                    TotalForward += 1
                    Direction = "+"
                elif Start != PutativeStart and len(Start) == 3:
                    Direction = "-"
                    TotalReverse += 1
                else:
                    Direction = "NaN"
                    TotalOther += 1
            Total+=1
        else:
            Direction = "NaN"
        NewLigne = ligne + "," + Direction
        ListeFinale.append(NewLigne)
    print Total
    print ("TotalForward = " + str(float(100)*float(TotalForward)/float(Total)))
    print ("TotalReverse = " + str(float(100)*float(TotalReverse)/float(Total)))
    print ("TotalOther = " + str(float(100)*float(TotalOther)/float(Total)))
    return ListeFinale


# Step 4 MakeListe Finale with transcription information
def access_transcription_status(Start, Stop, Direction, ListeTranscript):
    TranscriptFull = False
    TranscriptHalf = False
    TranscriptReverse = False
    TranscriptHalfReverse = False
    NoTranscript = True
    for i in ListeTranscript:
        StartTranscript = i[0]
        EndTranscript = i[1]
        DirectionTranscript = i[2]
        if StartTranscript<=Start and EndTranscript>Start and EndTranscript<Stop:
            if Direction == DirectionTranscript:
                TranscriptHalf = True
                NoTranscript = False
            else:
                TranscriptHalfReverse = True
                NoTranscript = False
        elif StartTranscript>Start and StartTranscript<Stop and EndTranscript>=Stop:
            if Direction == DirectionTranscript:
                TranscriptHalf = True
                NoTranscript = False
            else:
                TranscriptHalfReverse = True
                NoTranscript = False
        elif StartTranscript>Start and EndTranscript<Stop:
            if Direction == DirectionTranscript:
                TranscriptHalf = True
                NoTranscript = False
            else:
                TranscriptHalfReverse = True
                NoTranscript = False
        elif StartTranscript<=Start and EndTranscript>=Stop:
            if Direction == DirectionTranscript:
                TranscriptFull = True
                NoTranscript = False
                break
            else:
                TranscriptReverse = True
                NoTranscript = False
    if TranscriptFull == True:
        return "TranscriptFull"
    elif TranscriptFull == False and TranscriptHalf == True:
        return "TranscriptHalf"
    elif TranscriptFull == False and TranscriptHalf == False and TranscriptReverse == True:
        return "TranscriptReverse"
    elif TranscriptFull == False and TranscriptHalf == False and TranscriptReverse == False and TranscriptHalfReverse == True:
        return "TranscriptHalfReverse"
    else:
        return "NoTranscript"


def retrieve_transcription_event(ListeFinaleStep1, DicoTranscriptFI, DicoTranscriptDK, DicoTranscriptES, DicoTranscriptSE, DicoTranscriptUA, DicoTranscriptTR, DicoTranscriptZI):
    FinalListe = []
    NbTranscriptFull = 0
    NbTranscriptHalf = 0
    NbTranscriptReverse = 0
    NbTranscriptHalfReverse = 0
    NbNoTranscript = 0
    Total = 0
    for i in ListeFinaleStep1:
        Data = i.split(",")
        PopTarget = Data[11]
        Chromosome = Data[10]
        TargetStart = Data[8]
        TargetStop = Data[9]
        ORFdirection = Data[21]
        if TargetStart!= "NaN" and TargetStop!= "NaN":
            TargetStart = int(TargetStart)
            TargetStop = int(TargetStop)
            if PopTarget == "FI":
                ListeTranscripts = DicoTranscriptFI[Chromosome]
                TranscriptionStatus = access_transcription_status(TargetStart, TargetStop, ORFdirection, ListeTranscripts)
                if TranscriptionStatus == "TranscriptFull":
                    NbTranscriptFull+=1
                elif TranscriptionStatus == "TranscriptHalf":
                    NbTranscriptHalf+=1
                elif TranscriptionStatus == "TranscriptReverse":
                    NbTranscriptReverse+=1
                elif TranscriptionStatus == "TranscriptHalfReverse":
                    NbTranscriptHalfReverse+=1
                else:
                    NbNoTranscript+=1
                Total+=1
            elif PopTarget == "DK":
                ListeTranscripts = DicoTranscriptDK[Chromosome]
                TranscriptionStatus = access_transcription_status(TargetStart, TargetStop, ORFdirection, ListeTranscripts)
                if TranscriptionStatus == "TranscriptFull":
                    NbTranscriptFull+=1
                elif TranscriptionStatus == "TranscriptHalf":
                    NbTranscriptHalf+=1
                elif TranscriptionStatus == "TranscriptReverse":
                    NbTranscriptReverse+=1
                elif TranscriptionStatus == "TranscriptHalfReverse":
                    NbTranscriptHalfReverse+=1
                else:
                    NbNoTranscript+=1
                Total+=1
            elif PopTarget == "ES":
                ListeTranscripts = DicoTranscriptES[Chromosome]
                TranscriptionStatus = access_transcription_status(TargetStart, TargetStop, ORFdirection, ListeTranscripts)
                if TranscriptionStatus == "TranscriptFull":
                    NbTranscriptFull+=1
                elif TranscriptionStatus == "TranscriptHalf":
                    NbTranscriptHalf+=1
                elif TranscriptionStatus == "TranscriptReverse":
                    NbTranscriptReverse+=1
                elif TranscriptionStatus == "TranscriptHalfReverse":
                    NbTranscriptHalfReverse+=1
                else:
                    NbNoTranscript+=1
                Total+=1
            elif PopTarget == "SE":
                ListeTranscripts = DicoTranscriptSE[Chromosome]
                TranscriptionStatus = access_transcription_status(TargetStart, TargetStop, ORFdirection, ListeTranscripts)
                if TranscriptionStatus == "TranscriptFull":
                    NbTranscriptFull+=1
                elif TranscriptionStatus == "TranscriptHalf":
                    NbTranscriptHalf+=1
                elif TranscriptionStatus == "TranscriptReverse":
                    NbTranscriptReverse+=1
                elif TranscriptionStatus == "TranscriptHalfReverse":
                    NbTranscriptHalfReverse+=1
                else:
                    NbNoTranscript+=1
                Total+=1
            elif PopTarget == "TR":
                ListeTranscripts = DicoTranscriptTR[Chromosome]
                TranscriptionStatus = access_transcription_status(TargetStart, TargetStop, ORFdirection, ListeTranscripts)
                if TranscriptionStatus == "TranscriptFull":
                    NbTranscriptFull+=1
                elif TranscriptionStatus == "TranscriptHalf":
                    NbTranscriptHalf+=1
                elif TranscriptionStatus == "TranscriptReverse":
                    NbTranscriptReverse+=1
                elif TranscriptionStatus == "TranscriptHalfReverse":
                    NbTranscriptHalfReverse+=1
                else:
                    NbNoTranscript+=1
                Total+=1
            elif PopTarget == "UA":
                ListeTranscripts = DicoTranscriptUA[Chromosome]
                TranscriptionStatus = access_transcription_status(TargetStart, TargetStop, ORFdirection, ListeTranscripts)
                if TranscriptionStatus == "TranscriptFull":
                    NbTranscriptFull+=1
                elif TranscriptionStatus == "TranscriptHalf":
                    NbTranscriptHalf+=1
                elif TranscriptionStatus == "TranscriptReverse":
                    NbTranscriptReverse+=1
                elif TranscriptionStatus == "TranscriptHalfReverse":
                    NbTranscriptHalfReverse+=1
                else:
                    NbNoTranscript+=1
                Total+=1
            elif PopTarget == "ZI":
                ListeTranscripts = DicoTranscriptZI[Chromosome]
                TranscriptionStatus = access_transcription_status(TargetStart, TargetStop, ORFdirection, ListeTranscripts)
                if TranscriptionStatus == "TranscriptFull":
                    NbTranscriptFull+=1
                elif TranscriptionStatus == "TranscriptHalf":
                    NbTranscriptHalf+=1
                elif TranscriptionStatus == "TranscriptReverse":
                    NbTranscriptReverse+=1
                elif TranscriptionStatus == "TranscriptHalfReverse":
                    NbTranscriptHalfReverse+=1
                else:
                    NbNoTranscript+=1
                Total+=1
        else:
            TranscriptionStatus = "NaN"
        TheNewLigne = i+","+TranscriptionStatus
        FinalListe.append(TheNewLigne)
    print ("Nb Transcript full = "+str(NbTranscriptFull))
    print ("Nb transcripts half = "+str(NbTranscriptHalf))
    print ("Nb Transcripts reverse = "+ str(NbTranscriptReverse))
    print ("Nb Transcripts half reverse = "+str(NbTranscriptHalfReverse))
    print ("Nb no transcript  ="+str(NbNoTranscript))
    print ("Nb total seqs = "+str(Total))
    return FinalListe


# Step 5 Add info about premature Stop codons
def find_premature_stops(ListeFinaleStep2):
    ListeFinalStep3 = []
    ListeStops = ["TAG", "TAA", "TGA"]
    TotalPrematureStopBefore75 = 0
    TotalPrematureStopAfter75 = 0
    TotalNoPrematureStop = 0
    c = 0
    for i in ListeFinaleStep2:
        PrematureStopBefore75 = False
        PrematureStopAfter75 = False
        Data = i.split(",")
        SequenceTarget = Data[12][3:len(Data[12])-3]
        SeqToAnalyse = SequenceTarget[0:3]
        Threshold = float(75)*float(len(SequenceTarget))/float(100)
        PosCurseur = 0
        while len(SeqToAnalyse) == 3:
            if SeqToAnalyse in ListeStops:
                if PosCurseur<Threshold:
                    PrematureStopBefore75 = True
                else:
                    PrematureStopAfter75 = True
                break
            else:
                SequenceTarget = SequenceTarget[3:]
                if len(SequenceTarget)>=3:
                    SeqToAnalyse = SequenceTarget[0:3]
                else:
                    SeqToAnalyse = (SequenceTarget)
            PosCurseur+=3
        if PrematureStopBefore75 == True:
            TheNewLigne = i+","+"PrematureStopBefore75"
            TotalPrematureStopBefore75+=1
        elif PrematureStopAfter75 == True:
            TheNewLigne = i+","+"PrematureStopAfter75"
            TotalPrematureStopAfter75 +=1
        else:
            TheNewLigne = i+","+"NoPrematureStop"
            TotalNoPrematureStop+=1
        ListeFinalStep3.append(TheNewLigne)
    print ("NUAber of Premature Stop before 75 : "+str(TotalPrematureStopBefore75))
    print ("NUAber of Premature Stop After 75 : "+str(TotalPrematureStopAfter75))
    print ("NUAber of No Premature Stop : "+str(TotalNoPrematureStop))
    return ListeFinalStep3
            

# Step 6 MakeListe Finale with transcription information
def MakeFinalFile(ListeFinaleStep3, InnitialResultFile):
    F = open("FinalBLAST.rst", "w")
    L1 = InnitialResultFile[0].split("\n")[0]
    NewL1 = L1+",TargetORFdirection,Taraccess_transcription_status,PrematureStop"+"\n"
    F.write(NewL1)
    for i in ListeFinaleStep3:
        F.write(i+"\n")
    F.close()
    

def main_function():
    InnitialResultFile = open_file('211214_results.csv')
    GTF_FI = open_file("FItranscriptome.gtf")
    GTF_DK = open_file("DKtranscriptome.gtf")
    GTF_ES = open_file("EStranscriptome.gtf")
    GTF_SE = open_file("SEtranscriptome.gtf")
    GTF_UA = open_file("UAtranscriptome.gtf")
    GTF_TR = open_file("TRtranscriptome.gtf")
    GTF_ZI = open_file("ZItranscriptome.gtf")
    FIgenome = open_file("FIfinalGenome.masked.fa")
    DKgenome = open_file("DKfinalGenome.masked.fa")
    ESgenome = open_file("ESfinalGenome.masked.fa")
    SEgenome = open_file("SEfinalGenome.masked.fa")
    UAgenome = open_file("UAfinalGenome.masked.fa")
    TRgenome = open_file("TRfinalGenome.masked.fa")
    ZIgenome = open_file("ZIfinalGenome.masked.fa")
    DicoTranscriptFI = {}
    DicoTranscriptDK = {}
    DicoTranscriptES = {}
    DicoTranscriptSE = {}
    DicoTranscriptUA = {}
    DicoTranscriptTR = {}
    DicoTranscriptZI = {}
    fill_transcript_dic(DicoTranscriptFI, GTF_FI)
    fill_transcript_dic(DicoTranscriptDK, GTF_DK)
    fill_transcript_dic(DicoTranscriptES, GTF_ES)
    fill_transcript_dic(DicoTranscriptSE, GTF_SE)
    fill_transcript_dic(DicoTranscriptUA, GTF_UA)
    fill_transcript_dic(DicoTranscriptTR, GTF_TR)
    fill_transcript_dic(DicoTranscriptZI, GTF_ZI)
    DicoGenomeFI = {}
    DicoGenomeDK = {}
    DicoGenomeES = {}
    DicoGenomeSE = {}
    DicoGenomeUA = {}
    DicoGenomeTR = {}
    DicoGenomeZI = {}
    fill_genomes_dic(DicoGenomeFI, FIgenome)
    fill_genomes_dic(DicoGenomeDK, DKgenome)
    fill_genomes_dic(DicoGenomeES, ESgenome)
    fill_genomes_dic(DicoGenomeSE, SEgenome)
    fill_genomes_dic(DicoGenomeUA, UAgenome)
    fill_genomes_dic(DicoGenomeTR, TRgenome)
    fill_genomes_dic(DicoGenomeZI, ZIgenome)
    ListeFinaleStep1 = build_list_best_blast_hit(InnitialResultFile, DicoGenomeFI, DicoGenomeDK, DicoGenomeES, DicoGenomeSE, DicoGenomeUA, DicoGenomeTR, DicoGenomeZI)
    ListeFinaleStep2 = retrieve_transcription_event(ListeFinaleStep1, DicoTranscriptFI, DicoTranscriptDK, DicoTranscriptES, DicoTranscriptSE, DicoTranscriptUA, DicoTranscriptTR, DicoTranscriptZI)
    ListeFinaleStep3 = find_premature_stops(ListeFinaleStep2)
    MakeFinalFile(ListeFinaleStep3, InnitialResultFile)




















































































