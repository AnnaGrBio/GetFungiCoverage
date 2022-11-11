import cv2
import os
os.chdir("path to folder")
import random


""" docstring

    This function builds the file containing main informations about the sequences that can be accessed from GTF file and previous scripts
    
    input
    ------------
    - transcriptome
    - unspliced ORFs positions
    - DNA seq of de novo ORFs
    
    output
    ------------
    - information file
"""


def open_file(NameFile):
    F = open(NameFile, "r")
    L = F.readlines()
    return L  


def retrieve_information_GTF(File):
    # This function retrieves and store information from the gtf file for each transcript:
    # Direction, Coverage, FPKM, TPM, start, end, size unspliced
    dic_transcript = {}
    NbExon = 0
    transcript_properties = []
    for i in File[2:]:
        ligne = i.split("	")
        if ligne[2] == "transcript":
            if NbExon!=0:
                transcript_properties.append(NbExon)
                dic_transcript[TranscriptName] = transcript_properties
                transcript_properties = []
                NbExon = 0
            Chrom = ligne[0]
            transcript_start = ligne[3]
            transcript_end = ligne[4]
            Direction = ligne[6]
            lala = ligne[8].split(" ")
            TranscriptName = lala[3][1:len(lala[3])-2]
            Cov = lala[5][1:len(lala[5])-2]
            FPKM = lala[7][1:len(lala[7])-2]
            TPM = lala[9][1:len(lala[9])-3]
            transcript_properties = [Chrom,transcript_start,transcript_end,Direction,Cov,FPKM,TPM]
        else:
            NbExon += 1
    transcript_properties.append(NbExon)
    dic_transcript[TranscriptName] = transcript_properties
    return dic_transcript


def build_dic_ORF(File):
    # This function builds a dictionary of ORFs sequences
    D = {}
    for i in File[1:]:
        ligne = i.split(",")
        Name = ligne[0]
        Start = ligne[1]
        End = ligne[2]
        D[Name] = [Start, End]
    return D


def Access_all_de_novo_genes(File):
    # This function retrieves the de novo ORFs
    L = []
    for i in File:start
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            L.append(ligne[1:])
    return L


def build_final_dic_info_and_file(DicoElementsGTF,DicoStartStop,ListAllDeNovoGenes, Name):
    # This function access final informations from de novo genes and builds a final info file
    DicoGenes = {}
    for i in DicoElementsGTF.keys():
        TranscriptName = i
        GeneName = TranscriptName.split(".")[0]+"."+TranscriptName.split(".")[1]
        ListElements = DicoElementsGTF[TranscriptName]
        Chrom = ListElements[0]
        transcript_start = ListElements[1]
        transcript_end = ListElements[2]
        Direction = ListElements[3]
        Cov = ListElements[4]
        FPKM = ListElements[5]
        TPM = ListElements[6]
        NbExon = ListElements[7]
        ORF = False
        ORFname = ""
        for ORF in DicoStartStop.keys():
            TranscriptORF = ORF.split("_")[1]
            if TranscriptName == TranscriptORF:
                DeNovoORF = "DeNovoORF"
                unspliced_ORF_start = DicoStartStop[ORF][0]
                unspliced_ORF_end = DicoStartStop[ORF][1]
                ORFname = ORF
                ORF = True
                break
        if ORF == True:
            BigListe = [TranscriptName,NbExon,Chrom,transcript_start,transcript_end,Direction,Cov,FPKM,TPM,DeNovoORF,ORFname,unspliced_ORF_start,unspliced_ORF_end]
        else:
            BigListe = [TranscriptName,NbExon,Chrom,transcript_start,transcript_end,Direction,Cov,FPKM,TPM,"NotDeNovoORF","-","-","-"]
        if GeneName in DicoGenes.keys():
            DicoGenes[GeneName].append(BigListe)
        else:
            DicoGenes[GeneName] = [BigListe]
    F = open(Name + "_InformationFile", "w")
    F.write("GeneName,NbTranscripts,TranscriptName,NbExon,Chrom,transcript_start,transcript_end,Direction,Cov,FPKM,TPM,DeNovoORF,ORFname,unspliced_ORF_start,unspliced_ORF_end"+"\n")
    for Gene in DicoGenes.keys():
        NbTranscripts = len(DicoGenes[Gene])
        ListeTranscripts = DicoGenes[Gene]
        for sublists in ListeTranscripts:
            F.write(Gene + "," + str(NbTranscripts))
            for info in sublists:
                F.write("," + str(info))
            F.write("\n")
    F.close()
    

def main_function():
    GTFfile = open_file("FItranscriptome.gtf")
    DicoElementsGTF = retrieve_information_GTF(GTFfile)
    ORFstartAndStop = open_file("FI_UnsplicedORFpositions")
    DicoStartStop = build_dic_ORF(ORFstartAndStop)
    ORF = open_file("FI_DeNovoORFwithStop_DNA.fa")
    ListAllDeNovoGenes = Access_all_de_novo_genes(ORF)
    build_final_dic_info_and_file(DicoElementsGTF,DicoStartStop,ListAllDeNovoGenes, "AK5")



