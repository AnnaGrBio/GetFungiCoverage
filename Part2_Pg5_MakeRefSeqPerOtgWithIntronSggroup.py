import cv2
import os
os.chdir("path to folder")
import random
from collections import Counter


""" docstring

    This program select a referent sequences for orthogroups in which the ORF is not present in all lines and contain at least one intron
    It onlu search orthogroups where one sequence per line is present

    input
    ------------
    - orthogroup file from previous step

    output
    ------------
    - file with referent seq per orthogropus + some info about the orthogroup/sequence

"""

def open_file(NameFile):
    F = open(NameFile, "r")
    L = F.readlines()
    return L  


def retrieve_missing_lines(ListePresentes):
    ListeAllPops = ["AK5", "DK5", "GI5", "SW5", "UM", "YE", "Zamb"]
    ListeMissing = []
    for pop in ListeAllPops:
        if pop not in ListePresentes:
            ListeMissing.append(pop)
    return ListeMissing


def search_introns(Sequence):
	Presence = False
	for i in Sequence:
		if i == "a" or i =="t" or i =="g" or i == "c":
			Presence = True
	return Presence


def build_final_file(File, DicoSeqsDeNovoWithIntrons, ListeOrthogroupNotToTake):
    FinalFile = open("InfoFile_DeNovo_Enabling_WithIntrons.txt", "w")
    FinalFile.write("Orthogroup,NameDeNovo,NamePop,Exons,Chrom,Start,End,PopsToSearchIn"+"\n")
    PopsPresentes = []
    RefSeq = ""
    for i in File[1:]:
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            if len(PopsPresentes)>0:
                if len(RefSeq)>0 and len(PopsPresentes)<7:
                    ListeMissingPops = retrieve_missing_lines(PopsPresentes)
                    FinalFile.write(RefSeq)
                    for j in ListeMissingPops:
                        FinalFile.write(","+j)
                    FinalFile.write("\n")
            PopsPresentes = []
            RefSeq = ""
	    PresenceIntronsInitial = True
	    SizeSeqPrec = 10000000000
            Orthogroup = ligne[1:len(ligne)-1]
        else:
            Data = ligne.split(",")
            NameDeNovo = Data[0]
	    #print Data
            NamePop = NameDeNovo.split("_")[0]
            if NamePop not in PopsPresentes:
                PopsPresentes.append(NamePop)
            Chrom = Data[2]
            Exons = Data[1]
            Start = Data[5]
            End = Data[6]
	    Sequence = DicoSeqsDeNovoWithIntrons[NameDeNovo]
	    PresenceIntrons = search_introns(Sequence)
	    SizeSeq = len(Sequence)
            if Orthogroup not in ListeOrthogroupNotToTake:
		if RefSeq=="":
                	RefSeq = Orthogroup+","+NameDeNovo+","+NamePop+","+Exons+","+Chrom+","+Start+","+End
		else:
			if PresenceIntrons == False:
				RefSeq = Orthogroup+","+NameDeNovo+","+NamePop+","+Exons+","+Chrom+","+Start+","+End
				PresenceIntronsInitial = False
			else:
				if PresenceIntronsInitial == True:
					if SizeSeq<SizeSeqPrec:
						SizeSeqPrec = SizeSeq
						RefSeq = Orthogroup+","+NameDeNovo+","+NamePop+","+Exons+","+Chrom+","+Start+","+End
    if len(RefSeq)>0 and len(PopsPresentes)<7:
        ListeMissingPops = retrieve_missing_lines(PopsPresentes)
        FinalFile.write(RefSeq)
        for j in ListeMissingPops:
            FinalFile.write(","+j)
    FinalFile.close()
        

def build_list_presents(File):
	Liste = []
	for i in File[1:]:
		OrthogroupName = i.split(",")[0]
		Liste.append(OrthogroupName)
	print len(Liste)
	return Liste


def create_dic_seq(DataAK5, DataDK5, DataGI5, DataSW5, DataUM, DataYE, DataZamb):
    Dico = {}
    for i in DataAK5:
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            Name = ligne[1:]
        else:
            Seq = ligne
            Dico[Name] = Seq
    for i in DataDK5:
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            Name = ligne[1:]
        else:
            Seq = ligne
            Dico[Name] = Seq
    for i in DataGI5:
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            Name = ligne[1:]
        else:
            Seq = ligne
            Dico[Name] = Seq
    for i in DataSW5:
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            Name = ligne[1:]
        else:
            Seq = ligne
            Dico[Name] = Seq
    for i in DataUM:
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            Name = ligne[1:]
        else:
            Seq = ligne
            Dico[Name] = Seq
    for i in DataYE:
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            Name = ligne[1:]
        else:
            Seq = ligne
            Dico[Name] = Seq
    for i in DataZamb:
        ligne = i.split("\n")[0]
        if ligne[0] == ">":
            Name = ligne[1:]
        else:
            Seq = ligne
            Dico[Name] = Seq
    return Dico


def main_function():
    DataAK5 = open_file("FI_DeNovoGene_Intron_Exon.fa")
    DataDK5 = open_file("DK_DeNovoGene_Intron_Exon.fa")
    DataGI5 = open_file("ES_DeNovoGene_Intron_Exon.fa")
    DataSW5 = open_file("SE_DeNovoGene_Intron_Exon.fa")
    DataUM = open_file("UA_DeNovoGene_Intron_Exon.fa")
    DataYE = open_file("TR_DeNovoGene_Intron_Exon.fa")
    DataZamb = open_file("ZI_DeNovoGene_Intron_Exon.fa")
    DicoAllSeqsWithIntrons = create_dic_seq(DataAK5, DataDK5, DataGI5, DataSW5, DataUM, DataYE, DataZamb)
    InfoFileWithoutIntron = open_file("InfoFile_DeNovo_Enabling_WithoutIntrons.txt")
    ListeOrthogroupsWithoutIntrons = build_list_presents(InfoFileWithoutIntron)
    Orthogroups = open_file("Orthogroups_info.txt")
    build_final_file(Orthogroups, DicoAllSeqsWithIntrons, ListeOrthogroupsWithoutIntrons)






