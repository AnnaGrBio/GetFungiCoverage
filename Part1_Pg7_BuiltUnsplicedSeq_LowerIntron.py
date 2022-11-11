import cv2
import os
os.chdir("path to folder")
import random


""" docstring

	This program builds a new file containing the de novo ORFs, in which the introns have been lowered
	
	input
	------------
	- fasta genome
	- gtf file of the transcriptome assembly
	- info file with position of unspliced sequences
	
	output
	------------
	fasta file of de novo ORFs with introns lowered in the sequence
"""


def openFile(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  


########################### Make Dico Chrom
def build_dic_chromosomes(File):
	# This function builds a dictionary with all chromosomes fasta seqs inside
	Dico = {}
	Seq = ""
	for i in File:
		ligne = i.split("\n")[0]
		if ligne[0] == ">":
			if len(Seq)>0:
				Dico[NameChrom] = Seq
				Seq = ""
			NameChrom = ligne[1:]
		else:
			Seq += ligne
	Dico[NameChrom] = Seq
	return Dico


##############" Get ExonIntron Strcture
def build_dic_exon_intron_structure(Liste):
	# This function builds a dictionary with exon intron according to the position
	D = {}
	for i in range(len(Liste)):
		liste_exon = Liste[i]
		for j in range(liste_exon[0],liste_exon[1]+1):
			D[j] = "E"
		if i!=len(Liste)-1:
			for k in range(liste_exon[1]+1,Liste[i+1][0]):
				D[k]="I"
	return D


def get_structure(File):
	# This function determines exon intron position and indeed structure according to the GTF informations
	Dico = {}
	dico_chrom_transcript = {}
	for i in File[2:]:
		ligne = i.split()
		Chrom = ligne[0]
		name_transcript = (ligne[11][1:len(ligne[11])-2])
		dico_chrom_transcript[name_transcript] = Chrom
		Element = ligne[2]
		Start = int(ligne[3])
		End = int(ligne[4])
		CoordExon = [Start, End]
		if Element == "exon":
			if name_transcript in Dico.keys():
				Dico[name_transcript].append(CoordExon)
			else:
				Dico[name_transcript] = [CoordExon]
	DicoFinal = {}
	for transcript in Dico:
		CoordExons = Dico[transcript]
		dico_intermediate = {}
		dico_intermediate = build_dic_exon_intron_structure(CoordExons)
		DicoFinal[transcript] = dico_intermediate
	return DicoFinal, dico_chrom_transcript


################### Retrieve Pos unspliced ORFs
def get_unspliced_position(File, dico_chrom_transcript):
	# This function retrieves unspliced position of transcripts and ORFs
	D = {}
	for i in File[1:]:
		ligne = i.split(",")
		DeNovo = ligne[0]
		transcript_name = DeNovo.split("_")[1]
		Chrom = dico_chrom_transcript[transcript_name]
		Start = int(ligne[1])
		End = int(ligne[2])
		D[DeNovo] = [transcript_name, Chrom, Start, End]
	return D


### Make Final de novo
def reverse_seq(Seq):
	# THis function reverses the input sequence
	D = {"A":"T", 
	  "T":"A", 
	  "G":"C", 
	  "C":"G", 
	  "N":"N", 
	  "n":"n", 
	  "a":"t", 
	  "t":"a", 
	  "g":"c", 
	  "c":"g"}
	NewSeq1 = ''.join(reverse_seqd(Seq))
	NewSeq2 = ""
	for i in NewSeq1:
		New = D[i]
		NewSeq2+=New
	return NewSeq2


def final(DicoChrom,dic_structure_transcript, dico_spliced_ORF):
	# this function creates a final dictionary in wich all de novo genes have their exonic sequence in upper letter and intronic in lower letter
	c = 0
	final_dico_de_novo_gene = {}
	for DeNovo in dico_spliced_ORF.keys():
		transcript_name = dico_spliced_ORF[DeNovo][0]
		Chrom = dico_spliced_ORF[DeNovo][1]
		Start = dico_spliced_ORF[DeNovo][2]
		End = dico_spliced_ORF[DeNovo][3]
		DicoTranscriptSplicing = dic_structure_transcript[transcript_name]
		print DeNovo
		#print Start
		#print End
		if End>Start:
			de_novo_sequence = ""
			for i in range(Start,End+1):
				Nuc = DicoChrom[Chrom][i]
				if i in DicoTranscriptSplicing.keys():
					if DicoTranscriptSplicing[i] == "E":
						FinalNuc = Nuc.upper()
					else:
						FinalNuc = Nuc.lower()
				else:
					FinalNuc = Nuc.upper()
				de_novo_sequence+=FinalNuc
			final_dico_de_novo_gene[DeNovo] = de_novo_sequence
		else:
                        de_novo_sequence = ""
                        for i in range(End,Start+1):
                                Nuc = DicoChrom[Chrom][i]
                                if i in DicoTranscriptSplicing.keys():
                                        if DicoTranscriptSplicing[i] == "E":
                                                FinalNuc = Nuc.upper()
                                        else:
                                                FinalNuc = Nuc.lower()
                                else:
                                        FinalNuc = Nuc.upper()
                                de_novo_sequence+=FinalNuc
			Finalde_novo_sequence = reverse_seq(de_novo_sequence)
                        final_dico_de_novo_gene[DeNovo] = Finalde_novo_sequence
		c+=1
		print c
	return final_dico_de_novo_gene


def build_final_file(Dico, Name):
	# this function builds the final file with de novo ORF with introns lowered in the sequence
	F = open(Name+"_"+"DeNovoGene_Intron_Exon.fa", "w")
	for i in Dico.keys():
		F.write(">"+i+"\n")
		F.write(Dico[i]+"\n")
	F.close()


def main_function():
	fasta_genome = openFile("FIfinalGenome.masked.fa")
	DicoChrom = build_dic_chromosomes(fasta_genome)
	GTF = openFile("FItranscriptome.gtf")
	dic_structure_transcript, dico_chrom_transcript = get_structure(GTF)
	pos_unspliced_ORF = openFile("FI_UnsplicedORFpositions")
	dico_spliced_ORF = get_unspliced_position(pos_unspliced_ORF, dico_chrom_transcript)
	final_dic_de_novo = final(DicoChrom,dic_structure_transcript, dico_spliced_ORF)
	build_final_file(final_dic_de_novo, "FI")
