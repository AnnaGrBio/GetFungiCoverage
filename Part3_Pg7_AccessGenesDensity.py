import os
os.chdir("path to folder")
import random


""" docstring

    This program Determine the density of established genes in genomes on chromosomic segments of 160.000 bp

    input
    ------------
    - GFF files of the established genes in lines
    
    output
    ------------  
    - Density of established genes in line s chromosomic segments 
    
"""


chrom_size_FI = {"X_Chromosome":22379092, "Y_Chromosome":1150320, "mitochondrion_Chromosome":28136, "2L_Chromosome":23724664, "2R_Chromosome":24122957, "3L_Chromosome":28474363, "3R_Chromosome":30734454, "4_Chromosome":1351024}
chrom_size_DK = {"X_Chromosome":23873564, "Y_Chromosome":1316362, "mitochondrion_Chromosome":33222, "2L_Chromosome":26865308, "2R_Chromosome":24018429, "3L_Chromosome":21699023, "3R_Chromosome":29845977, "4_Chromosome":1305808}
chrom_size_ES = {"X_Chromosome":24198538, "Y_Chromosome":2215418, "2L_Chromosome":23416247, "2R_Chromosome":23996592, "3L_Chromosome":28025302, "3R_Chromosome":32392276, "4_Chromosome":1366154}
chrom_size_SE = {"X_Chromosome":22394049, "Y_Chromosome":2753288, "2L_Chromosome":24092724, "2R_Chromosome":24334582, "3L_Chromosome":27038293, "3R_Chromosome":30544244, "4_Chromosome":248932}
chrom_size_UA = {"X_Chromosome":22476591, "Y_Chromosome":709849, "2L_Chromosome":22785099, "2R_Chromosome":24114283, "3L_Chromosome":27885699, "3R_Chromosome":31625211, "4_Chromosome":1335514}
chrom_size_TR = {"X_Chromosome":23856858, "Y_Chromosome":2111415, "2L_Chromosome":22886371, "2R_Chromosome":24687807, "3L_Chromosome":27928844, "3R_Chromosome":30901240, "4_Chromosome":1331028}
chrom_size_ZI = {"X_Chromosome":22531871, "Y_Chromosome":2279389, "2L_Chromosome":22767714, "2R_Chromosome":23394925, "3L_Chromosome":26770529, "3R_Chromosome":29115850, "4_Chromosome":1263532}


def open_file(NameFile):
    F = open(NameFile, "r")
    L = F.readlines()
    return L  


def sort_liste(Liste):
    NewListe = []
    for i in Liste:
        Place = False
        for j in range(len(NewListe)):
            Other = NewListe[j]
            if Other >= i:
                NewListe.insert(j,i)
                Place = True
                break
        if Place == False:
            NewListe.append(i)
    return NewListe


def make_list_chromosomes(File):
    Liste2L = []
    Liste2R = []
    Liste3L = []
    Liste3R = []
    Liste4 = []
    ListeX = []
    ListeY = []
    for i in File[2:]:
        ligne = i.split("	")
        if "2L_Chromosome" in ligne:
            if len(ligne)>8:
                if "ID=gene" in ligne[8]:
                    Start = int(ligne[3])
                    Liste2L.append(Start)    
        if "2R_Chromosome" in ligne:
            if len(ligne)>8:
                if "ID=gene" in ligne[8]:
                    Start = int(ligne[3])
                    Liste2R.append(Start)
        if "3L_Chromosome" in ligne:
            if len(ligne)>8:
                if "ID=gene" in ligne[8]:
                    Start = int(ligne[3])
                    Liste3L.append(Start)
        if "3R_Chromosome" in ligne:
            if len(ligne)>8:
                if "ID=gene" in ligne[8]:
                    Start = int(ligne[3])
                    Liste3R.append(Start)
        if "4_Chromosome" in ligne:
            if len(ligne)>8:
                if "ID=gene" in ligne[8]:
                    Start = int(ligne[3])
                    Liste4.append(Start)
        if "X_Chromosome" in ligne:
            if len(ligne)>8:
                if "ID=gene" in ligne[8]:
                    Start = int(ligne[3])
                    ListeX.append(Start)
        if "Y_Chromosome" in ligne:
            if len(ligne)>8:
                if "ID=gene" in ligne[8]:
                    Start = int(ligne[3])
                    ListeY.append(Start)
    NewListe2L = sort_liste(Liste2L)
    NewListe2R = sort_liste(Liste2R)
    NewListe3L = sort_liste(Liste3L)
    NewListe3R = sort_liste(Liste3R)
    NewListe4 = sort_liste(Liste4)
    NewListeX = sort_liste(ListeX)
    NewListeY = sort_liste(ListeY)
    Dico = {"2L_Chromosome":NewListe2L,"2R_Chromosome":NewListe2R,"3L_Chromosome":NewListe3L,"3R_Chromosome":NewListe3R,"4_Chromosome":NewListe4,"X_Chromosome":NewListeX,"Y_Chromosome":NewListeY}
    return Dico
    
            
def determine_genomic_segment(DicoGene, DicoChromSize, Ecart):
    FinalDicoGene = {}
    for Chrom in DicoGene.keys():
        ListeGene = DicoGene[Chrom]
        SizeChrom = DicoChromSize[Chrom]
        NewDico =  {}
        Start = 0
        End = Ecart
        NbSegment = 1
        while End<SizeChrom:
            NbGene = 0
            for Gene in ListeGene:
                if Gene >= Start and Gene<=End:
                    NbGene += 1
                    ListeGene.remove(Gene)
            NewDico[NbSegment] = NbGene
            Start += Ecart
            End += Ecart
            NbSegment += 1
        FinalDicoGene[Chrom] = NewDico
    return FinalDicoGene


def make_final_file(PopName, DicoFinalGene):
    F = open("GeneDensity"+PopName, "w")
    F.write("Chromosome" + "," + "Segment" + "," + "NbGenes" + "\n")
    for i in DicoFinalGene:
        Chrom = i
        DicoChrom = DicoFinalGene[i]
        for j in range(1,100000):
            if j in DicoChrom.keys():
                F.write(Chrom+"," + str(j) + "," + str(DicoChrom[j]) + "\n")
    F.close()
            

def main_function():
	Gene_FI = open_file("FIdeNovoGenomeAnnotation.gff")
	DicoGeneFI = make_list_chromosomes(Gene_FI)   
	DicoFinalGene_FI = determine_genomic_segment(DicoGeneFI, chrom_size_FI, 160000)
	make_final_file("FI",DicoFinalGene_FI)
	Gene_DK = open_file("DKdeNovoGenomeAnnotation.gff")
	DicoGeneDK = make_list_chromosomes(Gene_DK)   
	DicoFinalGene_DK = determine_genomic_segment(DicoGeneDK, chrom_size_DK, 160000)
	make_final_file("DK",DicoFinalGene_DK)
	Gene_ES = open_file("ESdeNovoGenomeAnnotation.gff")
	DicoGeneES = make_list_chromosomes(Gene_ES)   
	DicoFinalGene_ES = determine_genomic_segment(DicoGeneES, chrom_size_ES, 160000)
	make_final_file("ES",DicoFinalGene_ES)
	Gene_SE = open_file("SEdeNovoGenomeAnnotation.gff")
	DicoGeneSE = make_list_chromosomes(Gene_SE)   
	DicoFinalGene_SE = determine_genomic_segment(DicoGeneSE, chrom_size_SE, 160000)
	make_final_file("SE",DicoFinalGene_SE)
	Gene_UA = open_file("UAdeNovoGenomeAnnotation.gff")
	DicoGeneUA = make_list_chromosomes(Gene_UA)   
	DicoFinalGene_UA = determine_genomic_segment(DicoGeneUA, chrom_size_UA, 160000)
	make_final_file("UA",DicoFinalGene_UA)
	Gene_TR = open_file("TRdeNovoGenomeAnnotation.gff")
	DicoGeneTR = make_list_chromosomes(Gene_TR)   
	DicoFinalGene_TR = determine_genomic_segment(DicoGeneTR, chrom_size_TR, 160000)
	make_final_file("TR",DicoFinalGene_TR)
	Gene_ZI = open_file("ZIdeNovoGenomeAnnotation.gff")
	DicoGeneZI = make_list_chromosomes(Gene_ZI)   
	DicoFinalGene_ZI = determine_genomic_segment(DicoGeneZI, chrom_size_ZI, 160000)
	make_final_file("ZI",DicoFinalGene_ZI)


main_function()













