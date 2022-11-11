import os
os.chdir("path to folder")
import random


""" docstring

    This program Determine the density of established genes in genomes on chromosomic segments of 160.000 bp, plus calculate tge Ds average value of the genes in the considered genomic regions

    input
    ------------
    - GFF files of the established genes in lines
    - Ds values of establishd genes per lines
    
    output
    ------------  
    - Density of established genes in line s chromosomic segments and DS value of the associated region
    
"""


ChromSizeFI = {"X_Chromosome":22379092, "Y_Chromosome":1150320, "mitochondrion_Chromosome":28136, "2L_Chromosome":23724664, "2R_Chromosome":24122957, "3L_Chromosome":28474363, "3R_Chromosome":30734454, "4_Chromosome":1351024}
ChromSizeDK = {"X_Chromosome":23873564, "Y_Chromosome":1316362, "mitochondrion_Chromosome":33222, "2L_Chromosome":26865308, "2R_Chromosome":24018429, "3L_Chromosome":21699023, "3R_Chromosome":29845977, "4_Chromosome":1305808}
ChromSizeES = {"X_Chromosome":24198538, "Y_Chromosome":2215418, "2L_Chromosome":23416247, "2R_Chromosome":23996592, "3L_Chromosome":28025302, "3R_Chromosome":32392276, "4_Chromosome":1366154}
ChromSizeSE = {"X_Chromosome":22394049, "Y_Chromosome":2753288, "2L_Chromosome":24092724, "2R_Chromosome":24334582, "3L_Chromosome":27038293, "3R_Chromosome":30544244, "4_Chromosome":248932}
ChromSizeUA = {"X_Chromosome":22476591, "Y_Chromosome":709849, "2L_Chromosome":22785099, "2R_Chromosome":24114283, "3L_Chromosome":27885699, "3R_Chromosome":31625211, "4_Chromosome":1335514}
ChromSizeTR = {"X_Chromosome":23856858, "Y_Chromosome":2111415, "2L_Chromosome":22886371, "2R_Chromosome":24687807, "3L_Chromosome":27928844, "3R_Chromosome":30901240, "4_Chromosome":1331028}
ChromSizeZI = {"X_Chromosome":22531871, "Y_Chromosome":2279389, "2L_Chromosome":22767714, "2R_Chromosome":23394925, "3L_Chromosome":26770529, "3R_Chromosome":29115850, "4_Chromosome":1263532}


def open_file(NameFile):
    F = open(NameFile, "r")
    L = F.readlines()
    return L  


def sort_file(Liste):
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


def make_list_chromosomes(File,DicoDnDs):
    Liste2L = []
    Dico2L = {}
    Liste2R = []
    Dico2R = {}
    Liste3L = []
    Dico3L = {}
    Liste3R = []
    Dico3R = {}
    Liste4 = []
    Dico4 = {}
    ListeX = []
    DicoX = {}
    ListeY = []
    DicoY = {}
    for i in File[2:]:
        ligne = i.split("	")
        if "2L_Chromosome" in ligne:
            if len(ligne)>8:
                if "ID=rna" in ligne[8]:
                    NomGene = ligne[8].split("=")[1]
                    NomGene = NomGene.split(";")[0]
                    Start = int(ligne[3])
                    Liste2L.append(Start)
                    if NomGene in DicoDnDs.keys():
                        Dico2L[Start] = DicoDnDs[NomGene]
        if "2R_Chromosome" in ligne:
            if len(ligne)>8:
                if "ID=rna" in ligne[8]:
                    NomGene = ligne[8].split("=")[1]
                    NomGene = NomGene.split(";")[0]
                    Start = int(ligne[3])
                    Liste2R.append(Start)
                    if NomGene in DicoDnDs.keys():
                        Dico2R[Start] = DicoDnDs[NomGene]
        if "3L_Chromosome" in ligne:
            if len(ligne)>8:
                if "ID=rna" in ligne[8]:
                    NomGene = ligne[8].split("=")[1]
                    NomGene = NomGene.split(";")[0]
                    Start = int(ligne[3])
                    Liste3L.append(Start)
                    if NomGene in DicoDnDs.keys():
                        Dico3L[Start] = DicoDnDs[NomGene]
        if "3R_Chromosome" in ligne:
            if len(ligne)>8:
                if "ID=rna" in ligne[8]:
                    NomGene = ligne[8].split("=")[1]
                    NomGene = NomGene.split(";")[0]
                    Start = int(ligne[3])
                    Liste3R.append(Start)
                    if NomGene in DicoDnDs.keys():
                        Dico3R[Start] = DicoDnDs[NomGene]
        if "4_Chromosome" in ligne:
            if len(ligne)>8:
                if "ID=rna" in ligne[8]:
                    NomGene = ligne[8].split("=")[1]
                    NomGene = NomGene.split(";")[0]
                    Start = int(ligne[3])
                    Liste4.append(Start)
                    if NomGene in DicoDnDs.keys():
                        Dico4[Start] = DicoDnDs[NomGene]
        if "X_Chromosome" in ligne:
            if len(ligne)>8:
                if "ID=rna" in ligne[8]:
                    NomGene = ligne[8].split("=")[1]
                    NomGene = NomGene.split(";")[0]
                    Start = int(ligne[3])
                    ListeX.append(Start)
                    if NomGene in DicoDnDs.keys():
                        DicoX[Start] = DicoDnDs[NomGene]
        if "Y_Chromosome" in ligne:
            if len(ligne)>8:
                if "ID=rna" in ligne[8]:
                    NomGene = ligne[8].split("=")[1]
                    NomGene = NomGene.split(";")[0]
                    Start = int(ligne[3])
                    ListeY.append(Start)
                    if NomGene in DicoDnDs.keys():
                        DicoY[Start] = DicoDnDs[NomGene]
    NewListe2L = sort_file(Liste2L)
    NewListe2R = sort_file(Liste2R)
    NewListe3L = sort_file(Liste3L)
    NewListe3R = sort_file(Liste3R)
    NewListe4 = sort_file(Liste4)
    NewListeX = sort_file(ListeX)
    NewListeY = sort_file(ListeY)
    Dico = {"2L_Chromosome":[NewListe2L,Dico2L],"2R_Chromosome":[NewListe2R,Dico2R],"3L_Chromosome":[NewListe3L,Dico3L],"3R_Chromosome":[NewListe3R,Dico3R],"4_Chromosome":[NewListe4,Dico4],"X_Chromosome":[NewListeX,DicoX],"Y_Chromosome":[NewListeY,DicoY]}
    return Dico
    
            
def build_segments(DicoGene, DicoChromSize, Ecart):
    FinalDicoGene = {}
    for Chrom in DicoGene.keys():
        ListeGene = DicoGene[Chrom][0]
        DicoSelection = DicoGene[Chrom][1]
        SizeChrom = DicoChromSize[Chrom]
        NewDico =  {}
        Start = 0
        End = Ecart
        NbSegment = 1
        while End<SizeChrom:
            NbGene = 0
            SUADs = 0.0
            NbValuesDs = 0.0
            FinalDsSegment = 0.0
            for Gene in ListeGene:
                if Gene >= Start and Gene<=End:
                    NbGene += 1
                    if Gene in DicoSelection.keys():
                        SUADs += DicoSelection[Gene]
                        NbValuesDs += 1
                    ListeGene.remove(Gene)
            if NbValuesDs != 0:
                FinalDsSegment = float(SUADs)/float(NbValuesDs)
            else:
                FinalDsSegment = "NA"
            NewDico[NbSegment] = [NbGene,FinalDsSegment]
            Start += Ecart
            End += Ecart
            NbSegment += 1
        FinalDicoGene[Chrom] = NewDico
    return FinalDicoGene


def get_max_Ds(DicoChrom):
    Max = 0.0
    for i in DicoChrom.keys():
        Liste = DicoChrom[i]
        ValueDs = Liste[1]
        if ValueDs!="NA":
            if ValueDs>Max:
                Max = ValueDs
    return Max


def build_file(PopName, DicoFinalGene):
    F = open("GeneDensity"+PopName, "w")
    F.write("Chromosome"+","+"Segment"+","+"NbGenes"+","+"Ds"+","+"MaxDs"+"\n")
    for i in DicoFinalGene.keys():
        Chrom = i
        DicoChrom = DicoFinalGene[i]
        MaxDs = get_max_Ds(DicoChrom)
        for j in range(1,100000):
            if j in DicoChrom.keys():
                F.write(Chrom+","+str(j)+","+str(DicoChrom[j][0])+","+str(DicoChrom[j][1])+","+str(MaxDs)+"\n")
    F.close()
            
            
def build_dico(DnDs):
    D = {}
    for i in DnDs[1:]:
        ligne = i.split("\n")[0]
        Values = ligne.split(",")
        GeneName = Values[0]
        Ds = Values[2]
        if float(Ds)<float(2):
            D[GeneName] = float(Ds)
    return D
    

def main_function():
	DnDsFI = open_file("z_FI_FinalDnDs.rst")
	DicoDnDsFI = build_dico(DnDsFI)
	Gene_FI = open_file("FIdeNovoGenomeAnnotation.gff")
	DicoGeneFI = make_list_chromosomes(Gene_FI,DicoDnDsFI)   
	DicoFinalGene_FI = build_segments(DicoGeneFI, ChromSizeFI, 160000)
	build_file("FI",DicoFinalGene_FI)
	DnDsDK = open_file("z_DK_FinalDnDs.rst")
	DicoDnDsDK = build_dico(DnDsDK)
	Gene_DK = open_file("DKdeNovoGenomeAnnotation.gff")
	DicoGeneDK = make_list_chromosomes(Gene_DK,DicoDnDsDK)   
	DicoFinalGene_DK = build_segments(DicoGeneDK, ChromSizeDK, 160000)
	build_file("DK",DicoFinalGene_DK)
	DnDsES = open_file("z_ES_FinalDnDs.rst")
	DicoDnDsES = build_dico(DnDsES)
	Gene_ES = open_file("ESdeNovoGenomeAnnotation.gff")
	DicoGeneES = make_list_chromosomes(Gene_ES,DicoDnDsES)   
	DicoFinalGene_ES = build_segments(DicoGeneES, ChromSizeES, 160000)
	build_file("ES",DicoFinalGene_ES)
	DnDsSE = open_file("z_SE_FinalDnDs.rst")
	DicoDnDsSE = build_dico(DnDsSE)
	Gene_SE = open_file("SEdeNovoGenomeAnnotation.gff")
	DicoGeneSE = make_list_chromosomes(Gene_SE,DicoDnDsSE)   
	DicoFinalGene_SE = build_segments(DicoGeneSE, ChromSizeSE, 160000)
	build_file("SE",DicoFinalGene_SE)
	DnDsUA = open_file("z_UA_FinalDnDs.rst")
	DicoDnDsUA = build_dico(DnDsUA)
	Gene_UA = open_file("UAdeNovoGenomeAnnotation.gff")
	DicoGeneUA = make_list_chromosomes(Gene_UA,DicoDnDsUA)   
	DicoFinalGene_UA = build_segments(DicoGeneUA, ChromSizeUA, 160000)
	build_file("UA",DicoFinalGene_UA)
	DnDsTR = open_file("z_TR_FinalDnDs.rst")
	DicoDnDsTR = build_dico(DnDsTR)
	Gene_TR = open_file("TRdeNovoGenomeAnnotation.gff")
	DicoGeneTR = make_list_chromosomes(Gene_TR,DicoDnDsTR)   
	DicoFinalGene_TR = build_segments(DicoGeneTR, ChromSizeTR, 160000)
	build_file("TR",DicoFinalGene_TR)
	DnDsZI = open_file("z_ZI_FinalDnDs.rst")
	DicoDnDsZI = build_dico(DnDsZI)
	Gene_ZI = open_file("ZIdeNovoGenomeAnnotation.gff")
	DicoGeneZI = make_list_chromosomes(Gene_ZI,DicoDnDsZI)   
	DicoFinalGene_ZI = build_segments(DicoGeneZI, ChromSizeZI, 160000)
	build_file("ZI",DicoFinalGene_ZI)


main_function()










