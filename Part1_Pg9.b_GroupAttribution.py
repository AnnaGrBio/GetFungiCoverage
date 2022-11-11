import cv2
import os
os.chdir("path to folder")
import random


""" docstring

    This programs search the position in genome of each proto-gene, and attribute a genomic position
    
    input
    ------------
    - bed file of de novo genes
    
    output
    ------------
    - genomic position of de novo genes
"""

def open_file(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  


def work_on_exonic_category(NewName, Signe):
    # this function search in details the position of proto-genes that overlap to an existing exon
    F = open(NewName, "r")
    L = F.readlines()
    CoordBegTransc = int(L[0].split()[1])
    CoordEndTranscr = int(L[0].split()[2])
    ListeAllPosBegin = []
    ListeAllPosEnd = []
    for i in L:
        SigneOfOfficialGene = i.split()[10]
        if SigneOfOfficialGene == Signe:
            RefPosBeg = int(i.split()[6])
            RefPosEnd = int(i.split()[7])
            ListeAllPosBegin.append(RefPosBeg)
            ListeAllPosEnd.append(RefPosEnd)
    BeginEarlier = True
    for i in ListeAllPosBegin:
        if i<CoordBegTransc:
            BeginEarlier = False
            break
    EndLater = True
    for i in ListeAllPosEnd:
        if i>CoordEndTranscr:
            EndLater = False
            break
    if BeginEarlier == True or EndLater == True:
        return "ExonLonger"
    else:
        return "ExonInside"
        
        
def work_category_with_point_direction(NewName, Signe):
    # This function work with the direction "." which later will be changed in "+"
    F = open(NewName, "r")
    L = F.readlines()
    CoordBegTransc = int(L[0].split()[1])
    CoordEndTranscr = int(L[0].split()[2])
    ListeAllPosBegin = []
    ListeAllPosEnd = []
    for i in L:
        SigneOfOfficialGene = i.split()[10]
        RefPosBeg = int(i.split()[6])
        RefPosEnd = int(i.split()[7])
        ListeAllPosBegin.append(RefPosBeg)
        ListeAllPosEnd.append(RefPosEnd)
    BeginEarlier = True
    for i in ListeAllPosBegin:
        if i<CoordBegTransc:
            BeginEarlier = False
            break
    EndLater = True
    for i in ListeAllPosEnd:
        if i>CoordEndTranscr:
            EndLater = False
            break
    if BeginEarlier == True or EndLater == True:

        return "UndirectionalExonLonger"
    else:
        return "UndirectionalExonInside"


def define_position(NewName,Signe):
    # THis function attribute a position for all categories but exonic
    F = open(NewName, "r")
    L = F.readlines()
    FinalCategory = ""
    if len(L) == 0:
        FinalCategory = "Intergenic"
    else:
        Elements = []
        TheSignesAreTheSame = False
        for i in L:
            SigneOfOfficialGene = i.split()[10]
            if SigneOfOfficialGene == Signe:
                TheSignesAreTheSame = True
                TypeOfEltRecovered = i.split()[12]
                if TypeOfEltRecovered not in Elements:
                    Elements.append(TypeOfEltRecovered)
        if TheSignesAreTheSame == False and Signe!=".":
            FinalCategory = "ReverseGenic"
        elif TheSignesAreTheSame == False and Signe==".":
            FinalCategory = work_category_with_point_direction(NewName, Signe)
        else:
            if "pseudogene" in Elements: 
                FinalCategory = "Pseudogene"
            elif "ncRNA_gene" in Elements:
                FinalCategory = "NcRNA"
            elif "gene" in Elements and "CDS" not in Elements:
                FinalCategory = "Intronic"
            elif "gene" in Elements and "CDS" in Elements:
                FinalCategory = work_on_exonic_category(NewName, Signe)
            elif "gene" not in Elements and "CDS" in Elements:
                FinalCategory = work_on_exonic_category(NewName, Signe)
            else:
                FinalCategory = "Weird"
    return FinalCategory


def Main(File, NameFinalFile, Reference):
    # This function creates the main output file
    c = 0
    DicoFinal = {}
    for i in File:
        Signe = i.split("	")[3]
        NameORF = i.split("	")[4].split("\n")[0]
        Intermediar = open("Intermediar", "w")
        Intermediar.write(i)
        Intermediar.close()
        print NameORF
        NewName = "Test"   #+str(c)
        os.system('bedtools intersect -wa -wb -a Intermediar -b ' +Reference+' -sorted  > ' +NewName)
        Category = define_position(NewName, Signe)
        DicoFinal[NameORF] = Category
        c+=1
    F = open(NameFinalFile, "w")
    for i in DicoFinal.keys():
        F.write(i+","+DicoFinal[i]+"\n")
    F.close()


def main_function():
    Data = open_file("FIdenovo_bed")
    Main(Data, "FI_PositionsInGenome", "FIgenome.bed")
