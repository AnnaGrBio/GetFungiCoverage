#import cv2
import os
os.chdir("path to folder")
import random
import math

""" docstring

    This program align the established genes that are found in all populations to their ortholog in the referent genome, and calculate the Dn Ds ratio between each gene and the orthologous gene of the reference genome

    input
    ------------
    - FASTA file of all established genes and the reference gene
    
    outpur
    ------------  
    - Average Dn Ds ration of the genes compared to their reference sequence in genomic regions 
    
"""


DNA_Codons = {
    # 'M' - START, '_' - STOP
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATA": "I", "ATT": "I", "ATC": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "_", "TAG": "_", "TGA": "_"
}


# GetAllFile Names
def open_file(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  


def build_list_names(File):
    L = []
    for i in File:
        Name = i.split("\n")[0]
        L.append(Name)
    #print (L)
    return L


def convert_proteins(Seq):
    Seq = Seq.split("\n")[0]
    FinalProt = ""
    Trio = ""
    #ompteur = 0
    for Nuc in Seq:
        if len(Trio)!=3:
            Trio += Nuc
        else:
            Trio = Trio.upper()
            if Trio not in DNA_Codons.keys():
                Aa = "N"
            else:
                Aa = DNA_Codons[Trio]
            FinalProt+=Aa
            #print (Trio)
            #print (Aa)
            Trio = Nuc
            #compteur+=1
            #if compteur>5:
                #break
    FinalProt = FinalProt+"_"+"\n"
    return FinalProt


def build_sub_file(File):
    ListSequences = open_file(File)
    NomRef = ""
    SeqRef = ""
    NomPop = ""
    SeqPop = ""
    IDPop = ""
    for i in ListSequences:
        if i[0] == ">":
            ligne = i.split("\n")[0]
            Nom = ligne.split("__")[1]
            #print (Nom)
            if Nom == "REF":
                NomRef = i
            else:
                NomPop = i
                IDPop = Nom 
        else:
            if SeqRef == "" and NomRef!= "":
                SeqRef = i
            else:
                SeqPop = i
                NameFile = IDPop+"_Nuc.fa"
                F = open(NameFile, "w")
                F.write(NomRef)
                F.write(SeqRef)
                F.write(NomPop)
                F.write(SeqPop)
                F.close()
                
                SeqRefProt = convert_proteins(SeqRef)
                SeqPopProt = convert_proteins(SeqPop)
                
                NameFile = IDPop+"_Prot.fa"
                F = open(NameFile, "w")
                F.write(NomRef)
                F.write(SeqRefProt)
                F.write(NomPop)
                F.write(SeqPopProt)
                F.close()
                
                
# Get Dn Ds 
#DicoCodeGenetique = {"TTT": "Phe", "TTC":"Phe", "TTA":"Leu", "TTG":"Leu", "CTT":"Leu", "CTC":"Leu", "CTA":"Leu", "CTG":"Leu", "ATT":"Ile", "ATC":"Ile", "ATA":"Ile", "ATG":"Met", "GTT":"Val", "GTC":"Val", "GTA":"Val", "GTG":"Val", 
def define_trios(TrioRef,TrioPop):
    NSitesyn = 0
    NSitenonsyn = 0
    NSubsyn = 0
    NSubnonsyn = 0
    AaRef = DNA_Codons[TrioRef]
    ListeNuc = ["A", "T", "G", "C"]
    for i in range(len(TrioRef)):
        NucRef = TrioRef[i]
        NucPop = TrioPop[i]
        NSitesynInter = 0
        NSitenonsynInter = 0
        NSubsynInter = 0
        NSubnonsynInter = 0
        for NewNuc in ListeNuc:
            if NewNuc!= NucRef:
                if i == 0:
                    putativeCodon = NewNuc+TrioRef[1:]
                elif i == 1:
                    putativeCodon = TrioRef[0]+NewNuc+TrioRef[2]
                else:
                    putativeCodon = TrioRef[0]+TrioRef[1]+NewNuc
                PutativeAa = DNA_Codons[putativeCodon]
                if PutativeAa == AaRef:
                    NSitesynInter+=1
                else:
                    NSitenonsynInter+=1
        if NucPop!=NucRef:
            if i == 0:
                putativeCodon = NucPop+TrioRef[1:]
            elif i == 1:
                putativeCodon = TrioRef[0]+NucPop+TrioRef[2]
            else:
                putativeCodon = TrioRef[0]+TrioRef[1]+NucPop
            PutativeAa = DNA_Codons[putativeCodon]
            if PutativeAa == AaRef:
                NSubsynInter+=1
            else:
                NSubnonsynInter+=1
        NSitesyn+=(float(NSitesynInter)/float(3))
        #print (NSitesyn)
        NSitenonsyn+=(float(NSitenonsynInter)/float(3))
        #print (NSitenonsyn)
        NSubsyn+=(float(NSubsynInter))
        #print (NSubsyn)
        NSubnonsyn+=(float(NSubnonsynInter))
        #print (NSubnonsyn)
        #print ("************")
    return NSitesyn,NSitenonsyn,NSubsyn,NSubnonsyn
                    
                
def calculate_dn_ds(Sitesyn, Sitenonsyn, Subsyn, Subnonsyn):
    if float(Sitenonsyn)!=0 or float(Sitenonsyn)!=-0:
        Pn = float(Subnonsyn)/float(Sitenonsyn)
    else:
        Pn = 100000
    if float(Sitesyn)!=0 or float(Sitesyn)!=-0:
        Ps = float(Subsyn)/float(Sitesyn)
    else:
        Ps = 100000
    if float(1)-((float(4)*Pn)/float(3))>0:
        Dn = -(float(3)/float(4))*math.log(float(1)-((float(4)*Pn)/float(3)))
    else:
        Dn = 100000
    if float(1)-((float(4)*Ps)/float(3))>0:
        Ds = -(float(3)/float(4))*math.log(float(1)-((float(4)*Ps)/float(3)))
    else:
        Ds = 100000
    return Dn, Ds


def calculate_dn_ds_internal(SeqRef,SeqPop):
    TrioRef = ""
    TrioPop = ""
    FinalNSitesyn = 0
    FinalNSitenonsyn = 0
    FinalNSubsyn = 0
    FinalNSubnonsyn = 0
    for i in range(len(SeqRef)):
        NucRef = SeqRef[i]
        NucPop = SeqPop[i]
        if len(TrioRef)!=3:
            TrioRef += NucRef
            TrioPop += NucPop
        else:
            TrioRef = TrioRef.upper()
            TrioPop = TrioPop.upper()
            if "-" not in TrioRef and "-" not in TrioPop and "N" not in TrioRef and "N" not in TrioPop and "n" not in TrioRef and "n" not in TrioPop:
                NSitesyn,NSitenonsyn,NSubsyn,NSubnonsyn = define_trios(TrioRef,TrioPop)
                FinalNSitesyn += NSitesyn
                FinalNSitenonsyn += NSitenonsyn
                FinalNSubsyn += NSubsyn
                FinalNSubnonsyn += NSubnonsyn
                #break
            TrioRef = NucRef
            TrioPop = NucPop
    #print (FinalNSitesyn)
    #print (FinalNSitenonsyn)
    #print (FinalNSubsyn)
    #print (FinalNSubnonsyn)
    Dn, Ds = calculate_dn_ds(FinalNSitesyn, FinalNSitenonsyn, FinalNSubsyn, FinalNSubnonsyn)
    return Dn, Ds


def main_get_dn_ds(Ali):
    NomSeqRef = ""
    NomSeqPop = ""
    NUASeq = 0
    for i in Ali:
        ligne = i.split("\n")[0]
        if ligne[0]==">":
            NUASeq+=1
        else:
            if NUASeq == 1:
                NomSeqRef+=ligne
            else:
                NomSeqPop+=ligne
    #print (len(NomSeqRef))
    #print (len(NomSeqPop))
    
    Dn, Ds = calculate_dn_ds_internal(NomSeqRef,NomSeqPop)
    ListeDnDs = [Dn,Ds]
    return ListeDnDs


# Make Alignments
def make_sequence_alignments(Liste):
    DicoFI = {}
    DicoDK = {}
    DicoES = {}
    DicoSE = {}
    DicoUA = {}
    DicoTR = {}
    DicoZI = {}
    for Names in Liste:
        build_sub_file(Names)
        Command = "/global/projects/programs/bin/mafft "+"FI_Prot.fa"+" > Ali_FI_Prot.fa"
        os.system(Command)
        os.system('perl /global/projects/programs/bin/pal2nal.pl Ali_FI_Prot.fa  FI_Nuc.fa  -output fasta > Ali_FI_NucPal2nal.fa')
        Ali = open_file("Ali_FI_NucPal2nal.fa")
        ListeDnDsFI = main_get_dn_ds(Ali)
        DicoFI[Names]= ListeDnDsFI
        os.system('rm Ali_FI_Prot.fa Ali_FI_NucPal2nal.fa FI_Nuc.fa FI_Prot.fa')
        
        Command = "/global/projects/programs/bin/mafft "+"DK_Prot.fa"+" > Ali_DK_Prot.fa"
        os.system(Command)
        os.system('perl /global/projects/programs/bin/pal2nal.pl Ali_DK_Prot.fa  DK_Nuc.fa  -output fasta > Ali_DK_NucPal2nal.fa')
        Ali = open_file("Ali_DK_NucPal2nal.fa")
        ListeDnDsDK = main_get_dn_ds(Ali)
        DicoDK[Names]= ListeDnDsDK
        os.system('rm Ali_DK_Prot.fa Ali_DK_NucPal2nal.fa DK_Nuc.fa DK_Prot.fa')
        
        Command = "/global/projects/programs/bin/mafft "+"ES_Prot.fa"+" > Ali_ES_Prot.fa"
        os.system(Command)
        os.system('perl /global/projects/programs/bin/pal2nal.pl Ali_ES_Prot.fa  ES_Nuc.fa  -output fasta > Ali_ES_NucPal2nal.fa')
        Ali = open_file("Ali_ES_NucPal2nal.fa")
        ListeDnDsES = main_get_dn_ds(Ali)
        DicoES[Names]= ListeDnDsES
        os.system('rm Ali_ES_Prot.fa Ali_ES_NucPal2nal.fa ES_Nuc.fa ES_Prot.fa')
        
        Command = "/global/projects/programs/bin/mafft "+"SE_Prot.fa"+" > Ali_SE_Prot.fa"
        os.system(Command)
        os.system('perl /global/projects/programs/bin/pal2nal.pl Ali_SE_Prot.fa  SE_Nuc.fa  -output fasta > Ali_SE_NucPal2nal.fa')
        Ali = open_file("Ali_SE_NucPal2nal.fa")
        ListeDnDsSE = main_get_dn_ds(Ali)
        DicoSE[Names]= ListeDnDsSE
        os.system('rm Ali_SE_Prot.fa Ali_SE_NucPal2nal.fa SE_Nuc.fa SE_Prot.fa')
        
        Command = "/global/projects/programs/bin/mafft "+"UA_Prot.fa"+" > Ali_UA_Prot.fa"
        os.system(Command)
        os.system('perl /global/projects/programs/bin/pal2nal.pl Ali_UA_Prot.fa  UA_Nuc.fa  -output fasta > Ali_UA_NucPal2nal.fa')
        Ali = open_file("Ali_UA_NucPal2nal.fa")
        ListeDnDsUA = main_get_dn_ds(Ali)
        DicoUA[Names]= ListeDnDsUA
        os.system('rm Ali_UA_Prot.fa Ali_UA_NucPal2nal.fa UA_Nuc.fa UA_Prot.fa')
        
        Command = "/global/projects/programs/bin/mafft "+"TR_Prot.fa"+" > Ali_TR_Prot.fa"
        os.system(Command)
        os.system('perl /global/projects/programs/bin/pal2nal.pl Ali_TR_Prot.fa  TR_Nuc.fa  -output fasta > Ali_TR_NucPal2nal.fa')
        Ali = open_file("Ali_TR_NucPal2nal.fa")
        ListeDnDsTR = main_get_dn_ds(Ali)
        DicoTR[Names]= ListeDnDsTR
        os.system('rm Ali_TR_Prot.fa Ali_TR_NucPal2nal.fa TR_Nuc.fa TR_Prot.fa')
        
        Command = "/global/projects/programs/bin/mafft "+"ZI_Prot.fa"+" > Ali_ZI_Prot.fa"
        os.system(Command)
        os.system('perl /global/projects/programs/bin/pal2nal.pl Ali_ZI_Prot.fa  ZI_Nuc.fa  -output fasta > Ali_ZI_NucPal2nal.fa')
        Ali = open_file("Ali_ZI_NucPal2nal.fa")
        ListeDnDsZI = main_get_dn_ds(Ali)
        DicoZI[Names]= ListeDnDsZI
        os.system('rm Ali_ZI_Prot.fa Ali_ZI_NucPal2nal.fa ZI_Nuc.fa ZI_Prot.fa')

    FFI = open("z_FI_FinalDnDs.rst", "w")
    FFI.write("GeneName"+","+"Dn"+","+"Ds"+"\n")
    for Genes in DicoFI.keys():
        NameGeneDecomposed = Genes.split(".")
        NameGeneDecomposed = NameGeneDecomposed[0]+"."+NameGeneDecomposed[1]
        NameGeneDecomposed = NameGeneDecomposed.split("_")
        NameGene = NameGeneDecomposed[1]+"_"+NameGeneDecomposed[2]+"_"+NameGeneDecomposed[3]
        FFI.write(NameGene+","+str(DicoFI[Genes][0])+","+str(DicoFI[Genes][1])+"\n")
    FFI.close()
    
    FDK = open("z_DK_FinalDnDs.rst", "w")
    FDK.write("GeneName"+","+"Dn"+","+"Ds"+"\n")
    for Genes in DicoDK.keys():
        NameGeneDecomposed = Genes.split(".")
        NameGeneDecomposed = NameGeneDecomposed[0]+"."+NameGeneDecomposed[1]
        NameGeneDecomposed = NameGeneDecomposed.split("_")
        NameGene = NameGeneDecomposed[1]+"_"+NameGeneDecomposed[2]+"_"+NameGeneDecomposed[3]
        FDK.write(NameGene+","+str(DicoDK[Genes][0])+","+str(DicoDK[Genes][1])+"\n")
    FDK.close()
    
    FES = open("z_ES_FinalDnDs.rst", "w")
    FES.write("GeneName"+","+"Dn"+","+"Ds"+"\n")
    for Genes in DicoES.keys():
        NameGeneDecomposed = Genes.split(".")
        NameGeneDecomposed = NameGeneDecomposed[0]+"."+NameGeneDecomposed[1]
        NameGeneDecomposed = NameGeneDecomposed.split("_")
        NameGene = NameGeneDecomposed[1]+"_"+NameGeneDecomposed[2]+"_"+NameGeneDecomposed[3]
        FES.write(NameGene+","+str(DicoES[Genes][0])+","+str(DicoES[Genes][1])+"\n")
    FES.close()
    
    FSE = open("z_SE_FinalDnDs.rst", "w")
    FSE.write("GeneName"+","+"Dn"+","+"Ds"+"\n")
    for Genes in DicoSE.keys():
        NameGeneDecomposed = Genes.split(".")
        NameGeneDecomposed = NameGeneDecomposed[0]+"."+NameGeneDecomposed[1]
        NameGeneDecomposed = NameGeneDecomposed.split("_")
        NameGene = NameGeneDecomposed[1]+"_"+NameGeneDecomposed[2]+"_"+NameGeneDecomposed[3]
        FSE.write(NameGene+","+str(DicoSE[Genes][0])+","+str(DicoSE[Genes][1])+"\n")
    FSE.close()
    
    FUA = open("z_UA_FinalDnDs.rst", "w")
    FUA.write("GeneName"+","+"Dn"+","+"Ds"+"\n")
    for Genes in DicoUA.keys():
        NameGeneDecomposed = Genes.split(".")
        NameGeneDecomposed = NameGeneDecomposed[0]+"."+NameGeneDecomposed[1]
        NameGeneDecomposed = NameGeneDecomposed.split("_")
        NameGene = NameGeneDecomposed[1]+"_"+NameGeneDecomposed[2]+"_"+NameGeneDecomposed[3]
        FUA.write(NameGene+","+str(DicoUA[Genes][0])+","+str(DicoUA[Genes][1])+"\n")
    FUA.close()
    
    FTR = open("z_TR_FinalDnDs.rst", "w")
    FTR.write("GeneName"+","+"Dn"+","+"Ds"+"\n")
    for Genes in DicoTR.keys():
        NameGeneDecomposed = Genes.split(".")
        NameGeneDecomposed = NameGeneDecomposed[0]+"."+NameGeneDecomposed[1]
        NameGeneDecomposed = NameGeneDecomposed.split("_")
        NameGene = NameGeneDecomposed[1]+"_"+NameGeneDecomposed[2]+"_"+NameGeneDecomposed[3]
        FTR.write(NameGene+","+str(DicoTR[Genes][0])+","+str(DicoTR[Genes][1])+"\n")
    FTR.close()
    
    FZI = open("z_ZI_FinalDnDs.rst", "w")
    FZI.write("GeneName"+","+"Dn"+","+"Ds"+"\n")
    for Genes in DicoZI.keys():
        NameGeneDecomposed = Genes.split(".")
        NameGeneDecomposed = NameGeneDecomposed[0]+"."+NameGeneDecomposed[1]
        NameGeneDecomposed = NameGeneDecomposed.split("_")
        NameGene = NameGeneDecomposed[1]+"_"+NameGeneDecomposed[2]+"_"+NameGeneDecomposed[3]
        FZI.write(NameGene+","+str(DicoZI[Genes][0])+","+str(DicoZI[Genes][1])+"\n")
    FZI.close()
    

def main_function():       
    FileName = open_file("AllFileNames.txt")
    ListeNames = build_list_names(FileName)
    make_sequence_alignments(ListeNames)

