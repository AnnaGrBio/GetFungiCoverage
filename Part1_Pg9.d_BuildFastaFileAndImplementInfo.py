import cv2
import os
os.chdir("path to folder")
import random

def openFile(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L  



def MakeListe(Names):
    L = []
    for i in Names:
        ligne = i.split(",")
        #print ligne
        Name = ligne[0]
        Cat = ligne[1][0:len(ligne[1])-1]
        L.append(Name)
    #print len(L)
    return L



def MakeFinalFile(Lgoods, FASTA, NameFile):
    F = open(NameFile, "w")
    for i in FASTA:
        if i[0] == ">":
            Name = i[1:len(i)-1]
            LongName = i
        else:
            Seq = i
            if Name in Lgoods:
                F.write(LongName)
                F.write(Seq)
    F.close()



def MakeDicoPos(File):
    Dico = {}
    for i in Names:
        ligne = i.split("\n")[0]
        Name = ligne.split(",")[0]
        Pos = ligne.split(",")[1]
        Dico[Name] = Pos
    print len(Dico)
    return Dico


def RedoInfoFile(Dico, File, NewName):
    Liste = []
    Title = File[0].split("\n")[0]
    for i in File[1:]:
        ligne = i.split("\n")[0]
        lala = ligne.split(",")
        if lala[12] in Dico.keys():
            ligne+=","+Dico[lala[12]]+"\n"
        else:
            ligne+=","+"-"+"\n"
        Liste.append(ligne)
    F = open(NewName, "w")
    F.write(Title+",GenomicPosition"+"\n")
    for i in Liste:
        F.write(i)
    F.close()
    
        


Names = openFile("FI_PositionsInGenome")
Nuc = openFile("FI_DeNovoORFwithStop_DNA.fa")
Prot = openFile("FI_DeNovoORFwithStop_Protein.fa")

Lgoods = MakeListe(Names)
MakeFinalFile(Lgoods, Nuc, "FI_DeNovoGenesNuc.fa")
MakeFinalFile(Lgoods, Prot, "FI_DeNovoGenesProt.fa")


DicoPos = MakeDicoPos(Names)
Info = openFile("FI_InformationFile")
RedoInfoFile(DicoPos, Info, "FI_Final_InformationFile")
