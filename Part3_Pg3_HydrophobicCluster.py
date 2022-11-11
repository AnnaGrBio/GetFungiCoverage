import cv2
import os
os.chdir("path to folder")
from matplotlib import patches
from os import listdir
from os.path import isfile, join


""" docstring

    This program retrieve individuals results from seg HCA and format the output to make stats with R

    input
    ------------

    output
    ------------
    - file with formated HCA per proto-genes and orthogroups

"""


def build_list_HCA():
    # This function retrieve all HCA results and store them in a list
    mypath = "path to results of seghca"
    ListeFile = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    return ListeFile


def open_file(NameFile):
    F = open(NameFile, "r")
    L = F.readlines()
    return L


def get_dico_hca(ListeFile):
    # This funtion associate the set of HCA result to their de novo gene
    DicoData = {}
    c = 0
    for i in ListeFile:
        if ".domain" in i or ".amas" in i or "_profil17.res" in i:
            ligne = i.split("\n")
            if len(ligne) == 1:
                print ligne
            else:
                NameDeNovo = ligne[0][1:]
                FileCategory = ligne[1]
                if FileCategory == ".domain":
                    #print NameDeNovo
                    Lala = open_file(i)
                    ListeInter = []
                    for j in Lala:
                        if j[0] != "#":
                            L2 = j.split("\n")
                            ListeInter.append(L2[0])
                    #print ListeInter
                    DicoData[NameDeNovo] = ListeInter
                    #c +=1
                    #if c>10:
                        #break
    return DicoData


def build_two_file_rsts(Dico):
    # This funtion format the results in 2 files
    F = open("ZnbHCAdomainsPerDeNovo", "w")
    F.write("ORFname"+","+"NbHCAdomains"+","+"NbLongHCA"+","+"NbShortHCA"+"\n")
    for i in Dico.keys():
        Name = i
        NbDomains = len(Dico[i])
        NbSmall = 0
        NbBig = 0
        if len(Dico[i])>0:
            for j in Dico[i]:
                ligne = j.split("-")
                Beg = int(ligne[0])
                End = int(ligne[1])
                Taille = End-Beg
                if Taille <50:
                    NbSmall+=1
                else:
                    NbBig +=1
        F.write(Name+","+str(NbDomains)+","+str(NbBig)+","+str(NbSmall)+"\n")
    F.close()
    F = open("ZsizeHCAdomains", "w")
    F.write("ORFname"+","+"HCAdomainsSize"+"\n")
    for i in Dico.keys():
        Name = i
        if len(Dico[i])>0:
            for j in Dico[i]:
                ligne = j.split("-")
                Beg = int(ligne[0])
                End = int(ligne[1])
                Taille = End-Beg
                F.write(Name+","+str(Taille)+"\n")
    F.close()
                    

def main_function():
    liste_HCA = build_list_HCA
    Dico = get_dico_hca(liste_HCA)
    build_two_file_rsts(Dico)


