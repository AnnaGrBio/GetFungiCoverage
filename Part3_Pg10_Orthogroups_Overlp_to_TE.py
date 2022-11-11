import os
os.chdir("path to folder")
import random
from collections import Counter


""" docstring

    This program calculate the proportion of proto-genes that are found in several copies in one single line and are located in a TE or overlapping to it

    input
    ------------
    - orthogroups information file
    - TE overlap statis of each proto-gene
    
    output
    ------------  
    - proportion of proto-genes that are found in several copies in one single line and are located in a TE or overlapping to it
    
"""


def open_file(NameFile):
    F = open(NameFile, "r")
    L = F.readlines()
    return L  


def build_list_orthogroups(File):
    liste_all_orthogroups = []
    name_orthogroup = ""
    liste_proto_gene = []
    for ligne in File:
        if ligne[0] == ">":
            if name_orthogroup == "":
                name_orthogroup = ligne
            else:
                liste_all_orthogroups.append(liste_proto_gene)
                liste_proto_gene = []
        else:
            name_proto_gene = ligne.split(",")[0]
            liste_proto_gene.append(name_proto_gene)
    return liste_all_orthogroups


def build_dic_proto_gene_status(my_file):
    dico = {}
    for ligne in my_file[1:]:
        liste_elements = ligne.split("\n")[0].split(",")
        dico[liste_elements[0]] = liste_elements[1]
    return dico


def calculate_percentage(liste_orthogroups,dico_status,name_pop_query,name_file):
    f = open(name_file,"w")
    nb_orthologs = 0
    nb_orthologs_inside_TE = 0
    nb_orthologs_overlap_TE = 0
    for orthogroups in liste_orthogroups:
        if len(orthogroups)>1:
            liste_orthologs = []
            for gene in orthogroups:
                name_pop_gene = gene.split("_")[0]
                if name_pop_gene == name_pop_query:
                    liste_orthologs.append(gene)
            if len(liste_orthologs)>1:
                nb_orthologs+=len(liste_orthologs)
                for orthologs in liste_orthologs:
                    status = dico_status[orthologs]
                    if status == "InsideTE":
                        nb_orthologs_inside_TE+=1
                    elif status == "OverlapTE":
                        nb_orthologs_overlap_TE+=1
                    f.write(orthologs+","+status+"\n")
                f.write("\n"+"\n")
    print (name_pop_query)
    print ("percentage orthologs inside TE : "+str(100*nb_orthologs_inside_TE/nb_orthologs))
    print ("percentage orthologs overlap TE : "+str(100*nb_orthologs_overlap_TE/nb_orthologs))
    print ("")
    f.close()
                

def main_function():
	file_orthogroups = open_file("Orthogroups_info.txt")
	liste_orthogroups = build_list_orthogroups(file_orthogroups)
	file_proto_gene_status = open_file("Name_proto_gene_status.txt")
	dico_status = build_dic_proto_gene_status(file_proto_gene_status)
	calculate_percentage(liste_orthogroups,dico_status,"FI","AK5_ortholog_TE_status.txt")
	calculate_percentage(liste_orthogroups,dico_status,"DK","DK5_ortholog_TE_status.txt")
	calculate_percentage(liste_orthogroups,dico_status,"ES","GI5_ortholog_TE_status.txt")
	calculate_percentage(liste_orthogroups,dico_status,"SE","SW5_ortholog_TE_status.txt")
	calculate_percentage(liste_orthogroups,dico_status,"UA","UM_ortholog_TE_status.txt")
	calculate_percentage(liste_orthogroups,dico_status,"TR","YE_ortholog_TE_status.txt")
	calculate_percentage(liste_orthogroups,dico_status,"ZI","Zamb_ortholog_TE_status.txt")








