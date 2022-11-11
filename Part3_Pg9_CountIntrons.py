import os
os.chdir("path to folder")
import random
from collections import Counter


""" docstring

    This program count the average number of proto-genes with intron in the lines

    input
    ------------
    - info file of the lines
    
    output
    ------------  
    - average number of proto-genes with intron in the lines
    
"""


def open_file(NameFile):
    F = open(NameFile, "r")
    L = F.readlines()
    return L  


def count_intron(File, pop):
    nb_gene_with_intron = 0
    total_nb_gene = 0
    for line in File:
        elmts_in_line = line.split(",")
        name_pop = elmts_in_line[12].split("_")[0]
        if name_pop == pop:
            total_nb_gene += 1
            nb_exon = int(elmts_in_line[3])
            if nb_exon>1:
                nb_gene_with_intron+=1
    print ("Pop "+pop+" : "+str(nb_gene_with_intron)+" genes with introns ")
    print ("Pop "+pop+" : "+str(total_nb_gene)+" genes in total ")
    print ("Pop "+pop+" : "+str(nb_gene_with_intron*100/total_nb_gene)+" percent of genes with introns ")
    print ("")


def main_function():
InfoFile = open_file("DeNovo_InfoFile")
	count_intron(InfoFile, "AK5")
	count_intron(InfoFile, "DK5")
	count_intron(InfoFile, "GI5")
	count_intron(InfoFile, "SW5")
	count_intron(InfoFile, "UM")
	count_intron(InfoFile, "YE")
	count_intron(InfoFile, "Zamb")
