import subprocess
import pandas as pd
from Bio import SeqIO
import os


""" docstring

    This program BLAST proto-genes against a target region or a target genome, and store the results

    input
    ------------
    - query list CVS from program 8
    - proto-genes

    output
    ------------
    - BLAST results sequences
"""


# Gets synteny sequence based on position provided through synteny_locate().
def target_region(query_id, search_population_list):
    # Path to BLAST default database dir
    target_regions = r"c:\PROGRA~1\NCBI\blast-2.12.0+\db\{}_region".format(query_id)
    target_output = open(target_regions, "a")
    for population in search_population_list:
        print("Finding target region for: {} in {} ".format(query_id, population))
        start = int(query.get(population + '_Start'))
        end = int(query.get(population + '_End'))
        if population == "FI":
            seq_dict = FI_seq_dict
        elif population == "DK":
            seq_dict = DK_seq_dict
        elif population == "SE":
            seq_dict = SE_seq_dict
        elif population == "ES":
            seq_dict = ES_seq_dict
        elif population == "UA":
            seq_dict = UA_seq_dict
        elif population == "TR":
            seq_dict = TR_seq_dict
        elif population == "ZI":
            seq_dict = ZI_seq_dict
        else:
            print(population)
            continue
        seq_chr = list(seq_dict.get(chrom))
        # If left surrounding gene is after right gene located
        if start > end:
            target_list = seq_chr[end:start]
        else:
            target_list = seq_chr[start:end]
        target_seq = "".join(target_list)
        print(">{}_{}_{}_region".format(population, start, end), file=target_output)
        print(target_seq, file=target_output)
    target_output.close()
    print("Done")
    return target_regions


# blast function. Query is direction to a temporary fasta file. Stores blast results as .xml file in output_dir
def blast(query_id, query_path, db_path, output_dir):
    global result
    if 'finalGenome' in db_path:
        result = str(output_dir) + "{}.xml".format("\\" + str(query_id) + '_' + pop + '_Genome')
    elif '_region' in db_path:
        result = str(output_dir) + "{}.xml".format("\\" + str(query_id) + '_region')
    cmd = ['blastn', '-db', db_path, '-query', query_path, '-outfmt', '5', '-out', result]
    print("Start blasting:  ", query_id)
    subprocess.call(cmd, shell=True)
    print("Done")
    return result, db_path


# Constructs blast database based on provided fasta file
def blast_db_construct(blast_db_path):
    subprocess.call("makeblastdb -in " + blast_db_path + " -dbtype nucl", shell=True)


# Path to constructed csv file containing the synteny position information
df = pd.read_csv(r'path to folder\20211210_query_dict_list.csv')
queries = df.to_dict('records')

# Checks if temp fasta file is already constructed
try:
    fasta_temp = r"path to folder\blast_query.fa"
    # creates blast fasta file
    with open(fasta_temp, 'x') as f:
        pass
except FileExistsError:
    print('Temporary BLAST fasta file already exists, starting BLAST operation')
    print("")


# Opens each population genome sequence as dictionary
print("Open population sequences..")
FI_seq_dict = {rec.id: rec.seq for rec in
                SeqIO.parse(r"path to folder\FIfinalGenome.masked.fa",
                            "fasta")}
DK_seq_dict = {rec.id: rec.seq for rec in
                SeqIO.parse(r"path to folder\DKfinalGenome.masked.fa",
                            "fasta")}
SE_seq_dict = {rec.id: rec.seq for rec in
                SeqIO.parse(r"path to folder\SEfinalGenome.masked.fa",
                            "fasta")}
ES_seq_dict = {rec.id: rec.seq for rec in
                SeqIO.parse(r"path to folder\ESfinalGenome.masked.fa",
                            "fasta")}
UA_seq_dict = {rec.id: rec.seq for rec in
               SeqIO.parse(r"path to folder\UAfinalGenome.masked.fa",
                           "fasta")}
TR_seq_dict = {rec.id: rec.seq for rec in
               SeqIO.parse(r"path to folder\TRfinalGenome.masked.fa",
                           "fasta")}
ZI_seq_dict = {rec.id: rec.seq for rec in SeqIO.parse(
    r"path to folder\ZIfinalGenome.masked.fa", "fasta")}
print("Done")
print("")
# Constructing BLAST database for all genomic sequences in genome_seq_path
genome_seq_path = r'path to folder\Fasta'
for filename in os.listdir(genome_seq_path):
    blast_db_construct(genome_seq_path+'\\'+filename)

total = len(queries)  # Total nUAber of queries
for n, query in enUAerate(queries):
    print("Start Query No. {} of {} {} %".format(n + 1, total, round(float(n / total) * 100, 2)))
    ID = query.get('ID')
    seq = query.get('Seq')
    chrom = query.get('Chrom')
    try:
        search_pops = query.get('SearchPop')
        search_pops = search_pops.strip('[', )
        search_pops = search_pops.strip(']', )
        search_pops = search_pops.replace("', '", ",")
        search_pops = search_pops.replace("'", "")
        search_pops = search_pops.split(",")
    except AttributeError:
        continue
    # Writes fasta query to temp blast fasta file
    fasta_query = """>{}_{}\n{}""".format(ID, chrom, seq)
    fasta = open(r"path to folder\blast_query.fa", "w")
    fasta.write(fasta_query)
    fasta.close()

    target_dir = target_region(query_id=ID, search_population_list=search_pops)
    blast_db_construct(target_dir)

    # Blast against whole genome of search populations
    for pop in search_pops:
        result, db = blast(query_id=ID,
                           query_path=r"path to folder\blast_query.fa",
                           db_path=pop + "finalGenome.masked.fa",
                           output_dir=r"path to folder\211210_blast_results")

    # Blast against located syntenic target regions of search populations
    result, db = blast(query_id=ID,
                       query_path=r"path to folder\blast_query.fa",
                       db_path=ID + "_region",
                       output_dir=r"path to folder\211210_blast_results")
    print("No more files in directory found, finished!")

