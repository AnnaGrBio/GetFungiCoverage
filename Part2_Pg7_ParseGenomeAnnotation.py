import gffutils
import json
import pandas as pd
import os


""" docstring

    This functions parse and store informations in genome annotation for search of syntenw betzeen genomes

    input
    ------------
    - gtf file of de novo genomes annotations

    output
    ------------   
    - cvs file with genes informations in genomes
"""


# Creates SQL database based on .gff file
def db_create(gff_path, output_path):
    output_path = output_path[:-4] + '.db'
    gffutils.create_db(gff_path, output_path)
    return output_path


# Reads annotation SQL database
def db_read(db_path):
    db = gffutils.FeatureDB(db_path, keep_order=True)
    return db


# Parses SQL database to csv file with ID, Chromosome, Start, Stop and Gene name
def db_parse(database, csv_output_path):
    gene_dict = {"ID": [], "Chrom": [], "start": [], "stop": [], "Ref-gene": []}
    gen_list = []
    ref_list = []
    gene_features = database.features_of_type("mRNA")
    genes = list(gene_features)
    # Grabs gene on position n and ID, Chrom, Position and Name
    for n, gene in enumerate(genes):
        gene = genes[n]
        id = gene.astuple()[0]
        chrom = gene.astuple()[1]
        start_pos = gene.astuple()[4]
        stop_pos = gene.astuple()[5]

        # Modifies ref_gene to fit str format
        ref_gene = str(json.loads(gene.astuple()[9])["ref-gene"])
        for char in "'[]":
            ref_gene = ref_gene.replace(char, "")

        # Creates list of dictionaries for each ID
        gen_list.append(
            {"ID": id, "Chrom": chrom, "Start": int(start_pos), "Stop": int(stop_pos), "Ref-Gene": str(ref_gene)})
        ref_list.append(ref_gene)

    # Create pandas Dataframe and adds ref_gen column
    df = pd.DataFrame(gen_list, columns=["ID", "Chrom", "Start", "Stop", "Gene"])
    df['Gene'] = ref_list
    df.to_csv(csv_output_path, index=False)
    return df


# Path to .gff annotation directory
path = r'path to folder'
def main_function():
    for filename in os.listdir(path):
        gff_file = os.path.join(path, filename)
        print("Database of file {} will be constructed...".format(filename))
        db_path = db_create(gff_file,
                            r'path to folder\{}'.format(filename))
        print("Database constructed at: ", db_path)
        print("Reading Database")
        db = db_read(db_path)
        print("Database read, start parsing to csv file")
        db_parse(db, r'path to folder\CSV\{}'.format(filename[:-4] + '.csv'))
        print("Database parsed to .csv file")
    print("Done")

