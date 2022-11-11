from bisect import bisect_left, bisect_right
import pandas as pd
from Bio import SeqIO


""" docstring

    This program search for synteny between genomes. For each proto-genes, the surrounding established genes are retrieved if presents in the target genome.
    The syntenic region between homologous established genes is retrieved

    input
    ------------
    - File from genome parsing
    - ref seq per orthogroups

    output
    ------------
    - csv file with syntenic regions and infos

"""

# Path to de novo candidate fasta file
de_novo_dict = {rec.id: rec.seq for rec in SeqIO.parse(
    r"path to folder\OrthogroupsRefSeq.fa", "fasta")}
info_dict = {
    'ID': '',
    'Seq': '',
    'Chrom': '',
    'NamePop': '',
    'Start': '',
    'End': '',
    'Orthogroup': '',
    'SearchPop': ''
}
info_list = []
info_list_neigh = []
query_dict_list = []
ID = list(de_novo_dict.keys())
Seq = list(de_novo_dict.values())


# Parses the de novo info and search population file into a list of dictionaries
def query_parser():
    # Path to de novo info file
    info_df = pd.read_csv(r"path to folder\InfoFile.csv")
    # Path to search population csv file
    pop_search_df = pd.read_csv(
        r"path to folder\Popstosearchin.csv", header=None)
    print("Getting info from fasta header...")

    for n, i in enumerate(ID):
        row = info_df.iloc[n]
        info_dict['ID'] = ID[n]
        info_dict['Seq'] = str(Seq[n])
        info_dict['Chrom'] = row["Chrom"]
        info_dict['NamePop'] = row["NamePop"]
        info_dict['Start'] = row["Start"]
        info_dict['End'] = row["End"]
        info_dict['Orthogroup'] = row["Orthogroup"]
        info_dict['Start'] = row["Start"]

        row_search = pop_search_df.iloc[n]
        search_str = row_search.to_string(index=False)
        # Converts search_str into a list with each population as separate entry
        search_list = search_str.splitlines()
        search_list = [x.strip('  ') for x in search_list]
        search_list = [x for x in search_list if str(x) != 'NaN']
        info_dict['SearchPop'] = str(search_list)
        info_list.append(info_dict.copy())
    print("Done")


# Finds surrounding gene in query genome based on positional data in infofile and genome annotation csv.
# Neighboring_degree argument default is set to 1 to find the next surrounding gene.
def neighboring_gene(neighboring_degree=1, query=None):
    if query is None:
        query = info_list
    if type(query) == dict:
        pop = query.get('NamePop')
        chrom = query.get('Chrom')
        start = query.get('Start')
        end = query.get('End')

        # Path to genome annotation directory
        annotated = pd.read_csv(
            r"path to folder\{}_genes.csv".format(pop), delimiter=",")
        chrom_filter = annotated[annotated["Chrom"] == str(chrom)]
        array_stop = pd.array(chrom_filter["Stop"])
        array_start = pd.array(chrom_filter["Start"])
        # If start position is greater end position 3'-5' orientation is assumed, switches start and end information
        if start >= end:
            start = query.get('End')
            end = query.get('Start')
        # Neighboring gene is determined by bisection of the arrays
        try:
            lower = array_stop[bisect_left(array_stop, start) - neighboring_degree]
            above = array_start[bisect_right(array_start, end) + (neighboring_degree - 1)]

            lower_neigh = chrom_filter.loc[chrom_filter["Stop"] == lower]
            above_neigh = chrom_filter.loc[chrom_filter["Start"] == above]

            # No gene left, RoiStart = Chromosome start
            if lower > above:
                query['GeneLeft'] = 'NaN'

            query['GeneRight'] = above_neigh['Gene'].to_string(index=False)
            query['GeneLeft'] = lower_neigh['Gene'].to_string(index=False)

            # update info_list_neigh if entry found
            if query in info_list_neigh:
                idx = info_list_neigh.index(query)
                info_list_neigh.remove(query)
                info_list_neigh.insert(idx, query.copy())

        except IndexError:
            lower_neigh = chrom_filter.loc[chrom_filter["Stop"] == lower]
            query['GeneLeft'] = lower_neigh['Gene'].to_string(index=False)
            query['GeneRight'] = 'NaN'
            # update info_list_neigh if entry found
            if query in info_list_neigh:
                idx = info_list_neigh.index(query)
                info_list_neigh.remove(query)
                info_list_neigh.insert(idx, query.copy())
        return query

    else:
        for dictionary in info_list:
            pop = dictionary.get('NamePop')
            chrom = dictionary.get('Chrom')
            start = dictionary.get('Start')
            end = dictionary.get('End')
            annotated = pd.read_csv(
                r"path to folder/{}_genes.csv".format(pop),
                delimiter=",")
            chrom_filter = annotated[annotated["Chrom"] == str(chrom)]
            array_stop = pd.array(chrom_filter["Stop"])
            array_start = pd.array(chrom_filter["Start"])

            # DeNovo candiate is 3'-5', swap start and end to get Roi
            if start >= end:
                start = dictionary.get('End')
                end = dictionary.get('Start')

            try:
                lower = array_stop[bisect_left(array_stop, start) - neighboring_degree]
                above = array_start[bisect_right(array_start, end) + (neighboring_degree - 1)]

                lower_neigh = chrom_filter.loc[chrom_filter["Stop"] == lower]
                above_neigh = chrom_filter.loc[chrom_filter["Start"] == above]

                # No gene left, RoiStart = Chromosome start
                if lower > above:
                    dictionary['GeneLeft'] = 'NaN'

                dictionary['GeneRight'] = above_neigh['Gene'].to_string(index=False)
                dictionary['GeneLeft'] = lower_neigh['Gene'].to_string(index=False)
                info_list_neigh.append(dictionary.copy())

            # No gene right, RoiEnd should be chromosome end
            except IndexError:
                lower_neigh = chrom_filter.loc[chrom_filter["Stop"] == lower]
                dictionary['GeneLeft'] = lower_neigh['Gene'].to_string(index=False)
                dictionary['GeneRight'] = 'NaN'
                info_list_neigh.append(dictionary.copy())
                continue


# Finds the position of the previously found surrounding gene by name in the annotation of the search population
# and stores it in a list of dictionaries
def synteny_locate():
    for count, dictionary in enumerate(info_list_neigh):
        # For each population in SearchPop, finds gene in corresponding pop. file and gets positional information
        for pop in dictionary.get('SearchPop'):
            # Path to genome annotation csv files
            pop_df = pd.read_csv(
                r"path to folder/{}_genes.csv".format(pop),
                delimiter=",")
            roi_start = pop_df.loc[pop_df["Gene"] == dictionary['GeneLeft']]['Start'].to_string(index=False)
            if dictionary['GeneRight'] == "NaN":
                roi_end = 'NaN'
            else:
                roi_end = pop_df.loc[pop_df["Gene"] == dictionary['GeneRight']]['Stop'].to_string(index=False)
            try:
                # If no gene with this name found in genome, search for second neighbor
                if roi_start == 'Series([], )':
                    n = 2
                    while roi_start == 'Series([], )':
                        dictionary = neighboring_gene(neighboring_degree=n, query=dictionary)
                        roi_start = pop_df.loc[pop_df["Gene"] == dictionary['GeneLeft']]['Start'].to_string(index=False)
                        n += 1
                    if roi_end == 'Series([], )':
                        n = 2
                        while roi_end == 'Series([], )':
                            dictionary = neighboring_gene(neighboring_degree=n, query=dictionary)
                            roi_end = pop_df.loc[pop_df["Gene"] == dictionary['GeneRight']]['Stop'].to_string(
                                index=False)
                            n += 1
                elif roi_end == 'Series([], )':
                    n = 2
                    while roi_end == 'Series([], )':
                        dictionary = neighboring_gene(neighboring_degree=n, query=dictionary)
                        roi_end = pop_df.loc[pop_df["Gene"] == dictionary['GeneRight']]['Stop'].to_string(index=False)
                        n += 1
                    if roi_start == 'Series([], )':
                        n = 2
                        while roi_start == 'Series([], )':
                            dictionary = neighboring_gene(neighboring_degree=n, query=dictionary)
                            roi_start = pop_df.loc[pop_df["Gene"] == dictionary['GeneLeft']]['Start'].to_string(
                                index=False)
                            n += 1

                dictionary[str(pop) + '_Start'] = roi_start.split('\n')[0]
                dictionary[str(pop) + '_End'] = roi_end.split('\n')[-1]

            except:
                dictionary['ERROR'] = 'ERROR'

                continue
        query_dict_list.append(dictionary.copy())
        print("No.{} of {} {}%".format(count + 1, len(info_list),
                                       round(float((count + 1) / int(len(info_list))) * int(100)), 2))
        print(dictionary.items())
        print("")


query_parser()
neighboring_gene()
synteny_locate()
# Converts list of dicts in dataframe and saves dataframe as csv file.
df = pd.DataFrame(query_dict_list)
df.to_csv(r'path to folder\20211210_query_dict_list.csv', index=False)

