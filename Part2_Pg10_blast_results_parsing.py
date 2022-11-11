import pandas as pd
from Bio import SeqIO, SearchIO
import re
import os
import Levenshtein

""" docstring
    
    This program work on BLAST result, align the proto-gene to non-coding homolog, and remove putative introns from the alignment

    input
    ------------
    - BLAST file csv from previous steps
    - referent proto-gene per orthogroup

    output
    ------------
    - aligned non coding homolog without intron
    - information about BLAST, alignments and mutations
"""


# Locates intronic region in blast alignments based on raw denovo sequence letter case. (lower letters = intron)
def intron_locater(denovo_seq, query_seq, hit_seq):
    region_dict = {}  # Key: Denovo str index, Value: E for exons or I for introns
    index_denovo = {}  # Key: Query index, Value: Query index with hyphen shift
    query_seq_constructer = []
    hit_seq_constructer = []
    intron_indicator = 0
    # find start index of query sequence in denovo sequence based on a search sequence of 30 nucleotides long
    denovo_align = re.sub(r'[-]', "", str(query_seq))[:30].lower()
    start = denovo_seq.lower().find(denovo_align)
    if start < 0:
        start = 0

    # Label initial seq with E for exon and I for intronic nucleotides
    for idx, nt in enumerate(denovo_seq):
        if nt.islower():
            region_dict[idx] = 'I'
            intron_indicator = 1
        else:
            region_dict[idx] = "E"

    # integrates offset if hyphen appears
    offset = 0
    for idx, nt in enumerate(query_seq):
        if nt == '-':
            offset += 1
        index_denovo[idx] = idx - offset

    # Enumerates over query sequence and gets the corr. denovo position and region info from dict.
    for idx, nt in enumerate(query_seq):
        denovo_pos = index_denovo.get(idx)
        region = region_dict.get(start + denovo_pos)
        if region == "E":
            query_seq_constructer.append(query_seq[idx])
        elif region == "I":
            lower_nt = query_seq[idx].lower()
            query_seq_constructer.append(lower_nt)

    # Enumerates over hit sequence and gets the corr. denovo position and region info from dict.
    for idx, nt in enumerate(hit_seq):
        denovo_pos = index_denovo.get(idx)
        region = region_dict.get(start + denovo_pos)
        if region == "E":
            hit_seq_constructer.append(hit_seq[idx])
        elif region == "I":
            lower_nt = hit_seq[idx].lower()
            hit_seq_constructer.append(lower_nt)

    # Joins lists to string
    query_seq_lower = "".join(query_seq_constructer)
    hit_seq_lower = "".join(hit_seq_constructer)

    # Removes lower case letters in sequence through regex function
    rm_lower = lambda text: re.sub('[a-z]', '', text)
    denovo_seq_without_intron = rm_lower(denovo_seq)
    query_seq = rm_lower(query_seq_lower)
    hit_seq = rm_lower(hit_seq_lower)

    # returns following information without introns: query sequence, target sequence and original denovo sequence.
    # intron_indicator bool if intron found
    return query_seq, hit_seq, intron_indicator, denovo_seq_without_intron


# Blast parser for result
def blast_parser(result_file_path, denovo_sequence_path, n):
    # Dict for genomic matches
    output_dict = {
        'ID': '',
        'Hit': '',
        'Synteny': '',
        'ATG': '',
        'Stop Codon': '',
        'Triplet': '',
        'Query Size': '',
        'Target Size': '',
        'Target Start Pos.': '',
        'Target Stop Pos.': '',
        'Target Chromosome': '',
        'Population': '',
        'Target Sequence': '',
        'BitScore': '',
        'Insertion': '',
        'Deletion': '',
        'Identical Length': '',
        'evalue': '',
        'Query start': '',
        'Query end': '',
        'Query stop codon': ''
    }
    # Dict for syntenic matches
    output_dict_syn = {
        'ID': '',
        'Hit': '',
        'Synteny': '',
        'ATG': '',
        'Stop Codon': '',
        'Triplet': '',
        'Query Size': '',
        'Target Size': '',
        'Target Start Pos.': '',
        'Target Stop Pos.': '',
        'Target Chromosome': '',
        'Population': '',
        'Target Sequence': '',
        'BitScore': '',
        'Insertion': '',
        'Deletion': '',
        'Identical Length': '',
        'evalue': '',
        'Query start': '',
        'Query end': '',
        'Query stop codon': ''
    }

    if 'Genome' in result_file_path:
        print('Genome Database')
        records = SearchIO.parse(result_file_path, 'blast-xml')
        for record in records:
            if 'mitochondrion' in record.id:
                query_chr = 'mitochondrion_Chromosome'
            else:
                query_chr = re.search('[2,3,4,X,Y][L,R]?_Chromosome', record.id)
                query_chr = query_chr.group()
            # Filters Hits on different Chromosomes not in query name
            if query_chr in record.hit_keys:
                # Constructs new class based on filter without hits on different Chromosomes
                filter_chr = lambda hit: hit.id in hit.query_id
                filtered_chr = record.hit_filter(filter_chr)
                for hits in filtered_chr:
                    for hit in hits:
                        # removes Chrom from query ID
                        if 'mitochondrion' in hit.query.id:
                            query_id = hit.query.id[:-25]
                        else:
                            query_id = re.sub('_[2,3,4,X,Y][L,R]?_Chromosome', "", hit.query.id)

                        # Original query sequence
                        denovo_seq = str(denovo_sequence_path.get(query_id))
                        # Query sequence of blast alignment
                        query_seq_align = hit.query.seq
                        # Target sequence of blast alignment
                        hit_seq_align = hit.hit.seq
                        # calls intron_locater
                        seqs_nointron = intron_locater(denovo_seq=denovo_seq, query_seq=query_seq_align,
                                                       hit_seq=hit_seq_align)

                        # Finds gaps based on hyphen in sequence after removing potential intron
                        gap_query = len(re.findall('[-]', seqs_nointron[0]))
                        gap_hit = len(re.findall('[-]', seqs_nointron[1]))

                        # Removes gaps of sequence after removing potential intron
                        query_seq = re.sub('[-]', "", str(seqs_nointron[0]))
                        hit_seq = re.sub('[-]', "", str(seqs_nointron[1]))

                        # Calculates Levenshtein ratio of original sequence (without intron) to hit_seq
                        perc_ident = Levenshtein.ratio(seqs_nointron[3], hit_seq)

                        # Gets ATG based on hit_seq
                        if hit_seq[0:3] == "ATG":
                            output_dict['ATG'] = 1
                        else:
                            output_dict['ATG'] = 0

                        # Gets Stop codon based on hit_seq
                        if hit_seq[len(hit_seq) - 3:len(hit_seq)] == "TAG":
                            output_dict['Stop Codon'] = 1
                        elif hit_seq[len(hit_seq) - 3:len(hit_seq)] == "TAA":
                            output_dict['Stop Codon'] = 1
                        elif hit_seq[len(hit_seq) - 3:len(hit_seq)] == "TGA":
                            output_dict['Stop Codon'] = 1
                        else:
                            output_dict['Stop Codon'] = 0

                        if len(hit_seq) % 3 == 0:
                            output_dict['Triplet'] = 1
                        else:
                            output_dict['Triplet'] = 0
                        if len(hit_seq) > 0:
                            if len(seqs_nointron[3]) / len(hit_seq) > 0.9:
                                output_dict['Identical Length'] = 1
                        else:
                            output_dict['Identical Length'] = 0
                        if len(denovo_seq) == hit.aln_span:
                            output_dict['Query end'] = 1
                        else:
                            output_dict['Query end'] = 0

                        # Adds output to dictionary
                        output_dict['ID'] = query_id
                        output_dict['Hit'] = 1
                        output_dict['Synteny'] = 0
                        output_dict['Query Size'] = len(seqs_nointron[3])
                        output_dict['Target Size'] = len(hit_seq)
                        output_dict['Intron'] = seqs_nointron[2]
                        output_dict['Target Start Pos.'] = hit.hit_start
                        output_dict['Target Stop Pos.'] = hit.hit_end
                        output_dict['Target Chromosome'] = hit.hit.id
                        if record.target[2] == "f":
                            output_dict['Population'] = record.target[0:2]
                        elif record.target[0] == 'Z':
                            output_dict['Population'] = 'Zamb'
                        else:
                            output_dict['Population'] = record.target[0:3]
                        output_dict['Target Sequence'] = hit_seq
                        output_dict['BitScore'] = hit.bitscore
                        output_dict['Identical'] = round(perc_ident, 4)
                        output_dict['Insertion'] = gap_query
                        output_dict['Deletion'] = gap_hit
                        output_dict['evalue'] = hit.evalue
                        output_dict['Query start'] = hit.query_start + 1
                        output_dict['Query stop codon'] = denovo_seq[len(denovo_seq) - 3:len(denovo_seq)]

                        # Copies dictionary to list
                        output_list.append(output_dict.copy())

            # No hit for Query found on same Chromosome, adds output to dictionary
            else:
                # removes Chrom from query ID
                if 'mitochondrion' in record.id:
                    query_id = record.id[:-25]
                else:
                    query_id = re.sub('_[2,3,4,X,Y][L,R]?_Chromosome', "", record.id)

                # Updates output dict with information for missing matches
                output_dict['ID'] = query_id
                output_dict['Hit'] = 0
                output_dict['Synteny'] = 0
                output_dict['ATG'] = 0
                output_dict['Stop Codon'] = 0
                output_dict['Triplet'] = 0
                output_dict['Query Size'] = record.seq_len
                output_dict['Target Size'] = 'NaN'
                output_dict['Intron'] = 'NaN'
                output_dict['Target Start Pos.'] = 'NaN'
                output_dict['Target Stop Pos.'] = 'NaN'
                output_dict['Target Chromosome'] = 'NaN'
                if record.target[2] == "f":
                    output_dict['Population'] = record.target[0:2]
                elif record.target[0] == 'Z':
                    output_dict['Population'] = 'Zamb'
                else:
                    output_dict['Population'] = record.target[0:3]
                output_dict['Target Sequence'] = 'NaN'
                output_dict['BitScore'] = 'NaN'
                output_dict['Identical'] = 'NaN'
                output_dict['Insertion'] = 'NaN'
                output_dict['Deletion'] = 'NaN'
                output_dict['Identical Length'] = 'NaN'
                output_dict['evalue'] = 'NaN'
                output_dict['Query start'] = 'NaN'
                output_dict['Query end'] = 'NaN'
                output_dict['Query stop codon'] = 'NaN'
                output_list.append(output_dict.copy())

            return n + 1

    # Blast results are based on target region
    else:
        print('Target Database')
        records = SearchIO.parse(result_file_path, 'blast-xml')
        for record in records:
            if 'mitochondrion' in record.id:
                query_chr = 'mitochondrion_Chromosome'
            else:
                query_chr = re.search('[2,3,4,X,Y][L,R]?_Chromosome', record.id)
                query_chr = query_chr.group()
            # Constructs new class based on filter without hits on different Chromosomes
            for hits in record:
                for hit in hits:
                    try:
                        query_id = re.sub('_[2,3,4,X,Y][L,R]?_Chromosome', "", record.id)
                        chrom = re.search('_[2,3,4,X,Y][L,R]?_Chromosome', record.id)

                        denovo_seq = str(denovo_sequence_path.get(query_id))

                        query_seq_align = hit.query.seq
                        hit_seq_align = hit.hit.seq
                        denovo_align = re.sub(r'[-]', "", str(query_seq_align))[:30].lower()
                        denovo_start = denovo_seq.lower().find(denovo_align)
                        if denovo_start < 0:
                            denovo_start = 0

                        seqs_nointron = intron_locater(denovo_seq=denovo_seq, query_seq=query_seq_align,
                                                       hit_seq=hit_seq_align)

                        # Finds gaps based on hypen in sequence after removing potential intron
                        gap_query = len(re.findall('[-]', seqs_nointron[0]))
                        gap_hit = len(re.findall('[-]', seqs_nointron[1]))

                        # Removes gaps of sequence after removing potential intron
                        query_seq = re.sub('[-]', "", str(seqs_nointron[0]))
                        hit_seq = re.sub('[-]', "", str(seqs_nointron[1]))

                        hit_start_codon = hit_seq[0:3]
                        hit_stop_codon = hit_seq[len(hit_seq) - 3:len(hit_seq)]
                        perc_ident = Levenshtein.ratio(seqs_nointron[3], hit_seq)

                        if hit_start_codon == "ATG":
                            output_dict_syn['ATG'] = 1
                        else:
                            output_dict_syn['ATG'] = 0

                        if hit_stop_codon == "TAG":
                            output_dict_syn['Stop Codon'] = 1
                        elif hit_stop_codon == "TAA":
                            output_dict_syn['Stop Codon'] = 1
                        elif hit_stop_codon == "TGA":
                            output_dict_syn['Stop Codon'] = 1
                        else:
                            output_dict_syn['Stop Codon'] = 0

                        if len(hit_seq) % 3 == 0:
                            output_dict_syn['Triplet'] = 1
                        else:
                            output_dict_syn['Triplet'] = 0
                        if len(hit_seq) > 0:
                            if len(seqs_nointron[3]) / len(hit_seq) > 0.9:
                                output_dict_syn['Identical Length'] = 1
                        else:
                            output_dict_syn['Identical Length'] = 0

                        if len(denovo_seq) == hit.aln_span:
                            output_dict_syn['Query end'] = 1
                        else:
                            output_dict_syn['Query end'] = 0

                        if 'mitochondrion' in hit.query.id:
                            query_id = hit.query.id[:-25]
                        else:
                            query_id = re.sub('_[2,3,4,X,Y][L,R]?_Chromosome', "", hit.query.id)

                        # Adds output to dictionary
                        output_dict_syn['ID'] = query_id
                        output_dict_syn['Hit'] = 1
                        output_dict_syn['Synteny'] = 1
                        output_dict_syn['Query Size'] = len(seqs_nointron[3])
                        output_dict_syn['Target Size'] = len(hit_seq)
                        output_dict_syn['Intron'] = seqs_nointron[2]
                        output_dict_syn['Target Start Pos.'] = hit.hit_start
                        output_dict_syn['Target Stop Pos.'] = hit.hit_end
                        output_dict_syn['Target Chromosome'] = chrom.group()[1:]
                        if hit.hit.id[3] == "_":
                            output_dict_syn['Population'] = hit.hit.id[0:3]
                        elif hit.hit.id[0] == 'Z':
                            output_dict_syn['Population'] = 'Zamb'
                        else:
                            output_dict_syn['Population'] = hit.hit.id[0:2]
                        output_dict_syn['Target Sequence'] = hit_seq
                        output_dict_syn['BitScore'] = hit.bitscore
                        output_dict_syn['Identical'] = round(perc_ident, 4)
                        output_dict_syn['Insertion'] = gap_query
                        output_dict_syn['Deletion'] = gap_hit
                        output_dict_syn['evalue'] = hit.evalue
                        output_dict_syn['Query start'] = hit.query_start + 1
                        output_dict_syn['Query stop codon'] = denovo_seq[len(denovo_seq) - 3:len(denovo_seq)]
                        output_list_target.append(output_dict_syn.copy())

                    except IndexError:
                        print("IndexError")
                        pass
        return n + 1


# Function fuses
def synteny_fuse():
    synteny_list = target_df.to_dict('records')
    genome_list = genome_df.to_dict('records')

    for dic in synteny_list:
        # Finds match by ID, Bitscore, Population and Sequence
        match = [item for item in genome_list if (item["ID"] == dic['ID']) & (item["BitScore"] == dic['BitScore']) & (
                item["Population"] == dic['Population']) & (item["Target Sequence"] == dic['Target Sequence'])]
        for entry in match:
            idx = genome_list.index(entry)
            print("{} gets updated with Synteny information. Index: {}".format(entry['ID'], idx))

            entry['Synteny'] = 1
            entry['ATG'] = dic['ATG']
            entry['Stop Codon'] = dic['Stop Codon']
            entry['Triplet'] = dic['Triplet']
            entry['Query Size'] = dic['Query Size']
            entry['Target Size'] = dic['Target Size']
            entry['Intron'] = dic['Intron']
            entry['Target Start Pos.'] = dic['Target Start Pos.']
            entry['Target Stop Pos.'] = dic['Target Stop Pos.']
            entry['Identical'] = dic['Identical']
            entry['Insertion'] = dic['Insertion']
            entry['Deletion'] = dic['Deletion']
            entry['Identical Length'] = dic['Identical Length']
            entry['evalue'] = dic['evalue']
            entry['Query start'] = dic['Query start']
            entry['Query end'] = dic['Query end']
            entry['Query stop codon'] = dic['Query stop codon']
            genome_list.pop(idx)
            genome_list.insert(idx, entry)

    fused_df = pd.DataFrame(genome_list)
    return fused_df


# Lists to store created dict by blast_parser
output_list = []
output_list_target = []
# Path to blast result dict
path = r'path to folder\211210_blast_results'

# Sequences of corresponding denovo genes
denovo_seqs = {rec.id: rec.seq for rec in
               SeqIO.parse(r"path to folder\OrthogroupsRefSeq.fa",
                           "fasta")}

# Number of files in Path
num_files = len(os.listdir(path))
n = 1
# Loop over files in Path
for filename in os.listdir(path):
    # Path to result file
    result = os.path.join(path, filename)
    print("Parsing: {} {} of {} {} % ".format(filename, n, num_files, round(float(n / int(num_files)) * 100, 2)))
    # Checks if result file is empty or has no entry for syntenic blasts
    if "_region.xml" in filename:
        # Checks if result file is empty
        if os.path.getsize(result) == 0:
            continue
        # Checks if result file has no match
        elif os.path.getsize(result) < 2500:
            continue
        else:
            n = blast_parser(result, denovo_seqs, n)
    else:
        n = blast_parser(result, denovo_seqs, n)
    print("")
# Constructs dataframes out of list of dicts
genome_df = pd.DataFrame(output_list)
target_df = pd.DataFrame(output_list_target)
# Fuses genome_df and target_df
fused_df = synteny_fuse()
# Saves returned fused_df as csv
fused_df.to_csv(r'path to folder\211214_results.csv', index=False)

