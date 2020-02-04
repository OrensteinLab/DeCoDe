import itertools as it
from tqdm import *
from collections import Counter
import json
from Bio.Data.CodonTable import standard_dna_table as sdt
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from scipy.stats import wasserstein_distance
import numpy as np

def generate_deg_codon_table():
    # Set up the codon table
    codon_table_r = {}
    for codon, aa in sdt.forward_table.items():
        if aa not in codon_table_r.keys():
            codon_table_r[aa] = [codon]
        else:
            codon_table_r[aa].append(codon)

    codon_table_r['*'] = sdt.stop_codons

    codon_table_f = {}
    for aa, codons in codon_table_r.items():
        for codon in codons:
            codon_table_f[codon] = aa

    # Set up the nucleotide table
    deg_nuc_f = {
        'A': 'A',
        'C': 'C',
        'G': 'G',
        'T': 'T',
        'W': 'AT',
        'S': 'CG',
        'M': 'AC',
        'K': 'GT',
        'R': 'AG',
        'Y': 'CT',
        'B': 'CGT',
        'D': 'AGT',
        'H': 'ACT',
        'V': 'ACG',
        'N': 'ACGT'
    }

    deg_nuc_r = {value: key for key, value in deg_nuc_f.items()}

    all_codons = list(it.product(deg_nuc_f, repeat=3))

    codon_aa_mapping = {}

    for codon in all_codons:
        single_codons = it.product(deg_nuc_f[codon[0]], deg_nuc_f[codon[1]], deg_nuc_f[codon[2]])

        aas = []

        for s in single_codons:
            aas.append(codon_table_f[''.join(s)])

        codon_aa_mapping[''.join(codon)] = Counter(aas)

    codon_aa_mapping['___'] = Counter(['_'])
    
    return codon_aa_mapping


# Generate the degenerate codon table 
codon_table = generate_deg_codon_table()


def extract_aa_combos_per_position(lib_dict):
    # Get amino acid combinations at each position
    libs = []
    for lib in lib_dict['parsed_lib']:
        sublib = []
        for pos in range(len(lib)):
            if isinstance(lib[pos], dict):
                counts = [codon_table[lib[pos][key][0]] for key in lib[pos]]
                sublib.extend(counts)
            else:
                sublib.append(Counter([lib[pos]]))
        libs.append(sublib)
    
    return libs


def extract_aa_pos_combos(libs, pos_combos):
    # Get all of the position/amino acid combinations
    lib_mers = []
    for lib in libs:
        for combo in pos_combos:
            position_ids = [''.join([key * lib[i][key] for key in lib[i]]) for i in combo]
            for aas in it.product(*position_ids):
                lib_mers.append((combo, aas))
                
    return lib_mers
            
    
def extract_targets(lib_dict):
    # Extract the target sequences
    target_seqs = []
    for seq_num in range(len(lib_dict['sequences'])):
        target_vars = list(lib_dict['sequences'][seq_num])
        seq = []
        for i in range(len(lib_dict['parsed_lib'][0])):
            if i in lib_dict['variable_positions']:
                seq.append(target_vars.pop(0))
            else:
                seq.append(lib_dict['fixed_positions'][str(i)])

        target_seqs.append(''.join(seq))
    
    return target_seqs
    
    
def extract_target_kmers(target_seqs, pos_combos):
    # Get the kmer combinations from the target sequences
    target_mers = []
    for seq in target_seqs:
        for combo in pos_combos:
            ids = [seq[pos] for pos in combo]
            target_mers.append((combo, tuple(ids)))
            
    return target_mers


def proc_gfp_comp(dc_filename, sl_filename, lib_size, sublibs, kmers, target_kmer_dict={},
                  target_kmer_dict_all={}):
    with open(dc_filename, 'r') as jsonfile:
        dc_lib_dict = json.load(jsonfile)
        
    with open(sl_filename, 'r') as jsonfile:
        sl_lib_dict = json.load(jsonfile)
    
    # Process DeCoDe kmers
    dc_sublibs = extract_aa_combos_per_position(dc_lib_dict)
    sl_sublibs = extract_aa_combos_per_position(sl_lib_dict)
    
    # Get targets
    targets = extract_targets(dc_lib_dict)
    
    # Get and return the lines (one per kmer) with the kmer counts
    # and set intersections
    base_line = '{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n'
    lines = []
    
    # Loop through analysis for each file pair/kmer combo
    for k in tqdm(kmers):
        pos_combos = list(it.combinations(dc_lib_dict['variable_positions'], k))
        
        dc_kmers_all = extract_aa_pos_combos(dc_sublibs, pos_combos)
        sl_kmers_all = extract_aa_pos_combos(sl_sublibs, pos_combos)
        
        dc_kmers = set(dc_kmers_all)
        sl_kmers = set(sl_kmers_all)
        
        # Let's not repeat ourselves
        if (k not in target_kmer_dict) or (k not in target_kmer_dict_all):
            target_kmer_dict_all[k] = extract_target_kmers(targets, pos_combos)
            target_kmer_dict[k] = set(target_kmer_dict_all[k])
        
        target_kmers_all = target_kmer_dict_all[k]
        target_kmers = target_kmer_dict[k]
        
        # Get the size of the set intersections for all combinations
        DCiSL = len(dc_kmers.intersection(sl_kmers))
        DCiT = len(dc_kmers.intersection(target_kmers))
        SLiT = len(sl_kmers.intersection(target_kmers))
        
        # Count weighted scores
        target_counts = Counter(target_kmers_all)
        
        DCiT_all = sum([target_counts[kmer] for kmer in dc_kmers])
        SLiT_all = sum([target_counts[kmer] for kmer in sl_kmers])
        
        # Wasserstein distances
#         dc_wass_kmers = list(target_kmers.union(dc_kmers))
#         sl_wass_kmers = list(target_kmers.union(sl_kmers))
        
#         DC_kmer_counts = Counter(dc_kmers_all)
#         SL_kmer_counts = Counter(sl_kmers_all)
        
#         target_dc_kmer_spectrum = np.array([target_counts[k] for k in dc_wass_kmers])
#         target_dc_kmer_spectrum = target_dc_kmer_spectrum / np.sum(target_dc_kmer_spectrum)
#         dc_kmer_spectrum = np.array([DC_kmer_counts[k] for k in dc_wass_kmers])
#         dc_kmer_spectrum = dc_kmer_spectrum / np.sum(dc_kmer_spectrum)
        
#         target_sl_kmer_spectrum = np.array([target_counts[k] for k in sl_wass_kmers])
#         target_sl_kmer_spectrum = target_sl_kmer_spectrum / np.sum(target_sl_kmer_spectrum)
#         sl_kmer_spectrum = np.array([SL_kmer_counts[k] for k in sl_wass_kmers])
#         sl_kmer_spectrum = sl_kmer_spectrum / np.sum(sl_kmer_spectrum)
        
#         wass_DC_T = wasserstein_distance(dc_kmer_spectrum, target_dc_kmer_spectrum)
#         wass_SL_T = wasserstein_distance(sl_kmer_spectrum, target_sl_kmer_spectrum)
        
        # Add output line to list
        lines.append(base_line.format('gfp', lib_size, sublibs, k,
                                      len(dc_kmers), len(sl_kmers), len(target_kmers),
                                      DCiSL, DCiT, SLiT, len(dc_kmers_all), len(sl_kmers_all),
                                      len(target_kmers_all), DCiT_all, SLiT_all))
        
    return target_kmer_dict, target_kmer_dict_all, lines


def extract_sl_p1sol2(k):
    # Solution from SwiftLib Problem 1 DP Sol. 2    
    sl_p1_sol2 = [
        ['VNS'],
        ['DBG', 'RAM'],
        ['DYA'],
        ['KCA'],
        ['DYG'],
        ['RSC'],
        ['RVM'],
        ['VAM', 'WBG'],
        ['VNS']
    ]

    pos_combos = list(it.combinations([187, 188, 190, 191, 192, 196, 257, 258, 259], k))
    pos_combos_renum = list(it.combinations(list(range(9)), k))

    sl_p1_sol2 = [list(i) for i in sl_p1_sol2]
    
    solution = []
    
    for pos in sl_p1_sol2:
        aas = Counter()
        for codon in pos:
            aas += codon_table[codon]
            
        solution.append(''.join([key * value for key, value in aas.items()]))

    sl_mers = []

    for num in range(len(pos_combos)):
        combo = pos_combos_renum[num]
        position_ids = [solution[i] for i in combo]
        for aas in it.product(*position_ids):
            sl_mers.append((pos_combos[num], aas))

    return sl_mers


def proc_rosetta_comp(dc_filename, kmers):
    with open(dc_filename, 'r') as jsonfile:
        dc_lib_dict = json.load(jsonfile)

    # Process DeCoDe kmers
    dc_sublibs = extract_aa_combos_per_position(dc_lib_dict)

    # Get targets
    targets = extract_targets(dc_lib_dict)

    # Get and return the lines (one per kmer) with the kmer counts
    # and set intersections
    base_line = '{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n'
    lines = []

    # Loop through analysis for each file pair/kmer combo
    for k in tqdm(kmers):
        pos_combos = list(it.combinations(dc_lib_dict['variable_positions'], k))
        
        dc_kmers_all = extract_aa_pos_combos(dc_sublibs, pos_combos)
        sl_kmers_all = extract_sl_p1sol2(k)
        
        dc_kmers = set(dc_kmers_all)
        sl_kmers = set(sl_kmers_all)
        
        target_kmers_all = extract_target_kmers(targets, pos_combos)
        target_kmers = set(target_kmers_all)
        
        # Get the size of the set intersections for all combinations
        DCiSL = len(dc_kmers.intersection(sl_kmers))
        DCiT = len(dc_kmers.intersection(target_kmers))
        SLiT = len(sl_kmers.intersection(target_kmers))
        
        # Count weighted scores
        target_counts = Counter(target_kmers_all)
        
        DCiT_all = sum([target_counts[kmer] for kmer in dc_kmers])
        SLiT_all = sum([target_counts[kmer] for kmer in sl_kmers])
        
        # Wasserstein distances
#         dc_wass_kmers = list(target_kmers.union(dc_kmers))
#         sl_wass_kmers = list(target_kmers.union(sl_kmers))
        
#         DC_kmer_counts = Counter(dc_kmers_all)
#         SL_kmer_counts = Counter(sl_kmers_all)
        
#         target_dc_kmer_spectrum = np.array([target_counts[k] for k in dc_wass_kmers])
#         target_dc_kmer_spectrum = target_dc_kmer_spectrum / np.sum(target_dc_kmer_spectrum)
#         dc_kmer_spectrum = np.array([DC_kmer_counts[k] for k in dc_wass_kmers])
#         dc_kmer_spectrum = dc_kmer_spectrum / np.sum(dc_kmer_spectrum)
        
#         target_sl_kmer_spectrum = np.array([target_counts[k] for k in sl_wass_kmers])
#         target_sl_kmer_spectrum = target_sl_kmer_spectrum / np.sum(target_sl_kmer_spectrum)
#         sl_kmer_spectrum = np.array([SL_kmer_counts[k] for k in sl_wass_kmers])
#         sl_kmer_spectrum = sl_kmer_spectrum / np.sum(sl_kmer_spectrum)
        
#         wass_DC_T = wasserstein_distance(dc_kmer_spectrum, target_dc_kmer_spectrum)
#         wass_SL_T = wasserstein_distance(sl_kmer_spectrum, target_sl_kmer_spectrum)
        
        # Add output line to list
        lines.append(base_line.format('1xbi', 320000000, 4, k,
                                      len(dc_kmers), len(sl_kmers), len(target_kmers),
                                      DCiSL, DCiT, SLiT, len(dc_kmers_all), len(sl_kmers_all),
                                      len(target_kmers_all), DCiT_all, SLiT_all))

    return lines


if __name__ == '__main__':
    
    # Set up the list for the output file
    all_lines = []
    all_lines.append('task,lib_limit,sublibs,k,dc_kmers,sl_kmers,target_kmers,DCiSL,DCiT,SLiT,dc_kmers_all,sl_kmers_all,target_kmers_all,DCiT_all,SLiT_all\n')
    
    # List the kmers to analyze
    kmers = [2, 3, 4]
    
    # Process 1xbi Rosetta experiment
    rosetta_file = 'results/sl_comparison/ilp/1xbi_320000000_4.json'

    all_lines.extend(proc_rosetta_comp(rosetta_file, kmers))
    
    # Process GFP experiments
    filename_blanks = [
        [100000, 1000000, 10000000, 100000000, 1000000000],
        [1, 2]
    ]

    filename_combos = list(it.product(*filename_blanks))
    gfp_files_base = 'results/sl_comparison/{}/gfp_239_{}_{}.json'

    target_kmer_dict = {}
    target_kmer_dict_all = {}

    for fn in tqdm(filename_combos):
        lib_size = fn[0]
        lib_count = fn[1]

        dc_filename = gfp_files_base.format('ilp', *fn)
        sl_filename = gfp_files_base.format('sl', *fn)

        target_kmer_dict, target_kmer_dict_all, comp_lines = proc_gfp_comp(dc_filename, sl_filename, lib_size,
                                                 lib_count, kmers, target_kmer_dict, target_kmer_dict_all)
        all_lines.extend(comp_lines)
        
    # Write output to file
    with open('results/sl_comparison/kmers.csv', 'w') as f:
        for item in all_lines:
            f.write(item)
