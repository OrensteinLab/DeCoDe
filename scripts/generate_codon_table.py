from Bio.Data.CodonTable import standard_dna_table as sdt
import itertools as it
from collections import Counter
import numpy as np
import pickle


# Define a function to extract the least degenerate codon
# for each distribution of amino acids
def extract_smallest_codon(codon_map, keys):
    keys = list(keys)
    unique_rows = np.unique(codon_map, axis=0)
    unique_row_keys = []
    
    for row in unique_rows:
        key_inds = np.where((codon_map == row).all(axis=1))[0]
        unique_row_keys.append([keys[i] for i in key_inds])
        
    smallest_key_set = []
        
    for key_set in unique_row_keys:
        min_deg_value = np.inf
        min_keys = []
        for key in key_set:
            if key == '___':
                deg_value = 1
                
            else:
                deg_value = np.product([len(deg_nuc_f[nt]) for nt in key])
            
            if deg_value == min_deg_value:
                min_keys.append(key)
                
            elif deg_value < min_deg_value:
                min_deg_value = deg_value
                min_keys = [key]
                
        smallest_key_set.append(min_keys)
        
    return unique_rows, smallest_key_set


if __name__ == '__main__':
    
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

    # Get a list of all possible degenerate codons
    all_codons = list(it.product(deg_nuc_f, repeat=3))

    codon_aa_mapping = {}

    for codon in all_codons:
        single_codons = it.product(deg_nuc_f[codon[0]], deg_nuc_f[codon[1]], deg_nuc_f[codon[2]])

        aas = []

        for s in single_codons:
            aas.append(codon_table_f[''.join(s)])

        codon_aa_mapping[''.join(codon)] = Counter(aas)
        
    # Add in the gap codon

    codon_aa_mapping['___'] = Counter(['-'])
    
    # Generate the binary version of D first to get coverage
    D_hat = np.zeros((len(codon_aa_mapping), 22))

    aa_idx = 'ACDEFGHIKLMNPQRSTVWY*-'

    for idx, key in enumerate(codon_aa_mapping.keys()):
        mapping_keys = codon_aa_mapping[key].keys()
        for aa in mapping_keys:
            D_hat[idx, aa_idx.index(aa)] = 1

    D_hat, D_keys = extract_smallest_codon(D_hat, codon_aa_mapping.keys())

    # Then generate D with just the reduced codon set
    D = np.zeros_like(D_hat)

    for idx, key_set in enumerate(D_keys):
        key = key_set[0]
        mapping_keys = codon_aa_mapping[key].keys()
        for aa in mapping_keys:
            D[idx, aa_idx.index(aa)] = codon_aa_mapping[key][aa]
    
    # Save the D and D_hat matrices for use in the ILP
    np.save('../data/D', D)
    np.save('../data/D_hat', D_hat)

    # Save the codon sets and mapping for use in parsing the final library
    with open('../data/codons.pickle', 'wb') as codon_file:
        pickle.dump(D_keys, codon_file)

    with open('../data/codon_aa_mapping.pickle', 'wb') as codon_file:
        pickle.dump(codon_aa_mapping, codon_file)
    