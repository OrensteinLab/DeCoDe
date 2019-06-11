import numpy as np
import pickle
from Bio import AlignIO as aln

#################################################
#### Datasets and  data generation functions ####
#################################################

with open('codons.pickle', 'rb') as codon_file:
    all_codons = pickle.load(codon_file)

    
with open('codon_aa_mapping.pickle', 'rb') as codon_file:
    codon_aa_mapping = pickle.load(codon_file)
    

def generate_seqs(count, length, dist=None):
    if dist is None:
        dist = np.random.exponential(size=length)
        dist /= dist.sum()
    
    seqs = []
    aa_choices = np.array(list('ACDEFGHIKLMNPQRSTVWY'))
    wt = ''.join([i for i in np.random.choice(aa_choices, length)])
    seqs.append(wt)
    
    while len(seqs) < count:
        seq = np.random.choice(seqs)
        pos = np.random.choice(list(range(length)), p=dist)
        mutant = np.random.choice(aa_choices, 1)[0]
        seq_list = list(seq)
        seq_list[pos] = mutant
        new_seq = ''.join(seq_list)
        if new_seq not in seqs:
            seqs.append(new_seq)
        
    return seqs


def gen_fasta(sequences):
    fasta = []
    i=0
    for seq in sequences:
        fasta.append('>Seq{}'.format(i))
        fasta.append(seq)
        i+=1
    fasta = '\n'.join(fasta)
    
    return fasta

#############################################
#### Input sequence processing functions ####
#############################################

def process_msa(filename, filetype):
    with open(filename, 'r') as handle:
        sequences = aln.read(handle, filetype)
    
    fixed_positions = {}
    variable_positions = []

    for position in range(len(sequences[0])):
        residues = list(set([seq[position] for seq in sequences]))
        if len(residues) == 1:
            fixed_positions[position] = residues[0]
        
        else:
            variable_positions.append(position)
        
    out_seqs = [''.join([seq[pos] for pos in variable_positions]) for seq in sequences]
    
    return (fixed_positions, variable_positions, out_seqs)

def create_O(sequences):
    n_targets = len(sequences)
    n_var_pos = len(sequences[0])
    
    O = {i: np.zeros((len(sequences[0]), 22)) for i in range(len(sequences))}
    
    aa_idx = 'ACDEFGHIKLMNPQRSTVWY*-'
    for i in O.keys():
        seq = sequences[i]
        for j, aa in enumerate(seq):
            O[i][j, aa_idx.index(aa)] = 1
            
    return (n_targets, n_var_pos, O)

#######################################
#### Solution extraction functions ####
#######################################

def retrieve_codons(codon_selection, codons=all_codons):
    keys = [[all_codons[i] for i in np.argmax(codon_selection, axis=2)[j]] for j in range(len(codon_selection))]
    return keys


def calc_num_seqs(keys, codon_aa_mapping=codon_aa_mapping):
    oligo_products = []
    
    for key_set in keys:
        oligo_products.append(np.product([len(codon_aa_mapping[key[0]]) for key in key_set]))
    
    return np.sum(oligo_products)


def calc_seq_prob(sequence, codon_keys, codon_aa_mapping=codon_aa_mapping):

    in_sl = []
    total_seqs = []
    
    for key_set in codon_keys:
        
        seq_total = 1
        p_i = 1
        
        for i, aa in enumerate(sequence):
            aa_dist = codon_aa_mapping[key_set[i][0]]
            
            total_aa_count = np.sum([aa_dist[key] for key in aa_dist])
            
            seq_total *= total_aa_count
            p_i *= aa_dist[aa]
            
        in_sl.append(p_i)
        total_seqs.append(seq_total)
        
    return np.sum(in_sl) / np.sum(total_seqs)


def calc_prob_on_target(target_seqs, codon_keys, codon_aa_mapping=codon_aa_mapping):
    p_on_target = 0
    
    for seq in target_seqs:
        p_on_target += calc_seq_prob(seq, codon_keys, codon_aa_mapping=codon_aa_mapping)
        
    return p_on_target


def check_if_seqs_in_soln(sequences, codon_keys, codon_aa_mapping=codon_aa_mapping):
    in_lib = []
    
    for seq in sequences:
        p = calc_seq_prob(seq, codon_keys, codon_aa_mapping)
        if p > 0:
            in_lib.append(True)
        else:
            in_lib.append(False)
            
    return in_lib
    
    
def parse_lib(seq_length, fixed_positions, codon_keys):
    
    parsed_lib = []
    
    for i, key_set in enumerate(codon_keys):
        
        parsed_sublib = []
        
        k = 0
        
        for j in range(seq_length):
            
            if j in fixed_positions:
                
                res = fixed_positions[j]
                
            else:
                aa_dist = codon_aa_mapping[key_set[k][0]]
                
                if np.sum([aa_dist[key] for key in aa_dist]) == 1:
                    res = max(aa_dist)
                    
                else:
                    res = {'/'.join(aa_dist.keys()): key_set[k]}
                    
                k += 1
                    
            if res != '-':
                    parsed_sublib.append(res)
                    
        parsed_lib.append(parsed_sublib)
    
    return parsed_lib
    