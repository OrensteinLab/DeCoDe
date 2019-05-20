import numpy as np
import pickle

with open('codons.txt', 'rb') as codon_file:
    all_codons = pickle.load(codon_file)

    
with open('codon_aa_mapping.txt', 'rb') as codon_file:
    codon_aa_mapping = pickle.load(codon_file)
    

def generate_seqs(count, length, dist=None):
    if dist is None:
        dist = np.random.exponential(size=length)
        dist /= dist.sum()
    
    seqs = []
    aa_choices = np.array(list('ACDEFGHIKLMNPQRSTVWY'))
    wt = ''.join([i for i in np.random.choice(aa_choices, length)])
    seqs.append(wt)
    
    for i in range(count - 1):
        seq = np.random.choice(seqs)
        pos = np.random.choice(list(range(length)), p=dist)
        mutant = np.random.choice(aa_choices, 1)[0]
        seq_list = list(seq)
        seq_list[pos] = mutant
        seqs.append(''.join(seq_list))
        
    return seqs


def create_S(sequences):
    n_targets = len(sequences)
    n_var_pos = len(sequences[0])
    
    S = {i: np.zeros((len(sequences[0]), 22)) for i in range(len(sequences))}
    
    aa_idx = 'ACDEFGHIKLMNPQRSTVWY*_'
    for i in S.keys():
        seq = sequences[i]
        for j, aa in enumerate(seq):
            S[i][j, aa_idx.index(aa)] = 1
            
    return (n_targets, n_var_pos, S)


def extract_solution(problem):
    final_objective = problem.value
    t, X, Z, B = problem.variables()
    return (t, X, Z, B)


def retrieve_codons(X, codons=all_codons):
    keys = [all_codons[i] for i in np.argmax(X.value, axis=1)]
    return keys


def calc_num_seqs(keys, codon_aa_mapping=codon_aa_mapping):
    return np.product([len(codon_aa_mapping[key]) for key in keys])


def calc_seq_prob(sequence, codon_keys, codon_aa_mapping=codon_aa_mapping):
    p = 0
    
    for i, aa in enumerate(sequence):
        aa_dist = codon_aa_mapping[codon_keys[i]]
        p_i = aa_dist[aa] / np.sum([aa_dist[key] for key in aa_dist])
        if p_i == 0:
            return p_i
        
        p += np.log(p_i)
        
    return np.exp(p)


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
    