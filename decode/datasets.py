import pickle
import numpy as np
from importlib.resources import open_binary


with open_binary('decode.data', 'codons.pickle') as codon_file:
    all_codons = pickle.load(codon_file)

with open_binary('decode.data', 'codon_aa_mapping.pickle') as codon_aa_mapping_file:
    codon_aa_mapping = pickle.load(codon_aa_mapping_file)

with open_binary('decode.data', 'D.npy') as D_file:
    D = np.load(D_file)
    
with open_binary('decode.data', 'D_hat.npy') as D_hat_file:
    D_hat = np.load(D_hat_file)
