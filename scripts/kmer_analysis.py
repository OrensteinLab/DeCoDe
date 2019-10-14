import itertools as it
from tqdm import *
import json

def extract_aa_combos_per_position(lib_dict):
    # Get amino acid combinations at each position
    libs = []
    for lib in lib_dict['parsed_lib']:
        sublib = []
        for pos in range(len(lib)):
            if isinstance(lib[pos], dict):
                sublib.append([aa for aa in list(lib[pos].keys())[0].split('/') if aa != '*'])
            else:
                sublib.append([lib[pos]])
        libs.append(sublib)
    
    return libs


def extract_aa_pos_combos(libs, pos_combos):
    # Get all of the position/amino acid combinations
    lib_mers = []
    for lib in libs:
        for combo in pos_combos:
            position_ids = [lib[i] for i in combo]
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

def proc_gfp_comp(dc_filename, sl_filename, lib_size, sublibs, kmers, target_kmer_dict={}):
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
    base_line = '{},{},{},{},{},{},{},{},{},{}\n'
    lines = []
    
    # Loop through analysis for each file pair/kmer combo
    for k in tqdm(kmers):
        pos_combos = list(it.combinations(dc_lib_dict['variable_positions'], k))
        
        dc_kmers = set(extract_aa_pos_combos(dc_sublibs, pos_combos))
        sl_kmers = set(extract_aa_pos_combos(sl_sublibs, pos_combos))
        
        # Let's not repeat ourselves
        if k not in target_kmer_dict:
            target_kmer_dict[k] = set(extract_target_kmers(targets, pos_combos))
            
        target_kmers = target_kmer_dict[k]
        
        # Get the size of the set intersections for all combinations
        DCiSL = len(dc_kmers.intersection(sl_kmers))
        DCiT = len(dc_kmers.intersection(target_kmers))
        SLiT = len(sl_kmers.intersection(target_kmers))
        
        # Add output line to list
        lines.append(base_line.format('gfp', lib_size, sublibs, k,
                                      len(dc_kmers), len(sl_kmers), len(target_kmers),
                                      DCiSL, DCiT, SLiT))
        
    return target_kmer_dict, lines


def extract_sl_p1sol2(k):
    # Solution from SwiftLib Problem 1 DP Sol. 2
    sl_p1_sol2 = [
        'ADEGHIKLMNPQRSTV',
        'ADEGKLMNRSTVW',
        'AILSTV',
        'AS',
        'ALMSTV',
        'AGST',
        'ADEGKNRST',
        'DEHKLMNQRSTW',
        'ADEGHIKLMNPQRSTV'
    ]

    pos_combos = list(it.combinations([187, 188, 190, 191, 192, 196, 257, 258, 259], k))
    pos_combos_renum = list(it.combinations(list(range(9)), k))

    sl_p1_sol2 = [list(i) for i in sl_p1_sol2]

    sl_mers = []

    for num in range(len(pos_combos)):
        combo = pos_combos_renum[num]
        position_ids = [sl_p1_sol2[i] for i in combo]
        for aas in it.product(*position_ids):
            sl_mers.append((pos_combos[num], aas))

    return set(sl_mers)


def proc_rosetta_comp(dc_filename, kmers):
    with open(dc_filename, 'r') as jsonfile:
        dc_lib_dict = json.load(jsonfile)

    # Process DeCoDe kmers
    dc_sublibs = extract_aa_combos_per_position(dc_lib_dict)

    # Get targets
    targets = extract_targets(dc_lib_dict)

    # Get and return the lines (one per kmer) with the kmer counts
    # and set intersections
    base_line = '{},{},{},{},{},{},{},{},{},{}\n'
    lines = []

    # Loop through analysis for each file pair/kmer combo
    for k in tqdm(kmers, leave=False):
        pos_combos = list(it.combinations(dc_lib_dict['variable_positions'], k))

        dc_kmers = set(extract_aa_pos_combos(dc_sublibs, pos_combos))
        sl_kmers = extract_sl_p1sol2(k)

        target_kmers = set(extract_target_kmers(targets, pos_combos))

        # Get the size of the set intersections for all combinations
        DCiSL = len(dc_kmers.intersection(sl_kmers))
        DCiT = len(dc_kmers.intersection(target_kmers))
        SLiT = len(sl_kmers.intersection(target_kmers))

        # Add output line to list
        lines.append(base_line.format('1xbi', 320000000, 2, k,
                                      len(dc_kmers), len(sl_kmers), len(target_kmers),
                                      DCiSL, DCiT, SLiT))

    return lines


if __name__ == '__main__':
    
    # Set up the list for the output file
    all_lines = []
    all_lines.append('task,lib_limit,sublibs,k,dc_kmers,sl_kmers,target_kmers,DCiSL,DCiT,SLiT\n')
    
    # List the kmers to analyze
    kmers = [2, 3, 4]
    
    # Process 1xbi Rosetta experiment
    rosetta_file = 'results/sl_comparison/ilp/1xbi_320000000_2.json'

    all_lines.extend(proc_rosetta_comp(rosetta_file, kmers))
    
    # Process GFP experiments
    filename_blanks = [
        [100000, 1000000, 10000000, 100000000, 1000000000],
        [1, 2]
    ]

    filename_combos = list(it.product(*filename_blanks))
    gfp_files_base = 'results/sl_comparison/{}/gfp_239_{}_{}.json'

    target_kmer_dict = {}

    for fn in tqdm(filename_combos):
        lib_size = fn[0]
        lib_count = fn[1]

        dc_filename = gfp_files_base.format('ilp', *fn)
        sl_filename = gfp_files_base.format('sl', *fn)

        target_kmer_dict, comp_lines = proc_gfp_comp(dc_filename, sl_filename, lib_size,
                                                 lib_count, kmers, target_kmer_dict)
        all_lines.extend(comp_lines)
        
    # Write output to file
    with open('results/sl_comparison/kmers.csv', 'w') as f:
        for item in all_lines:
            f.write(item)
