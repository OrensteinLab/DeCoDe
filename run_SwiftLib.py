import click
import numpy as np
from decode.swiftlib_reader import run_swift_lib
from decode.seq_utils import *
import time
import json

@click.command()
@click.option('--limit', required=True, type=int, prompt='Total lib size limit', help='Total library size limit.')
@click.option('--sublib', required=True, type=int, prompt='Sublibrary limit', help='Total sublibrary limit')
@click.option('-q', '--quiet', is_flag=True, default=False, show_default=True, help='Run quietly.')

# solve_library(S, lib_lim, n_oligos, A=codon_table, bins=bins, solver='gurobi', verbose=True, parallel=True, approximate=False)

def optimize_lib(limit, sublib, quiet):
    """Optimize a redundent codon library given a set of target sequences, a total size limit, and a sublibrary count limit."""
    
    if not quiet:
        click.echo('Reading MSA file...')
    
    fixed_positions, variable_positions, sequences = process_msa('examples/gfp/gfp_239.aln', 'clustal')
    
    if not quiet:
        click.echo('')
        click.echo('Library size limit:\t\t{}'.format(limit))
        click.echo('Sublibrary count limit:\t\t{}'.format(sublib))
        click.echo('Number of targets:\t\t{}'.format(len(sequences)))
        click.echo('Number of variable positions:\t{}'.format(len(variable_positions)))
        click.echo('')
    
    if not quiet:
        click.echo('Creating CSV input string...')
                   
    with open('aux_data/gfp_239_input.csv', 'r') as csv_file:
        csv_string = csv_file.read()
        
    split_string = csv_string.split('\n')
    split_string[2] = split_string[2].replace('1', str(sublib))
    input_string = '\n'.join(split_string)
        
    if not quiet:
        click.echo('Solving for degenerate codon library with SwiftLib...\n')

    codons = run_swift_lib(input_string, 'csv', limit, sublib)
    
    n_var_pos = len(variable_positions)
    coverage = check_if_seqs_in_soln(sequences, codons)
    n_covered = int(sum(coverage))
    total_lib_size = int(calc_num_seqs(codons))
    on_target_p = calc_prob_on_target(sequences, codons)
    
    if not quiet:
        click.echo('')
        click.echo('Number of covered targets:\t{:d}'.format(n_covered))
        click.echo('Total library size:\t\t{:d}'.format(total_lib_size))
        click.echo('Probability on target:\t\t{:0.5f}'.format(on_target_p))
        click.echo('')
    
    # Parse the library
    total_seq_length = len(fixed_positions) + len(variable_positions)
    parsed_lib = parse_lib(total_seq_length, fixed_positions, codons)
    
    if not quiet:
        click.echo('Writing output...')
    
    data = {
        'fixed_positions': fixed_positions,
        'variable_positions': variable_positions,
        'sequences': sequences,
        'coverage': coverage,
        'n_var_pos': n_var_pos,
        'n_covered': n_covered,
        'total_lib_size': total_lib_size,
        'on_target_p': on_target_p,
        'parsed_lib': parsed_lib
    }
    
    with open('results/sl_comparison/sl/gfp_239_{:d}_{:d}.json'.format(limit, sublib), 'w') as fp:
        json.dump(data, fp)
    
if __name__ == '__main__':
    optimize_lib()
    