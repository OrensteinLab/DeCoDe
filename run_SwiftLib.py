import click
import numpy as np
from swiftlib_reader import run_swift_lib
from seq_utils import *
import time
import json

@click.command()
@click.argument('alignment_file', required=True, type=click.Path(exists=True, dir_okay=False, readable=True))
@click.option('--limit', required=True, type=int, prompt='Total lib size limit', help='Total library size limit.')
@click.option('--sublib', required=True, type=int, prompt='Sublibrary limit', help='Total sublibrary limit')
@click.option('-q', '--quiet', is_flag=True, default=False, show_default=True, help='Run quietly.')

# solve_library(S, lib_lim, n_oligos, A=codon_table, bins=bins, solver='gurobi', verbose=True, parallel=True, approximate=False)

def optimize_lib(alignment_file, limit, sublib, codon_table, bins, solver, verbose, parallel, approx, quiet):
    """Optimize a redundent codon library given a set of target sequences, a total size limit, and a sublibrary count limit."""
    
    if not quiet:
        click.echo('Reading MSA file...')
    
    fixed_positions, variable_positions, sequences = process_msa(alignment_file, 'clustal')
    
    if not quiet:
        click.echo('')
        click.echo('Library size limit:\t\t{}'.format(limit))
        click.echo('Sublibrary count limit:\t\t{}'.format(sublib))
        click.echo('Number of targets:\t\t{}'.format(len(sequences)))
        click.echo('Number of variable positions:\t{}'.format(len(variable_positions))
        click.echo('')
    
    if not quiet:
        click.echo('Creating FASTA string...')
                   
    fasta_string = gen_fasta(sequences)
                   
    if not quiet:
        click.echo('Solving for degenerate codon library with SwiftLib...\n')
    
    start = time.time()
    
    lib = run_swift_lib(input_string, 'fasta', limit, oligo_limit)
    
    end = time.time()
    
    total_time = end - start
    solve_time = solution['problem'].solver_stats.solve_time
    construct_time = total_time - solve_time
    
    n_covered = int(np.sum(solution['binary_coverage']))
    codons = retrieve_codons(solution['codon_selection'])
    total_lib_size = int(calc_num_seqs(codons))
    on_target_p = calc_prob_on_target(sequences, codons)
    
    if not quiet:
        click.echo('')
        click.echo('Number of covered targets:\t{:d}'.format(n_covered))
        click.echo('Total library size:\t\t{:d}'.format(total_lib_size))
        click.echo('Probability on target:\t\t{:0.5f}'.format(on_target_p))
        click.echo('')
    
    parsed_lib = parse_lib(fixed_positions, codons)
    
    if not quiet:
        click.echo('Writing output...')
    
    data = {
        'fixed_positions': fixed_positions,
        'variable_positions': variable_positions,
        'sequences': sequences,
        'coverage': list(solution['binary_coverage']),
        'n_var_pos': n_var_pos,
        'n_covered': n_covered,
        'total_lib_size': total_lib_size,
        'on_target_p': on_target_p,
        'parsed_lib': parsed_lib,
        'construct_time': construct_time,
        'solve_time': solve_time,
        'total_time': total_time
    }
    
    with open('ILP_lib_opt_{:d}_{:d}.json'.format(limit, sublib), 'w') as fp:
        json.dump(data, fp)
    
if __name__ == '__main__':
    optimize_lib()
    