import click
import numpy as np
from decode.ilp import solve_library
from decode.seq_utils import *
import time
import json

@click.command()
@click.argument('alignment_file', required=True, type=click.Path(exists=True, dir_okay=False, readable=True))
@click.argument('output_file', required=True, type=click.Path(exists=False, dir_okay=False, readable=True))
@click.option('--limit', required=True, type=int, prompt='Total lib size limit', help='Total library size limit.')
@click.option('--sublib', required=True, type=int, prompt='Sublibrary limit', help='Total sublibrary limit')
@click.option('--bins', default=100, show_default=True, help='Specify the number of bins for approximation of multi-sublibrary optimizations.')
@click.option('--time-limit', default=0, show_default=True, help='Time limit in seconds for the ILP solver.')
@click.option('--threads', default=0, show_default=True, help='Time limit in seconds for the ILP solver.')
@click.option('-q', '--quiet', is_flag=True, default=False, show_default=True, help='Run quietly.')


def optimize_lib(alignment_file, output_file, limit, sublib, bins, time_limit, threads, quiet):
    """Optimize a degenerate codon library given a set of target sequences, a total size limit, and a sublibrary count limit."""
    
    if quiet:
        verbose = False
    else: verbose=True
    
    if not quiet:
        click.echo('Reading MSA file...')
    
    # Read the MSA file
    fixed_positions, variable_positions, sequences = process_msa(alignment_file, 'clustal')
    
    if not quiet:
        click.echo('Removing fixed positions...')
    
    # Generate the one-hot encoded input
    n_targets, n_var_pos, O = create_O(sequences)
    
    if not quiet:
        click.echo('')
        click.echo('Library size limit:\t\t{}'.format(limit))
        click.echo('Sublibrary count limit:\t\t{}'.format(sublib))
        click.echo('Number of targets:\t\t{}'.format(n_targets))
        click.echo('Number of variable positions:\t{}'.format(n_var_pos))
        click.echo('')
    
    # Start timing
    start = time.time()
    
    # Get the ILP solver solution
    solution = solve_library(O, limit, sublib, bins=bins, verbose=verbose, time_limit=time_limit, threads=threads)
    
    # End timing
    end = time.time()
    
    # Extract setup and solve times
    total_time = end - start
    solve_time = solution['problem'].solver_stats.solve_time
    construct_time = total_time - solve_time
    
    # Get library stats
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
    
    # Parse the library
    total_seq_length = len(fixed_positions) + len(variable_positions)
    parsed_lib = parse_lib(total_seq_length, fixed_positions, codons)

    if not quiet:
        click.echo('Writing output...')

    if solution['problem'].status == 'optimal':
        solution_optimal = True

    elif solution['problem'].status == 'optimal_inaccurate':
        solution_optimal = False

    # Generate and write the output
    data = {
        'solution_optimal': solution_optimal,
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
    
    with open(output_file, 'w') as fp:
        json.dump(data, fp)
    
if __name__ == '__main__':
    optimize_lib()
    