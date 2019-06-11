import cvxpy
import numpy as np

def solve_library(O, lib_lim, n_templates, bins=1e3, verbose=True, parallel=True,
                  approximate=False, time_limit=0, threads=0):
    ####################
    # i = n_targets    #
    # s = n_templates  #
    # p = n_var_pos    #
    # a = n_aas        #
    ####################

    # D.shape = (n_codons, n_aas)
    D = np.load('D.npy')
    D_hat = np.load('D_hat.npy')
    
    # Extract number of codons and number of amino acid options
    n_codons, n_aas = D.shape
    n_targets = len(O)
    n_var_pos = len(O[0])
    
    # Set up variables
    t = cvxpy.Variable(n_targets, boolean=True)
    G = {s: cvxpy.Variable((n_var_pos, n_codons), boolean=True) for s in range(n_templates)}
    C = []
    B = cvxpy.Variable((n_targets, n_templates), boolean=True)

    # Set up constraints
    constraints = []

    # Define relationship between C, G, and D
    for s in range(n_templates):
        C.append(G[s] * D_hat)

    # Constrain only one deg. codon can be used
    for s in range(n_templates):
        constraints.append(cvxpy.sum(G[s], axis=1) == 1)

    # Constrain "and" for all positions (check for cover of target)
    for i in range(n_targets):
        for s in range(n_templates):
            expression = cvxpy.sum(cvxpy.multiply(O[i], C[s])) - n_var_pos + n_var_pos * (1 - B[i, s])
            constraints.append(expression >= 0)
            constraints.append(expression <= n_var_pos)

    # Constrain "or" for all oligos (check whether target is covered by at least one oligo)
    expression =  - cvxpy.sum(B, axis=1) + (n_templates + 1) * t
    constraints.append(expression >= 0)
    constraints.append(expression <= n_templates)

    # Constrain library size for a single oligo
    if n_templates == 1 and not approximate:
        print('Using exact library size.\n')
        # Constrain library size
        lib_size = cvxpy.sum(G[0] * cvxpy.log(cvxpy.sum(D, axis=1)))
        constraints.append(lib_size <= cvxpy.log(lib_lim))
        
    else:
        print('Using approximate library size.\n')
        if bins > lib_lim:
            bins = lib_lim
        bin_constraints, lib_size = bin_oligo_count(lib_lim, n_templates, G, D, bins)
        constraints.extend(bin_constraints)

    # Define the objective            
    objective = cvxpy.sum(t)
    
    # Maximize covered sequences subject to constraints
    problem = cvxpy.Problem(cvxpy.Maximize(objective), constraints)
    
    aux_params = {}
    
    if time_limit > 0:
        aux_params['TimeLimit'] = time_limit
        
    if time_limit > 0:
        aux_params['Threads'] = threads

    # Solving the problem
    problem.solve(solver=cvxpy.GUROBI, verbose=verbose, parallel=parallel, **aux_params)
    
    # Make all variables available within a dictionary
    solution = {
        'binary_coverage': t.value,
        'coverage_count': np.sum(B.value, axis=1),
        'codon_selection': np.stack([G[x].value for x in G]),
        'problem': problem
    }
    
    return solution


def bin_oligo_count(lib_lim, n_templates, G, D, n_bins=1e3):
    if n_bins < lib_lim:
        n_bins = int(n_bins)
    else:
        n_bins = int(lib_lim)
        
    first_lower_bin = 0
    
    bin_upper_lim = np.log(np.linspace(1, lib_lim, num=n_bins))
    bin_lower_lim = np.array([bin_upper_lim[i-1] if i > 0 else first_lower_bin for i in range(n_bins)])
    
    bins = cvxpy.Variable((n_templates, n_bins), boolean=True)
    
    constraints = []
    
    for s in range(len(G)):
        log_n_seq = cvxpy.sum(G[s] * cvxpy.log(cvxpy.sum(D, axis=1)))
        
        constraints.append(log_n_seq <= cvxpy.log(lib_lim))
        
        constraints.append(bins[s, :] * bin_lower_lim <= log_n_seq)
        constraints.append(bins[s, :] * bin_upper_lim >= log_n_seq)
        constraints.append(cvxpy.sum(bins[s, :]) == 1)
        
    lib_size = cvxpy.sum(bins * np.exp(bin_upper_lim))
    constraints.append(lib_size <= lib_lim)
        
    return constraints, lib_size
