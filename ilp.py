import cvxpy
import numpy as np

def solve_library(S, lib_lim, n_oligos, A=None, bins=1e3, solver='gurobi', verbose=True, parallel=True, approximate=False, time_limit=0):
    #################
    # i = n_targets #
    # j = s = n_oligos  #
    # k = p = n_var_pos #
    # l = d = n_codons  #
    # m = a = n_aas     #
    #################

    # A.shape = (n_codons, n_aas)
    if A is None:
        A = np.load('A.npy')
    
    # Extract number of codons and number of amino acid options
    n_codons, n_aas = A.shape
    n_targets = len(S)
    n_var_pos = len(S[0])
    
    # Set up variables
    t = cvxpy.Variable(n_targets, boolean=True)
    X = {j: cvxpy.Variable((n_var_pos, n_codons), boolean=True) for j in range(n_oligos)}
    
#     Z = {j: cvxpy.Variable((n_var_pos, n_aas), boolean=True) for j in range(n_oligos)}
    Z = []
    
    B = cvxpy.Variable((n_targets, n_oligos), boolean=True)

    # Set up constraints
    constraints = []

    # Define relationship between Z, X, and A
    for j in range(n_oligos):
#         constraints.append(Z[j] == X[j] * A)
        Z.append(X[j] * A)

    # Constrain only one deg. codon can be used
    for j in range(n_oligos):
        constraints.append(cvxpy.sum(X[j], axis=1) == 1)

    # Constrain "and" for all positions (check for cover of target)
    for i in range(n_targets):
        for j in range(n_oligos):
            expression = cvxpy.sum(cvxpy.multiply(S[i], Z[j])) - n_var_pos + n_var_pos * (1 - B[i, j])
            constraints.append(expression >= 0)
            constraints.append(expression <= n_var_pos)

    # Constrain "or" for all oligos (check whether target is covered by at least one oligo)
    expression =  - cvxpy.sum(B, axis=1) + (n_oligos + 1) * t
    constraints.append(expression >= 0)
    constraints.append(expression <= n_oligos)

    # Constrain library size for a single oligo
    if n_oligos == 1 and not approximate:
        print('Using exact library size.\n')
        # Constrain library size
        lib_size = cvxpy.sum(X[0] * cvxpy.log(cvxpy.sum(A, axis=1)))
        constraints.append(lib_size <= cvxpy.log(lib_lim))
        
    else:
        print('Using approximate library size.\n')
        if bins > lib_lim:
            bins = lib_lim
        bin_constraints, lib_size = bin_oligo_count(lib_lim, n_oligos, X, A, bins)
        constraints.extend(bin_constraints)

    # Define the objective            
    objective = cvxpy.sum(t)
    
    # We tell cvxpy that we want to maximize sequence count 
    # subject to constraints. All constraints in 
    # cvxpy must be passed as a list.
    problem = cvxpy.Problem(cvxpy.Maximize(objective), constraints)

    # Solving the problem
    problem.solve(solver=cvxpy.GUROBI, verbose=verbose, parallel=parallel)
    
    # Make all variables available within a dictionary
    solution = {
        'binary_coverage': t.value,
        'coverage_count': np.sum(B.value, axis=1),
        'codon_selection': np.stack([X[x].value for x in X]),
        'problem': problem
    }
    
    return solution


def bin_oligo_count(lib_lim, n_oligos, X, A, n_bins=1e3):
    if n_bins < lib_lim:
        n_bins = int(n_bins)
    else:
        n_bins = int(lib_lim)
        
    first_lower_bin = -1 #(np.log(lib_lim) - np.log(lib_lim + 1)) / 10
    
    bin_upper_lim = np.log(np.linspace(1, lib_lim, num=n_bins))
    bin_lower_lim = np.array([bin_upper_lim[i-1] if i > 0 else first_lower_bin for i in range(n_bins)])
    
#     bin_lims = {j: cvxpy.Variable((bins, 2), boolean=True) for j in range(n_oligos)}
    bins = cvxpy.Variable((n_oligos, n_bins), boolean=True)
    
    constraints = []
    
    for j in range(len(X)):
        log_n_seq = cvxpy.sum(X[j] * cvxpy.log(cvxpy.sum(A, axis=1)))
        
        constraints.append(log_n_seq <= cvxpy.log(lib_lim))
        
        constraints.append(log_n_seq - bins[j, :] * bin_lower_lim >= 0)
        constraints.append(bins[j, :] * bin_upper_lim >= log_n_seq)
        constraints.append(cvxpy.sum(bins[j, :]) == 1)
        
#         upper_lim_check = bin_upper_lim - log_n_seq + np.log(lib_lim) * (1 - bin_lims[j][:, 0])
#         constraints.append(upper_lim_check >= 0)
#         constraints.append(upper_lim_check <= np.log(lib_lim))
        
#         lower_lim_check = bin_lower_lim - log_n_seq + np.log(lib_lim) * bin_lims[j][:, 1]
#         constraints.append(lower_lim_check >= 0)
#         constraints.append(lower_lim_check <= np.log(lib_lim))
        
#         bin_assignment = bin_lims[j][:, 0] + bin_lims[j][:, 1] - 2 + 3 * (1 - bins[j, :])
#         constraints.append(bin_assignment >= 0)
#         constraints.append(bin_assignment <= 2)
        
    lib_size = cvxpy.sum(bins * np.exp(bin_upper_lim))
    constraints.append(lib_size <= lib_lim)
        
    return constraints, lib_size
