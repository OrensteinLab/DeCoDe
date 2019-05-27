import cvxpy
import numpy as np

def solve_library(S, lib_lim, n_oligos, A=None, solver='gurobi', verbose=True, parallel=True, approximate=False):
    #################
    # i = n_targets #
    # j = n_oligos  #
    # k = n_var_pos #
    # l = n_codons  #
    # m = n_aas     #
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
    Z = {j: cvxpy.Variable((n_var_pos, n_aas), boolean=True) for j in range(n_oligos)}
    B = cvxpy.Variable((n_targets, n_oligos), boolean=True)

    # Set up constraints
    constraints = []

    # Define relationship between Z, X, and A
    for j in range(n_oligos):
        constraints.append(Z[j] == X[j] * A)

    # Constrain only one deg. codon can be used
    for j in range(n_oligos):
        for k in range(n_var_pos):
            constraints.append(cvxpy.sum(X[j][k, :]) == 1)

    # Constrain "and" for all positions (check for cover of target)
    for i in range(n_targets):
        for j in range(n_oligos):
            expression = cvxpy.sum(cvxpy.multiply(S[i], Z[j])) - n_var_pos + n_var_pos * (1 - B[i, j])
            constraints.append(expression >= 0)
            constraints.append(expression <= n_var_pos)

        # Constrain "or" for all oligos (check whether target is covered by at least one oligo)
        expression =  - cvxpy.sum(B[i, :]) + (n_oligos + 1) * t[i]
        constraints.append(expression >= 0)
        constraints.append(expression <= n_oligos)

    # Constrain library size for a single oligo
    if n_oligos == 1 and not approximate:
        print('Using exact library size')
        # Constrain library size
        log_n_seq = cvxpy.sum(X[0] * cvxpy.log(cvxpy.sum(A, axis=1)))
        constraints.append(log_n_seq <= cvxpy.log(lib_lim))
        
    else:
        print('Using approximate library size')
        bin_constraints = bin_oligo_count(lib_lim, n_oligos, X, A, bins=1e3)
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
    var_dict = {
        't': t,
        'X': X,
        'Z': Z,
        'B': B,
    }
    
    return problem


def bin_oligo_count(lib_lim, n_oligos, X, A, bins=1e3):
    bins = int(bins)
    first_lower_bin = (np.log(lib_lim) - np.log(lib_lim + 1)) / 10
    
    bin_upper_lim = np.log(np.linspace(1, lib_lim, num=bins))
    bin_lower_lim = np.array([bin_upper_lim[i-1] if i > 0 else first_lower_bin for i in range(bins)])
    
    bin_lims = {j: cvxpy.Variable((bins, 2), boolean=True) for j in range(n_oligos)}
    bins = cvxpy.Variable((n_oligos, bins), boolean=True)
    
    constraints = []
    
    for j in range(len(X)):
        log_n_seq = cvxpy.sum(X[j] * cvxpy.log(cvxpy.sum(A, axis=1)))
        
        constraints.append(log_n_seq <= np.log(lib_lim))
        
        upper_lim_check = bin_upper_lim - log_n_seq + np.log(lib_lim + 1) * (1 - bin_lims[j][:, 0])
        constraints.append(upper_lim_check >= 0)
        constraints.append(upper_lim_check <= np.log(lib_lim))
        
        lower_lim_check = bin_lower_lim - log_n_seq + np.log(lib_lim + 1) * bin_lims[j][:, 1]
        constraints.append(lower_lim_check >= 0)
        constraints.append(lower_lim_check <= np.log(lib_lim))
        
        bin_assignment = bin_lims[j][:, 0] + bin_lims[j][:, 1] - 2 + 3 * (1 - bins[j, :])
        constraints.append(bin_assignment >= 0)
        constraints.append(bin_assignment <= 2)
        
    constraints.append(cvxpy.sum(bins * np.exp(bin_upper_lim)) <= lib_lim)
        
    return constraints
