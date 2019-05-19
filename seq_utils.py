import cvxpy
import numpy as np

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


def retrieve_codons(all_codons, X):
    keys = [all_codons[i] for i in np.argmax(X, axis=1)]
    return keys


def calc_num_seqs(codon_aa_mapping, keys):
    return np.product([len(codon_aa_mapping[key]) for key in keys])


def extract_solution(problem):
    final_objective = problem.value
    t, X, Z, B = problem.variables()
    return (t, X, Z, B)
    
    
def solve_library(lib_lim, n_oligos, A=None, solver=, verbose=True, parallel=True):
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

    # Constrain library size
    for j in range(n_oligos):
        # Constrain library size
        log_n_seq = cvxpy.sum(X[j] * cvxpy.log(cvxpy.sum(A, axis=1)))
        constraints.append(log_n_seq <= cvxpy.log(lib_lim))

    # Define the objective            
    objective = cvxpy.sum(t)

    # We tell cvxpy that we want to maximize sequence count 
    # subject to constraints. All constraints in 
    # cvxpy must be passed as a list.
    problem = cvxpy.Problem(cvxpy.Maximize(objective), constraints)

    # Solving the problem
    problem.solve(solver=cvxpy.GUROBI, verbose=verbose, parallel=parallel)
    
    return problem
    