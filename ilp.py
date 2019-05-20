import cvxpy

def solve_library(S, lib_lim, n_oligos, A=None, solver='gurobi', verbose=True, parallel=True):
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