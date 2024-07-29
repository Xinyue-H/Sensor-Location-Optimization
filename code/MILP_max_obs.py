# sensor location optimization
import numpy as np
import pandas as pd
import cvxpy as cp

def maximize_observability_MILP(A_list, A_dict, n, p):
    '''
    ex.     
    K = 585  # number of modes
    p = 19  # number of sensors to place
    n = 22  # number of total states
    '''
    K = len(A_list)

    # define variables
    z = cp.Variable(n, boolean=True)
    x = cp.Variable((K, n))

    A = []
    for matrix in A_list:
        temp_matrix = np.zeros((n, n))
        temp_matrix[matrix > 0] = 1
        A.append(temp_matrix - np.diag(np.diagonal(temp_matrix)))
    A = np.array(A)
    w = np.fromiter(A_dict.values(), dtype=float) / sum(np.fromiter(A_dict.values(), dtype=float))

    # define constraints & objectives
    constraints = []
    ones = np.ones(n).reshape(n, 1)

    # Ensure exactly p sensors
    constraints += [cp.sum(z) == p]

    for k in range(K):
        constraints += [x[k] >= 0]
        constraints += [x[k] <= 1]
        constraints += [x[k] <= (A[k].T) @ x[k] + z]

    objective = cp.Maximize(w.T @ (x @ ones))

    prob = cp.Problem(objective, constraints)
    prob.solve()

    # Print the results
    # print(f"Optimal sensor placement (z): {z.value}")
    # print(f"Objective value: {prob.value}")
    return z.value
