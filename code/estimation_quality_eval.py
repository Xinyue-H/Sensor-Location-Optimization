import numpy as np
import matplotlib.pyplot as plt

def compute_error(x_list, density_list):
    n = len(density_list[0])  # number of links
    Absolute_percentage_error = 0  
    # add up error for all time steps
    for i in range(len(x_list)):
        # err is a vector of estimation errors (with length n)
        err = np.array(x_list[i]) - np.array(density_list[i])
        Absolute_percentage_error += np.sum(abs(err)/np.array(density_list[i]))
    MAPE = Absolute_percentage_error/len(x_list)/n
    return MAPE

def plot_mode_prob(mu_list, mode_list):
    true_mode_prob = list()
    T = len(mu_list)
    for i in range(T):
        true_mode_prob.append(mu_list[i][mode_list[i]])
    t_list = np.arange(T)
    plt.figure(figsize=(10, 5)) 
    plt.title("True mode probability in IMM filtering")
    plt.plot(t_list, true_mode_prob)



def observability_score(A_list, C, w):
    from numpy.linalg import matrix_power
    obs_num_list = []
    n = C.shape[0]
    # find number of reachable nodes in each graph
    obs_num_list=[((C@matrix_power(temp_matrix,n-1)).sum(axis = 0)>0).sum() for temp_matrix in A_list]
    return (np.array(obs_num_list)*w).sum()

