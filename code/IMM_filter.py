import numpy as np

def compute_PI(A_list, mode_transition_num):
    N_m = len(A_list) # number modes
    # initialization
    PI = (mode_transition_num+np.eye(N_m)) / np.sum(mode_transition_num+np.eye(N_m), axis=1, keepdims=True)
    for i in range(N_m):
        for j in range(N_m):
            if PI[i, j]<0.005:
                PI[j, j] = PI[j, j] - (0.005-PI[i, j]) 
                PI[i, j] = 0.005
    return PI


def IMM_filter(A_list, u_list, boundary_list, PI, C, Q, y_k):
    '''
    This function will implement the Interacting Multiple Model to estimate the traffic state
    It takes in:
        A_list: a list of A 
    '''
    from scipy.stats import multivariate_normal
    T = len(boundary_list)
    x_list = []; mu_list = []; Cov_list = []
    # calculate the mixing probabilities 
    N_m = len(A_list) # number of modes
    # intialize same probability for all modes
    mu_k = np.ones(N_m)/N_m  

    p, n = C.shape # number of sensors, number of states

    initial_mean = np.zeros(n).reshape(n, 1)
    initial_cov = np.eye(n)

    # get matrix B, all set to identity 
    B = np.eye(n)
    D = np.eye(n)
    R = np.eye(p)*0

    x_k = list(); Cov_k = list()
    y_prediction = [None] * N_m
    x_0i = [None] * N_m; Cov_0i = [None] * N_m
    x_prediction = [None] * N_m
    Cov_prediction = [None] * N_m

    for i in range(N_m):
        x_k.append(initial_mean)
        Cov_k.append(initial_cov)
        
    # filtering starts
    for k in range(T):
        # mixing step
        ### calculate the mixing probabilities mu_mixing
        mu_mixing = np.zeros((N_m, N_m)) 
        for i in range(N_m):
            mu_i_sum = sum(PI[:, i] * mu_k)
            for j in range(N_m):
                mu_mixing[j, i] = PI[j, i] * mu_k[j] / mu_i_sum
        ### calculate the mixed estimates and covariances

        for i in range(N_m):
            x_0i[i] = 0; Cov_0i[i] = 0
            for j in range(N_m):
                x_0i[i] += mu_mixing[j, i] * x_k[j]
            for j in range(N_m):
                Cov_0i[i] += mu_mixing[j, i] * (Cov_k[j] + (x_k[j] - x_0i[i])@(x_k[j] - x_0i[i]).T)

        # updating model likelihood 
        likelihood_list = []
        temp_sum = []
        # Mode matched prediction update
        for i in range(N_m):
            x_prediction[i] = A_list[i] @ x_0i[i] + u_list[i] + boundary_list[k]

            #Cov_prediction[i] = A_list[i] @ Cov_0i[i] @ A_list[i].T + B @ Q @ B.T
            Cov_prediction[i] = A_list[i] @ Cov_0i[i] @ A_list[i].T + Q
            
            S = C@Cov_prediction[i]@C.T + R
            K = Cov_prediction[i] @ C.T @ np.linalg.inv(S)
            x_k[i] = x_prediction[i] + K@(y_k[k] - C@x_prediction[i])
            Cov_k[i] = Cov_prediction[i] - K@C@Cov_prediction[i]
            y_prediction[i] = C@x_prediction[i]
            
            try:    
                mvn = multivariate_normal(mean=y_prediction[i].reshape(p,), cov=S)
            except:
                min_eigenval = min(np.linalg.eigvals(S))
                S = S-min_eigenval*np.eye(p)+np.eye(p)
                mvn = multivariate_normal(mean=y_prediction[i].reshape(p,), cov=S)
                
#             mvn = multivariate_normal(mean=y_prediction[i].reshape(p,), cov=S)     
            likelihood_list.append(mvn.pdf(y_k[k].reshape(p,)))


        for i in range(N_m):
#             mvn = multivariate_normal(mean=y_prediction[i].reshape(p,), cov=S)
#             likelihood_list.append(mvn.pdf(y_k[k].reshape(p,)))
            temp_s = 0
            for j in range(N_m):
                temp_s += PI[j, i] * mu_k[j]
            temp_sum.append(temp_s)

        likelihood_list = np.array(likelihood_list); temp_sum = np.array(temp_sum)
        mu_k = likelihood_list*temp_sum / sum(likelihood_list*temp_sum)
        
        # record
        x_posterior = 0; Cov_posterior = 0
        for i in range(N_m):
            x_posterior += mu_k[i] * x_k[i]
        x_list.append(x_posterior)
        mu_list.append(mu_k)
        
    return x_list, mu_list