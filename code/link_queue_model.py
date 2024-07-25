import numpy as np

class fundamental_diagram:
    # define the fundamental diagram for a link
    '''
    The fundamental diagram is a triangular fundamental diagram which takes in 3 parameters:
        1. vf (float, positive): free-flow travel speed (mph)
        2. w (float, positive): congestion wave speed (mph) (should be lower than vf, typically 1/6 - 1/3 of vf)
        3. kj (float, positive): jam density
    '''
    def __init__(self, vf, w, kj):
        self.vf = vf
        self.w = w
        self.kj = kj
        self.kc = kj*w/(vf+w)
        self.C = vf * self.kc
    def demand(self, k):
        # defines the demand function using LQM (CTM)
        return (min(self.vf * k, self.C)/3600)      # (veh/s)
    def supply(self, k):
        # defines the supply function using LQM (CTM)
        return (min(self.C, self.w*(self.kj-k))/3600)   #(veh/s)




class link:
    def __init__(self, density, length, FD):
        self.up = []
        self.k = density
        self.down = []
        self.length = length
        self.FD = FD
        self.C = FD.C
    
    def demand(self):
        # defines the demand function using LQM (CTM)
        k = self.k
        return (min(self.FD.vf * k, self.FD.C)/3600)      # (veh/s)
        
    def supply(self):
        # defines the supply function using LQM (CTM)
        return (min(self.FD.C, self.FD.w*(self.FD.kj-self.k))/3600)   #(veh/s)   
    
    def demand_supply_values(self):
        '''
        Computes the supply and demand of a link
        Note: the supply and demand is in veh/s
        '''
        return self.demand(), self.supply()




def compute_boundary_flow(i, inflow_list, outflow_list, link_list):
    '''
        x = Ax + b + u; this part computes u, the boundary flow
    '''
    # get the boundary flow first
    n = len(link_list)
    boundary_flow = np.zeros((n,1))
    for j in range(n):
        # find the boundary links
        if inflow_list[i][j] > 0: 
            d_j, s_j = link_list[j].demand_supply_values()
            boundary_flow[j][0] = min(s_j, inflow_list[i][j])
        if outflow_list[i][j] > 0: 
            d_j, s_j = link_list[j].demand_supply_values()
            boundary_flow[j][0] = -min(d_j, outflow_list[i][j])
    return boundary_flow




def check_occurrence(A, u, A_list, A_dict, u_list, u_dict, mode_list, mode_transition_num):
    '''
    Record distinct A, u, and record the occurrence for mode frequency tracking purpose
    '''
    check_flag = 0
    for j in range(len(A_list)):
        A_ = A_list[j]
        if np.all(A_ == A):
            check_flag = 1
            continue
    if check_flag == 0:
        A_list.append(A)
        u_list.append(u)
        A_dict[len(A_list)-1] = 1
        u_dict[len(u_list)-1] = 1
    else:
        A_dict[j] += 1
        u_dict[j] += 1
    # Find the index of A in A_list
    mode = next(i for i, x in enumerate(A_list) if np.array_equal(x, A))
    mode_list.append(mode)
    if mode >= mode_transition_num.shape[0]:
        # Add a new row of zeros
        new_row = np.zeros((1, mode_transition_num.shape[1]))
        mode_transition_num = np.append(mode_transition_num, new_row, axis=0)

        # Add a new column of zeros
        new_column = np.zeros((mode_transition_num.shape[0], 1))
        mode_transition_num = np.append(mode_transition_num, new_column, axis=1)
    
    if len(mode_list) > 1: 
        mode_transition_num[mode_list[-2], mode] += 1
    
    return A_list, A_dict, u_list, u_dict, mode_list, mode_transition_num




def LQM_matrix(inflow_list, outflow_list, network_geometry, link_list, Q, WITH_NOISE=False, alpha_merge=0.5, split_ratio=0.5, T=3600):

    def get_ordinary(A, u, l):
        '''
        From the supply and demand, compute the system update equation (matrix)
        x_(k+1) = A x_k + u_k
        Compute A and u_k
        l: indices of the 2 ordinary links ex. [1, 2]
        '''
        l1 = l[0]; l2 = l[1]
        link_1 = link_list[l1]
        link_2 = link_list[l2]  
        d1,s1 = link_1.demand_supply_values()
        d2,s2 = link_2.demand_supply_values()
        g1 = min(d1, s2)
        if g1 == d1:
            # demand of link1 is less than supply of link2
            if d1<link_1.C/3600:
                # demand of link1 is less than capacity (v*k)
                A[l1][l1] += -link_1.FD.vf/3600
                A[l2][l1] += link_1.FD.vf/3600
            else:
                # f2 = g1 = C
                u[l1, 0] += -link_1.C/3600
                u[l2, 0] += link_1.C/3600
        elif s2 < link_2.C/3600:
            A[l2][l2] += -link_2.FD.w /3600
            A[l1][l2] += link_2.FD.w/3600
            u[l1, 0] += -link_2.FD.w * link_2.FD.kj / 3600
            u[l2, 0] += link_2.FD.w * link_2.FD.kj / 3600     
        else:
            u[l1, 0] = u[l1, 0] -link_2.C/3600
            u[l2, 0] += link_2.C/3600
        return A, u


    def get_merge(A, u, l, alpha_merge = 0.5):
        '''
        From the supply and demand, compute the system update equation (matrix)
        x_(k+1) = A x_k + u_k
        Compute A and u_k
        l: indices of the merging junction links ex. [link1_id, link2_id, link3_id] >-  
        alpha_merge: merging priority, fair merge (0.5) by default
        '''
        l1 = l[0]; l2 = l[1]; l3 = l[2]
        link_1 = link_list[l1]
        link_2 = link_list[l2]
        link_3 = link_list[l3]
        
        d1,s1 = link_1.demand_supply_values()
        d2,s2 = link_2.demand_supply_values()
        d3,s3 = link_3.demand_supply_values()
        f3 = min(d1+d2, s3)
        g1 = min(d1, max(s3-d2, alpha_merge *s3))
        g2 = min(d2, max(s3-d1, (1-alpha_merge) *s3))
        
        # if 1
        if f3 == d1+d2:  # demand is less than supply   
            # if 2
            if d1<link_1.C/3600:  # d1 = k1 v1 (will not involve u)
                A[l3][l1] += link_1.FD.vf/3600
                A[l1][l1] += -link_1.FD.vf/3600
            else: 
                # d1 = C1, need to update u 
                u[l1, 0] += -link_1.C/3600
                u[l3, 0] += link_1.C/3600
            if d2<link_2.C/3600: 
                A[l3][l2] += link_2.FD.vf/3600
                A[l2][l2] += -link_2.FD.vf/3600
            else:
                # d2 = C2, need to update u
                u[l2, 0] += -link_2.C/3600
                u[l3, 0] += link_2.C/3600
        else:  # d1+d2>s3
            # if 1
            if s3<link_3.C/3600:
                # s3 = w3*(kj - k3)
                A[l3][l3] += -link_3.FD.w/3600
                u[l3, 0] += link_3.FD.w*link_3.FD.kj/3600
                
                # if 2
                if g1 == d1:  
                    if d1<link_1.C/3600: # g1=k1*v1, g2=s3-d1
                        A[l2][l1] += link_1.FD.vf/3600
                        A[l1][l1] += -link_1.FD.vf/3600
                        A[l2][l3] += link_3.FD.w/3600
                        u[l2, 0] += -link_3.FD.w * link_3.FD.kj /3600
                    else:   #g1 = C1, g2 = s3-C1
                        u[l1, 0] += -link_1.C/3600
                        u[l2, 0] += link_1.C/3600 - link_3.FD.w*link_3.FD.kj/3600
                        A[l2][l3] += link_3.FD.w/3600
                        
                elif g2 == d2:  
                    if d2<link_2.C/3600: # g2=k2*v2, g1=s3-d2
                        A[l1][l2] += link_2.FD.vf/3600
                        A[l2][l2] += -link_2.FD.vf/3600
                        A[l1][l3] += link_3.FD.w/3600
                        u[l1, 0] += -link_3.FD.w * link_3.FD.kj /3600
                    else:   #g1 = C1, g2 = s3-C1
                        u[l2, 0] += -link_2.C/3600
                        u[l1, 0] += link_2.C/3600 - link_3.FD.w*link_3.FD.kj/3600
                        A[l1][l3] += link_3.FD.w/3600    
                
                else:  # g1=g2=1/2 s3 for fair merge
                    A[l1][l3] += alpha_merge * link_3.FD.w/3600
                    A[l2][l3] += (1-alpha_merge)*link_3.FD.w/3600
                    u[l1, 0] += -alpha_merge * link_3.FD.w * link_3.FD.kj / 3600
                    u[l2, 0] += -(1-alpha_merge) * link_3.FD.w * link_3.FD.kj / 3600    

            else: 
                # s3>=link_3.C/3600
                u[l3, 0] += link_3.C/3600
                
                if g1 == d1:  
                    if d1<link_1.C/3600: # g1=k1*v1, g2=s3-d1
                        A[l2][l1] += link_1.FD.vf/3600
                        A[l1][l1] += -link_1.FD.vf/3600
                        u[l2, 0] += -link_3.C/3600
                    else:   #g1 = C1, g2 = s3-C1
                        u[l1, 0] += -link_1.C/3600
                        u[l2, 0] += link_1.C/3600 - link_3.C/3600
                        
                elif g2 == d2:  
                    if d2<link_2.C/3600: # g2=k2*v2, g1=s3-d2
                        A[l1][l2] += link_2.FD.vf/3600
                        A[l2][l2] += -link_2.FD.vf/3600
                        u[l1, 0] += -link_3.C/3600
                    else:   #g1 = C1, g2 = s3-C1
                        u[l2, 0] += -link_2.C/3600
                        u[l1, 0] += link_2.C/3600 -  link_3.C/3600
                
                else:  # g1=g2=1/2 s3 for fair merge
                    u[l1, 0] += -alpha_merge *  link_3.C / 3600
                    u[l2, 0] += -(1-alpha_merge) *  link_3.C / 3600 
                
        return A, u      


    def get_diverge(A,u, l, r=0.5):
        # r is the split ratio from link1 to link2, 0.5: 50% of vehicles travel to link2 from link1
        l1 = l[0]; l2 = l[1]; l3 = l[2]; r12 = r; r13 = 1-r
        link_1 = link_list[l1]
        link_2 = link_list[l2]
        link_3 = link_list[l3]
        d1,s1 = link_1.demand_supply_values()
        d2,s2 = link_2.demand_supply_values()
        d3,s3 = link_3.demand_supply_values()
        g1 = min(d1, s2/r12, s3/r13)
        f2 = r12*g1
        f3 = r13*g1
        
        if g1 == d1:  #all flow can get to downstreams
            if d1<link_1.C/3600:
                A[l1][l1] += -link_1.FD.vf /3600
                A[l2][l1] += link_1.FD.vf * r12 /3600
                A[l3][l1] += link_1.FD.vf * r13 /3600
            else:
                u[l1, 0] += -link_1.C/3600
                u[l2, 0] += link_1.C*r12/3600
                u[l3, 0] += link_1.C*r13/3600
        elif g1 == s2/r12:  #link2 is the bottleneck
            if s2<link_2.C/3600:
                A[l1][l2] += 1/r12 * link_2.FD.w /3600; u[l1, 0] += -1/r12 * link_2.FD.w * link_2.FD.kj  /3600
                A[l2][l2] += -link_2.FD.w/3600;  u[l2, 0] += link_2.FD.w * link_2.FD.kj  /3600
                A[l3][l2] += -r13/r12*link_2.FD.w/3600; u[l3, 0] +=r13/r12 * link_2.FD.w * link_2.FD.kj  /3600
            else:
                u[l1, 0] += -1/r12 * link_2.C/3600
                u[l2, 0] += link_2.C/3600
                u[l3, 0] += r13/r12*link_2.C/3600
                
        else: #link3 is the bottleneck
            if s3<link_3.C/3600:
                A[l1][l3] += 1/r13 * link_3.FD.w /3600; u[l1, 0] += -1/r13 * link_3.FD.w * link_3.FD.kj /3600
                A[l3][l3] += -link_3.FD.w/3600;  u[l3, 0] += link_3.FD.w * link_3.FD.kj  /3600
                A[l2][l3] += -r12/r13*link_3.FD.w/3600; u[l2, 0] +=r13/r12 * link_3.FD.w * link_3.FD.kj /3600
            else: 
                u[l1, 0] += -1/r13 * link_3.C/3600
                u[l2, 0] += link_3.C/3600
                u[l3, 0] += r12/r13*link_3.C/3600
        return A, u


    '''
    boundary_flow: np.array (dim = n) a vector of n dimensions with non-zero elements at the boundary locations
                  defines the flow (veh/s) at each boundary link
    network_geometry: list of list of list, defines the network geometry including ordinary, merge and diverge junctions
                    ex. network >--< [[[3,4]], [[1,2,3]], [[4,5,6]]]
    '''

    n = len(link_list)
    link_density = np.array([link_list[i].k for i in range(n)]).reshape(n, 1)
    density_list = list(); boundary_list = list()
    A_list = list(); u_list = list()
    A_dict = {}; u_dict = {}
    mode_list = list()
    mode_transition_num = np.zeros((1, 1)).astype(int)
    
    for i in range(T):    
        # update mode matrix 
        A = np.eye(n)
        # update the Bu term 
        u = np.zeros((n,1))
        # compute noise
        
        noise = (np.random.multivariate_normal(np.zeros(n), Q)).reshape(n,1)
        
        # compute boundary flow 
        boundary_flow = compute_boundary_flow(i, inflow_list, outflow_list, link_list)
        boundary_list.append(boundary_flow)       
        
        # update the matrix elements based on all junction types
        for ordinary_junction_index in network_geometry[0]:
            A, u = get_ordinary(A, u, ordinary_junction_index)  
            
        for merge_junction_index in network_geometry[1]:
            A, u = get_merge(A, u, merge_junction_index, alpha_merge)  
        
        for diverge_junction_index in network_geometry[2]:
            A, u = get_diverge(A, u, diverge_junction_index, split_ratio)  
            
        
        if WITH_NOISE:
            link_density = A @ link_density + u + boundary_flow + noise
        else:
            link_density = A @ link_density + u + boundary_flow
        # update the link_list
        for j in range(n):
            link_list[j].k = link_density[j, 0]
        # record the results
        density_list.append(link_density)
        
        A_list, A_dict, u_list, u_dict, mode_list, mode_transition_num = check_occurrence(A, u, A_list, A_dict, u_list, u_dict, mode_list, mode_transition_num)
    return density_list, boundary_list, A_list, u_list, mode_list, mode_transition_num, A_dict