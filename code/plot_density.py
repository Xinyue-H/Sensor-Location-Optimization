
import matplotlib.pyplot as plt
import numpy as np

def plot_density(density_list, mode_list):
    plt.rc('font', size=20)
    T = len(density_list)
    n = len(density_list[0])
    t_list = list(range(T))
    
#     for j in range(n):
#         k_list = [i[j][0] for i in density_list]
#         plt.plot(t_list, k_list, label ="k" +str(j)+"(veh/mi)")
    
    # Create the plot

    fig, ax = plt.subplots(figsize=(10, 8)) 
    
    for j in range(n):
        k_list = [i[j][0] for i in density_list]   
        ax.plot(t_list, k_list, label='k'+str(j+1)+" (veh/h)")
    

#     # Fill the background between x=0 and x=1 with red
#     ax.axvspan(0, 1, alpha=0.1, color='red')

#     # Fill the background between x=1 and x=3 with blue
#     ax.axvspan(1, 3, alpha=0.3, color='blue')

    # Define a colormap with different colors for each index
    cmap = plt.cm.get_cmap('viridis', max(mode_list)+1)
    
    start_index = 0
    for end_index in range(1, T):
        if mode_list[end_index]!=mode_list[start_index]:
            ax.axvspan(start_index, end_index-1, alpha=0.3, color=cmap(mode_list[start_index]))
            start_index = end_index
    ax.axvspan(start_index, end_index, alpha=0.3, color=cmap(mode_list[start_index]))         
    
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    #plt.legend(loc="upper left")
    plt.show()




def plot_2densities(density_list, x_list, vect = None):
    plt.rc('font', size=20)
    T = len(density_list)
    t_list = list(range(T))
    
#     for j in range(n):
#         k_list = [i[j][0] for i in density_list]
#         plt.plot(t_list, k_list, label ="k" +str(j)+"(veh/mi)")
    
    # Create the plot

    fig, ax = plt.subplots(figsize=(10, 8)) 
    
    color_ind = 0
    n = len(density_list[0])
    for j in range(n):
        k1_list = [i[j][0] for i in density_list]   
        k2_list = [i[j][0] for i in x_list]   
        
        if vect is None:
            color = plt.cm.plasma(j / n)
            ax.plot(t_list, k1_list, label='k'+str(j+1)+" (veh/h)", linestyle='-', color=color)
            ax.plot(t_list, k2_list, label='k'+str(j+1)+" (veh/h) predicted", linestyle='--', color=color)
        elif vect[j]==0:
            color = plt.cm.plasma(color_ind / len(np.where(vect == 0)[0]))
            color_ind += 1
            ax.plot(t_list, k1_list, label='k'+str(j+1)+" (veh/h)", linestyle='-', color=color)
            ax.plot(t_list, k2_list, label='k'+str(j+1)+" (veh/h) predicted", linestyle='--', color=color)
            
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Density (veh/mi)')
    # Put the legend outside the plot
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))

    # Define a colormap with different colors for each index
    cmap = plt.cm.get_cmap('viridis', max(mode_list)+1)
    
    start_index = 0
    for end_index in range(1, T):
        if mode_list[end_index]!=mode_list[start_index]:
            ax.axvspan(start_index, end_index-1, alpha=0.3, color=cmap(mode_list[start_index]))
            start_index = end_index
    ax.axvspan(start_index, end_index, alpha=0.3, color=cmap(mode_list[start_index]))         

    plt.show()




def plot_Cov(Cov_list, mode_list):
    plt.rc('font', size=20)
    T = len(mode_list)
    t_list = list(range(T))
    
#     for j in range(n):
#         k_list = [i[j][0] for i in density_list]
#         plt.plot(t_list, k_list, label ="k" +str(j)+"(veh/mi)")
    
    # Create the plot

    fig, ax = plt.subplots(figsize=(10, 8)) 
    Cov_norms = []
    for i in range(T):
        Cov_norms.append(np.linalg.norm(Cov_list[i]))
    
    ax.plot(t_list, Cov_norms, label='norm')
    

#     # Fill the background between x=0 and x=1 with red
#     ax.axvspan(0, 1, alpha=0.1, color='red')

#     # Fill the background between x=1 and x=3 with blue
#     ax.axvspan(1, 3, alpha=0.3, color='blue')

    # Define a colormap with different colors for each index
    cmap = plt.cm.get_cmap('viridis', max(mode_list)+1)
    
    start_index = 0
    for end_index in range(1, T):
        if mode_list[end_index]!=mode_list[start_index]:
            ax.axvspan(start_index, end_index-1, alpha=0.3, color=cmap(mode_list[start_index]))
            start_index = end_index
    ax.axvspan(start_index, end_index, alpha=0.3, color=cmap(mode_list[start_index]))         

    plt.legend(loc="upper left")
    plt.show()