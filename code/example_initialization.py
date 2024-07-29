# initialize 22 link example

import numpy as np
from link_queue_model import fundamental_diagram, link


def initialize_22link():
    np.random.seed(10)
    T = 3600 # time steps
    l = 1  # length of road

    FD1 = fundamental_diagram(108, 18, 150)
    FD2 = fundamental_diagram(90, 12, 150)
    FD3 = fundamental_diagram(60, 10, 150)

    density_list = []

    inflow_list = []; outflow_list = []
    inflow1_list = np.ones(3600)
    inflow2_list = np.ones(3600)
    inflow3_list = np.ones(3600)

    n = 22
    Q = np.eye(n)*0.1
    R = np.eye(n)*0
    inflow = np.zeros(n)
    inflow[0:3] = 1
    for i in range(T):      
        inflow_list.append(inflow)

    outflow1 = np.zeros(n)
    outflow2 = np.zeros(n)
    outflow1[-3:] = 10
    outflow2[-3:] = 10
    outflow2[-1] = 0.2


    for i in range(T):      
        inflow_list.append(inflow)
        if i > 1000 and i<2000:
            outflow_list.append(outflow2)
        else:
            outflow_list.append(outflow1)

    link1 = link(0, l, FD1)
    link2 = link(0, l, FD1)
    link3 = link(0, l, FD1)
    link4 = link(0, l, FD2)
    link5 = link(0, l, FD1)
    link6 = link(0, l, FD2)
    link7 = link(0, l, FD1)
    link8 = link(0, l, FD2)
    link9 = link(0, l, FD1)
    link10 = link(0, l, FD1)
    link11 = link(0, l, FD1)
    link12 = link(0, l, FD1)
    link13 = link(0, l, FD1)
    link14 = link(0, l, FD1)
    link15 = link(0, l, FD2)
    link16 = link(0, l, FD1)
    link17 = link(0, l, FD2)
    link18 = link(0, l, FD1)
    link19 = link(0, l, FD2)
    link20 = link(0, l, FD1)
    link21 = link(0, l, FD1)
    link22 = link(0, l, FD1)

    link_list = [link1,link2,link3,link4, link5, link6,link7,link8,link9,link10,
                    link11,link12,link13,link14, link15, link16,link17,link18,link19,link20,link21,link22]


    network_geometry = ([[3, 9], [9, 13], [8, 12], [12, 18]],
    [[4, 5, 10], [6, 7, 11], [13, 14, 19], [15, 16, 20], [17, 18, 21]],
    [[0, 3, 4], [1, 5, 6], [2, 7, 8], [10, 14, 15], [11, 16, 17]])

    return inflow_list, outflow_list, link_list, network_geometry, Q



def initialize_22link_new():
    np.random.seed(10)
    T = 3600 # time steps
    l = 1  # length of road

    FD_list = []
    for i in range(22):
        FD_list.append(fundamental_diagram(70+2*i, (70+2*i)/4, 150))

    density_list = []

    inflow_list = []; outflow_list = []
    inflow1_list = np.ones(3600)
    inflow2_list = np.ones(3600)
    inflow3_list = np.ones(3600)

    n = 22
    Q = np.eye(n)*0.1
    R = np.eye(n)*0
    inflow = np.zeros(n)
    inflow[0:3] = 1
    for i in range(T):      
        inflow_list.append(inflow)

    outflow1 = np.zeros(n)
    outflow2 = np.zeros(n)
    outflow1[-3:] = 10
    outflow2[-3:] = 10
    outflow2[-1] = 0.2


    for i in range(T):      
        inflow_list.append(inflow)
        if i > 1000 and i<2000:
            outflow_list.append(outflow2)
        else:
            outflow_list.append(outflow1)

    link1 = link(0, l, FD_list[1])
    link2 = link(0, l, FD_list[2])
    link3 = link(0, l, FD_list[3])
    link4 = link(0, l, FD_list[4])
    link5 = link(0, l, FD_list[5])
    link6 = link(0, l, FD_list[6])
    link7 = link(0, l, FD_list[7])
    link8 = link(0, l, FD_list[8])
    link9 = link(0, l, FD_list[9])
    link10 = link(0, l, FD_list[10])
    link11 = link(0, l, FD_list[11])
    link12 = link(0, l, FD_list[12])
    link13 = link(0, l, FD_list[13])
    link14 = link(0, l, FD_list[14])
    link15 = link(0, l, FD_list[15])
    link16 = link(0, l, FD_list[16])
    link17 = link(0, l, FD_list[17])
    link18 = link(0, l, FD_list[18])
    link19 = link(0, l, FD_list[19])
    link20 = link(0, l, FD_list[20])
    link21 = link(0, l, FD_list[21])
    link22 = link(0, l, FD_list[22])

    link_list = [link1,link2,link3,link4, link5, link6,link7,link8,link9,link10,
                    link11,link12,link13,link14, link15, link16,link17,link18,link19,link20,link21,link22]


    network_geometry = ([[3, 9], [9, 13], [8, 12], [12, 18]],
    [[4, 5, 10], [6, 7, 11], [13, 14, 19], [15, 16, 20], [17, 18, 21]],
    [[0, 3, 4], [1, 5, 6], [2, 7, 8], [10, 14, 15], [11, 16, 17]])

    return inflow_list, outflow_list, link_list, network_geometry, Q





def initialize_22link_new():
    np.random.seed(10)
    T = 3600 # time steps
    l = 1  # length of road

    FD_list = []
    for i in range(22):
        FD_list.append(fundamental_diagram(70+2*i, (70+2*i)/4, 150))

    density_list = []

    inflow_list = []; outflow_list = []
    inflow1_list = np.ones(3600)
    inflow2_list = np.ones(3600)
    inflow3_list = np.ones(3600)

    n = 22
    Q = np.eye(n)*0.1
    R = np.eye(n)*0
    inflow = np.zeros(n)
    inflow[0:3] = 2
    for i in range(T):      
        inflow_list.append(inflow)

    outflow1 = np.zeros(n)
    outflow2 = np.zeros(n)
    outflow1[-3:] = 10
    outflow2[-3:] = 10
    outflow2[-1:] = 0.2


    for i in range(T):      
        inflow_list.append(inflow)
        if i > 1500 and i<2000:
            outflow_list.append(outflow2)
        else:
            outflow_list.append(outflow1)

    link1 = link(0, l, FD_list[21])
    link2 = link(0, l, FD_list[20])
    link3 = link(0, l, FD_list[2])
    link4 = link(0, l, FD_list[4])
    link5 = link(0, l, FD_list[5])
    link6 = link(0, l, FD_list[6])
    link7 = link(0, l, FD_list[7])
    link8 = link(0, l, FD_list[8])
    link9 = link(0, l, FD_list[9])
    link10 = link(0, l, FD_list[10])
    link11 = link(0, l, FD_list[11])
    link12 = link(0, l, FD_list[12])
    link13 = link(0, l, FD_list[13])
    link14 = link(0, l, FD_list[14])
    link15 = link(0, l, FD_list[15])
    link16 = link(0, l, FD_list[16])
    link17 = link(0, l, FD_list[17])
    link18 = link(0, l, FD_list[18])
    link19 = link(0, l, FD_list[19])
    link20 = link(0, l, FD_list[20])
    link21 = link(0, l, FD_list[21])
    link22 = link(0, l, FD_list[19])

    link_list = [link1,link2,link3,link4, link5, link6,link7,link8,link9,link10,
                    link11,link12,link13,link14, link15, link16,link17,link18,link19,link20,link21,link22]


    network_geometry = ([[3, 9], [9, 13], [8, 12], [12, 18]],
    [[4, 5, 10], [6, 7, 11], [13, 14, 19], [15, 16, 20], [17, 18, 21]],
    [[0, 3, 4], [1, 5, 6], [2, 7, 8], [10, 14, 15], [11, 16, 17]])

    return inflow_list, outflow_list, link_list, network_geometry, Q