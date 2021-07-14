# Copyright 1996-2015 PSERC. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

"""Power flow data for 6 bus, 3 gen case from Wood & Wollenberg.
"""

from numpy import int64, array, zeros, ones, arange, ix_, r_, flatnonzero as find, dot as mult
from scipy.sparse import csr_matrix as sparse
from numpy.linalg import inv

def case6ww_GCEP():
    """Power flow data for 6 bus, 3 gen case from Wood & Wollenberg.
    Please see L{caseformat} for details on the case file format.

    This is the 6 bus example from pp. 104, 112, 119, 123-124, 549 of
    I{"Power Generation, Operation, and Control, 2nd Edition"},
    by Allen. J. Wood and Bruce F. Wollenberg, John Wiley & Sons, NY, Jan 1996.

    @return: Power flow data for 6 bus, 3 gen case from Wood & Wollenberg.
    """
    ppc = {"version": '2'}

    ##-----  Power Flow Data  -----##
    ## system MVA base
    ppc["baseMVA"] = 100.0
    
    ## bus and load data
    
    ppc["lf"] = array([
        [450,  500,  550,  600],    
        [550,  600,  650,  700],
        [650,  700,  750,  800],
        [750,  800,  850,  900],
        [850,  900,  950,  1000]]);
    
    ppc["num_h"] = array([2920, 1460, 2555, 1825])  
                       
    # bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin
    ppc["bus"] = array([
        [1, 3,  0,  0, 0, 0, 1, 1.05, 0, 230, 1, 1.05, 1.05],
        [2, 2,  0,  0, 0, 0, 1, 1.05, 0, 230, 1, 1.05, 1.05],
        [3, 2,  0,  0, 0, 0, 1, 1.07, 0, 230, 1, 1.07, 1.07],
        [4, 1, (600/210)*70, 70, 0, 0, 1, 1,    0, 230, 1, 1.05, 0.95],
        [5, 1, (600/210)*70, 70, 0, 0, 1, 1,    0, 230, 1, 1.05, 0.95],
        [6, 1, (600/210)*70, 70, 0, 0, 1, 1,    0, 230, 1, 1.05, 0.95]
    ])

    ## generator data
    # bus, Pg, Qg, Qmax, Qmin, Vg, mBase, status, Pmax, Pmin, Pc1, Pc2,
    # Qc1min, Qc1max, Qc2min, Qc2max, ramp_agc, ramp_10, ramp_30, ramp_q, apf
    # ppc["gen"] = array([
    #     [1, 0, 0, 100, -100, 1.05, 100, 1, 200, 50,   0, 0, 0, 0, 0, 0, 0, 0, 0, 9*10,  8.5*10],
    #     [2, 0, 0, 100, -100, 1.05, 100, 1, 150, 37.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12*10, 12*10],
    #     [3, 0, 0, 100, -100, 1.07, 100, 1, 180, 45,   0, 0, 0, 0, 0, 0, 0, 0, 0, 11*10, 10.1*10]
    # ])

    ## branch data
    # fbus, tbus, r, x, b, rateA, rateB, rateC, ratio, angle, status, angmin, angmax
    ppc["branch"] = array([                       #TS           #SC
        [1, 2, 0.1,  0.2,  0.04, 300, 40, 40, 0, 0, 0, -360, 360, 1],
        [1, 4, 0.05, 0.2,  0.04, 300, 60, 60, 0, 0, 0, -360, 360, 1], #50
        [1, 5, 0.08, 0.3,  0.06, 300, 40, 40, 0, 0, 0, -360, 360, 1],
        [2, 3, 0.05, 0.25, 0.06, 300, 40, 40, 0, 0, 0, -360, 360, 1],
        [2, 4, 0.05, 0.1,  0.02, 300, 60, 60, 0, 0, 0, -360, 360, 1], 
        [2, 5, 0.1,  0.3,  0.04, 300, 30, 30, 0, 0, 0, -360, 360, 1],
        [2, 6, 0.07, 0.2,  0.05, 300, 90, 90, 0, 0, 0, -360, 360, 1], 
        [3, 5, 0.12, 0.26, 0.05, 300, 70, 70, 0, 0, 0, -360, 360, 1], #50
        [3, 6, 0.02, 0.1,  0.02, 300, 80, 80, 0, 0, 0, -360, 360, 1], 
        [4, 5, 0.2,  0.4,  0.08, 300, 20, 20, 0, 0, 0, -360, 360, 1],
        [5, 6, 0.1,  0.3,  0.06, 150, 40, 40, 0, 0, 0, -360, 360, 1]
    ])

    ##-----  OPF Data  -----##
    ## generator cost data
    # 1 startup shutdown n x1 y1 ... xn yn
    # 2 startup shutdown n c(n-1) ... c0
    # ppc["gencost"] = array([
    #     [2, 0, 0, 3, 0.00533, 11.669, 213.1],
    #     [2, 0, 0, 3, 0.00889, 10.333, 200],
    #     [2, 0, 0, 3, 0.00741, 10.833, 240]
    # ])      

    # Power generation data
             # bus c b a p_max     p_tin  Min_up  Min_do  H_cost C_cost C_time I_state Ramp-up and -down limit C_inv C_O&M/ENS Inv_dec Pg_0
    ppc["units"] = array([
            [1, 1000, 18.50, 0.00048, 250,  30, 8, 8, 4500,  9000, 5,  8, 225,     0,     0, 0, 0],   # Gen1 1
            [4,  970, 14.36, 0.00031, 100,  20, 8, 8, 5000, 10000, 5,  8, 225,     0,     0, 0, 0],   # Gen4 2
            [4,  700, 22.11, 0.002  ,  50,   5, 5, 5,  550,  1100, 4, -5, 50,      0,     0, 0, 0],   # Gen4 3
            [1,  680, 20.41, 0.00211, 150,  20, 5, 5,  560,  1120, 4, -5, 50, 300000, 12000, 1, 0],   # Gen1 1
            [1,  680, 20.41, 0.00211, 150,  20, 5, 5,  560,  1120, 4, -5, 50, 300000, 12000, 1, 0],   # Gen1 2
            [1,  680, 20.41, 0.00211, 150,  20, 5, 5,  560,  1120, 4, -5, 50, 300000, 12000, 1, 0],   # Gen1 3
            [3,  450, 25.95, 0.00398, 100,  15, 6, 6,  900,  1800, 4, -6, 60, 250000, 30000, 1, 0],   # Gen3 4
            [3,  450, 25.95, 0.00398, 100,  15, 6, 6,  900,  1800, 4, -6, 60, 250000, 30000, 1, 0],   # Gen3 5
            [3,  450, 25.95, 0.00398, 100,  15, 6, 6,  900,  1800, 4, -6, 60, 250000, 30000, 1, 0],   # Gen3 6
            [5,  370, 14.08, 0.00712, 250,  10, 3, 3,  170,   340, 2, -3, 60, 350000, 36000, 1, 0],   # Gen5 7
            [5,  370, 14.08, 0.00712, 250,  10, 3, 3,  170,   340, 2, -3, 60, 350000, 36000, 1, 0],   # Gen5 8
            [5,  370, 14.08, 0.00712, 250,  10, 3, 3,  170,   340, 2, -3, 60, 350000, 36000, 1, 0],   # Gen5 9
            [5,  370, 14.08, 0.00712, 250,  10, 3, 3,  170,   340, 2, -3, 60, 350000, 36000, 1, 0],   # Gen5 10
            [4,  480, 10000, 0.00079, 999,   0, 3, 3,  260,   520, 2, -3, 60,      0,     0, 2, 0],   # Gen1 (Lost of Load)
            [5,  660, 10000, 0.00413, 999,   0, 1, 1,   30,    60, 0, -1, 135,     0,     0, 2, 0],   # Gen2
            [6,  665, 10000, 0.00222, 999,   0, 1, 1,   30,    60, 0, -1, 135,     0,     0, 2, 0]
            ])   # Gen3   
           # 0     1      2        3    4    5  6  7     8      9 10  11   12 13 14 15 16

    ng = len(ppc["units"])
    nb = len(ppc["bus"])
    nl = len(ppc["branch"])
    # Load
    # aux = sum(ppc["bus"][:,2]) # Load is divided proportionally
    for b in range(nb):
        if ppc["bus"][b,2] < 0:
            ppc["bus"][b,2] = 0
        # if ppc["bus"][b,2] > 0:  
        #     ppc["bus"][b,2]  = ppc["bus"][b,2] / aux
    
    # Generation
    pos_g = (ppc["units"][:,0]-1).astype(int64) # gen location (astype... for integer positions in the SF matrix)
    ppc["pos_g"] = pos_g
    Cg = sparse((ones(ng), (pos_g, range(ng))), (nb, ng)) #conection gen matrix
    ppc["Cg"] = Cg
    
    # Transmission
    b = 1 / ppc["branch"][:,3]
    f = ppc["branch"][:, 0]-1
    t = ppc["branch"][:, 1]-1
    I = r_[range(nl), range(nl)]
    S = sparse((r_[ones(nl), -ones(nl)], (I, r_[f, t])), (nl, nb))#total
    Bf = sparse((r_[b, -b], (I, r_[f, t])), (nl,nb))
    Bbus = S.T * Bf
    slack_bus=find(ppc["bus"][:,1]==3)
    ppc["SL"] = slack_bus
    buses = arange(1, nb)
    noslack = find(arange(nb) != slack_bus)
    SF = zeros((nl, nb))
    SF[:,noslack]=Bf[:, buses].todense()*inv(Bbus[ix_(noslack, buses)].todense()) 
    ppc["SF"] = SF

    # new transmission candidates
    pos_new_L = find(ppc["branch"][:,13]==0) # new
    f = ppc["branch"][pos_new_L, 0]-1
    t = ppc["branch"][pos_new_L, 1]-1
    In = r_[range(len(pos_new_L)), range(len(pos_new_L))]
    ppc["Cf"] = sparse((r_[ones(len(pos_new_L)), -ones(len(pos_new_L))], (In, r_[f, t])), (len(pos_new_L), nb)) # S new lines
    
    return ppc
