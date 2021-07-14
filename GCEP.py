# V. Hinojosa. Date: 15/01/2021.
from gurobipy import *
from numpy import zeros, savetxt, absolute, dot as mult, size, flatnonzero as find
import time

t00 = time.time() #formulation time   

import case6ww_GCEP as mpc
sep = mpc.case6ww_GCEP()

# Power system data & parameters
int_rate = 0.1 # interest rate
R = 0.2 # Reserve criterium

# Load modeling  
T = 5 # long-term horizon
factor = 0.1666 # load increasing (1/3)
Load = sep["bus"][:,2] # initial load data
ldc = [1.0, 0.92, 0.83, 0.75] #sep["p_load"] # LDC curve (peak 1, peak 2, mean, rest, valley)
hs = [1825, 2555, 1460, 2920] #sep["num_h"] # LDC time
Load_M = find(max(ldc))[0] # maximum load
nh = len(hs) # number of time-steps per year

# tranmission power flows
transmission = 0
nb = len(sep['bus']); # number of buses
nl = len(sep['branch']); # number of transmission lines

# Generation planning information
Gen = sep["units"] # generation data
ng = len(Gen) # number of total power units
PM = Gen[:,4] # maximum power unit
# Pm = Gen[:,5] # minimum power unit
pos_exi_G = find(Gen[:,15]==0) # existing power units position
pos_new_G = find(Gen[:,15]==1) # new power units position
pos_ens = find(Gen[:,15]==2) # new power units position
N_e_G = len(pos_exi_G) # power generation existing units
N_x_G = len(pos_new_G) # power generation candidates
CI_G = Gen[pos_new_G,13] * PM[pos_new_G] # Investment costs
CM_G = Gen[pos_new_G,14] * PM[pos_new_G] # O&M costs


""" #########################     Mathematical formulation     ########################################################################################## """
# Initializing m
model = Model('GCEP')
model.Params.MIPGap = 1e-6
model.Params.OutputFlag = 0 #m.setParam('OutputFlag', False) 

# VARIABLE DEFINITIONS
p = model.addVars(ng, nh, T, vtype=GRB.CONTINUOUS, lb=0, name='p') # Power output of unit j in period t. # these variables are used like p[tt,ii,jj]
n_G = model.addVars(N_x_G, T, vtype=GRB.BINARY, name='new_p') # Power investment decision

# OPTIMIZATION PROBLEM
""" Investment and operational costs """  
CInv_G = 0; COM_G = 0; COp = 0;
for t in range(T):
    f_anual = 1/pow(1+int_rate,t)
    if t == 0:
        CInv_G += f_anual * (quicksum(CI_G * n_G.select('*',t)))
        COM_G += f_anual * (quicksum(CM_G * n_G.select('*',t)))
    else:
        CInv_G += f_anual * (quicksum(CI_G * n_G.select('*',t)) - quicksum(CI_G * n_G.select('*',t-1)))
        COM_G += f_anual * (quicksum(CM_G * n_G.select('*',t)))       
    for h in range(nh):
        COp += f_anual * quicksum(Gen[:,2] * hs[h] * p.select('*',h,t))

""" Objective function """
model.setObjective((CInv_G + COM_G + COp)/1e6, GRB.MINIMIZE) 
 
""" Constraints (s.t.) """   

# Reserve (1)
# for t in range(T):
#     model.addConstr(quicksum(PM[pos_new_G] * n_G.select('*',t)) + sum(PM[pos_exi_G]) >= (1+R) * (1+t*factor) * ldc[Load_M] * sum(Load), 'Res[T=%d]' % (t+1)) 

# Balance (2)
for t in range(T):
    for h in range(nh):
        model.addConstr(quicksum(p.select('*',h,t)) == (1+t*factor) * ldc[h] * sum(Load), 'Bal[h=%d,T=%d]' % (h+1,t+1)) 
       
# New power units disjunctive constraints (3)
for t in range(T):
    for h in range(nh):
        for g in range(N_x_G):
            model.addConstr(p[pos_new_G[g],h,t] <= PM[pos_new_G[g]] * n_G[g,t], 'p_new[n_G=%d,t=%d,T=%d]' % (pos_new_G[g],h+1,t+1))       
       
# Sequential instalation (4)
for t in range(T):
    if t>0:
        for g in range(N_x_G):
            model.addConstr(n_G[g,t] >= n_G[g,t-1], 'Seq[%d,%d]' % (pos_new_G[g],t+1))       
   
# Maximum power generation (5) => existing
for t in range(T):
    for h in range(nh):
        for g in range(N_e_G):
            model.addConstr(p[g,h,t] <= PM[g], 'p_M[%d,%d]' % (g,t+1))
            
# Transmission limits (6)
if transmission:    
    for t in range(T):
        for h in range(nh):
            sep["bus"][:,2] = (1+t*factor) * ldc[h] * sep["bus"][:,2]
            PF_sl=mult(sep['SF'],sep["bus"][:,2]) #Slack PF
            rhs_1=-sep["branch"][:,5]-PF_sl #FM - Slack PF
            rhs_2=-sep["branch"][:,5]+PF_sl        
            for l in range(nl):
                expr = quicksum(p.select('*',h,t) * sep['SF'][l,sep['pos_g']])
                txt = str(int(sep["branch"][l,0])) + str('-') + str(int(sep["branch"][l,1])) + str(',') + str(int(h)) + str(',') + str(int(t))
                model.addConstr(-expr >= rhs_1[l], 'FM[%s]' % txt) 
                model.addConstr(expr >= rhs_2[l], 'Fm[%s]' % txt) 
            sep["bus"][:,2] = sep["bus"][:,2] / ((1+t*factor) * ldc[h])      

t11 = time.time() #formulation time

if 1: # write LP
    model.write('GCEP.lp'); 
    # sys.exit(1);  

# SOLVER & INFO
model.optimize()
  
status = model.Status
if status == GRB.Status.OPTIMAL:
    print ('\nCost = %.2f (MM$) => CI_G = %.2f (MM$) + COM_G = %.2f (MM$) + COp = %.2f (MM$)' % (model.objVal,CInv_G.getValue()/1e6,COM_G.getValue()/1e6,COp.getValue()/1e6))
    print('num_Vars =  %d / num_Const =  %d / num_NonZeros =  %d' % (model.NumVars,model.NumConstrs,model.DNumNZs)) 
    if 1:
        print ('\nGeneration planning solution:')
        sol_n_G = model.getAttr('x',n_G); 
        for t in range(T):
            for l in range(N_x_G): 
                if t == 0:
                    aux = sol_n_G[l,t]
                else:
                    aux = sol_n_G[l,t] - sol_n_G[l,t-1]
                if aux > 0:
                    print('T%d: G_bus_%d' % (t+1,Gen[pos_new_G[l],0])) 

        if 0:
            vX = zeros([T+1,N_x_G]) # investment decisions          
            vX[0] = Gen[pos_new_G,0]
            for t in range(T):
                vX[t+1] = sol_v.select('*',t)        
            # savetxt('GCEP.txt', vX, delimiter=",", fmt='%d')
        if 0:
            print ('\nPower generation solution:')
            sol_p = model.getAttr('x',p); 
            pX = zeros([T,nh,ng]) 
            for t in range(T):
                for h in range(nh):
                    print('\nyear = %d / period = %d' % (t+1,h+1))
                    pX[t,h] = sol_p.select('*',h,t)
                    for g in range(ng):
                        if pX[t,h,g] > 0:
                            if Gen[g,15] == 2:
                                print('\nENS[bus_%.0f] = %.2f (MW)' % (Gen[g,0], pX[t,h,g]))
                            else:
                                if Gen[g,15] == 0:
                                    print('Pg[bus_%.0f] = %.2f (MW)' % (Gen[g,0], pX[t,h,g]))
                                else:
                                    print('Pg_new[bus_%.0f] = %.2f (MW)' % (Gen[g,0], pX[t,h,g]))
        if transmission * 0:    
            print ('\nTransmission power flows:')
            F = zeros([T,nh,nl]) 
            for t in range(T):
                for h in range(nh):
                    print('\nyear = %d / period = %d' % (t+1,h+1))
                    sep["bus"][:,2] = (1+t*factor) * ldc[h] * sep["bus"][:,2]
                    F[t,h]=mult(sep['SF'], sep["Cg"]*sol_p.select('*',h,t)-sep['bus'][:,2])
                    for l in range(nl):
                        print('f[%.0f-%.0f] = %.3f (MW)' % (sep["branch"][l,0], sep["branch"][l,1], F[t,h,l]))
                    sep["bus"][:,2] = sep["bus"][:,2] / ((1+t*factor) * ldc[h])   
    print('=> Formulation time: %.4f (s)'% (t11-t00))
    print('=> Solver time: %.4f (s)' % (model.Runtime))  
elif status == GRB.Status.INF_OR_UNBD or \
   status == GRB.Status.INFEASIBLE  or \
   status == GRB.Status.UNBOUNDED:
   print('The m cannot be solved because it is infeasible or unbounded => status "%d"' % status)