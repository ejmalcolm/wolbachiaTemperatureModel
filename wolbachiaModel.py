from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import PercentFormatter
from fractions import Fraction
import math
import numpy as np

def infection(x,t,v,da,dj,u,r,tau,D,b,q):
    M_i = x[0]
    M_u = x[1]
    F_i = x[2]
    F_u = x[3]
    J_i = x[4]
    J_u = x[5]
    
    M_u0 = M_u #used to calculate weekly release

    A = .0002
    beta = 3.4

    # v = juvenile -> adulthood rate
    # da = adult death rate
    # dj = juvenile death rate
    # u = weekly released mosquitos, expressed as percentage of uninfected equilibrium
    # r = fraction of wolbachia-infected mosquitos that are male
    # tau = vertical transmission rate
    # D = additional fitness cost
    # b = birth and survival to juvenile status rate
    # q = rate of cytoplasmic incompatibility

    temp = 25.1030 - 4.4346*math.cos( (2*math.pi/365) * (t + 338.4344) )
    Lv = -0.0449*(temp**2) + 2.1474*temp - 0.9722

    da = 1 / Lv

    M_i = (.5 * v * J_i) - (da*M_i) + r*(u*M_u0/7)
    M_u = (.5 * v * J_u) - (da*M_u)
    F_i = (.5 * v * J_i) - (da*F_i) + (1-r)*(u*M_u0/7)
    F_u = (.5 * v * J_u) - (da*F_u)
    J_i = (tau*b*F_i) - ( ( A * (J_i+J_u) ) ** (beta-1) )*J_i - (v*J_i)
    J_u = (1-tau)*b*F_i + b*( (1-q) * (M_i / (M_i+M_u) ) )*F_u - ( ( A * (J_i+J_u) ) ** (beta-1) )*J_u - (v*J_u)
    
    return [M_i, M_u, F_i, F_u, J_i, J_u]

def graph_for_ratios(ratios,params):
    t = np.linspace(0,52,1000)
    ###ADD BACK IN for ratio in ratios:
    y0 = [0,100,0,100,0,100] # mi, mu, fi, fu, ji, ju
    solution = odeint(infection,y0,t,params)
    plt.xlabel("Time")
    plt.ylabel("Number in Group")
    plt.plot(t, solution[:,0], label="Mi")
    plt.plot(t, solution[:,1], label="Mu")
    plt.plot(t, solution[:,2], label="Fi")
    plt.plot(t, solution[:,3], label="Fu")
    plt.plot(t, solution[:,4], label="Ji")
    plt.plot(t, solution[:,5], label="Ju")
    plt.legend()

ratio_to_test = (.5,1,2,3) #I:U ratios to simulate
params = ( 
    (1/15,1/15,1/15,1/2,1,1,.05,9,1),
        )  # v, da, dj, u, r, tau, D, b, q
n = 0 # used for looping through the coordinates of the grid
for set_params in params:
    # ax = plt.subplot(gs[0,n],title="v:%s,da:%s,dj:%s,u:%s,r:%s, tau:%s, D:%s, b:%s, q:%s" % 
    # tuple([Fraction(set_params[i]).limit_denominator() for i in range(len(set_params))]))
    # ^this is hideous but it's the best i got right now
    graph_for_ratios(ratio_to_test,set_params)
    n += 1

plt.show()