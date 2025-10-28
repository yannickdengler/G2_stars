import numpy as np
# import EoS


from EoS import eps_OM as eps_OM
from EoS import n_OM as n_OM
from EoS import c_s2_OM as c_s2_OM
from EoS import rescale_OM as rescale_OM

from EoS import eps_DM as eps_DM
from EoS import n_of_P_DM as n_DM
from EoS import c_s2_DM as c_s2_DM

def d_nu_r(r, P, M):
    return (M+4*np.pi*r**3*P)/(r*(r-2*M))

def d_P_r_OM(r, vals):
    P_OM = vals[0]
    P_DM = vals[1]
    P = P_OM + P_DM
    e_OM = eps_OM(P_OM)
    M = vals[2]+vals[3]
    return -(P_OM+e_OM)*d_nu_r(r, P, M)

def d_M_r_P_OM(r, vals):
    return 4*np.pi*r**2*eps_OM(vals[0])

def d_P_r_DM(r, vals):
    P_OM = vals[0]
    P_DM = vals[1]
    P = P_OM + P_DM
    e_DM = eps_DM(P_DM)
    M = vals[2]+vals[3]
    return -(P_DM+e_DM)*d_nu_r(r, P, M)

def d_M_r_P_DM(r, vals):
    return 4*np.pi*r**2*eps_DM(vals[1])

def Q_P_2_fluid(r, P_OM, P_DM, M):
    e_OM = eps_OM(P_OM)
    e_DM = eps_DM(P_DM)
    c_s_OM = c_s2_OM(P_OM)
    c_s_DM = c_s2_DM(P_DM)
    if c_s_OM == 0:
        OM_term = 0
    else:
        OM_term = (e_OM+P_OM)/c_s_OM
    if c_s_DM == 0:
        DM_term = 0
    else:
        DM_term = (e_DM+P_DM)/c_s_DM
    P = P_OM + P_DM
    eps = e_OM + e_DM
    return (4*np.pi*r)*(5*eps+9*P+OM_term+DM_term-6/(4*np.pi*r**2))/(r-2*M)-4*((M+4*np.pi*r**3*P)/(r**2*(1-2*M/r)))**2
   
def F_P_2_fluid(r, P_OM, P_DM, M):
    P = P_OM + P_DM
    eps = eps_OM(P_OM)+eps_DM(P_DM)
    return (r-4*np.pi*r**3*(eps-P))/(r-2*M)

def d_y_r_P_2(r, vals):
    P_OM = vals[0]
    P_DM = vals[1]
    M = vals[2]+vals[3]
    # y = vals[6]
    y = min([max([vals[6],-3]),2])
    # print(vals,r)
    return -(y**2 + y*F_P_2_fluid(r, P_OM, P_DM, M)+r**2*Q_P_2_fluid(r, P_OM, P_DM, M))/r

def k_2_y(y, C):
    return 8*C**5/5*(1-2*C)**2*(2+2*C*(y-1)-y)/(2*C*(6-3*y+3*C*(5*y-8))+4*C**3*(13-11*y+C*(3*y-2)+2*C**2*(1+y))+3*(1-2*C)**2*(2-y+2*C*(y-1))*np.log(1-2*C))

def d_N_r_OM(r, vals):
    P_OM = vals[0]
    M = vals[2]+vals[3]
    return 4*np.pi*r**2*n_OM(P_OM)/np.sqrt(1-2*M/r)

def d_N_r_DM(r, vals):
    P_DM = vals[1]
    M = vals[2]+vals[3]
    return 4*np.pi*r**2*n_DM(P_DM)/np.sqrt(1-2*M/r)