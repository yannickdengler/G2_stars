import numpy as np
import global_vars as v
import sys
sys.path.insert(0, '/home/dengler_yannick/Documents/G2_stars/EoS/scripts')
import EoS

#                   HIERACHY OF VARIABLES!!
#
#           r, mu, P, eps, M, y         (OM,DM)
#
#           vals: mu_OM, mu_DM, M_OM, M_DM, y
#
#


######################################################### mu TOV

def d_nu_r(r, P, M):
    return (M+4*v.pi*r**3*P)/(r*(r-2*M))

########################### 1 fluid: vals: (mu, M)

def d_mu_r(r, vals):
    P = EoS.P_mu(vals[0])
    M = vals[1]
    return -vals[0]*d_nu_r(r, P, M)

def d_M_r_mu(r, vals):
    return 4*v.pi*r**2*EoS.eps_mu(vals[0])

########################### 1 fluid: vals: (P, M)

def d_P_r(r, vals):
    P = vals[0]
    M = vals[1]
    eps = EoS.eps_P(P)
    return -(P+eps)*d_nu_r(r, P, M)

def d_M_r_P(r, vals):
    return 4*v.pi*r**2*EoS.eps_P(vals[0])

########################### 2 fluid: vals: (mu_OM, mu_DM, M_OM, M_DM)

def d_mu_r_OM(r, vals):
    P = EoS.P_OM_mu(vals[0])+EoS.P_DM_mu(vals[1])
    M = vals[2]+vals[3]
    return -vals[0]*d_nu_r(r, P, M)

def d_M_r_mu_OM(r, vals):
    return 4*v.pi*r**2*EoS.eps_OM_mu(vals[0])

def d_mu_r_DM(r, vals):
    P = EoS.P_OM_mu(vals[0])+EoS.P_DM_mu(vals[1])
    M = vals[2]+vals[3]
    return -vals[1]*d_nu_r(r, P, M)

def d_M_r_mu_DM(r, vals):
    return 4*v.pi*r**2*EoS.eps_DM_mu(vals[1])

########################### 2 fluid: vals: (P_OM, P_DM, M_OM, M_DM)

def d_P_r_OM(r, vals):
    P_OM = vals[0]
    P_DM = vals[1]
    P = P_OM + P_DM
    eps_OM = EoS.eps_OM_P(P_OM)
    M = vals[2]+vals[3]
    return -(P_OM+eps_OM)*d_nu_r(r, P, M)

def d_M_r_P_OM(r, vals):
    return 4*v.pi*r**2*EoS.eps_OM_P(vals[0])

def d_P_r_DM(r, vals):
    P_OM = vals[0]
    P_DM = vals[1]
    P = P_OM + P_DM
    eps_DM = EoS.eps_DM_P(P_DM)
    M = vals[2]+vals[3]
    return -(P_DM+eps_DM)*d_nu_r(r, P, M)

def d_M_r_P_DM(r, vals):
    return 4*v.pi*r**2*EoS.eps_DM_P(vals[1])

########################## TIDAL 1 fluid vals[2] = y

def Q_mu(r, mu, M):
    P = EoS.P_mu(mu)
    eps = EoS.eps_mu(mu)
    c_s = EoS.c_s_mu(mu)
    term = 0
    if c_s != 0:
        term = (eps+P)/c_s
    return (4*v.pi*r)*(5*eps+9*P+term-6/(4*v.pi*r**2))/(r-2*M)-4*((M+4*v.pi*r**3*P)/(r**2*(1-2*M/r)))**2

def F_mu(r, mu, M):
    P = EoS.P_mu(mu)
    eps = EoS.eps_mu(mu)
    return (r-4*v.pi*r**3*(eps-P))/(r-2*M)

def d_y_r_mu(r, vals):
    mu = vals[0]
    M = vals[1]
    y = vals[2]
    return -(y**2 + y*F_mu(r, mu, M)+r**2*Q_mu(r, mu, M))/r

def Q_P(r, P, M):
    eps = EoS.eps_P(P)
    c_s = EoS.c_s_P(P)
    term = 0
    if c_s != 0:
        term = (eps+P)/c_s
    return (4*v.pi*r)*(5*eps+9*P+term-6/(4*v.pi*r**2))/(r-2*M)-4*((M+4*v.pi*r**3*P)/(r**2*(1-2*M/r)))**2

def F_P(r, P, M):
    eps = EoS.eps_P(P)
    return (r-4*v.pi*r**3*(eps-P))/(r-2*M)

def d_y_r_P(r, vals):
    P = vals[0]
    M = vals[1]
    y = vals[2]
    return -(y**2 + y*F_P(r, P, M)+r**2*Q_P(r, P, M))/r

########################## TIDAL 2 fluid vals[4] = y

def Q_mu_2_fluid(r, mu_OM, mu_DM, M):
    P_OM = EoS.P_OM_mu(mu_OM)
    P_DM = EoS.P_DM_mu(mu_DM)
    eps_OM = EoS.eps_OM_mu(mu_OM)
    eps_DM = EoS.eps_DM_mu(mu_DM)
    c_s_OM = EoS.c_s_OM_mu(mu_OM)
    c_s_DM = EoS.c_s_DM_mu(mu_DM)
    if c_s_OM == 0:
        OM_term = 0
    else:
        OM_term = (eps_OM+P_OM)/c_s_OM
    if c_s_DM == 0:
        DM_term = 0
    else:
        DM_term = (eps_DM+P_DM)/c_s_DM
    P = P_OM + P_DM
    eps = eps_OM + eps_DM
    return (4*v.pi*r)*(5*eps+9*P+OM_term+DM_term-6/(4*v.pi*r**2))/(r-2*M)-4*((M+4*v.pi*r**3*P)/(r**2*(1-2*M/r)))**2
   
def F_mu_2_fluid(r, mu_OM, mu_DM, M):
    P = EoS.P_OM_mu(mu_OM)+EoS.P_DM_mu(mu_DM)
    eps = EoS.eps_OM_mu(mu_OM)+EoS.eps_DM_mu(mu_DM)
    return (r-4*v.pi*r**3*(eps-P))/(r-2*M)

def d_y_r_mu_2(r, vals):
    mu_OM = vals[0]
    mu_DM = vals[1]
    M = vals[2]+vals[3]
    y = vals[4]
    return -(y**2 + y*F_mu_2_fluid(r, mu_OM, mu_DM, M)+r**2*Q_mu_2_fluid(r, mu_OM, mu_DM, M))/r


def Q_P_2_fluid(r, P_OM, P_DM, M):
    eps_OM = EoS.eps_OM_P(P_OM)
    eps_DM = EoS.eps_DM_P(P_DM)
    c_s_OM = EoS.c_s_OM_P(P_OM)
    c_s_DM = EoS.c_s_DM_P(P_DM)
    if c_s_OM == 0:
        OM_term = 0
    else:
        OM_term = (eps_OM+P_OM)/c_s_OM
    if c_s_DM == 0:
        DM_term = 0
    else:
        DM_term = (eps_DM+P_DM)/c_s_DM
    P = P_OM + P_DM
    eps = eps_OM + eps_DM
    return (4*v.pi*r)*(5*eps+9*P+OM_term+DM_term-6/(4*v.pi*r**2))/(r-2*M)-4*((M+4*v.pi*r**3*P)/(r**2*(1-2*M/r)))**2
   
def F_P_2_fluid(r, P_OM, P_DM, M):
    P = P_OM + P_DM
    eps = EoS.eps_OM_P(P_OM)+EoS.eps_DM_P(P_DM)
    return (r-4*v.pi*r**3*(eps-P))/(r-2*M)

def d_y_r_P_2(r, vals):
    P_OM = vals[0]
    P_DM = vals[1]
    M = vals[2]+vals[3]
    y = vals[4]
    return -(y**2 + y*F_P_2_fluid(r, P_OM, P_DM, M)+r**2*Q_P_2_fluid(r, P_OM, P_DM, M))/r


def k_2_y(y, C):
    return 8*C**5/5*(1-2*C)**2*(2+2*C*(y-1)-y)/(2*C*(6-3*y+3*C*(5*y-8))+4*C**3*(13-11*y+C*(3*y-2)+2*C**2*(1+y))+3*(1-2*C)**2*(2-y+2*C*(y-1))*np.log(1-2*C))

# def get_maximum_value_of_array_inter(array):








