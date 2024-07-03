import numpy as np
from scipy.interpolate import interp1d as inter
import sys

MeV_fm = 197.32698  
mass_DM = float(sys.argv[1])   # in MeV
conv_P = mass_DM**4/MeV_fm**3

EoS_III = np.transpose(np.genfromtxt("EoS/data/EoS_III_inter_P_eps"))
G2_EoS = np.transpose(np.genfromtxt("EoS/data/EoS_G2_PP_inter_P_eps"))

P_OM = EoS_III[0]/conv_P
eps_OM = EoS_III[1]/conv_P
P_DM = G2_EoS[0]
eps_DM = G2_EoS[1]

OM_inter = inter(P_OM,eps_OM)
DM_inter = inter(P_DM,eps_DM)

def eps_OM_P(P):
    if P <= 0:
        return 0
    else:
        return OM_inter(P)
def eps_DM_P(P):
    if P <= 0:
        return 0
    else:
        return DM_inter(P)

def c_s_OM_P(P):
    if P <= 0:
        return 0
    else:
        small = 1e-3*P
        return 2*small/(eps_OM_P(P+small)-eps_OM_P(P-small))
def c_s_DM_P(P):
    if P <= 0:
        return 0
    else:
        small = 1e-3*P
        return 2*small/(eps_DM_P(P+small)-eps_DM_P(P-small))



def P_mu(mu):        # not needed
    return 0
def eps_mu(mu):        # not needed
    return 0
def c_s_mu(mu):        # not needed
    return 0

def mu_P(P):        # not needed
    return 0
def eps_P(P):        # not needed
    return 0
def c_s_P(P):        # not needed
    return 0

##########################################################

def P_OM_mu(mu):        # not needed
    return 0
def eps_OM_mu(mu):        # not needed
    return 0
def c_s_OM_mu(mu):        # not needed
    return 0

def mu_OM_P(P):        # not needed
    return 0


def P_DM_mu(mu):        # not needed
    return 0
def eps_DM_mu(mu):        # not needed
    return 0
def c_s_DM_mu(mu):        # not needed
    return 0

def mu_DM_P(P):        # not needed
    return 0