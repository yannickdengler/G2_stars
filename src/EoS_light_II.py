import numpy as np
import matplotlib.pyplot as plt

almost_zero = 1e-50

m_neutron = 939.5654205    

c_DM, K_DM, g_DM, P_arr_DM = np.transpose(np.genfromtxt("EoS_calc/output/cKgP_DM_light.csv", delimiter=","))

def get_index_DM(P):
    for i in range(len(P_arr_DM)):
        if P <= P_arr_DM[i]:
            return i
    return 999
    
def n_of_P_DM(P):
    if P > 0:
        ind = get_index_DM(P)
        return (P/K_DM[ind])**(1/g_DM[ind])
    else:
        return 0

def eps_DM(P):
    if P > 0:
        ind = get_index_DM(P)
        n = (P/K_DM[ind])**(1/g_DM[ind])
        return c_DM[ind]*n + K_DM[ind]*n**g_DM[ind]/(g_DM[ind]-1)
    else:
        return 0

def c_s2_DM(P):
    if P > 0:
        ind = get_index_DM(P)
        n = (P/K_DM[ind])**(1/g_DM[ind])
        return g_DM[ind]*P/(P+c_DM[ind]*n + K_DM[ind]*n**g_DM[ind]/(g_DM[ind]-1))
    else:
        return 0
    
c_OMII, K_OMII, g_OMII, P_arr_OMII = np.transpose(np.genfromtxt("EoS_calc/output/cKgP_EoSII.csv", delimiter=","))

def rescale_OMII(m_DM):
    for i in range(len(P_arr_OMII)):
        c_OMII[i] = c_OMII[i]*m_neutron/m_DM
        P_arr_OMII[i] = P_arr_OMII[i]*m_neutron**4/m_DM**4
        K_OMII[i] = K_OMII[i]*(m_neutron**(4-3*g_OMII[i]))/(m_DM**(4-3*g_OMII[i]))

def get_index_OMII(P):
    for i in range(len(P_arr_OMII)):
        if P <= P_arr_OMII[i]:
            return i
    return 999
    
def n_of_P_OMII(P):
    if P > 0:
        ind = get_index_OMII(P)
        return (P/K_OMII[ind])**(1/g_OMII[ind])
    else:
        return 0

def eps_OMII(P):
    if P > 0:
        ind = get_index_OMII(P)
        n = (P/K_OMII[ind])**(1/g_OMII[ind])
        return c_OMII[ind]*n + K_OMII[ind]*n**g_OMII[ind]/(g_OMII[ind]-1)
    else:
        return 0

def c_s2_OMII(P):
    if P > 0:
        ind = get_index_OMII(P)
        n = (P/K_OMII[ind])**(1/g_OMII[ind])
        return g_OMII[ind]*P/(P+c_OMII[ind]*n + K_OMII[ind]*n**g_OMII[ind]/(g_OMII[ind]-1))
    else:
        return 0    

n_OM = n_of_P_OMII
eps_OM = eps_OMII
c_s2_OM = c_s2_OMII
rescale_OM = rescale_OMII


if __name__ == "__main__":
    Parr = np.logspace(np.log10(0.1*min(P_arr_DM)),np.log10(0.999*max(P_arr_DM)))
    plt.plot(Parr, [n_of_P_DM(P) for P in Parr], label = "n DM", color = "red")
    plt.plot(Parr, [eps_DM(P) for P in Parr], label = "eps DM", ls = "dashed", color = "red")
    plt.plot(Parr, [c_s2_DM(P) for P in Parr], label = "cs2 DM", ls = "dotted", color = "red")

    Parr = np.logspace(np.log10(0.1*min(P_arr_OMII)),np.log10(0.999*max(P_arr_OMII)))
    plt.plot(Parr, [n_of_P_OMII(P) for P in Parr], label = "n OMII", color = "orange")
    plt.plot(Parr, [eps_OMII(P) for P in Parr], label = "eps OMII", ls = "dashed", color = "orange")
    plt.plot(Parr, [c_s2_OMII(P) for P in Parr], label = "cs2 OMII", ls = "dotted", color = "orange")

    plt.xscale("log")
    plt.yscale("log")
    plt.grid()
    plt.legend()
    plt.show()
