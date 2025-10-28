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
    
c_OMI, K_OMI, g_OMI, P_arr_OMI = np.transpose(np.genfromtxt("EoS_calc/output/cKgP_EoSI.csv", delimiter=","))

def rescale_OMI(m_DM):
    for i in range(len(P_arr_OMI)):
        c_OMI[i] = c_OMI[i]*m_neutron/m_DM
        P_arr_OMI[i] = P_arr_OMI[i]*m_neutron**4/m_DM**4
        K_OMI[i] = K_OMI[i]*(m_neutron**(4-3*g_OMI[i]))/(m_DM**(4-3*g_OMI[i]))

def get_index_OMI(P):
    for i in range(len(P_arr_OMI)):
        if P <= P_arr_OMI[i]:
            return i
    return 999
    
def n_of_P_OMI(P):
    if P > 0:
        ind = get_index_OMI(P)
        return (P/K_OMI[ind])**(1/g_OMI[ind])
    else:
        return 0

def eps_OMI(P):
    if P > 0:
        ind = get_index_OMI(P)
        n = (P/K_OMI[ind])**(1/g_OMI[ind])
        return c_OMI[ind]*n + K_OMI[ind]*n**g_OMI[ind]/(g_OMI[ind]-1)
    else:
        return 0

def c_s2_OMI(P):
    if P > 0:
        ind = get_index_OMI(P)
        n = (P/K_OMI[ind])**(1/g_OMI[ind])
        return g_OMI[ind]*P/(P+c_OMI[ind]*n + K_OMI[ind]*n**g_OMI[ind]/(g_OMI[ind]-1))
    else:
        return 0
    
n_OM = n_of_P_OMI
eps_OM = eps_OMI
c_s2_OM = c_s2_OMI
rescale_OM = rescale_OMI

if __name__ == "__main__":
    Parr = np.logspace(np.log10(0.1*min(P_arr_DM)),np.log10(0.999*max(P_arr_DM)))
    plt.plot(Parr, [n_of_P_DM(P) for P in Parr], label = "n DM", color = "red")
    plt.plot(Parr, [eps_DM(P) for P in Parr], label = "eps DM", ls = "dashed", color = "red")
    plt.plot(Parr, [c_s2_DM(P) for P in Parr], label = "cs2 DM", ls = "dotted", color = "red")

    Parr = np.logspace(np.log10(0.1*min(P_arr_OMI)),np.log10(0.999*max(P_arr_OMI)))
    plt.plot(Parr, [n_of_P_OMI(P) for P in Parr], label = "n OMI", color = "green")
    plt.plot(Parr, [eps_OMI(P) for P in Parr], label = "eps OMI", ls = "dashed", color = "green")
    plt.plot(Parr, [c_s2_OMI(P) for P in Parr], label = "cs2 OMI", ls = "dotted", color = "green")

    plt.xscale("log")
    plt.yscale("log")
    plt.grid()
    plt.legend()
    plt.show()
