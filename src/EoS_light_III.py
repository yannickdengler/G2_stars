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


c_OMIII, K_OMIII, g_OMIII, P_arr_OMIII = np.transpose(np.genfromtxt("EoS_calc/output/cKgP_EoSIII.csv", delimiter=","))

def rescale_OMIII(m_DM):
    for i in range(len(P_arr_OMIII)):
        c_OMIII[i] = c_OMIII[i]*m_neutron/m_DM
        P_arr_OMIII[i] = P_arr_OMIII[i]*m_neutron**4/m_DM**4
        K_OMIII[i] = K_OMIII[i]*(m_neutron**(4-3*g_OMIII[i]))/(m_DM**(4-3*g_OMIII[i]))

def get_index_OMIII(P):
    for i in range(len(P_arr_OMIII)):
        if P <= P_arr_OMIII[i]:
            return i
    return 999
    
def n_of_P_OMIII(P):
    if P > 0:
        ind = get_index_OMIII(P)
        return (P/K_OMIII[ind])**(1/g_OMIII[ind])
    else:
        return 0

def eps_OMIII(P):
    if P > 0:
        ind = get_index_OMIII(P)
        n = (P/K_OMIII[ind])**(1/g_OMIII[ind])
        return c_OMIII[ind]*n + K_OMIII[ind]*n**g_OMIII[ind]/(g_OMIII[ind]-1)
    else:
        return 0

def c_s2_OMIII(P):
    if P > 0:
        ind = get_index_OMIII(P)
        n = (P/K_OMIII[ind])**(1/g_OMIII[ind])
        return g_OMIII[ind]*P/(P+c_OMIII[ind]*n + K_OMIII[ind]*n**g_OMIII[ind]/(g_OMIII[ind]-1))
    else:
        return 0
    

n_OM = n_of_P_OMIII
eps_OM = eps_OMIII
c_s2_OM = c_s2_OMIII
rescale_OM = rescale_OMIII


if __name__ == "__main__":
    Parr = np.logspace(np.log10(0.1*min(P_arr_DM)),np.log10(0.999*max(P_arr_DM)))
    plt.plot(Parr, [n_of_P_DM(P) for P in Parr], label = "n DM", color = "red")
    plt.plot(Parr, [eps_DM(P) for P in Parr], label = "eps DM", ls = "dashed", color = "red")
    plt.plot(Parr, [c_s2_DM(P) for P in Parr], label = "cs2 DM", ls = "dotted", color = "red")

    Parr = np.logspace(np.log10(0.1*min(P_arr_OMIII)),np.log10(0.999*max(P_arr_OMIII)))
    plt.plot(Parr, [n_of_P_OMIII(P) for P in Parr], label = "n OMIII", color = "blue")
    plt.plot(Parr, [eps_OMIII(P) for P in Parr], label = "eps OMIII", ls = "dashed", color = "blue")
    plt.plot(Parr, [c_s2_OMIII(P) for P in Parr], label = "cs2 OMIII", ls = "dotted", color = "blue")

    plt.xscale("log")
    plt.yscale("log")
    plt.grid()
    plt.legend()
    plt.show()
