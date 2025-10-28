import numpy as np
import matplotlib.pyplot as plt

almost_zero = 1e-50


# c_DM, K_DM, g_DM, P_arr_DM = np.transpose(np.genfromtxt("cKgP_DM_light.csv", delimiter=","))
c_DM, K_DM, g_DM, P_arr_DM = np.transpose(np.genfromtxt("cKgP_DM_heavy.csv", delimiter=","))

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
    
m_neutron = 939.5654205                                                              # in MeV

# c_OMI, K_OMI, g_OMI, P_arr_OMI = np.transpose(np.genfromtxt("cKgP_EoSI.csv", delimiter=","))
# def rescale_OMI(m_DM):
#     for i in range(len(P_arr_OMI)):
#         c_OMI[i] = c_OMI[i]*m_neutron/m_DM
#         P_arr_OMI[i] = P_arr_OMI[i]*m_neutron**4/m_DM**4
#         K_OMI[i] = K_OMI[i]*(m_neutron**(4-3*g_OMI[i]))/(m_DM**(4-3*g_OMI[i]))

# def get_index_OMI(P):
#     for i in range(len(P_arr_OMI)):
#         if P <= P_arr_OMI[i]:
#             return i
#     return 999
    
# def n_of_P_OMI(P):
#     if P > 0:
#         ind = get_index_OMI(P)
#         return (P/K_OMI[ind])**(1/g_OMI[ind])
#     else:
#         return 0

# def eps_OMI(P):
#     if P > 0:
#         ind = get_index_OMI(P)
#         n = (P/K_OMI[ind])**(1/g_OMI[ind])
#         return c_OMI[ind]*n + K_OMI[ind]*n**g_OMI[ind]/(g_OMI[ind]-1)
#     else:
#         return 0

# def c_s2_OMI(P):
#     if P > 0:
#         ind = get_index_OMI(P)
#         n = (P/K_OMI[ind])**(1/g_OMI[ind])
#         return g_OMI[ind]*P/(P+c_OMI[ind]*n + K_OMI[ind]*n**g_OMI[ind]/(g_OMI[ind]-1))
#     else:
#         return 0

# c_OMII, K_OMII, g_OMII, P_arr_OMII = np.transpose(np.genfromtxt("cKgP_EoSII.csv", delimiter=","))

# def rescale_OMII(m_DM):
#     for i in range(len(P_arr_OMII)):
#         c_OMII[i] = c_OMII[i]*m_neutron/m_DM
#         P_arr_OMII[i] = P_arr_OMII[i]*m_neutron**4/m_DM**4
#         K_OMII[i] = K_OMII[i]*(m_neutron**(4-3*g_OMII[i]))/(m_DM**(4-3*g_OMII[i]))

# def get_index_OMII(P):
#     for i in range(len(P_arr_OMII)):
#         if P <= P_arr_OMII[i]:
#             return i
#     return 999
    
# def n_of_P_OMII(P):
#     if P > 0:
#         ind = get_index_OMII(P)
#         return (P/K_OMII[ind])**(1/g_OMII[ind])
#     else:
#         return 0

# def eps_OMII(P):
#     if P > 0:
#         ind = get_index_OMII(P)
#         n = (P/K_OMII[ind])**(1/g_OMII[ind])
#         return c_OMII[ind]*n + K_OMII[ind]*n**g_OMII[ind]/(g_OMII[ind]-1)
#     else:
#         return 0

# def c_s2_OMII(P):
#     if P > 0:
#         ind = get_index_OMII(P)
#         n = (P/K_OMII[ind])**(1/g_OMII[ind])
#         return g_OMII[ind]*P/(P+c_OMII[ind]*n + K_OMII[ind]*n**g_OMII[ind]/(g_OMII[ind]-1))
#     else:
#         return 0
    

c_OMIII, K_OMIII, g_OMIII, P_arr_OMIII = np.transpose(np.genfromtxt("cKgP_EoSIII.csv", delimiter=","))

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

    Parr = np.logspace(np.log10(0.1*min(P_arr_OMI)),np.log10(0.999*max(P_arr_OMI)))
    plt.plot(Parr, [n_of_P_OMI(P) for P in Parr], label = "n OMI", color = "green")
    plt.plot(Parr, [eps_OMI(P) for P in Parr], label = "eps OMI", ls = "dashed", color = "green")
    plt.plot(Parr, [c_s2_OMI(P) for P in Parr], label = "cs2 OMI", ls = "dotted", color = "green")

    Parr = np.logspace(np.log10(0.1*min(P_arr_OMII)),np.log10(0.999*max(P_arr_OMII)))
    plt.plot(Parr, [n_of_P_OMII(P) for P in Parr], label = "n OMII", color = "orange")
    plt.plot(Parr, [eps_OMII(P) for P in Parr], label = "eps OMII", ls = "dashed", color = "orange")
    plt.plot(Parr, [c_s2_OMII(P) for P in Parr], label = "cs2 OMII", ls = "dotted", color = "orange")

    Parr = np.logspace(np.log10(0.1*min(P_arr_OMIII)),np.log10(0.999*max(P_arr_OMIII)))
    plt.plot(Parr, [n_of_P_OMIII(P) for P in Parr], label = "n OMIII", color = "blue")
    plt.plot(Parr, [eps_OMIII(P) for P in Parr], label = "eps OMIII", ls = "dashed", color = "blue")
    plt.plot(Parr, [c_s2_OMIII(P) for P in Parr], label = "cs2 OMIII", ls = "dotted", color = "blue")

    plt.xscale("log")
    plt.yscale("log")
    plt.grid()
    plt.legend()
    plt.show()
