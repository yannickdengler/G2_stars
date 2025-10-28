import numpy as np
import matplotlib.pyplot as plt
import matplotlib

import styles

m_neutron = 939.5654205   

############ LIGHT #############

c_DM_l , K_DM_l , g_DM_l , P_arr_DM_l  = np.transpose(np.genfromtxt("EoS/cKgP_DM_light.csv", delimiter=","))

def get_index_DM_l(P):
    for i in range(len(P_arr_DM_l)):
        if P <= P_arr_DM_l[i]:
            return i
    return 999
    
def n_of_P_DM_l(P):
    if P > 0:
        ind = get_index_DM_l(P)
        return (P/K_DM_l[ind])**(1/g_DM_l[ind])
    else:
        return 0

def eps_DM_l(P):
    if P > 0:
        ind = get_index_DM_l(P)
        n = (P/K_DM_l[ind])**(1/g_DM_l[ind])
        return c_DM_l[ind]*n + K_DM_l[ind]*n**g_DM_l[ind]/(g_DM_l[ind]-1)
    else:
        return 0

def c_s2_DM_l(P):
    if P > 0:
        ind = get_index_DM_l(P)
        n = (P/K_DM_l[ind])**(1/g_DM_l[ind])
        return g_DM_l[ind]*P/(P+c_DM_l[ind]*n + K_DM_l[ind]*n**g_DM_l[ind]/(g_DM_l[ind]-1))
    else:
        return 0

############ HEAVY #############

c_DM_h , K_DM_h , g_DM_h , P_arr_DM_h  = np.transpose(np.genfromtxt("EoS/cKgP_DM_heavy.csv", delimiter=","))

def get_index_DM_h(P):
    for i in range(len(P_arr_DM_h)):
        if P <= P_arr_DM_h[i]:
            return i
    return 999
    
def n_of_P_DM_h(P):
    if P > 0:
        ind = get_index_DM_h(P)
        return (P/K_DM_h[ind])**(1/g_DM_h[ind])
    else:
        return 0

def eps_DM_h(P):
    if P > 0:
        ind = get_index_DM_h(P)
        n = (P/K_DM_h[ind])**(1/g_DM_h[ind])
        return c_DM_h[ind]*n + K_DM_h[ind]*n**g_DM_h[ind]/(g_DM_h[ind]-1)
    else:
        return 0

def c_s2_DM_h(P):
    if P > 0:
        ind = get_index_DM_h(P)
        n = (P/K_DM_h[ind])**(1/g_DM_h[ind])
        return g_DM_h[ind]*P/(P+c_DM_h[ind]*n + K_DM_h[ind]*n**g_DM_h[ind]/(g_DM_h[ind]-1))
    else:
        return 0

############ EoS I #############

c_OMI  , K_OMI  , g_OMI  , P_arr_OMI   = np.transpose(np.genfromtxt("EoS/cKgP_EoSI.csv", delimiter=","))

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

############ EoS II #############

c_OMII , K_OMII , g_OMII , P_arr_OMII  = np.transpose(np.genfromtxt("EoS/cKgP_EoSII.csv", delimiter=","))

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

############ EoS III #############

c_OMIII, K_OMIII, g_OMIII, P_arr_OMIII = np.transpose(np.genfromtxt("EoS/cKgP_EoSIII.csv", delimiter=","))

# def rescale_OMIII(m_DM):
#     for i in range(len(P_arr_OMIII)):
#         c_OMIII[i] = c_OMIII[i]*m_neutron/m_DM
#         P_arr_OMIII[i] = P_arr_OMIII[i]*m_neutron**4/m_DM**4
#         K_OMIII[i] = K_OMIII[i]*(m_neutron**(4-3*g_OMIII[i]))/(m_DM**(4-3*g_OMIII[i]))

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


if __name__ == "__main__":
    p_low_EoS_I    = 5.030643560643215e-6
    p_low_EoS_II   = 6.367580732002786e-6
    p_low_EoS_III  = 5.030643560643215e-6
    p_low_DM       = 0.00004

    p_high_EoS_I   = 0.01859823271207
    p_high_EoS_II  = 0.0171003357007
    p_high_EoS_III = 0.0218441625408
    p_high_DM      = 0.03

    p_arr_I   = np.logspace(np.log10(0.1*p_low_EoS_I),np.log10(p_high_EoS_I),num = 1000)
    p_arr_II  = np.logspace(np.log10(0.1*p_low_EoS_II),np.log10(p_high_EoS_II),num = 1000)
    p_arr_III = np.logspace(np.log10(0.1*p_low_EoS_III),np.log10(p_high_EoS_III),num = 1000)
    p_arr_DM  = np.logspace(np.log10(0.1*p_low_DM),np.log10(p_high_DM),num = 1000)

    plt.plot([101425.86 *(eps_OMI(P)   if eps_OMI(P)   <= 0.01859823271207 else np.nan ) for P in p_arr_I]  , p_arr_I*101425.86  , label = r"$\text{EoS I}$"    , color = "red")
    plt.plot([101425.86 *(eps_OMII(P)  if eps_OMII(P)  <= 0.0171003357007  else np.nan ) for P in p_arr_II] , p_arr_II*101425.86 , label = r"$\text{EoS II}$"   , color = "green")
    plt.plot([101425.86 *(eps_OMIII(P) if eps_OMIII(P) <= 0.0218441625408  else np.nan ) for P in p_arr_III], p_arr_III*101425.86, label = r"$\text{EoS III}$"  , color = "blue")
    plt.plot([101425.86 *(eps_DM_l(P)  if eps_DM_l(P)  <= 0.03             else np.nan ) for P in p_arr_DM] , p_arr_DM*101425.86 , label = r"$\text{light EoS}$", color = "purple", ls = "dashed")
    plt.plot([101425.86 *(eps_DM_h(P)  if eps_DM_h(P)  <= 0.03             else np.nan ) for P in p_arr_DM] , p_arr_DM*101425.86 , label = r"$\text{heavy EoS}$", color = "orange", ls = "dashed")

    plt.plot([1e-5,1e5],[1e-5,1e5], color = "black")

    light_points = np.genfromtxt("EoS_calc/light_points.dat")
    plt.scatter(light_points[0],light_points[1], color = "purple")
    heavy_points = np.genfromtxt("EoS_calc/heavy_points.dat")
    plt.scatter(heavy_points[0],heavy_points[1], color = "orange")
    I_points = np.genfromtxt("EoS_calc/EoS_I_points.dat")
    plt.scatter(I_points[0,(len(I_points)-5):],I_points[1,(len(I_points)-5):], color = "red")
    plt.scatter(I_points[0,:(len(I_points)-5)],I_points[1,:(len(I_points)-5)], color = "red", marker="*")
    II_points = np.genfromtxt("EoS_calc/EoS_II_points.dat")
    plt.scatter(II_points[0,(len(II_points)-5):],II_points[1,(len(II_points)-5):], color = "green")
    plt.scatter(II_points[0,:(len(II_points)-5)],II_points[1,:(len(II_points)-5)], color = "green", marker="*")
    III_points = np.genfromtxt("EoS_calc/EoS_III_points.dat")
    plt.scatter(III_points[0,(len(III_points)-5):],III_points[1,(len(III_points)-5):], color = "blue")
    plt.scatter(III_points[0,:(len(III_points)-5)],III_points[1,:(len(III_points)-5)], color = "blue", marker="*")

    plt.grid()

    plt.xscale("log")
    plt.yscale("log")

    xticks = [300,1000,3000]
    plt.xticks(xticks,[r"$%i$"%x for x in xticks])    
    yticks = [10,100]
    plt.yticks(yticks,[r"$%i$"%x for x in yticks])    

    plt.xlabel(r"$\varepsilon\,\text{ [MeV/fm}^3\text{]}$")
    plt.ylabel(r"$p\,\text{ [MeV/fm}^3\text{]}$")

    plt.xlim([1e-3*130149,4000])
    plt.ylim([1e-5*130149,950])

    plt.legend()
    plt.savefig("EoS.pdf", bbox_inches = "tight")
    # plt.show()
    plt.clf()

if __name__ == "__main__":
    p_low_EoS_I    = 5.030643560643215e-6
    p_low_EoS_II   = 6.367580732002786e-6
    p_low_EoS_III  = 5.030643560643215e-6
    p_low_DM       = 0.00004

    p_high_EoS_I   = 0.01859823271207
    p_high_EoS_II  = 0.0171003357007
    p_high_EoS_III = 0.0218441625408
    p_high_DM      = 0.03

    p_arr_I   = np.logspace(np.log10(0.1*p_low_EoS_I),np.log10(p_high_EoS_I),num = 1000)
    p_arr_II  = np.logspace(np.log10(0.1*p_low_EoS_II),np.log10(p_high_EoS_II),num = 1000)
    p_arr_III = np.logspace(np.log10(0.1*p_low_EoS_III),np.log10(p_high_EoS_III),num = 1000)
    p_arr_DM  = np.logspace(np.log10(0.1*p_low_DM),np.log10(p_high_DM),num = 1000)

    plt.plot([101425.86 *(eps_OMI(P)   if eps_OMI(P)   <= 0.01859823271207 else np.nan ) for P in p_arr_I]  , 
             [c_s2_OMI(P)   if eps_OMI(P)   <= 0.01859823271207 else np.nan for P in p_arr_I], label = r"$\text{EoS I}$"    , color = "red")
    plt.plot([101425.86 *(eps_OMII(P)  if eps_OMII(P)  <= 0.0171003357007  else np.nan ) for P in p_arr_II] , 
             [c_s2_OMII(P)  if eps_OMII(P)  <= 0.0171003357007  else np.nan for P in p_arr_II], label = r"$\text{EoS II}$"  , color = "green")
    plt.plot([101425.86 *(eps_OMIII(P) if eps_OMIII(P) <= 0.0218441625408  else np.nan ) for P in p_arr_III], 
             [c_s2_OMIII(P) if eps_OMIII(P) <= 0.0218441625408  else np.nan for P in p_arr_III], label = r"$\text{EoS III}$"  , color = "blue")
    plt.plot([101425.86 *(eps_DM_l(P)  if eps_DM_l(P)  <= 0.03             else np.nan ) for P in p_arr_DM] , 
             [c_s2_DM_l(P)  if eps_DM_l(P)  <= 0.03             else np.nan for P in p_arr_DM] , label = r"$\text{light EoS}$", color = "purple", ls = "dashed")
    plt.plot([101425.86 *(eps_DM_h(P)  if eps_DM_h(P)  <= 0.03             else np.nan ) for P in p_arr_DM] , 
             [c_s2_DM_h(P)  if eps_DM_h(P)  <= 0.03             else np.nan for P in p_arr_DM] , label = r"$\text{heavy EoS}$", color = "orange", ls = "dashed")

    # I_points = np.genfromtxt("EoS_calc/EoS_I_points.dat")
    # plt.scatter(I_points[0,(len(I_points)-5):],I_points[1,(len(I_points)-5):], color = "red")
    # plt.scatter(I_points[0,:(len(I_points)-5)],I_points[1,:(len(I_points)-5)], color = "red", marker="*")
    # II_points = np.genfromtxt("EoS_calc/EoS_II_points.dat")
    # plt.scatter(II_points[0,(len(II_points)-5):],II_points[1,(len(II_points)-5):], color = "green")
    # plt.scatter(II_points[0,:(len(II_points)-5)],II_points[1,:(len(II_points)-5)], color = "green", marker="*")

    I_points = np.genfromtxt("EoS_calc/EoS_I_points.dat")
    plt.scatter(I_points[0,(len(I_points)-5):],[c_s2_OMI(P) for P in I_points[1,(len(I_points)-5):]/101425.86], color = "red")
    plt.scatter(I_points[0,:(len(I_points)-5)],[c_s2_OMI(P) for P in I_points[1,:(len(I_points)-5)]/101425.86], color = "red", marker="*")
    II_points = np.genfromtxt("EoS_calc/EoS_II_points.dat")
    plt.scatter(II_points[0,(len(II_points)-5):],[c_s2_OMII(P) for P in II_points[1,(len(II_points)-5):]/101425.86], color = "green")
    plt.scatter(II_points[0,:(len(II_points)-5)],[c_s2_OMII(P) for P in II_points[1,:(len(II_points)-5)]/101425.86], color = "green", marker="*")
    III_points = np.genfromtxt("EoS_calc/EoS_III_points.dat")
    plt.scatter(III_points[0,(len(III_points)-5):],[c_s2_OMIII(P) for P in III_points[1,(len(III_points)-5):]/101425.86], color = "blue")
    plt.scatter(III_points[0,:(len(III_points)-5)],[c_s2_OMIII(P) for P in III_points[1,:(len(III_points)-5)]/101425.86], color = "blue", marker="*")
    light_points = np.genfromtxt("EoS_calc/light_points.dat")
    plt.scatter(light_points[0],[c_s2_DM_l(P) for P in light_points[1]/101425.86], color = "purple")
    heavy_points = np.genfromtxt("EoS_calc/heavy_points.dat")
    plt.scatter(heavy_points[0],[c_s2_DM_h(P) for P in heavy_points[1]/101425.86], color = "orange")

    plt.grid()

    plt.axhline(1, color = "black")

    plt.xscale("log")
    plt.yscale("log")

    # plt.xlabel("$c_s^2")
    # plt.ylabel("$p$ [MeV/fm$^3$]")
    plt.xlabel(r"$\varepsilon\,\text{ [MeV/fm}^3\text{]}$")
    plt.ylabel(r"$c_s^2$")

    plt.xlim([1e-3*130149,4000])
    plt.ylim([1e-2,1.1])

    xticks = [300,1000,3000]
    plt.xticks(xticks,[r"$%i$"%x for x in xticks])    
    yticks = np.logspace(-2,0,3)
    plt.yticks(yticks,[r"$10^{%i}$"%np.log10(x) for x in yticks])    

    plt.legend()
    plt.savefig("c_s_2.pdf", bbox_inches = "tight")
    # plt.show()
