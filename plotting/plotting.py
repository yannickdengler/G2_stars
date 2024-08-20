import sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

from matplotlib import cm
from matplotlib.ticker import LinearLocator

MeV_fm = 197.32698  

EoS = int(sys.argv[1])
EoS_str = ""
if EoS == 1:
    EoS_str = "I"
    def EoS_max(m_DM):
        return 582.1/m_DM**4*MeV_fm**3
elif EoS == 2:
    def EoS_max(m_DM):
        return 238.7/m_DM**4*MeV_fm**3
    EoS_str = "II"
elif EoS == 3:
    def EoS_max(m_DM):
        return 136.8/m_DM**4*MeV_fm**3
    EoS_str = "III"
else:
    print("Eos_%s is not supported"%EoS)
    exit()
MeV_fm = 197.32698  
MeV_solar = 8.97*1e-61                                                      # 1 MeV in solar masses
Planck_Mass = 1.095*1e-38                                                   # in solar masses
Planck_Length = 1.6163e-38                                                  # in km

dm_om_arr=(0.250000, 0.325367, 0.423453, 0.551110, 0.717251, 0.933478, 1.214891, 1.581139, 2.057799, 2.678155, 3.485528, 4.536296, 5.903836, 7.683642, 10.000000)
M_arr=(250.000000, 304.753414, 371.498572, 452.861832, 552.044757, 672.950096, 820.335356, 1000.000000, 1219.013654, 1485.994289, 1811.447329, 2208.179027, 2691.800385, 3281.341424, 4000.000000)

def get_M_R(EoS_str, dm_om, M_DM, norm=1):
    res = np.genfromtxt("results/full_EoS_%s_run/M_R_TIDAL_P_dm_om_%1.3e_M_%1.3e_2_fluid.out"%(EoS_str,dm_om, M_DM))
    con_mass = pow(Planck_Mass, 3)*pow(M_DM*MeV_solar, -2)
    con_radius = pow(Planck_Mass, 2)*pow(M_DM*MeV_solar, -2)*Planck_Length
    for ii in range(len(res)):
        res[ii][1] = con_radius*res[ii][1]*norm
        res[ii][2] = con_radius*res[ii][2]*norm
        res[ii][3] = con_mass*res[ii][3]*norm
        res[ii][4] = con_mass*res[ii][4]*norm
    return res

def get_stable(EoS_str, dm_om, M_DM, norm=1):
    M_R = get_M_R(EoS_str, dm_om, M_DM, norm)
    tmp = np.genfromtxt("results/full_EoS_%s_run/stability_P_dm_om_%1.3e_M_%1.3e_2_fluid"%(EoS_str,dm_om, M_DM))[1:]
    stability = []
    for k in range(len(tmp)):
        if (tmp[k][0] < EoS_max(M_DM)):
            if (tmp[k][1] <= EoS_max(M_DM)):
                stability.append(list(tmp[k]))
            else:
                stability.append([tmp[k][0],EoS_max(M_DM)])
    M_R_stable = []
    for m in range(len(stability)):               # num of stable in branch
        M_R_stable.append([])
        for k in range(len(M_R[0])):              # num of obs  
            M_R_stable[m].append([])
            for l in range(len(M_R)):             # num of P0
                if M_R[l][0] > stability[m][0] and M_R[l][0] < stability[m][1]:
                    M_R_stable[m][k].append(M_R[l][k])
    return M_R_stable

def get_stable_M_R_dm_om_M_arr():
    '''
        returns an array for the stable branches of all M_DM and dm_om of shape [dm_om][m_DM][stability][ops][p0].
    '''
    res = []
    for i in range(len(dm_om_arr)):
        res.append([])
        for j in range(len(M_arr)):
            res[i].append(get_stable(EoS_str, dm_om_arr[i], M_arr[j]))
    return res

def get_stable_M_R_pure(norm=1):
    '''
        returns two array (pure OM and pure DM) of shape [stability][ops][p0].
    '''
    return get_stable(EoS_str, 0, 1000), get_stable(EoS_str, -1, 1000, norm)

def combine_stables(M_R_stable):
    res = []
    for k in range(len(M_R_stable)):                  # stable
        for m in range(len(M_R_stable[k])):           # obs
            res.append([])
            for l in range(len(M_R_stable[k][m])):             # p0
                res[m].append(M_R_stable[k][m][l])
    return res

def combine_stables_dm_om_M_DM(M_R_stable):
    res = []
    for i in range(len(M_R_stable)):
        res.append([])
        for j in range(len(M_R_stable[i])):
            res[i].append([])
            combine = combine_stables(M_R_stable[i][j])
            for k in range(len(combine)):
                res[i][j].append(combine[k])
    return res

def get_data(M_R_stable):
    """
    Extracts data from M_R_stable of the form {obs}[p0 (not for all)] and translates it to a dicc
    """
    res = {}
    res["M_OM"]=M_R_stable[3]
    res["M_DM"]=M_R_stable[4]
    res["M_sum"]=[M_R_stable[3][i]+M_R_stable[4][i] for i in range(len(M_R_stable[3]))]
    ind_crit=np.argmax(res["M_sum"])
    res["ind_crit"]=ind_crit
    res["p_crit"]=M_R_stable[0][ind_crit]
    res["R_OM"]=M_R_stable[1]
    res["R_DM"]=M_R_stable[2]
    res["R_max"]=[max(M_R_stable[1][i],M_R_stable[2][i]) for i in range(len(M_R_stable[3]))]
    res["M_max"]=res["M_sum"][ind_crit]
    res["M_max_OM"]=M_R_stable[3][ind_crit]
    res["M_max_DM"]=M_R_stable[4][ind_crit]
    res["R_crit"]=res["R_max"][ind_crit]
    res["R_crit_OM"]=M_R_stable[1][ind_crit]
    res["R_crit_DM"]=M_R_stable[2][ind_crit]
    res["k2_max"]=max(M_R_stable[7])
    res["k2_crit"]=M_R_stable[7][ind_crit]
    res["C_max"]=max(M_R_stable[6])
    res["C_crit"]=M_R_stable[6][ind_crit]
    return res

def get_data_dm_om_M_arr(M_R_stable):
    """
    Calls get_data from M_R_stable of form [dm_om][M_DM] to give a dicc of shape {ops}[dm_om][M_DM][p0 (not for all)]
    """
    res = {}
    data = []
    for i in range(len(M_R_stable)):
        data.append([])
        for j in range(len(M_R_stable[i])):
            data[i].append(get_data(M_R_stable[i][j]))
    for key in data[0][0].keys():
        res[key] = []
        for i in range(len(M_R_stable)):
            res[key].append([])
            for j in range(len(M_R_stable[i])):
                res[key][i].append([])
                res[key][i][j] = data[i][j][key]
    return res

def plot_3D(arr, label):
    '''
        Creates a 3D plot in x-m_DM, y-dm_om, z-"arr" (with label=label)
    '''
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    M_plot, dm_om_plot = np.meshgrid(M_arr, dm_om_arr)

    surf = ax.plot_surface(np.log10(M_plot), np.log10(dm_om_plot), np.log10(arr), cmap=cm.coolwarm, linewidth=0, antialiased=False)

    ax.set_zlabel(label)
    ax.set_xlabel("$M_{DM}$")
    ax.set_ylabel("$p_{0,DM}/p_{0,OM}")
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter('{x:.02f}')

    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()

def plot_M_R_dm_om(data, M_ind, show=True, save = True):
    cmap = mpl.colormaps['viridis']
    colors = cmap(np.linspace(0, 1, len(data["R_max"])))
    plt.grid()
    for i in range(len(data["R_max"])):
        plt.plot(data["R_max"][i][M_ind], data["M_sum"][i][M_ind], color = colors[i])
    if save:
        plt.savefig("plots/full_EoS_%s_run/M_R_EoS_%s_dm_om_M_DM_%e.pdf"%(EoS_str,EoS_str,M_arr[M_ind]))
    if show:
        plt.show()
    # plt.clf()

def plot_M_R_dm_om_keep_stable(data, M_ind, show=True, save = True):
    cmap = mpl.colormaps['viridis']
    colors = cmap(np.linspace(0, 1, 15))
    plt.grid()
    for i in range(len(data)):
        for j in range(len(data[i][M_ind])):
            plt.plot(data[i][M_ind][j]["R_max"], data[i][M_ind][j]["M_sum"], color = colors[i])
    if save:
        plt.savefig("plots/full_EoS_%s_run/M_R_EoS_%s_dm_om_M_DM_%e.pdf"%(EoS_str,EoS_str,M_arr[M_ind]))
    if show:
        plt.show()

# def plot_M_R_dm_om_branches(data, M_ind, show=True, save = True):
#     cmap = mpl.colormaps['viridis']
#     colors = cmap(np.linspace(0, 1, len(data["R_max"])))
#     plt.grid()
#     for i in range(len(data["R_max"])):
#         plt.plot(data["R_max"][i][M_ind], data["M_sum"][i][M_ind], color = colors[i])
#     if save:
#         plt.savefig("plots/full_EoS_%s_run/M_R_EoS_%s_dm_om_M_DM_%e.pdf"%(EoS_str,EoS_str,M_arr[M_ind]))
#     if show:
#         plt.show()
    # plt.clf()

# def plot_M_R_M_DM(data, dm_om_ind, show=True, save = True):
#     cmap = mpl.colormaps['viridis']
#     colors = cmap(np.linspace(0, 1, len(data[0][0][0]["R_max"])))
#     plt.grid()
#     for dat in data:
#         for i in range(len(dat[0][0]["R_max"])):
#             plt.plot(dat[dm_om_ind][i]["R_max"], dat[dm_om_ind][i]["M_sum"], color = colors[i])
#     if save:
#         plt.savefig("plots/full_EoS_%s_run/M_R_EoS_%s_M_DM_dm_om_%e.pdf"%(EoS_str,EoS_str,dm_om_arr[dm_om_ind]))
#     if show:
#         plt.show()
#     plt.clf()

# def plot_M_R_M_DM_branches(data, dm_om_ind, show=True, save = True):
#     cmap = mpl.colormaps['viridis']
#     colors = cmap(np.linspace(0, 1, len(data["R_max"])))
#     plt.grid()
#     for dat in data:
#         for i in range(len(data["R_max"][0])):
#             plt.plot(dat["R_max"][dm_om_ind][i], dat["M_sum"][dm_om_ind][i], color = colors[i])
#     if save:
#         plt.savefig("plots/full_EoS_%s_run/M_R_EoS_%s_M_DM_dm_om_%e.pdf"%(EoS_str,EoS_str,dm_om_arr[dm_om_ind]))
#     if show:
#         plt.show()
#     plt.clf()


if __name__ == "__main__":
    stable = get_stable_M_R_dm_om_M_arr()
    # stable_c = combine_stables_dm_om_M_DM(stable)
    # data = get_data_dm_om_M_arr(stable_c)

    # data = []
    # for k in range(len(stable[0][0])):
    #     data.append(get_data_dm_om_M_arr([[stable[j][i][k] for i in range(len(stable[0]))] for j in range(len(stable))]))

    data = []
    for i in range(len(stable)):             # dm_om
        data.append([])
        for j in range(len(stable[i])):      # m_dm
            data[i].append([])
            for k in range(len(stable[i][j])):  # stable
                data[i][j].append(get_data(stable[i][j][k]))

    for i in range(15):
        stable_pure_OM, stable_pure_DM = get_stable_M_R_pure(norm=(M_arr[i]/1000)**(-2))
        stable_c_pure_OM = combine_stables(stable_pure_OM)
        data_OM = get_data(stable_c_pure_OM)
        stable_c_pure_DM = combine_stables(stable_pure_DM)
        data_DM = get_data(stable_c_pure_DM)

        plt.plot(data_OM["R_max"], data_OM["M_sum"], color = "black", label = "OM", ls = "dashed")
        plt.plot(data_DM["R_max"], data_DM["M_sum"], color = "black", label = "DM", ls = "dashdot")
        plt.legend()
        # for dat in data:
        # plot_M_R_dm_om(data, i, show = False, save=False)
        plot_M_R_dm_om_keep_stable(data, i, show = False, save=False)
        plt.savefig("plots/full_EoS_%s_run/M_R_EoS_%s_dm_om_M_DM_%e.pdf"%(EoS_str,EoS_str,M_arr[i]))
        plt.clf()

        







    # plot_3D(data["M_max"],"M_max")

    ##########################################################




# def get_data(M_R_stable):
#     res = {}
#     res["M_OM"]=[[[] for i in range(len(M_R_stable))] for j in range(len(M_R_stable[0]))]
#     res["M_DM"]=[[[] for i in range(len(M_R_stable))] for j in range(len(M_R_stable[0]))]
#     res["M_sum"]=[[[] for i in range(len(M_R_stable))] for j in range(len(M_R_stable[0]))]
#     res["ind_crit"]=[[[] for i in range(len(M_R_stable))] for j in range(len(M_R_stable[0]))]
#     res["p_crit"]=[[[] for i in range(len(M_R_stable))] for j in range(len(M_R_stable[0]))]
#     res["R_OM"]=[[[] for i in range(len(M_R_stable))] for j in range(len(M_R_stable[0]))]
#     res["R_DM"]=[[[] for i in range(len(M_R_stable))] for j in range(len(M_R_stable[0]))]
#     res["R_max"]=[[[] for i in range(len(M_R_stable))] for j in range(len(M_R_stable[0]))]
#     res["M_max"]=[[[] for i in range(len(M_R_stable))] for j in range(len(M_R_stable[0]))]
#     res["M_max_OM"]=[[[] for i in range(len(M_R_stable))] for j in range(len(M_R_stable[0]))]
#     res["M_max_DM"]=[[[] for i in range(len(M_R_stable))] for j in range(len(M_R_stable[0]))]
#     res["R_crit"]=[[[] for i in range(len(M_R_stable))] for j in range(len(M_R_stable[0]))]
#     res["R_crit_OM"]=[[[] for i in range(len(M_R_stable))] for j in range(len(M_R_stable[0]))]
#     res["R_crit_DM"]=[[[] for i in range(len(M_R_stable))] for j in range(len(M_R_stable[0]))]
#     res["k2_max"]=[[[] for i in range(len(M_R_stable))] for j in range(len(M_R_stable[0]))]
#     res["k2_crit"]=[[[] for i in range(len(M_R_stable))] for j in range(len(M_R_stable[0]))]
#     res["C_max"]=[[[] for i in range(len(M_R_stable))] for j in range(len(M_R_stable[0]))]
#     res["C_crit"]=[[[] for i in range(len(M_R_stable))] for j in range(len(M_R_stable[0]))]

#     for i in range(len(M_R_stable)):
#         for j in range(len(M_R_stable[i])):
#             print(i,j)
#             res["M_OM"][i][j] = M_R_stable[i][j][3]
#             res["M_DM"][i][j] = M_R_stable[i][j][4]
#             for k in range(len(M_R_stable[i][j][3])):
#                 res["M_sum"][i][j].append(M_R_stable[i][j][3][k]+M_R_stable[i][j][4][k])
#             ind_crit = np.argmax(res["M_sum"][i][j])
#             res["ind_crit"][i][j] = ind_crit
#             res["R_OM"][i][j] = M_R_stable[i][j][1]
#             res["R_DM"][i][j] = M_R_stable[i][j][2]
#             res["R_max"][i][j]=max(M_R_stable[i][j][1],M_R_stable[i][j][2])
#             res["M_max"][i][j] = res["M_sum"][i][j][ind_crit]
#             res["M_max_OM"][i][j] = M_R_stable[i][j][3][ind_crit]
#             res["M_max_DM"][i][j] = M_R_stable[i][j][4][ind_crit]
#             res["R_crit"][i][j] = res["R_max"][i][j][ind_crit]
#             res["R_crit_OM"][i][j] =  M_R_stable[i][j][1][ind_crit]
#             res["R_crit_DM"][i][j] =  M_R_stable[i][j][2][ind_crit]
#             res["k2_max"][i][j] = max(M_R_stable[i][j][7])
#             res["k2_crit"][i][j] = M_R_stable[i][j][7][ind_crit]
#             res["C_max"][i][j] = max(M_R_stable[i][j][6])
#             res["C_crit"][i][j] = M_R_stable[i][j][6][ind_crit]
#             res["p_crit"][i][j] = M_R_stable[i][j][0][ind_crit]
#     return res



# def combine_stables(M_R_stable):
#     '''
#         return an array which connects the stable branches of M_R_stable to give a shape of [dm_om][m_DM][ops][p0].
#     '''
#     res = []
#     for i in range(len(M_R_stable)):
#         res.append([])
#         for j in range(len(M_R_stable[i])):
#             res[i].append([])
#             for k in range(len(M_R_stable[i][j])):                  # stable
#                 for m in range(len(M_R_stable[i][j][k])):           # obs
#                     res[i][j].append([])
#                     for l in range(len(M_R_stable[i][j][k][m])):             # p0
#                         res[i][j][m].append(M_R_stable[i][j][k][m][l])
#     return res




# def get_stable_M_R_dm_om_M_arr():
#     '''
#         returns an array for the stable branches of all M_DM and dm_om of shape [dm_om][m_DM][stability][ops][p0].
#     '''
#     M_R_stable = []
#     for i in range(len(dm_om_arr)):
#         M_R_stable.append([])
#         for j in range(len(M_arr)):
#             M_R_stable[i].append([])
#             mass_DM = M_arr[j]
#             con_mass = pow(Planck_Mass, 3)*pow(mass_DM*MeV_solar, -2)
#             con_radius = pow(Planck_Mass, 2)*pow(mass_DM*MeV_solar, -2)*Planck_Length               #1.769 * 1e-76# *neutron_mass_fm**2
#             M_R = np.genfromtxt("results/full_EoS_%s_run/M_R_TIDAL_P_dm_om_%1.3e_M_%1.3e_2_fluid.out"%(EoS_str,dm_om_arr[i], M_arr[j]))
#             for ii in range(len(M_R)):
#                 M_R[ii][1] = con_radius*M_R[ii][1]
#                 M_R[ii][2] = con_radius*M_R[ii][2]
#                 M_R[ii][3] = con_mass*M_R[ii][3]
#                 M_R[ii][4] = con_mass*M_R[ii][4]
#             tmp = np.genfromtxt("results/full_EoS_%s_run/stability_P_dm_om_%1.3e_M_%1.3e_2_fluid"%(EoS_str,dm_om_arr[i], M_arr[j]))[1:]
#             stability = []
#             for k in range(len(tmp)):
#                 if (tmp[k][0] < EoS_max(mass_DM)):
#                     if (tmp[k][1] <= EoS_max(mass_DM)):
#                         stability.append(list(tmp[k]))
#                     else:
#                         stability.append([tmp[k][0],EoS_max(mass_DM)])
#             for m in range(len(stability)):               # num of stable in branch
#                 M_R_stable[i][j].append([])
#                 for k in range(len(M_R[0])):              # num of obs  
#                     M_R_stable[i][j][m].append([])
#                     for l in range(len(M_R)):             # num of P0
#                         if M_R[l][0] > stability[m][0] and M_R[l][0] < stability[m][1]:
#                             M_R_stable[i][j][m][k].append(M_R[l][k])
#     return M_R_stable
