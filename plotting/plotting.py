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
    def EoS_max(m_DM,dm_om):
        if dm_om == -1:
            return 1e50
        return 582.1/m_DM**4*MeV_fm**3
elif EoS == 2:
    def EoS_max(m_DM,dm_om):
        if dm_om == -1:
            return 1e50
        return 238.7/m_DM**4*MeV_fm**3
    EoS_str = "II"
elif EoS == 3:
    def EoS_max(m_DM,dm_om):
        if dm_om == -1:
            return 1e50
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
        if (tmp[k][0] < EoS_max(M_DM,dm_om)):
            if (tmp[k][1] <= EoS_max(M_DM,dm_om)):
                stability.append(list(tmp[k]))
            else:
                stability.append([tmp[k][0],EoS_max(M_DM,dm_om)])
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

# def combine_stables_dm_om_M_DM(M_R_stable):
#     res = []
#     for i in range(len(M_R_stable)):
#         res.append([])
#         for j in range(len(M_R_stable[i])):
#             res[i].append([])
#             combine = combine_stables(M_R_stable[i][j])
#             for k in range(len(combine)):
#                 res[i][j].append(combine[k])
#     return res

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
    res["M_crit_OM"]=M_R_stable[3][ind_crit]
    res["M_crit_DM"]=M_R_stable[4][ind_crit]
    res["R_crit"]=res["R_max"][ind_crit]
    res["R_crit_OM"]=M_R_stable[1][ind_crit]
    res["R_crit_DM"]=M_R_stable[2][ind_crit]
    res["k2"]=M_R_stable[7]
    res["C"]=M_R_stable[6]
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

def plot_3D(arr, label, save=False, show=False):
    '''
        Creates a 3D plot in x-m_DM, y-dm_om, z-"arr" (with label=label), arr as [dm_om][m_DM]
    '''
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    M_plot, dm_om_plot = np.meshgrid(M_arr, dm_om_arr)

    surf = ax.plot_surface(np.log10(M_plot), np.log10(dm_om_plot), np.log10(arr), cmap=cm.coolwarm, linewidth=0, antialiased=False)

    # ax.set_zlabel(label)
    ax.set_xlabel("$M_{DM}$")
    ax.set_ylabel("$p_{0,DM}/p_{0,OM}$")
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter('{x:.02f}')

    fig.colorbar(surf, shrink=0.5, aspect=5, label=label)
    if save:
        plt.savefig("plots/full_EoS_%s_run/%s_3D_EoS_%s.pdf"%(EoS_str,label,EoS_str))
    if show:
        plt.show()
    ax.cla()
    plt.close(fig)

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

def plot_M_R_dm_om_keep_stable(data, data_OM, data_DM, M_ind, show=True, save = True,what_R="max",what_M="sum"):
    if what_R == "max":
        R_key = "R_max"
    elif what_R == "OM":
        R_key = "R_OM"
    elif what_R == "DM":
        R_key = "R_DM"
    if what_M == "sum":
        M_key = "M_sum"
    elif what_M == "OM":
        M_key = "M_OM"
    elif what_M == "DM":
        M_key = "M_DM"

    fig, ax = plt.subplots()
    norm = mpl.colors.Normalize(vmin=np.log10(min(dm_om_arr)),vmax=np.log10(max(dm_om_arr)))
    cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.viridis)
    cmap.set_array([])
    colors = cmap.to_rgba(np.linspace(np.log10(min(dm_om_arr)), np.log10(max(dm_om_arr)), len(dm_om_arr)))
    plt.grid()
    x_max = 0
    y_max = 0
    for i in range(len(data)):
        for j in range(len(data[i][M_ind])):
            x_max = max([x_max, max(data[i][M_ind][j][R_key])])
            y_max = max([y_max, max(data[i][M_ind][j][M_key])])
            plt.plot(data[i][M_ind][j][R_key], data[i][M_ind][j][M_key], color = colors[i])
    ax.plot(data_OM[R_key], data_OM[M_key], color = "black", label = "OM", ls = "dashed")
    ax.plot(data_DM[R_key], data_DM[M_key], color = "black", label = "DM", ls = "solid")
    ax.set_xlabel("$R_{%s}$"%what_R)
    ax.set_ylabel("$M_{%s}$"%what_M)
    ax.set_xlim([0,1.1*x_max])
    ax.set_ylim([0,1.1*y_max])
    cb = plt.colorbar(cmap)
    cb.set_label("$p_{0,DM}/p_{0,OM}$")
    plt.title("$m_{DM}=%1.1e MeV$"%M_arr[M_ind])
    plt.legend()
    if save:
        plt.savefig("plots/full_EoS_%s_run/M_%s_R_%s_EoS_%s_dm_om_M_DM_%e.pdf"%(EoS_str,what_M,what_R,EoS_str,M_arr[M_ind]))
    if show:
        plt.show()
    ax.cla()
    plt.close(fig)


def plot_k_2_C_dm_om_keep_stable(data, data_OM, data_DM, M_ind, show=True, save = True):
    fig, ax = plt.subplots()
    norm = mpl.colors.Normalize(vmin=np.log10(min(dm_om_arr)),vmax=np.log10(max(dm_om_arr)))
    cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.viridis)
    cmap.set_array([])
    colors = cmap.to_rgba(np.linspace(np.log10(min(dm_om_arr)), np.log10(max(dm_om_arr)), len(dm_om_arr)))
    plt.grid()
    x_max = 0
    y_max = 0
    for i in range(len(data)):
        for j in range(len(data[i][M_ind])):
            x_max = max([x_max, max(data[i][M_ind][j]["C"])])
            y_max = max([y_max, max(data[i][M_ind][j]["k2"])])
            plt.plot(data[i][M_ind][j]["C"], data[i][M_ind][j]["k2"], color = colors[i])
    ax.plot(data_OM["C"], data_OM["k2"], color = "black", label = "OM", ls = "dashed")
    ax.plot(data_DM["C"], data_DM["k2"], color = "black", label = "DM", ls = "solid")
    ax.set_xlabel("C")
    ax.set_ylabel("$k2$")
    ax.set_xlim([0,1.1*x_max])
    ax.set_ylim([0,1.1*y_max])
    cb = plt.colorbar(cmap)
    cb.set_label("$p_{0,DM}/p_{0,OM}$")
    plt.title("$m_{DM}=%1.1e MeV$"%M_arr[M_ind])
    plt.legend()
    if save:
        plt.savefig("plots/full_EoS_%s_run/k_2_C_EoS_%s_dm_om_M_DM_%e.pdf"%(EoS_str,EoS_str,M_arr[M_ind]))
    if show:
        plt.show()
    ax.cla()
    plt.close(fig)

def get_data_stable():
    stable = get_stable_M_R_dm_om_M_arr()

    data = []
    for i in range(len(stable)):             # dm_om
        data.append([])
        for j in range(len(stable[i])):      # m_dm
            data[i].append([])
            for k in range(len(stable[i][j])):  # stable
                data[i][j].append(get_data(stable[i][j][k]))
    return data

def get_data_stable_OM_DM(norm):
    stable_pure_OM, stable_pure_DM = get_stable_M_R_pure(norm=norm)
    stable_c_pure_OM = combine_stables(stable_pure_OM)
    data_OM = get_data(stable_c_pure_OM)
    stable_c_pure_DM = combine_stables(stable_pure_DM)
    data_DM = get_data(stable_c_pure_DM)
    return data_OM, data_DM

def plot_all_M_tot_R_tot():
    data = get_data_stable()
    for i in range(15):
        data_OM, data_DM = get_data_stable_OM_DM((M_arr[i]/1000)**(-2))
        plot_M_R_dm_om_keep_stable(data,data_OM,data_DM, i, show = False, save=True)
        plot_M_R_dm_om_keep_stable(data,data_OM,data_DM, i, show = False, save=True,what_R="OM",what_M="OM")
        plot_M_R_dm_om_keep_stable(data,data_OM,data_DM, i, show = False, save=True,what_R="DM",what_M="DM")

def plot_all_k_2_C():
    data = get_data_stable()
    for i in range(15):
        data_OM, data_DM = get_data_stable_OM_DM((M_arr[i]/1000)**(-2))
        plot_k_2_C_dm_om_keep_stable(data,data_OM,data_DM, i, show = False, save=True)

def get_arr(data, label, stab_max_inds):
    tmp = [[[data[i][j][k][label] for k in range(len(data[i][j]))] for j in range(len(M_arr))] for i in range(len(dm_om_arr))]

    arr = []
    for i in range(len(dm_om_arr)):
        arr.append([])
        for j in range(len(M_arr)):
            arr[i].append(tmp[i][j][stab_max_inds[i][j]])
    return arr

def plot_all_3D():
    data = get_data_stable()

    tmp = [[[data[i][j][k]["M_max"] for k in range(len(data[i][j]))] for j in range(len(M_arr))] for i in range(len(dm_om_arr))]

    stab_max_inds = []
    M_max_arr = []
    for i in range(len(dm_om_arr)):
        stab_max_inds.append([])
        M_max_arr.append([])
        for j in range(len(M_arr)):
            stab_max_inds[i].append(np.argmax(tmp[i][j]))
            M_max_arr[i].append(tmp[i][j][stab_max_inds[i][j]])
    
    plot_3D(M_max_arr, label="$M_{max}$", save=True, show=False)
    for key in ["R_crit", "k2_crit", "C_crit", "M_crit_OM", "M_crit_DM", "R_crit_OM", "R_crit_DM"]:
        plot_3D(get_arr(data, key, stab_max_inds), key, save=True, show=False)
        
    tmp = [[[data[i][j][k]["k2_max"] for k in range(len(data[i][j]))] for j in range(len(M_arr))] for i in range(len(dm_om_arr))]
    k2_max_arr = []
    for i in range(len(dm_om_arr)):
        k2_max_arr.append([])
        for j in range(len(M_arr)):
            k2_max_arr[i].append(max(tmp[i][j]))
    plot_3D(k2_max_arr, label="$k_{2max}$", save=True, show=False)
    



if __name__ == "__main__":
    plot_all_M_tot_R_tot()
    plot_all_k_2_C()
    plot_all_3D()

        





