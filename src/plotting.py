import plot
import sys
import global_vars as g
import matplotlib.pyplot as plt
import numpy as np
import sys

from matplotlib import cm
from matplotlib.ticker import LinearLocator

MeV_fm = 197.32698  
def EoS_I_max(m_DM):
    return 582.1/m_DM**4*MeV_fm**3
def EoS_II_max(m_DM):
    return 238.7/m_DM**4*MeV_fm**3
def EoS_III_max(m_DM):
    return 136.8/m_DM**4*MeV_fm**3


EoS = int(sys.argv[1])
plot_max_mass=False
plot_k2_crit=False
plot_k2_max=False
plot_r_crit=False
plot_comp=False
if sys.argv[2] == "max_mass":
    plot_max_mass = True
elif sys.argv[2] == "k2_crit":
    plot_k2_crit = True
elif sys.argv[2] == "k2_max":
    plot_k2_max = True
elif sys.argv[2] == "r_crit":
    plot_r_crit = True
elif sys.argv[2] == "comp":
    plot_comp = True

MeV_fm = 197.32698  
MeV_solar = 8.97*1e-61                                                      # 1 MeV in solar masses
Planck_Mass = 1.095*1e-38                                                   # in solar masses
Planck_Length = 1.6163e-38                                                  # in km

dm_om_arr=(0.250000, 0.325367, 0.423453, 0.551110, 0.717251, 0.933478, 1.214891, 1.581139, 2.057799, 2.678155, 3.485528, 4.536296, 5.903836, 7.683642, 10.000000)
M_arr=(250.000000, 304.753414, 371.498572, 452.861832, 552.044757, 672.950096, 820.335356, 1000.000000, 1219.013654, 1485.994289, 1811.447329, 2208.179027, 2691.800385, 3281.341424, 4000.000000)

# M_max_arr = np.zeros((len(dm_om_arr), len(M_arr)))
# M_crit_OM_arr = np.zeros((len(dm_om_arr), len(M_arr)))
# M_crit_DM_arr = np.zeros((len(dm_om_arr), len(M_arr)))
# R_crit_max_arr = np.zeros((len(dm_om_arr), len(M_arr)))
# R_crit_OM = np.zeros((len(dm_om_arr), len(M_arr)))
# R_crit_DM = np.zeros((len(dm_om_arr), len(M_arr)))
# k_2_crit_arr = np.zeros((len(dm_om_arr), len(M_arr)))
# k_2_max_arr = np.zeros((len(dm_om_arr), len(M_arr)))
# C_crit_arr = np.zeros((len(dm_om_arr), len(M_arr)))
# comp = np.zeros((len(dm_om_arr), len(M_arr)))

def get_stable_M_R_dm_om_M_arr(EoS):
    '''
        return an array of shape [dm_om][m_DM][stability][ops][p0].
    '''
    M_R_stable = []
    for i in range(len(dm_om_arr)):
        M_R_stable.append([])
        for j in range(len(M_arr)):
            M_R_stable[i].append([])
            mass_DM = M_arr[j]
            con_mass = pow(Planck_Mass, 3)*pow(mass_DM*MeV_solar, -2)
            con_radius = pow(Planck_Mass, 2)*pow(mass_DM*MeV_solar, -2)*Planck_Length               #1.769 * 1e-76# *neutron_mass_fm**2
            if EoS == 1:
                M_R = np.genfromtxt("results/full_EoS_I_run/M_R_TIDAL_P_dm_om_%1.3e_M_%1.3e_2_fluid.out"%(dm_om_arr[i], M_arr[j]))
            elif EoS == 2:
                M_R = np.genfromtxt("results/full_EoS_II_run/M_R_TIDAL_P_dm_om_%1.3e_M_%1.3e_2_fluid.out"%(dm_om_arr[i], M_arr[j]))
            elif EoS == 3:
                M_R = np.genfromtxt("results/full_EoS_III_run/M_R_TIDAL_P_dm_om_%1.3e_M_%1.3e_2_fluid.out"%(dm_om_arr[i], M_arr[j]))
            for ii in range(len(M_R)):
                M_R[ii][1] = con_radius*M_R[ii][1]
                M_R[ii][2] = con_radius*M_R[ii][2]
                M_R[ii][3] = con_mass*M_R[ii][3]
                M_R[ii][4] = con_mass*M_R[ii][4]
            if EoS == 1:
                tmp = np.genfromtxt("results/full_EoS_I_run/stability_P_dm_om_%1.3e_M_%1.3e_2_fluid"%(dm_om_arr[i], M_arr[j]))[1:]
                stability = []
                EoS_max = EoS_I_max(mass_DM)
                for k in range(len(tmp)):
                    if (tmp[k][0] < EoS_max):
                        if (tmp[k][1] <= EoS_max):
                            stability.append(list(tmp[k]))
                        else:
                            stability.append([tmp[k][0],EoS_max])
            if EoS == 2:
                tmp = np.genfromtxt("results/full_EoS_II_run/stability_P_dm_om_%1.3e_M_%1.3e_2_fluid"%(dm_om_arr[i], M_arr[j]))[1:]
                stability = []
                EoS_max = EoS_II_max(mass_DM)
                for k in range(len(tmp)):
                    if (tmp[k][0] < EoS_max):
                        if (tmp[k][1] <= EoS_max):
                            stability.append(list(tmp[k]))
                        else:
                            stability.append([tmp[k][0],EoS_max])
            if EoS == 3:
                tmp = np.genfromtxt("results/full_EoS_III_run/stability_P_dm_om_%1.3e_M_%1.3e_2_fluid"%(dm_om_arr[i], M_arr[j]))[1:]
                stability = []
                EoS_max = EoS_III_max(mass_DM)
                for k in range(len(tmp)):
                    if (tmp[k][0] < EoS_max):
                        if (tmp[k][1] <= EoS_max):
                            stability.append(list(tmp[k]))
                        else:
                            stability.append([tmp[k][0],EoS_max])
            for m in range(len(stability)):               # num of stable in branch
                M_R_stable[i][j].append([])
                for k in range(len(M_R[0])):              # num of obs  
                    M_R_stable[i][j][m].append([])
                    for l in range(len(M_R)):             # num of P0
                        if M_R[l][0] > stability[m][0] and M_R[l][0] < stability[m][1]:
                            M_R_stable[i][j][m][k].append(M_R[l][k])
    return M_R_stable

def combine_stables(M_R_stable):
    '''
        return an array which connects the stable branches of M_R_stable to give a shape of [dm_om][m_DM][ops][p0].
    '''
    res = []
    # print(len(M_R_stable))
    # print(len(M_R_stable[0]))
    # print(len(M_R_stable[0][0]))
    # print(len(M_R_stable[0][0][0]))
    # print(len(M_R_stable[0][0][0][0]))
    # print(M_R_stable[0][0][0][0][0])
    # print()
    for i in range(len(M_R_stable)):
        res.append([])
        for j in range(len(M_R_stable[i])):
            res[i].append([])
            for k in range(len(M_R_stable[i][j])):                  # stable
                for m in range(len(M_R_stable[i][j][k])):           # obs
                    res[i][j].append([])
                    for l in range(len(M_R_stable[i][j][k][m])):             # p0
                        res[i][j][m].append([])
                        res[i][j][m][l].append(M_R_stable[i][j][k][m][l])
    print(len(res))
    print(len(res[0]))
    print(len(res[0][0]))
    print(len(res[0][0][0]))
    print()
    # exit()
    return res

# whatever = [[[] for i in range(3)] for j in range(3)]

# print(whatever)

# exit()

def get_data(M_R_stable):
    print(len(M_R_stable))
    print(len(M_R_stable[0]))
    print(len(M_R_stable[0][0]))
    res = {}
    res["M_OM"]=[[[] for i in range(len(M_R_stable))] for j in range(len(M_R_stable[0]))]
    res["M_DM"]=[ []*len(M_R_stable) for i in range(len(M_R_stable[0]))]
    res["M_sum"]=[ []*len(M_R_stable) for i in range(len(M_R_stable[0]))]
    res["ind_crit"]=[ []*len(M_R_stable) for i in range(len(M_R_stable[0]))]
    res["p_crit"]=[ []*len(M_R_stable) for i in range(len(M_R_stable[0]))]
    res["R_OM"]=[ []*len(M_R_stable) for i in range(len(M_R_stable[0]))]
    res["R_DM"]=[ []*len(M_R_stable) for i in range(len(M_R_stable[0]))]
    res["R_max"]=[ []*len(M_R_stable) for i in range(len(M_R_stable[0]))]
    res["M_max"]=[ []*len(M_R_stable) for i in range(len(M_R_stable[0]))]
    res["M_max_OM"]=[ []*len(M_R_stable) for i in range(len(M_R_stable[0]))]
    res["M_max_DM"]=[ []*len(M_R_stable) for i in range(len(M_R_stable[0]))]
    res["R_crit"]=[ []*len(M_R_stable) for i in range(len(M_R_stable[0]))]
    res["R_crit_OM"]=[ []*len(M_R_stable) for i in range(len(M_R_stable[0]))]
    res["R_crit_DM"]=[ []*len(M_R_stable) for i in range(len(M_R_stable[0]))]
    res["k2_max"]=[ []*len(M_R_stable) for i in range(len(M_R_stable[0]))]
    res["k2_crit"]=[ []*len(M_R_stable) for i in range(len(M_R_stable[0]))]
    res["C_max"]=[ []*len(M_R_stable) for i in range(len(M_R_stable[0]))]
    res["C_crit"]=[ []*len(M_R_stable) for i in range(len(M_R_stable[0]))]
    print(res["M_OM"])

    for i in range(len(M_R_stable)):
        for j in range(len(M_R_stable[i])):
            print(i,j)
            res["M_OM"][i][j] = M_R_stable[i][j][3]
            # res["M_DM"][i][j] = M_R_stable[i][j][4]
            # res["M_sum"][i][j] = [x[3]+x[4] for x in M_R_stable[i][j]]
            # ind_crit = np.argmax(res["M_sum"][i][j])
            # res["ind_crit"][i][j] = ind_crit
            # res["R_OM"][i][j] = M_R_stable[i][j][1]
            # res["R_DM"][i][j] = M_R_stable[i][j][2]
            # res["R_max"][i][j] = [max(x[1],x[2]) for x in M_R_stable[i][j]]
            # res["M_max"][i][j] = res["M_sum"][i][j][ind_crit]
            # res["M_max_OM"][i][j] = M_R_stable[i][j][3][ind_crit]
            # res["M_max_DM"][i][j] = M_R_stable[i][j][4][ind_crit]
            # res["R_crit"][i][j] = res["R_max"][i][j][3][ind_crit]
            # res["R_crit_OM"][i][j] =  M_R_stable[i][j][1][ind_crit]
            # res["R_crit_DM"][i][j] =  M_R_stable[i][j][2][ind_crit]
            # res["k2_max"][i][j] = max(M_R_stable[i][j][7])
            # res["k2_crit"][i][j] = M_R_stable[i][j][7][ind_crit]
            # res["C_max"][i][j] = max(M_R_stable[i][j][6])
            # res["C_crit"][i][j] = M_R_stable[i][j][6][ind_crit]
            # res["p_crit"][i][j] = M_R_stable[i][j][0][ind_crit]
    return res

        # M_sum = [x[3]+ x[4] for x in M_R]
        # M_sum_stable = [x[3]+x[4] for x in M_R_stable]
        # ind_crit = np.argmax(M_sum_stable)
        # R_max = [max([x[1],x[2]]) for x in M_R]
        # R_max_stable = [max([x[1],x[2]]) for x in M_R_stable]
        # M_max_arr[i][j] = M_sum_stable[ind_crit]
        # R_crit_max_arr[i][j] = max([M_R_stable[ind_crit][1],M_R_stable[ind_crit][2]])
        # M_crit_OM_arr[i][j] = M_R_stable[ind_crit][3]
        # M_crit_DM_arr[i][j] = M_R_stable[ind_crit][4]
        # R_crit_OM[i][j] = M_R_stable[ind_crit][1]
        # R_crit_DM[i][j] = M_R_stable[ind_crit][2]
        # C_crit_arr[i][j] = M_R_stable[ind_crit][6]
        # k_2_crit_arr[i][j] = M_R_stable[ind_crit][7]
        # k_2_max_arr[i][j] = max(np.transpose(M_R_stable)[7])
        # comp[i][j] = M_max_arr[i][j]/R_crit_max_arr[i][j]

if __name__ == "__main__":
    stable = get_stable_M_R_dm_om_M_arr(1)
    stable_c = combine_stables(stable)
    print("________-", stable_c[0][0])
    data = get_data(stable_c)
    # print(data)
    for key, val in data:
        print(key, val)    


    exit()

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Make data.
X = np.arange(-5, 5, 1)
Y = np.arange(-5, 5, .5)
X, Y = np.meshgrid(X, Y)
R = np.sqrt(X**2 + Y**2)
Z = np.sin(R)

# Make data.
M_plot, dm_om_plot = np.meshgrid(M_arr, dm_om_arr)


# for i in range(len(dm_om_plot)):
#     for j in range(len(dm_om_plot[i])):
#         print(M_plot[i][j], dm_om_plot[i][j], M_max_arr[i][j])

print(M_max_arr[0][7],np.log10(M_max_arr[0][13]))
# print(dm_om_plot)
# print(M_plot)
# print(M_max_arr)


# Plot the surface.
if plot_max_mass:
    surf = ax.plot_surface(np.log10(M_plot), np.log10(dm_om_plot), np.log10(M_max_arr), cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax.set_zlabel("$M_{max}$")
elif plot_k2_crit:
    surf = ax.plot_surface(np.log10(M_plot), np.log10(dm_om_plot), k_2_crit_arr, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax.set_zlabel("$k2_{crit}$")
elif plot_k2_max:
    surf = ax.plot_surface(np.log10(M_plot), np.log10(dm_om_plot), k_2_max_arr, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax.set_zlabel("$k2_{max}$")
elif plot_r_crit:
    surf = ax.plot_surface(np.log10(M_plot), np.log10(dm_om_plot), np.log10(R_crit_OM), cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax.set_zlabel("$R_{crit}$")
elif plot_comp:
    surf = ax.plot_surface(np.log10(M_plot), np.log10(dm_om_plot), np.log10(comp), cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax.set_zlabel("C")

ax.set_xlabel("M Dm")
ax.set_ylabel("dm_om")

# Customize the z axis.
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter('{x:.02f}')


# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

if plot_max_mass:
    plt.savefig("plots/max_mass_EoS_%s.pdf"%EoS)
elif plot_k2_crit:
    plt.savefig("plots/k2_crit_EoS_%s.pdf"%EoS)
elif plot_k2_max:
    plt.savefig("plots/k2_max_EoS_%s.pdf"%EoS)
elif plot_r_crit:
    plt.savefig("plots/r_crit_EoS_%s.pdf"%EoS)
elif plot_comp:
    plt.savefig("plots/comp_crit_EoS_%s.pdf"%EoS)

plt.show()