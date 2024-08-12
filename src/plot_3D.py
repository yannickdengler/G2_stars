import matplotlib.pyplot as plt
import numpy as np
import sys

from matplotlib import cm
from matplotlib.ticker import LinearLocator

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

dm_om_arr=(0.717251, 0.933478, 1.214891, 1.581139, 2.057799, 2.678155, 3.485528, 4.536296, 5.903836, 7.683642, 10.000000)
M_arr=(250.000000, 304.753414, 371.498572, 452.861832, 552.044757, 672.950096, 820.335356, 1000.000000, 1219.013654, 1485.994289, 1811.447329, 2208.179027, 2691.800385, 3281.341424, 4000.000000)

M_max_arr = np.zeros((len(dm_om_arr), len(M_arr)))
M_crit_OM_arr = np.zeros((len(dm_om_arr), len(M_arr)))
M_crit_DM_arr = np.zeros((len(dm_om_arr), len(M_arr)))
R_crit_max_arr = np.zeros((len(dm_om_arr), len(M_arr)))
R_crit_OM = np.zeros((len(dm_om_arr), len(M_arr)))
R_crit_DM = np.zeros((len(dm_om_arr), len(M_arr)))
k_2_crit_arr = np.zeros((len(dm_om_arr), len(M_arr)))
k_2_max_arr = np.zeros((len(dm_om_arr), len(M_arr)))
C_crit_arr = np.zeros((len(dm_om_arr), len(M_arr)))
comp = np.zeros((len(dm_om_arr), len(M_arr)))

for i in range(len(dm_om_arr)):
    for j in range(len(M_arr)):
        mass_DM = M_arr[j]
        # print(i,j)
        # print(mass_DM)
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
            stability = np.genfromtxt("results/full_EoS_I_run/stability_P_dm_om_%1.3e_M_%1.3e_2_fluid"%(dm_om_arr[i], M_arr[j]))
        elif EoS == 2:
            stability = np.genfromtxt("results/full_EoS_II_run/stability_P_dm_om_%1.3e_M_%1.3e_2_fluid"%(dm_om_arr[i], M_arr[j]))
        elif EoS == 3:
            stability = np.genfromtxt("results/full_EoS_III_run/stability_P_dm_om_%1.3e_M_%1.3e_2_fluid"%(dm_om_arr[i], M_arr[j]))
        M_R_stable = []
        for stables in stability:
            for k in range(len(M_R)):
                if M_R[k][0] > stables[0] and M_R[k][0] < stables[1]:
                    M_R_stable.append(M_R[k])
        M_sum = [x[3]+ x[4] for x in M_R]
        M_sum_stable = [x[3]+ x[4] for x in M_R_stable]
        ind_crit = np.argmax(M_sum_stable)
        R_max = [max([x[1],x[2]]) for x in M_R]
        R_max_stable = [max([x[1],x[2]]) for x in M_R_stable]
        M_max_arr[i][j] = M_sum_stable[ind_crit]
        R_crit_max_arr[i][j] = max([M_R_stable[ind_crit][1],M_R_stable[ind_crit][2]])
        M_crit_OM_arr[i][j] = M_R_stable[ind_crit][3]
        M_crit_DM_arr[i][j] = M_R_stable[ind_crit][4]
        R_crit_OM[i][j] = M_R_stable[ind_crit][1]
        R_crit_DM[i][j] = M_R_stable[ind_crit][2]
        C_crit_arr[i][j] = M_R_stable[ind_crit][6]
        k_2_crit_arr[i][j] = M_R_stable[ind_crit][7]
        k_2_max_arr[i][j] = max(np.transpose(M_R_stable)[7])
        comp[i][j] = M_max_arr[i][j]/R_crit_max_arr[i][j]


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