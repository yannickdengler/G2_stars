import matplotlib.pyplot as plt
import numpy as np

from matplotlib import cm
from matplotlib.ticker import LinearLocator

MeV_fm = 197.32698  
MeV_solar = 8.97*1e-61                                                      # 1 MeV in solar masses
Planck_Mass = 1.095*1e-38                                                   # in solar masses
Planck_Length = 1.6163e-38                                                  # in km

dm_om_arr=(0.25, 0.35355339059, 0.5, 0.70710678118, 1, 1.41421356237, 2, 2.82842712475, 4)
M_arr=(250, 353.55339059, 500, 707.10678118, 1000, 1414.21356237, 2000, 2828.42712475, 4000)

M_max_arr = np.zeros((len(dm_om_arr), len(M_arr)))
M_crit_OM_arr = np.zeros((len(dm_om_arr), len(M_arr)))
M_crit_DM_arr = np.zeros((len(dm_om_arr), len(M_arr)))
R_crit_max_arr = np.zeros((len(dm_om_arr), len(M_arr)))
R_crit_OM = np.zeros((len(dm_om_arr), len(M_arr)))
R_crit_DM = np.zeros((len(dm_om_arr), len(M_arr)))
k_2_crit_arr = np.zeros((len(dm_om_arr), len(M_arr)))
k_2_max_arr = np.zeros((len(dm_om_arr), len(M_arr)))
C_crit_arr = np.zeros((len(dm_om_arr), len(M_arr)))

# fig, ax = plt.subplots(9,9)

for i in range(len(dm_om_arr)):
    for j in range(len(M_arr)):
        mass_DM = M_arr[j]
        print(i,j)
        print(mass_DM)
        con_mass = pow(Planck_Mass, 3)*pow(mass_DM*MeV_solar, -2)
        con_radius = pow(Planck_Mass, 2)*pow(mass_DM*MeV_solar, -2)*Planck_Length               #1.769 * 1e-76# *neutron_mass_fm**2
        M_R = np.genfromtxt("results/M_R_TIDAL_P_dm_om_%1.3e_M_%1.3e_2_fluid.out"%(dm_om_arr[i], M_arr[j]))
        print(len(M_R))
        print("ääääääääääääääääääää")
        for ii in range(len(M_R)):
            M_R[ii][1] = con_radius*M_R[ii][1]
            M_R[ii][2] = con_radius*M_R[ii][2]
            M_R[ii][3] = con_mass*M_R[ii][3]
            M_R[ii][4] = con_mass*M_R[ii][4]
        stability = np.genfromtxt("results/stability_P_dm_om_%1.3e_M_%1.3e_2_fluid"%(dm_om_arr[i], M_arr[j]))
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
        # print(np.transpose(M_R_stable)[3])
        # print(np.transpose(M_R_stable)[4])

        # for m in range(len(M_sum_stable)):
        #     print(m,M_sum_stable[m])
        # # print(M_sum_stable)
        # print(ind_crit)
        M_max_arr[i][j] = M_sum_stable[ind_crit]
        R_crit_max_arr[i][j] = max([M_R_stable[ind_crit][1],M_R_stable[ind_crit][2]])
        M_crit_OM_arr[i][j] = M_R_stable[ind_crit][3]
        M_crit_DM_arr[i][j] = M_R_stable[ind_crit][4]
        R_crit_OM[i][j] = M_R_stable[ind_crit][1]
        R_crit_DM[i][j] = M_R_stable[ind_crit][2]
        C_crit_arr[i][j] = M_R_stable[ind_crit][6]
        k_2_crit_arr[i][j] = M_R_stable[ind_crit][7]
        k_2_max_arr[i][j] = max(np.transpose(M_R_stable)[7])

        # ax[i][j].plot(R_max_stable,np.transpose(M_R_stable)[3], label="M OM stable")
        # ax[i][j].plot(R_max_stable,np.transpose(M_R_stable)[4], label="M DM stable")
        # ax[i][j].plot(R_max_stable,M_sum_stable, label="M sum stable")
        # ax[i][j].plot(R_max,np.transpose(M_R)[3], label="M OM", ls = "dotted")
        # ax[i][j].plot(R_max,np.transpose(M_R)[4], label="M DM", ls = "dotted")
        # ax[i][j].plot(R_max,M_sum, label="M sum", ls = "dotted")
        # ax[i][j].set_title("%1.2e %1.2e"%(dm_om_arr[i], M_arr[j]))

# plt.show()
# # plt.legend()/

# plt.clf()




fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Make data.
X = np.arange(-5, 5, 1)
Y = np.arange(-5, 5, .5)
X, Y = np.meshgrid(X, Y)
R = np.sqrt(X**2 + Y**2)
Z = np.sin(R)
# print(R)
# print(Z)
# print(Z.shape)

# # Plot the surface.
# surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False)


# Make data.
M_plot, dm_om_plot = np.meshgrid(M_arr, dm_om_arr)

print(dm_om_plot)
print(M_plot)
print(M_max_arr)


# Plot the surface.
# surf = ax.plot_surface(np.log10(M_plot), np.log10(dm_om_plot), np.log10(M_max_arr), cmap=cm.coolwarm, linewidth=0, antialiased=False)
# surf = ax.plot_surface(np.log10(M_plot), np.log10(dm_om_plot), np.log10(R_crit_OM), cmap=cm.coolwarm, linewidth=0, antialiased=False)
# surf = ax.plot_surface(np.log10(M_plot), np.log10(dm_om_plot), k_2_crit_arr, cmap=cm.coolwarm, linewidth=0, antialiased=False)
surf = ax.plot_surface(np.log10(M_plot), np.log10(dm_om_plot), k_2_max_arr, cmap=cm.coolwarm, linewidth=0, antialiased=False)

ax.set_ylabel("dm_om")
ax.set_xlabel("M Dm")

# Customize the z axis.
ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')


# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()