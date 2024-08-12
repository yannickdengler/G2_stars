import routines
import sys
import global_vars as g

# routines.print_EoS_P_2_fluid(1e-2,1e8,5000, pref="_EoS_III")
# exit()

import matplotlib.pyplot as plt
import numpy as np

EoS_1 = np.transpose(np.genfromtxt("results/EoS_P_2_fluid_EoS_I.out"))
EoS_2 = np.transpose(np.genfromtxt("results/EoS_P_2_fluid_EoS_II.out"))
EoS_3 = np.transpose(np.genfromtxt("results/EoS_P_2_fluid_EoS_III.out"))

plt.plot(EoS_1[0], EoS_1[1], label = "1-OM")
# plt.plot(EoS_1[0], EoS_1[3], label = "1-DM", ls = "--")
plt.plot(EoS_2[0], EoS_2[1], label = "2-OM")
# plt.plot(EoS_2[0], EoS_2[3], label = "2-DM", ls = "--")
plt.plot(EoS_3[0], EoS_3[1], label = "3-OM")
# plt.plot(EoS_3[0], EoS_3[3], label = "3-DM", ls = "--")

# print(len(EoS_1))
# print(len(EoS_1[0]))

plt.xscale("log")
plt.yscale("log")

plt.legend()
plt.show()

# P_start = 1e-3/g.conv_P
# P_end = 1e3/g.conv_P
# dm_om = float(sys.argv[2])

# plot.plot_EoS_P_2_fluid("results/EoS_P_2_fluid.out")



# routines.M_R_TIDAL_2_fluid(var_start=P_start,var_end=P_end, steps=100, dm_om_ratio=dm_om, var="P")

# routines.single_P_2_fluid(P_start,0)
# import plot
# plot.plot_single_2_fluid("results/single_P_2.out")


# routines.print_EoS_P_2_fluid(1e-10,1e4,1000)

# plot.plot_EoS_P_2_fluid("results/EoS_P_2_fluid.out")