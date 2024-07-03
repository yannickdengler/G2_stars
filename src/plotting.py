import plot
import sys
import global_vars as g

dm_om = float(sys.argv[2])


# plot.plot_EoS_P_2_fluid("results/EoS_P_2_fluid.out", log=1)
plot.plot_M_R_2_fluid("results/M_R_TIDAL_P_dm_om_%1.3e_M_%1.3e_2_fluid.out"%(dm_om, g.mass_DM), log = 1)