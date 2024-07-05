import routines
import sys
import global_vars as g



P_start = 1e-3/g.conv_P
P_end = 1e3/g.conv_P
dm_om = float(sys.argv[2])



routines.M_R_TIDAL_2_fluid(var_start=P_start,var_end=P_end, steps=100, dm_om_ratio=dm_om, var="P")


# routines.print_EoS_P_2_fluid(1e-5,1e4,1000)

# import plot

# plot.plot_EoS_P_2_fluid("results/EoS_P_2_fluid.out")