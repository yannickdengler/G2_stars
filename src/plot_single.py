import plot
import sys
import global_vars as g

dm_om = float(sys.argv[2])

dm_om = 0
mass_DM = g.neutron_mass

# plot.plot_EoS_P_2_fluid("results/EoS_P_2_fluid.out", log=1)
# plot.plot_M_R_2_fluid("results/M_R_TIDAL_P_dm_om_%1.3e_M_%1.3e_2_fluid.out"%(dm_om, g.mass_DM), log = 1)
# plot.plot_M_R_2_fluid("results/M_R_TIDAL_P_dm_om_%1.3e_M_%1.3e_2_fluid.out"%(dm_om, mass_DM), log = 1, pref = "/full_EoS_II_run")

# plot.plot_M_R_2_fluid("results/full_EoS_II_run_new/M_R_TIDAL_P_dm_om_0.000e+00_M_1.000e+03_2_fluid.out", log = 1)

# plot.plot_EoS_P_2_fluid("results/EoS_P_2_fluid.out", log=1)
# plot.plot_EoS_P_2_fluid("results/EoS_P_2_fluid.out", log=1)

plot.plot_M_R_2_fluid("M_R_TIDAL_P_dm_om_0.000e+00_M_1.000e+03_2_fluid", log = 0, pref="full_EoS_I_run")
plot.plot_M_R_2_fluid("M_R_TIDAL_P_dm_om_0.000e+00_M_1.000e+03_2_fluid", log = 0, pref="full_EoS_II_run")
plot.plot_M_R_2_fluid("M_R_TIDAL_P_dm_om_0.000e+00_M_1.000e+03_2_fluid", log = 0, pref="full_EoS_III_run")