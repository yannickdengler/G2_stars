import routines

routines.write_stability_file("P_dm_om_%1.3e_M_%1.3e_2_fluid"%(0, 1000),"/full_EoS_I_run")
routines.write_stability_file("P_dm_om_%1.3e_M_%1.3e_2_fluid"%(0, 1000),"/full_EoS_II_run")
routines.write_stability_file("P_dm_om_%1.3e_M_%1.3e_2_fluid"%(0, 1000),"/full_EoS_III_run")
routines.write_stability_file("P_dm_om_%1.3e_M_%1.3e_2_fluid"%(-1, 1000),"/full_EoS_I_run")
routines.write_stability_file("P_dm_om_%1.3e_M_%1.3e_2_fluid"%(-1, 1000),"/full_EoS_II_run")
routines.write_stability_file("P_dm_om_%1.3e_M_%1.3e_2_fluid"%(-1, 1000),"/full_EoS_III_run")
routines.write_stability_file("P_dm_om_%1.3e_M_%1.3e_2_fluid"%(1, 1000),"/full_EoS_I_run")
routines.write_stability_file("P_dm_om_%1.3e_M_%1.3e_2_fluid"%(1, 1000),"/full_EoS_II_run")
routines.write_stability_file("P_dm_om_%1.3e_M_%1.3e_2_fluid"%(1, 1000),"/full_EoS_III_run")

dm_om_arr=(0.250000, 0.325367, 0.423453, 0.551110, 0.717251, 0.933478, 1.214891, 1.581139, 2.057799, 2.678155, 3.485528, 4.536296, 5.903836, 7.683642, 10.000000)
M_arr=(250.000000, 304.753414, 371.498572, 452.861832, 552.044757, 672.950096, 820.335356, 1000.000000, 1219.013654, 1485.994289, 1811.447329, 2208.179027, 2691.800385, 3281.341424, 4000.000000)

for dm_om in dm_om_arr:
   for M_DM in M_arr:
       routines.write_stability_file("P_dm_om_%1.3e_M_%1.3e_2_fluid"%(dm_om, M_DM),"/full_EoS_I_run")
       routines.write_stability_file("P_dm_om_%1.3e_M_%1.3e_2_fluid"%(dm_om, M_DM),"/full_EoS_II_run")
       routines.write_stability_file("P_dm_om_%1.3e_M_%1.3e_2_fluid"%(dm_om, M_DM),"/full_EoS_III_run")

