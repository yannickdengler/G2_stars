import routines

dm_om_arr=(0.25, 0.35355339059, 0.5, 0.70710678118, 1, 1.41421356237, 2, 2.82842712475, 4)
M_arr=(250, 353.55339059, 500, 707.10678118, 1000, 1414.21356237, 2000, 2828.42712475, 4000)
# dm_om_arr=(0.3, 0.4, 0.6, 0.85, 1.2, 1.7, 2.4, 3.4)
# M_arr=(300, 400, 600, 850, 1200, 1700, 2400, 3400)

for dm_om in dm_om_arr:
   for M_DM in M_arr:
       routines.write_stability_file("P_dm_om_%1.3e_M_%1.3e_2_fluid"%(dm_om, M_DM))