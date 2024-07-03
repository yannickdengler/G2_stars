dm_om_min = 0.25
dm_om_max = 10
dm_om_steps = 15

M_min = 250
M_max = 4000
M_steps = 15


with open("arr","w") as f:
    f.write("\ndm_om_arr=(")
    for i in range(dm_om_steps):
        f.write("%f"%(dm_om_min*(dm_om_max/dm_om_min)**(i/(dm_om_steps-1))))
        if i != dm_om_steps-1:
            f.write(" ")
    f.write(")\n")
    f.write("M_arr=(")
    for i in range(dm_om_steps):
        f.write("%f"%(M_min*(M_max/M_min)**(i/(M_steps-1))))
        if i != M_steps-1:
            f.write(" ")
    f.write(")\n")
    