import numpy as np
import matplotlib.pyplot as plt
import global_vars as v

# import EoS
# import run

########################################### plot EoS

def plot_EoS_mu_1_fluid(filename, log=1):
    EoS = np.transpose(np.genfromtxt(filename))
    plt.plot(EoS[0],EoS[1], label="P")
    plt.plot(EoS[0],EoS[2], label="eps")
    plt.plot(EoS[0],EoS[3], label="c_s")
    plt.title("EoS(mu) 1 fluid")
    if log == 1:
        plt.xscale("log")
        plt.yscale("log")
    plt.xlabel("mu")
    plt.legend()
    plt.savefig("results/EoS_mu_1_fluid.pdf")
    plt.show()

def plot_EoS_P_1_fluid(filename, log=1):
    EoS = np.transpose(np.genfromtxt(filename))
    plt.plot(EoS[0],EoS[1], label="mu")
    plt.plot(EoS[0],EoS[2], label="eps")
    plt.plot(EoS[0],EoS[3], label="c_s")
    plt.title("EoS(P) 1 fluid")
    if log == 1:
        plt.xscale("log")
        plt.yscale("log")
    plt.xlabel("P")
    plt.legend()
    plt.savefig("results/EoS_P_1_fluid.pdf")
    plt.show()

def plot_EoS_mu_2_fluid(filename, log=1):
    EoS = np.transpose(np.genfromtxt(filename))
    plt.plot(EoS[0],EoS[1], label="P_OM")
    plt.plot(EoS[0],EoS[2], label="eps_OM")
    plt.plot(EoS[0],EoS[3], label="c_s_OM")
    plt.plot(EoS[0],EoS[4], label="P_DM")
    plt.plot(EoS[0],EoS[5], label="eps_DM")
    plt.plot(EoS[0],EoS[6], label="c_s_DM")
    plt.title("EoS(mu) 1 fluid")
    if log == 1:
        plt.xscale("log")
        plt.yscale("log")
    plt.xlabel("mu")
    plt.legend()
    plt.savefig("results/EoS_mu_2_fluid.pdf")
    plt.show()

def plot_EoS_P_2_fluid(filename, log=1):
    EoS = np.transpose(np.genfromtxt(filename))
    print(v.conv_P)
    # # plt.plot(EoS[0]*conv_P,EoS[1]*v.con_chem_pot, label="mu_OM", ls="-")
    # plt.plot(EoS[0]*v.conv_P,EoS[2]*v.conv_P, label="eps_OM in MeV fm^-3", ls="-")
    # plt.plot(EoS[0]*v.conv_P,EoS[3], label="c_s_OM", ls="-")
    # # plt.plot(EoS[0]*conv_P,EoS[4]*v.con_chem_pot, label="mu_DM", ls="--")
    # plt.plot(EoS[0]*v.conv_P,EoS[5]*v.conv_P, label="eps_DM in MeV fm^-3", ls="--")
    # plt.plot(EoS[0]*v.conv_P,EoS[6], label="c_s_DM", ls="--")

    plt.plot(EoS[0]*v.conv_P,EoS[1]*v.conv_P, label="eps_OM in MeV fm^-3", ls="-")
    plt.plot(EoS[0]*v.conv_P,EoS[2], label="c_s_OM", ls="-")
    plt.plot(EoS[0]*v.conv_P,EoS[3]*v.conv_P, label="eps_DM in MeV fm^-3", ls="--")
    plt.plot(EoS[0]*v.conv_P,EoS[4], label="c_s_DM", ls="--")

    plt.plot((min(EoS[0])*v.conv_P,max(EoS[0])*v.conv_P),(min(EoS[0])*v.conv_P,max(EoS[0])*v.conv_P), label = "P=eps")
    plt.title("EoS(P) 2 fluid")
    if log == 1:
        plt.xscale("log")
        plt.yscale("log")
    plt.xlabel("P in MeV fm^-3")
    plt.legend()
    plt.grid()
    plt.savefig("results/EoS_P_2_fluid.pdf")
    plt.show()

########################################### single

def plot_single_1_fluid(filename):
    single = np.transpose(np.genfromtxt(filename))
    plt.plot(single[0],single[1], label="stepsize")
    plt.plot(single[0],single[2], label="M")
    plt.plot(single[0],single[3], label="mu")
    plt.plot(single[0],single[4], label="P")
    plt.plot(single[0],single[5], label="eps")
    plt.plot(single[0],single[6], label="c_s")
    plt.plot(single[0],single[7], label="y")
    plt.plot(single[0],single[8], label="Q")
    plt.plot(single[0],single[9], label="F")
    plt.plot(single[0],single[10], label="d_y_r")
    plt.title("all values")
    plt.xlabel("radius R")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    plt.grid()
    plt.savefig("results/single_2_fluid.pdf")
    plt.show()

def plot_single_2_fluid(filename):
    single = np.transpose(np.genfromtxt(filename))
    plt.plot(single[0],single[1], label="stepsize")
    plt.plot(single[0],single[2], label="M OM")
    plt.plot(single[0],single[3], label="M DM")
    plt.plot(single[0],single[4], label="mu OM")
    plt.plot(single[0],single[5], label="mu DM")
    plt.plot(single[0],single[6], label="P_OM")
    plt.plot(single[0],single[7], label="P_DM")
    plt.plot(single[0],single[8], label="eps OM")
    plt.plot(single[0],single[9], label="eps DM")
    plt.plot(single[0],single[10], label="c_s OM")
    plt.plot(single[0],single[11], label="c_s DM", ls = "--")
    plt.title("all values")
    plt.xlabel("radius R")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    plt.grid()
    plt.savefig("results/single_2_fluid.pdf")
    plt.show()

def plot_single_TIDAL_2_fluid(filename):
    single = np.transpose(np.genfromtxt(filename))
    plt.plot(single[0],single[1], label="stepsize")
    plt.plot(single[0],single[2], label="M OM")
    plt.plot(single[0],single[3], label="M DM")
    plt.plot(single[0],single[4], label="mu OM")
    plt.plot(single[0],single[5], label="mu DM")
    plt.plot(single[0],single[6], label="P_OM")
    plt.plot(single[0],single[7], label="P_DM")
    plt.plot(single[0],single[8], label="eps OM")
    plt.plot(single[0],single[9], label="eps DM")
    plt.plot(single[0],single[10], label="c_s OM")
    plt.plot(single[0],single[11], label="c_s DM", ls = "--")
    plt.plot(single[0],single[12], label="y", ls = "--")
    plt.plot(single[0],single[13], label="Q", ls = "--")
    plt.plot(single[0],single[14], label="F", ls = "--")
    plt.plot(single[0],single[15], label="d_y_r", ls = "--")
    plt.title("all values")
    plt.xlabel("radius R")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    plt.grid()
    plt.savefig("results/single_TIDAL_2_fluid.pdf")
    plt.show()



########################################### plot M_R

def plot_M_R_1_fluid(filename, log = 0):
    M_R = np.transpose(np.genfromtxt(filename))
    plt.plot(M_R[1]*v.con_radius,M_R[2]*v.con_mass, label="M")
    # plt.plot(M_R[1],M_R[2], label="M")
    plt.title("Mass Radius Relation")
    if log == 1:
        plt.xscale("log")
        plt.yscale("log")
    plt.xlabel("R in km")
    plt.ylabel("M in $M_0$")
    # plt.xlim([0,10])
    # plt.ylim([0,1.5])
    plt.legend()
    plt.savefig("results/M_R_1_fluid.pdf")
    plt.show()


def plot_M_R_2_fluid(filename, log = 0):
    M_R = np.transpose(np.genfromtxt(filename))
    R_tot = []
    for i in range(len(M_R[1])):
        R_tot.append(max(M_R[1][i]*v.con_radius,M_R[2][i]*v.con_radius))
    plt.plot(R_tot,M_R[3]*v.con_mass, label="M OM")
    plt.plot(R_tot,M_R[4]*v.con_mass, label="M DM")
    plt.title("Mass Radius Relation")
    if log == 1:
        plt.xscale("log")
        plt.yscale("log")
    plt.xlabel("R OM in km")
    plt.ylabel("M in solar masses")
    plt.grid()
    plt.legend()
    plt.title(filename)
    # print(filename)
    plt.savefig("results/M_R_2_fluid_%s.pdf"%filename[8:len(filename)-4])
    # plt.show()

def plot_k_2_C_2_fluid(filename, dm_om = 1):
    M_R = np.transpose(np.genfromtxt(filename))
    plt.plot(M_R[6],M_R[7], label="r=%1.1e"%dm_om)
    # plt.title("k_2 vs C\n"+EoS.K_gam_units()+", r=%1.1e"%dm_om)
    plt.xlabel("C")
    plt.ylabel("$k_2$")
    plt.legend()
    # plt.savefig("results/k_2_C_2_fluid.pdf")
    # plt.show()



########################################################################### TESTING

def plot_for_different_limits(log = 0):
    M_R_28_32 = np.transpose(np.genfromtxt("results/M_R_P_1_fluid_28_32.out"))
    M_R_28_33 = np.transpose(np.genfromtxt("results/M_R_P_1_fluid_28_33.out"))
    M_R_28_34 = np.transpose(np.genfromtxt("results/M_R_P_1_fluid_28_34.out"))
    M_R_29_34 = np.transpose(np.genfromtxt("results/M_R_P_1_fluid_29_34.out"))
    M_R_30_32 = np.transpose(np.genfromtxt("results/M_R_P_1_fluid_30_32.out"))
    M_R_30_34 = np.transpose(np.genfromtxt("results/M_R_P_1_fluid_30_34.out"))
    plt.plot(M_R_28_32[1]*v.con_radius,M_R_28_32[2]*v.con_mass, label="M low=1.04  high=1.251")
    plt.plot(M_R_28_33[1]*v.con_radius,M_R_28_33[2]*v.con_mass, label="M low=1.04  high=1.325")
    plt.plot(M_R_28_34[1]*v.con_radius,M_R_28_34[2]*v.con_mass, label="M low=1.04  high=1.398")
    plt.plot(M_R_29_34[1]*v.con_radius,M_R_29_34[2]*v.con_mass, label="M low=1.067  high=1.398")
    plt.plot(M_R_30_32[1]*v.con_radius,M_R_30_32[2]*v.con_mass, label="M low=1.104  high=1.251")
    plt.plot(M_R_30_34[1]*v.con_radius,M_R_30_34[2]*v.con_mass, label="M low=1.104  high=1.398")
    plt.title("Mass Radius Relation for different limits")
    if log == 1:
        plt.xscale("log")
        plt.yscale("log")
    plt.xlabel("R in km")
    plt.ylabel("M in $M_0$")
    plt.xlim([0,15])
    plt.ylim([0,2])
    plt.legend()
    plt.savefig("results/M_R_differnt_limits.pdf")
    plt.show()

def plot_stability(filename_template, log = 1):
    plt.figure(figsize=(5,15))
    M_R = np.transpose(np.genfromtxt("results/M_R_TIDAL_"+filename_template+".out"))
    stability = np.asarray(np.genfromtxt("results/stability_"+filename_template))
    # print(len(M_R))
    # print(len(stability), len(stability[0]))
    # M_R_stable = []
    # exit()
    # for i in range(len(stability)):
    #     M_R_stable.append([])
    #     for k in range(len(M_R[j])):
    #         M_R_stable[i].append([])
    #         if (M_R[0][k] > stability[i][0]) and (M_R[0][k] < stability[i][1]):
    #             for j in range(M_R_stable):
    #                 M_R_stable[j].append([])
    #                 M_R_stable[i][j].append(M_R[j][k])
    plt.subplot(311)
    plt.title(filename_template)
    plt.grid()
    if log == 1:
        plt.xscale("log")
        plt.yscale("log")
    plt.plot(M_R[1],M_R[3], label="M OM")
    plt.plot(M_R[1],M_R[4], label="M DM")
    plt.legend()
    plt.xlabel("R in km")
    plt.ylabel("M in solar masses")
    plt.subplot(312)
    plt.grid()
    if log == 1:
        plt.xscale("log")
        plt.yscale("log")
    # plt.plot(M_R[0]*v.conv_P,M_R[3]*v.con_mass, label="M OM")
    # plt.plot(M_R[0]*v.conv_P,M_R[4]*v.con_mass, label="M DM")
    plt.plot(M_R[0],M_R[3], label="M OM")
    plt.plot(M_R[0],M_R[4], label="M DM")
    for i in range(len(stability)):
        for j in range(len(stability[i])):
            plt.axvline(stability[i][j])
    plt.legend()
    plt.xlabel("P0")
    # plt.xlabel("P0 in MeV/fm^3")
    plt.ylabel("M in solar masses")
    plt.subplot(313)
    plt.grid()
    if log == 1:
        plt.xscale("log")
        plt.yscale("log")
    plt.plot(M_R[0],M_R[1], label="R OM")
    plt.plot(M_R[0],M_R[2], label="R DM")
    for i in range(len(stability)):
        for j in range(len(stability[i])):
            plt.axvline(stability[i][j])
    plt.xlabel("P0")
    # plt.plot(M_R[0]*v.conv_P,M_R[1]*v.con_mass, label="R OM")
    # plt.plot(M_R[0]*v.conv_P,M_R[2]*v.con_mass, label="R DM")
    # plt.xlabel("P0 in MeV/fm^3")
    plt.ylabel("R in km")

    # plt.xlabel("R in km")
    # plt.ylabel("M in solar masses")
    plt.legend()
    plt.savefig("results/Stability.pdf")
    plt.show()




























# #####################################################################

# def plot_M_R(filename):
#     output = np.transpose(np.genfromtxt(filename, "f"))

#     plt.plot(output[0]*v.con_radius, output[1]*v.con_mass, label = "M OM")
#     plt.plot(output[0]*v.con_radius, output[2]*v.con_mass, label = "M DM")
#     plt.title("MRR")
#     plt.legend()
#     plt.xlabel("R in km")
#     # plt.xscale("log")
#     # plt.yscale("log")
#     plt.ylabel("M in solar masses")
#     plt.legend()
#     plt.savefig("results/M_R_R.pdf")
#     plt.show()

def plot_k_2_C(filename):
    output = np.transpose(np.genfromtxt(filename, "f"))

    plt.plot(output[6], output[7], label = "k_2")
    plt.xlim((0,0.35))
    plt.ylim((0,0.35))
    # plt.title("k_2 vs C")
    plt.title(filename)
    plt.xlabel("C")
    plt.ylabel("k_2")
    plt.legend()
    plt.savefig("results/k_2_C.pdf")
    plt.show()

# def plot_EoS_mu():
#     EoS = np.transpose(np.genfromtxt("results/EoS.out"))
#     plt.plot(EoS[0],EoS[1], label="P OM")
#     plt.plot(EoS[0],EoS[2], label="eps OM")
#     plt.plot(EoS[0],EoS[3], label="c_s OM")
#     plt.plot(EoS[0],EoS[4], label="P DM")
#     plt.plot(EoS[0],EoS[5], label="eps DM")
#     plt.plot(EoS[0],EoS[6], label="c_s DM")
#     plt.title("EoS")
#     plt.xscale("log")
#     plt.yscale("log")
#     plt.xlabel("mu")
#     plt.legend()
#     plt.savefig("results/EoS.pdf")
#     plt.show()

# def plot_single():
#     single = np.transpose(np.genfromtxt("results/single.out"))
#     plt.plot(single[0],single[1], label="M OM")
#     plt.plot(single[0],single[2], label="M DM")
#     plt.plot(single[0],single[3], label="P OM")
#     plt.plot(single[0],single[4], label="P DM")
#     plt.plot(single[0],single[5], label="y")
#     plt.plot(single[0],single[6], label="stepsize")
#     plt.plot(single[0],single[7], label="eps OM")
#     plt.plot(single[0],single[8], label="eps DM")
#     plt.plot(single[0],single[9], label="c_s OM")
#     plt.plot(single[0],single[10], label="c_s DM")
#     plt.plot(single[0],single[11], label="Q",ls="--")
#     plt.plot(single[0],single[12], label="F",ls="--")
#     plt.plot(single[0],single[13], label="d_y_r",ls="--")
#     plt.title("all values")
#     plt.xlabel("radius R")
#     plt.xscale("log")
#     plt.yscale("log")
#     plt.legend()
#     plt.grid()
#     plt.savefig("results/single.pdf")
#     plt.show()

# def plot_single_short():
#     single = np.transpose(np.genfromtxt("results/single.out"))
#     plt.plot(single[0]*v.con_radius,single[1]*v.con_mass, label="M OM")
#     plt.plot(single[0]*v.con_radius,single[2]*v.con_mass, label="M DM")
#     plt.plot(single[0]*v.con_radius,single[3]*v.conv_P, label="P OM")
#     plt.plot(single[0]*v.con_radius,single[4]*v.conv_P, label="P DM")
#     plt.plot(single[0]*v.con_radius,single[5], label="y")
#     plt.plot(single[0]*v.con_radius,single[6]*v.con_radius, label="stepsize")
#     # plt.plot(single[0]*v.con_radius,single[7], label="eps OM")
#     # plt.plot(single[0]*v.con_radius,single[8], label="eps DM")
#     # plt.plot(single[0]*v.con_radius,single[9], label="c_s OM")
#     # plt.plot(single[0]*v.con_radius,single[10], label="c_s DM")
#     # plt.plot(single[0]*v.con_radius,single[11], label="Q")
#     # plt.plot(single[0]*v.con_radius,single[12], label="F")
#     # plt.plot(single[0]*v.con_radius,single[13], label="d_y_r")
#     plt.title("all values")
#     plt.xlabel("radius R")
#     plt.xscale("log")
#     plt.yscale("log")
#     plt.legend()
#     plt.grid()
#     plt.savefig("results/single.pdf")
#     plt.show()

# def plot_only_M_R_1_fluid(filename):
#     output = np.transpose(np.genfromtxt(filename, "f"))

#     plt.plot(output[1]*v.con_radius, output[2]*v.con_mass, label = "M")
#     plt.title("MRR")
#     plt.legend()
#     plt.xlabel("R in km")
#     # plt.xscale("log")
#     # plt.yscale("log")
#     plt.ylabel("M in solar masses")
#     plt.legend()
#     plt.savefig("results/M_R_R.pdf")
#     plt.show()

# def plot_EoS_mu_1_fluid():
#     EoS = np.transpose(np.genfromtxt("results/EoS.out"))
#     plt.plot(EoS[0],EoS[1], label="P")
#     plt.plot(EoS[0],EoS[2], label="eps")
#     plt.plot(EoS[0],EoS[3], label="c_s")
#     plt.title("EoS")
#     plt.xscale("log")
#     plt.yscale("log")
#     plt.xlabel("mu")
#     plt.legend()
#     plt.savefig("results/EoS.pdf")
#     plt.show()


# # exit()

# # f = open("EoS_Axel.out","w")
# # # f.write("%e\t%e\t%e\n"%(0, 0, 0))
# # # for init_mu in np.linspace(neutron_mass_fm, 10, 10000):
# # for init_mu in np.linspace(neutron_mass_fm, 7.7, 10000):
# #     f.write("%e\t%e\t%e\t%e\t%e\n"%(init_mu, n_mu(init_mu), P_mu(init_mu), eps_mu(init_mu), c_s_mu(init_mu)))
# # f.close()

# # EoS_plot = np.transpose(np.genfromtxt("EoS_Axel.out", "f"))

# # plt.plot(EoS_plot[0], EoS_plot[1], label = "n")
# # plt.plot(EoS_plot[0], EoS_plot[2], label = "P")
# # plt.plot(EoS_plot[0], EoS_plot[3], label = "eps")
# # plt.plot(EoS_plot[0], EoS_plot[4], label = "c_s")
# # plt.legend()
# # plt.show()

# # exit()

# # f = open("EoS_Axel_P_eps.out","w")
# # # f.write("%e\t%e\n"%(0, 0))
# # for init_mu in np.linspace(neutron_mass_fm, 10, 999):
# #     f.write("%e\t%e\n"%(P_mu(init_mu)*conv_P, eps_mu(init_mu)*conv_P))
# # f.write("%e\t%e"%(1e30, 1e30))
# # f.close()

# # EoS_plot_PE = np.transpose(np.genfromtxt("EoS_Axel_P_eps.out", "f"))
# # plt.xlim([EoS_plot_PE[0][1], EoS_plot_PE[0][len(EoS_plot_PE[0])-1]])
# # plt.ylim([EoS_plot_PE[1][1], EoS_plot_PE[1][len(EoS_plot_PE[1])-1]])
# # plt.xscale("log")
# # plt.yscale("log")
# # plt.plot(EoS_plot_PE[0], EoS_plot_PE[1], label = "eps")
# # plt.legend()
# # plt.show()









# # exit()









































