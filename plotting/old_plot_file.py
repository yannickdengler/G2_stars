#!/usr/bin/env python
# coding: utf-8

import os, sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
# from scipy.special import lambertw
# import subprocess
#plt.rcParams["font.size"] = 20
# import math
# plt.style.use('presentation.mplstyle')
from matplotlib import colors
from ellipses import draw_errors

import styles

import math

# def round_to_decades(lo, hi):
#     lo_exp = math.floor(math.log10(lo))
#     hi_exp = math.ceil(math.log10(hi))

#     lo_rounded = 10**lo_exp
#     hi_rounded = 10**hi_exp

#     return lo_rounded, hi_rounded

# font = {'size'   : 16}
# matplotlib.rc('font', **font)

mpl = 1.221E22
mplkg = 1.783E-30
msol = 1.988E30
lpl = 1.616E-35

hbarc=197.3             # in MeV fm

Planck_mass_solar = 1.095*1e-38             # in solar masses
Planck_length_km = 1.6163*1e-38             # in km
MeV_solar = 8.971e-61                       # in solar masses

N_OM = 60
N_DM = 60

scattersize = 10
scattersize_stab = 10

# stabilitycolor = 'cornflowerblue'
# consistencycolor = 'blue'
# consistency1color = 'crimson'

unstabilitycolor = "whitesmoke"#'gainsboro'
stabilitycolor = "lightgray"#'gainsboro'
consistencycolor = 'magenta'
consistency1color = 'magenta'

colormp = "viridis"

lower_frac_lim = 1
max_dm_frac = 10

def limit_func(what, EoS, m_DM):
    xmin,xmax,ymin,ymax = [0,0,0,0]

    if what == "M_tot_R_OM":
        xmin = 8
        xmax = 15
        ymin = 0.5
        ymax = 2.5
    elif what == "k2_C":
        xmin = 0
        xmax = 0.3
        ymin = 0
        ymax = 0.2
    elif what == "Lambda_M_tot":
        return None
    elif what == "M_tot_R_DM":
        return None
    elif what == "R_DM_R_OM":
        xmin = 0
        xmax = 15
        if m_DM == 371.498572:
            ymin = 0
            ymax = 90
        elif m_DM == 1000:
            ymin = 0
            ymax = 10
        elif m_DM == 4000:
            ymin = 0
            ymax = 0.65
    elif what == "M_DM_M_OM":
        xmin = 0
        xmax = 2.5
        if m_DM == 371.498572:
            ymin = 0
            ymax = 9
        elif m_DM == 1000:
            ymin = 0
            ymax = 1.2
        elif m_DM == 4000:
            ymin = 0
            ymax = 0.08
    elif what == "frac_OM":
        return None
    elif what == "frac_DM":
        return None
    elif what == "stability":
        return None
    return [xmin,xmax,ymin,ymax]

def EoS_name(EoS):
    if EoS == "light_EoS_I":
        return "EoS I light"
    elif EoS == "light_EoS_II":
        return "EoS II light"
    elif EoS == "light_EoS_III":
        return "EoS III light"
    elif EoS == "heavy_EoS_I":
        return "EoS I heavy"
    elif EoS == "heavy_EoS_II":
        return "EoS II heavy"
    elif EoS == "heavy_EoS_III":
        return "EoS III heavy"
    
def plot_M_R_constraints(ax):
    NICER = np.transpose(np.genfromtxt("constraints/NICER.csv",delimiter=","))
    GW170817 = np.transpose(np.genfromtxt("constraints/GW170817.csv",delimiter=","))

    ax.fill(NICER[0],NICER[1], alpha = 0.5)
    ax.fill(GW170817[0],GW170817[1], alpha = 0.5)
    PSR_low_up = [2.14-0.18,2.14+0.2]
    ax.axhspan(PSR_low_up[0],PSR_low_up[1], color = "grey", alpha = 0.3)

def p_c_lim_OM(EoS):
    # print(EoS[6:])
    if EoS[6:] == "EoS_I":
        return 538.7315
    elif EoS[6:] == "EoS_II":
        return 203.8201
    elif EoS[6:] == "EoS_III":
        return 112.3842
    else:
        raise ValueError
    
def p_c_lim_DM(EoS, mDM):
    # print(EoS[:7])
    if EoS[:5] == "light":
        if mDM == 500:
            return 25.49257
        elif mDM == 1000:
            return 407.88110
        elif mDM == 2000:
            return 6526.09800
        elif mDM == 4000:
            return 104417.60000
    elif EoS[:5] == "heavy":
        if mDM == 500:
            return 28.54073
        elif mDM == 1000:
            return 456.65170
        elif mDM == 2000:
            return 7306.42700
        elif mDM == 4000:
            return 116902.80000
    else:
        raise ValueError


def plot(EoS, mass, what="M_tot_R_OM", pref = "", show = False):
    file = "data/%s_data_m_DM_%e.dat"%(EoS,mass)
    OM_CP = np.genfromtxt(file, usecols = 0) # ordinary matter central pressure
    DM_CP = np.genfromtxt(file, usecols = 1) # DM matter central pressure
    OM_rad = np.genfromtxt(file, usecols = 2) # OM matter radius
    DM_rad = np.genfromtxt(file, usecols = 3) # DM matter radius
    max_rad = np.genfromtxt(file, usecols = 4) # max betweem OM and DM radius
    OM_mass = np.genfromtxt(file, usecols = 5) # OM mass
    DM_mass = np.genfromtxt(file, usecols = 6) # OM mass
    tot_mass = np.genfromtxt(file, usecols = 7) # total mass
    OM_npart = np.genfromtxt(file, usecols = 8) # OM number of particles
    DM_npart = np.genfromtxt(file, usecols = 9) # DM number of particles
    love_num_aux = np.genfromtxt(file, usecols = 10) # Auxillary parameter for love number 'y'
    compactness = np.genfromtxt(file, usecols = 11) # Compactness (M_OM+M_DM)/max(R_OM, R_DM)
    love_num_k2 = np.genfromtxt(file, usecols = 12) # Love number k_2
    OM_compactness = np.genfromtxt(file, usecols = 13) # compactness @R_OM
    DM_compactness = np.genfromtxt(file, usecols = 14) # compactness @R_DM
    OM_Ceps = np.genfromtxt(file, usecols = 15) # central energy density OM
    DM_Ceps = np.genfromtxt(file, usecols = 16) # central energy density DM
    stability =  np.genfromtxt(file, usecols = 17) # compactness @R_DM
    mass_DM =  np.genfromtxt(file, usecols = 18)[0] # compactness @R_DM
    Lambda =  np.genfromtxt(file, usecols = 19) # dimensionless Tidal Deform  
    file_OM = "data/%s_pure_OM_data_m_DM_%e.dat"%(EoS,mass)
    OM_CP_pure_OM = np.genfromtxt(file_OM, usecols = 0) # ordinary matter central pressure
    DM_CP_pure_OM = np.genfromtxt(file_OM, usecols = 1) # DM matter central pressure
    OM_rad_pure_OM = np.genfromtxt(file_OM, usecols = 2) # OM matter radius
    DM_rad_pure_OM = np.genfromtxt(file_OM, usecols = 3) # DM matter radius
    max_rad_pure_OM = np.genfromtxt(file_OM, usecols = 4) # max betweem OM and DM radius
    OM_mass_pure_OM = np.genfromtxt(file_OM, usecols = 5) # OM mass
    DM_mass_pure_OM = np.genfromtxt(file_OM, usecols = 6) # OM mass
    tot_mass_pure_OM = np.genfromtxt(file_OM, usecols = 7) # total mass
    OM_npart_pure_OM = np.genfromtxt(file_OM, usecols = 8) # OM number of particles
    DM_npart_pure_OM = np.genfromtxt(file_OM, usecols = 9) # DM number of particles
    love_num_aux_pure_OM = np.genfromtxt(file_OM, usecols = 10) # Auxillary parameter for love number 'y'
    compactness_pure_OM = np.genfromtxt(file_OM, usecols = 11) # Compactness (M_OM+M_DM)/max(R_OM, R_DM)
    love_num_k2_pure_OM = np.genfromtxt(file_OM, usecols = 12) # Love number k_2
    OM_compactness_pure_OM = np.genfromtxt(file_OM, usecols = 13) # compactness @R_OM
    DM_compactness_pure_OM = np.genfromtxt(file_OM, usecols = 14) # compactness @R_DM
    OM_Ceps_pure_OM = np.genfromtxt(file_OM, usecols = 15) # central energy density OM
    DM_Ceps_pure_OM = np.genfromtxt(file_OM, usecols = 16) # central energy density DM
    stability_pure_OM =  np.genfromtxt(file_OM, usecols = 17) # compactness @R_DM
    mass_DM_pure_OM =  np.genfromtxt(file_OM, usecols = 18)[0] # compactness @R_DM
    Lambda_pure_OM =  np.genfromtxt(file_OM, usecols = 19) # dimensionless Tidal Deform  
    file_DM = "data/%s_pure_DM_data_m_DM_%e.dat"%(EoS,mass)
    OM_CP_pure_DM = np.genfromtxt(file_DM, usecols = 0) # ordinary matter central pressure
    DM_CP_pure_DM = np.genfromtxt(file_DM, usecols = 1) # DM matter central pressure
    OM_rad_pure_DM = np.genfromtxt(file_DM, usecols = 2) # OM matter radius
    DM_rad_pure_DM = np.genfromtxt(file_DM, usecols = 3) # DM matter radius
    max_rad_pure_DM = np.genfromtxt(file_DM, usecols = 4) # max betweem OM and DM radius
    OM_mass_pure_DM = np.genfromtxt(file_DM, usecols = 5) # OM mass
    DM_mass_pure_DM = np.genfromtxt(file_DM, usecols = 6) # OM mass
    tot_mass_pure_DM = np.genfromtxt(file_DM, usecols = 7) # total mass
    OM_npart_pure_DM = np.genfromtxt(file_DM, usecols = 8) # OM number of particles
    DM_npart_pure_DM = np.genfromtxt(file_DM, usecols = 9) # DM number of particles
    love_num_aux_pure_DM = np.genfromtxt(file_DM, usecols = 10) # Auxillary parameter for love number 'y'
    compactness_pure_DM = np.genfromtxt(file_DM, usecols = 11) # Compactness (M_DM+M_DM)/max(R_DM, R_DM)
    love_num_k2_pure_DM = np.genfromtxt(file_DM, usecols = 12) # Love number k_2
    OM_compactness_pure_DM = np.genfromtxt(file_DM, usecols = 13) # compactness @R_DM
    DM_compactness_pure_DM = np.genfromtxt(file_DM, usecols = 14) # compactness @R_DM
    OM_Ceps_pure_DM = np.genfromtxt(file_DM, usecols = 15) # central energy density OM
    DM_Ceps_pure_DM = np.genfromtxt(file_DM, usecols = 16) # central energy density DM
    stability_pure_DM =  np.genfromtxt(file_DM, usecols = 17) # compactness @R_DM
    mass_DM_pure_DM =  np.genfromtxt(file_DM, usecols = 18)[0] # compactness @R_DM
    Lambda_pure_DM =  np.genfromtxt(file_DM, usecols = 19) # dimensionless Tidal Deform  
    dm_frac = np.array([100*DM_mass[x]/tot_mass[x] for x in range(len(tot_mass))])

    # print(stability_pure_OM)
    # print(stability_pure_DM)

    

    # condition_consistency = (stability == 1) & (tot_mass > 0.7) & (tot_mass < 2.2) & (OM_rad > 10) & (OM_rad < 16) & (Lambda < 2000)
    stable = (stability == 1)
    # stable_constraint = (stability == 1) & (tot_mass > 0.7) & (tot_mass < 2.2) & (OM_rad > 10) & (OM_rad < 16) & (Lambda < 2000) & (dm_frac > lower_frac_lim)
    # stable_constraint_1 = (stability == 1) & (tot_mass > 0.7) & (tot_mass < 2.2) & (OM_rad > 10) & (OM_rad < 16) & (Lambda < 2000) & (dm_frac <= lower_frac_lim)
    stable_constraint = (stability == 1) & (dm_frac > lower_frac_lim) & (dm_frac <= max_dm_frac)
    stable_constraint_1 = (stability == 1) & (dm_frac <= lower_frac_lim)

    # print(max(dm_frac[stable_constraint]))

    label_stable = r"$\text{stable}$"
    label_stable1 = r"$\text{stable}<1\%$"


    if what == "stability":
        fig, ax = plt.subplots(figsize = (8,7))
        ax.scatter(OM_CP[stability == 0],DM_CP[stability == 0],  s = scattersize_stab, color = "white",zorder=-2, marker="s")
        # ax.scatter(OM_CP[stability == 1],DM_CP[stability == 1],  s = scattersize_stab, label =  label_stable, color = "red")
        # ax.scatter(OM_CP[(stability == 0)],DM_CP[(stability == 0)],  s = scattersize_stab, label =   'Unstable', color = unstabilitycolor)
        ax.scatter(OM_CP[stable],DM_CP[stable],  s = scattersize_stab, label =  label_stable, color = stabilitycolor, marker="s")
        sc = ax.scatter(OM_CP[stable_constraint],DM_CP[stable_constraint],  c = dm_frac[stable_constraint],  s = scattersize_stab, cmap=plt.get_cmap(colormp), norm = colors.LogNorm(vmin=lower_frac_lim, vmax=max_dm_frac), marker="s")   #, norm = colors.Normalize(vmin=lower_frac_lim, vmax=100-lower_frac_lim))   
        ax.scatter(OM_CP[stable_constraint_1],DM_CP[stable_constraint_1],  s = scattersize_stab, color = consistency1color, label =  label_stable1, marker="s")
        # cb = plt.colorbar(sc, ticks=[5,10,30,90])
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r"$p_{0,O}\,\text{[MeV/fm}^3\text{]}$")
        ax.set_ylabel(r"$p_{0,D}\,\text{[MeV/fm}^3\text{]}$")
        
        scale = 0.7
        xmin = min(OM_CP)*scale
        xmax = max(OM_CP)/scale
        ymin = min(DM_CP)*scale
        ymax = max(DM_CP)/scale
        exp_x_min = math.ceil(math.log10(xmin))
        exp_x_max = math.floor(math.log10(xmax))
        exp_y_min = math.ceil(math.log10(ymin))
        exp_y_max = math.floor(math.log10(ymax))
        xticks = np.logspace(exp_x_min, exp_x_max, exp_x_max-exp_x_min)
        ax.set_xticks(xticks, [r"$10^{%i}$"%np.log10(x) for x in xticks])
        yticks = np.logspace(exp_y_min, exp_y_max, exp_y_max-exp_y_min)
        ax.set_yticks(yticks, [r"$10^{%i}$"%np.log10(x) for x in yticks])
        ax.set_xlim([xmin,xmax])
        ax.set_ylim([ymin,ymax])
    elif what == "M_tot_R_OM":    
        fig, ax = plt.subplots(figsize = (8,6))
        ax.scatter(OM_rad[stable],tot_mass[stable],  s = scattersize, label =  label_stable, color = stabilitycolor)
        sc = ax.scatter(OM_rad[stable_constraint],tot_mass[stable_constraint],  c = dm_frac[stable_constraint],  s = scattersize, cmap=plt.get_cmap(colormp), norm = colors.LogNorm(vmin=lower_frac_lim, vmax=max_dm_frac))   #, norm = colors.Normalize(vmin=lower_frac_lim, vmax=100-lower_frac_lim))   
        ax.scatter(OM_rad[stable_constraint_1],tot_mass[stable_constraint_1],  s = scattersize, color = consistency1color, label =  label_stable1)
        # ax.scatter(OM_rad_pure_OM[(stability_pure_OM == 1)],tot_mass_pure_OM[(stability_pure_OM == 1)],  s = scattersize, color = consistency1color)
        # ax.scatter(DM_rad_pure_DM[(stability_pure_DM == 1)],tot_mass_pure_DM[(stability_pure_DM == 1)],  s = scattersize, label =   'pure OM', color = "green")
        # plot_M_R_constraints(ax)
        draw_errors(ax)
        ax.set_xlabel(r"$R_{O}\,\text{[km]}$")
        ax.set_ylabel(r"$M_{tot}\,\text{[}M_\odot\text{]}$")
        xticks = np.linspace(8,14,4)
        ax.set_xticks(xticks, [r"$%i$"%x for x in xticks])
        yticks = np.linspace(0.5,2.5,5)
        ax.set_yticks(yticks, [r"$%1.1f$"%x for x in yticks])
    elif what == "M_tot_R_DM":    
        fig, ax = plt.subplots(figsize = (8,6))
        ax.scatter(DM_rad[stable],tot_mass[stable],  s = scattersize, label =  label_stable, color = stabilitycolor)
        sc = ax.scatter(DM_rad[stable_constraint],tot_mass[stable_constraint],  c = dm_frac[stable_constraint],  s = scattersize, cmap=plt.get_cmap(colormp), norm = colors.LogNorm(vmin=lower_frac_lim, vmax=max_dm_frac))   #, norm = colors.Normalize(vmin=lower_frac_lim, vmax=100-lower_frac_lim))   
        ax.scatter(DM_rad[stable_constraint_1],tot_mass[stable_constraint_1],  s = scattersize, color = consistency1color, label =  label_stable1)
        ax.scatter(OM_rad_pure_OM[(stability_pure_OM == 1)],tot_mass_pure_OM[(stability_pure_OM == 1)],  s = scattersize, label =   'pure OM', color = "red")
        ax.scatter(DM_rad_pure_DM[(stability_pure_DM == 1)],tot_mass_pure_DM[(stability_pure_DM == 1)],  s = scattersize, label =   'pure OM', color = "green")
        ax.set_xlabel(r"$R_{D}\,\text{[km]}$")
        ax.set_ylabel(r"$M_{tot}\,\text{[}M_\odot\text{]}$")
    elif what == "M_DM_M_OM":    
        fig, ax = plt.subplots(figsize = (8,6))
        ax.scatter(OM_mass[stable],DM_mass[stable],  s = scattersize, label =  label_stable, color = stabilitycolor)
        sc = ax.scatter(OM_mass[stable_constraint],DM_mass[stable_constraint],  c = dm_frac[stable_constraint],  s = scattersize, cmap=plt.get_cmap(colormp), norm = colors.LogNorm(vmin=lower_frac_lim, vmax=max_dm_frac))   #, norm = colors.Normalize(vmin=lower_frac_lim, vmax=100-lower_frac_lim))   
        ax.scatter(OM_mass[stable_constraint_1],DM_mass[stable_constraint_1],  s = scattersize, color = consistency1color, label =  label_stable1)
        # cb = plt.colorbar(sc, ticks=[5,10,30,90])
        # cb.set_label('Fraction of DM (%)', rotation=270, labelpad=15, size = 14)
        ax.set_xlabel(r"$R_{O}\,\text{[km]}$")
        ax.set_ylabel(r"$R_{D}\,\text{[km]}$")
    elif what == "R_DM_R_OM":    
        fig, ax = plt.subplots(figsize = (8,6))
        ax.scatter(OM_rad[stable],DM_rad[stable],  s = scattersize, label =  label_stable, color = stabilitycolor, zorder = 9)
        sc = ax.scatter(OM_rad[stable_constraint],DM_rad[stable_constraint],  c = dm_frac[stable_constraint],  s = scattersize, cmap=plt.get_cmap(colormp), norm = colors.LogNorm(vmin=lower_frac_lim, vmax=max_dm_frac), zorder = 10)   #, norm = colors.Normalize(vmin=lower_frac_lim, vmax=100-lower_frac_lim))   
        ax.scatter(OM_rad[stable_constraint_1],DM_rad[stable_constraint_1],  s = scattersize, color = consistency1color, label =  label_stable1, zorder = 11)
        plt.plot([0,100], [0,100], color = "black", ls = "dashed", label = "$R_{OM} = R_{DM}$", zorder = 5)
        ax.set_xlabel(r"$R_{O}\,\text{[km]}$")
        ax.set_ylabel(r"$R_{D}\,\text{[km]}$")
    elif what == "Lambda_M_tot":    
        fig, ax = plt.subplots(figsize = (8,6))
        # ax.set_xscale("log")
        ax.set_yscale("log")
        ax.scatter(tot_mass[stable],Lambda[stable],  s = scattersize, label =  label_stable, color = stabilitycolor, zorder = 9)
        sc = ax.scatter(tot_mass[stable_constraint],Lambda[stable_constraint],  c = dm_frac[stable_constraint],  s = scattersize, cmap=plt.get_cmap(colormp), norm = colors.LogNorm(vmin=lower_frac_lim, vmax=max_dm_frac), zorder = 10)   #, norm = colors.Normalize(vmin=lower_frac_lim, vmax=100-lower_frac_lim))   
        ax.scatter(tot_mass[stable_constraint_1],Lambda[stable_constraint_1],  s = scattersize, color = consistency1color, label =  label_stable1, zorder = 11)
        ax.set_xlabel(r"$M_{tot}\,\text{[}M_\odot\text{]}$")
        ax.set_ylabel(r"$\Lambda$")
        xticks = np.linspace(0,4,5)
        ax.set_xticks(xticks, [r"$%i$"%x for x in xticks])
        yticks = np.logspace(2,8,5)
        ax.set_yticks(yticks,  [r"$10^{%i}$"%np.log10(x) for x in yticks])
    elif what == "k2_C":    
        fig, ax = plt.subplots(figsize = (8,6))
        ax.scatter(compactness[stable],love_num_k2[stable],  s = scattersize, label =  label_stable, color = stabilitycolor, zorder = 9)
        sc = ax.scatter(compactness[stable_constraint],love_num_k2[stable_constraint],  c = dm_frac[stable_constraint],  s = scattersize, cmap=plt.get_cmap(colormp), norm = colors.LogNorm(vmin=lower_frac_lim, vmax=max_dm_frac), zorder = 10)   #, norm = colors.Normalize(vmin=lower_frac_lim, vmax=100-lower_frac_lim))   
        ax.scatter(compactness[stable_constraint_1],love_num_k2[stable_constraint_1],  s = scattersize, color = consistency1color, label =  label_stable1, zorder = 11)
        ax.set_xlabel(r"$\text{C}$")
        ax.set_ylabel(r"$k_2$")
        xticks = np.linspace(0,0.3,4)
        ax.set_xticks(xticks, [r"$%1.1f$"%x for x in xticks])
        yticks = np.linspace(00.05,0.2,4)
        ax.set_yticks(yticks, [r"$%1.2f$"%x for x in yticks])
    # elif what == "frac_OM":    
        # fig, ax = plt.subplots(figsize = (8,6))
    #     ax.scatter(OM_CP[stable],dm_frac[stable],  s = scattersize, label =  label_stable, color = stabilitycolor, zorder = 9)
    #     sc = ax.scatter(OM_CP[stable_constraint],dm_frac[stable_constraint],  c = dm_frac[stable_constraint],  s = scattersize, cmap=plt.get_cmap(colormp), norm = colors.LogNorm(vmin=lower_frac_lim, vmax=max_dm_frac), zorder = 10)   #, norm = colors.Normalize(vmin=lower_frac_lim, vmax=100-lower_frac_lim))   
    #     ax.scatter(OM_CP[stable_constraint_1],dm_frac[stable_constraint_1],  s = scattersize, color = consistency1color, label =  label_stable1, zorder = 11)
    #     ax.set_xlabel("$p_{OM}$ [km]")
    #     ax.set_ylabel("frac")
    # elif what == "frac_DM":    
        # fig, ax = plt.subplots(figsize = (8,6))
    #     ax.scatter(OM_CP[stable],dm_frac[stable],  s = scattersize, label =  label_stable, color = stabilitycolor, zorder = 9)
    #     sc = ax.scatter(OM_CP[stable_constraint],dm_frac[stable_constraint],  c = dm_frac[stable_constraint],  s = scattersize, cmap=plt.get_cmap(colormp), norm = colors.LogNorm(vmin=lower_frac_lim, vmax=max_dm_frac), zorder = 10)   #, norm = colors.Normalize(vmin=lower_frac_lim, vmax=100-lower_frac_lim))   
    #     ax.scatter(OM_CP[stable_constraint_1],dm_frac[stable_constraint_1],  s = scattersize, color = consistency1color, label =  label_stable1, zorder = 11)
    #     ax.set_xlabel("$p_{OM}$ [km]")
    #     ax.set_ylabel("frac")
    # plt.tight_layout()
    limits = limit_func(what, EoS, mass)
    # tick_arr = [lower_frac_lim,3,8,15,25,max_dm_frac]
    tick_arr = [lower_frac_lim,2,3,4,6,8,max_dm_frac]
    if limits != None:
        xmin,xmax,ymin,ymax = limits
        ax.set_xlim([xmin,xmax])
        ax.set_ylim([ymin,ymax])    
    cb = plt.colorbar(sc)#, ticks=tick_arr)
    cb.Ticks = []
    cb.set_ticks(tick_arr)
    cb.set_label(r'$\text{Fraction of DM (%)}$', rotation=270, labelpad=15)
    cb.ax.set_yticklabels([r"$\text{%i}$"%x for x in tick_arr])
    # cb.ax.set_yticklabels(["%i"%lower_frac_lim,"5","10","20","30","%i"%max_dm_frac])

    # ax.grid()
    
    plt.text(0.75, 1.05, r'$m_{\text{C}} = %0.2f\,\text{[GeV]}$' %(mass_DM/1000), ha='left', va='top', transform=ax.transAxes, fontsize=16)
    plt.text(0.01, 1.05, r"$\text{%s}$"%EoS_name(EoS), ha='left', va='top', transform=ax.transAxes, fontsize=16)

    if what=="M_tot_R_OM" or what == "stability":
        ax.legend(loc="upper left")
    else:
        ax.legend()
    fig.savefig('plots/%s_%s_mDM_%e.pdf' %(what, EoS, mass_DM), dpi  = 150, bbox_inches = "tight")
    if show:
        plt.show()
    plt.close(fig)

# mass_list=[125, 190,250, 304.753414,  371.498572, 452.861832, 552.044757, 672.950096, 820.335356, 1000.000000, 1219.013654, 1485.994289, 1811.447329, 2208.179027, 2691.800385, 3281.341424, 4000.000000]

EoS_list = ["light_EoS_I","light_EoS_II","light_EoS_III","heavy_EoS_I","heavy_EoS_II","heavy_EoS_III"]
# EoS_list = ["heavy_EoS_I","heavy_EoS_II","heavy_EoS_III"]
# EoS_list = ["light_EoS_II",]
mass_list=[500, 1000, 2000, 4000]
# mass_list=[500,]
for EoS in EoS_list:
    for mass in mass_list:
        plot(EoS, mass, what = "M_tot_R_OM")
for EoS in EoS_list:
    for mass in mass_list:
        plot(EoS, mass, what = "stability",show=False)
for EoS in EoS_list:
    for mass in mass_list:
        plot(EoS, mass, what = "Lambda_M_tot")
for EoS in EoS_list:
    for mass in mass_list:
        plot(EoS, mass, what = "k2_C")


# for EoS in EoS_list:
#     for mass in mass_list:
#         plot(EoS, mass, what = "M_DM_M_OM")
# for EoS in EoS_list:
#     for mass in mass_list:
#         plot(EoS, mass, what = "R_DM_R_OM")



# for EoS in EoS_list:
#     for mass in mass_list:
#         plot(EoS, mass, what = "frac_OM")
        # plot(EoS, mass, what = "frac_DM")

# for EoS in EoS_list:
#     for mass in mass_list:
#         plot(EoS, mass, what = "M_tot_R_OM")
#         plot(EoS, mass, what = "M_tot_R_DM")


























# filelist = []
# for i in range(len(EoS_list)):
#     for mass in mass_list:
#         filelist.append("data/%s_data_m_DM_%e.dat"%(EoS,mass))

# filelist.append("data/%s_data_m_DM_%e.dat"%(EoS_list[3],mass_list[0]))
# plot(filelist[0], EoS_list[3], what = "P_DM_P_OM", show = True)

    # elif what == "R_OM_M_DM":
    #     ax.scatter(OM_rad[stable],DM_mass[stable],  s = scattersize, label =  label_stable, color = stabilitycolor)
    #     sc = ax.scatter(OM_rad[stable_constraint],DM_mass[stable_constraint],  c = dm_frac[stable_constraint],  s = scattersize, cmap=plt.get_cmap(colormp), norm = colors.LogNorm(vmin=lower_frac_lim, vmax=100-lower_frac_lim))   #, norm = colors.Normalize(vmin=lower_frac_lim, vmax=100-lower_frac_lim))   
    #     ax.scatter(OM_rad[stable_constraint_1],DM_mass[stable_constraint_1],  s = scattersize, color = consistency1color)
    #     cb = plt.colorbar(sc, ticks=[5,10,30,90])
    #     cb.ax.set_yticklabels(["5","10","30","90"])
    #     cb.set_label('Fraction of DM (%)', rotation=270, labelpad=15, size = 14)
    #     ax.set_xlabel("$R_{OM}$")
    #     ax.set_ylabel("$M_{DM}$")
    #     plt.tight_layout()
    # elif what == "P_DM_M_DM":
    #     ax.scatter(DM_CP[stable],DM_mass[stable],  s = scattersize, label =  label_stable, color = stabilitycolor)
    #     sc = ax.scatter(DM_CP[stable_constraint],DM_mass[stable_constraint],  c = dm_frac[stable_constraint],  s = scattersize, cmap=plt.get_cmap(colormp), norm = colors.LogNorm(vmin=lower_frac_lim, vmax=100-lower_frac_lim))   #, norm = colors.Normalize(vmin=lower_frac_lim, vmax=100-lower_frac_lim))   
    #     ax.scatter(DM_CP[stable_constraint_1],DM_mass[stable_constraint_1],  s = scattersize, color = consistency1color)
    #     cb = plt.colorbar(sc, ticks=[5,10,30,90])
    #     cb.ax.set_yticklabels(["5","10","30","90"])
    #     cb.set_label('Fraction of DM (%)', rotation=270, labelpad=15, size = 14)
    #     ax.set_xlabel("$P_{DM}$")
    #     ax.set_ylabel("$M_{DM}$")
    #     # ax.set_xlim([lim[0],lim[1]])
    #     # ax.set_ylim([lim[2],lim[3]])
    #     plt.tight_layout()
    # elif what == "P_DM_R_DM":
    #     ax.scatter(DM_CP[stable],DM_rad[stable],  s = scattersize, label =  label_stable, color = stabilitycolor)
    #     sc = ax.scatter(DM_CP[stable_constraint],DM_rad[stable_constraint],  c = dm_frac[stable_constraint],  s = scattersize, cmap=plt.get_cmap(colormp), norm = colors.LogNorm(vmin=lower_frac_lim, vmax=100-lower_frac_lim))   #, norm = colors.Normalize(vmin=lower_frac_lim, vmax=100-lower_frac_lim))   
    #     ax.scatter(DM_CP[stable_constraint_1],DM_rad[stable_constraint_1],  s = scattersize, color = consistency1color)
    #     cb = plt.colorbar(sc, ticks=[5,10,30,90])
    #     cb.ax.set_yticklabels(["5","10","30","90"])
    #     cb.set_label('Fraction of DM (%)', rotation=270, labelpad=15, size = 14)
    #     ax.set_xlabel("$P_{DM}$")
    #     ax.set_ylabel("$R_{DM}$")
    #     # ax.set_xlim([lim[0],lim[1]])
    #     # ax.set_ylim([lim[2],lim[3]])
    #     plt.tight_layout()
    # elif what == "P_DM_P_OM":
    #     ax.scatter(DM_CP[stable],OM_CP[stable],  s = scattersize, label =  label_stable, color = stabilitycolor)
    #     sc = ax.scatter(DM_CP[stable_constraint],OM_CP[stable_constraint],  c = dm_frac[stable_constraint],  s = scattersize, cmap=plt.get_cmap(colormp), norm = colors.LogNorm(vmin=lower_frac_lim, vmax=100-lower_frac_lim))   #, norm = colors.Normalize(vmin=lower_frac_lim, vmax=100-lower_frac_lim))   
    #     ax.scatter(DM_CP[stable_constraint_1],OM_CP[stable_constraint_1],  s = scattersize, color = consistency1color)
    #     cb = plt.colorbar(sc, ticks=[5,10,30,90])
    #     cb.ax.set_yticklabels(["5","10","30","90"])
    #     cb.set_label('Fraction of DM (%)', rotation=270, labelpad=15, size = 14)
    #     ax.set_xlabel("$P_{DM}$")
    #     ax.set_ylabel("$P_{OM}$")
    #     # ax.set_xlim([lim[0],lim[1]])
    #     # ax.set_ylim([lim[2],lim[3]])
    #     plt.tight_layout()
