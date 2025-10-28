#!/usr/bin/env python
# coding: utf-8

import os, sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
from ellipses import draw_errors
from ellipses import draw_error_text

import styles

import math

scattersize = 10
scattersize_stab = 10

unstabilitycolor = "whitesmoke"#'gainsboro'
stabilitycolor = "lightgray"#'gainsboro'
# consistencycolor = 'magenta'
consistency1color = 'sandybrown'

colormp = "viridis"

lower_frac_lim = 1
max_dm_frac = 10

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
    
def get_data(file):
    OM_CP = np.genfromtxt(file, usecols = 0) # OM matter radius
    DM_CP = np.genfromtxt(file, usecols = 1) # OM matter radius
    tot_mass = np.genfromtxt(file, usecols = 7) # total mass
    stability =  np.genfromtxt(file, usecols = 17) # compactness @R_D
    DM_mass = np.genfromtxt(file, usecols = 6) # OM mass
    dm_frac = np.array([100*DM_mass[x]/tot_mass[x] for x in range(len(tot_mass))])
    return OM_CP, DM_CP, stability, dm_frac

label_stable = r"$\text{stable}$"
label_stable1 = r"$\text{stable}<1\%$"

scattersize_stab = 11.5

def plot_stability(data, ax):
    OM_CP, DM_CP, stability, dm_frac = data
    # unstable = (stability == 0)
    stable = (stability == 1)
    stable_constraint = (stability == 1) & (dm_frac > lower_frac_lim) & (dm_frac <= max_dm_frac)
    stable_constraint_1 = (stability == 1) & (dm_frac <= lower_frac_lim)
    # ax.scatter(OM_CP[unstable],DM_CP[unstable],  s = scattersize_stab, label =  label_stable, color = "white", marker="s", zorder = -4)
    ax.scatter(OM_CP[stable],DM_CP[stable],  s = scattersize_stab, label =  label_stable, color = stabilitycolor, marker="s", zorder = 2)
    sc = ax.scatter(OM_CP[stable_constraint],DM_CP[stable_constraint],
                       c = dm_frac[stable_constraint],  s = scattersize_stab, 
                       cmap=plt.get_cmap(colormp), norm = colors.LogNorm(vmin=lower_frac_lim, vmax=max_dm_frac), marker="s", zorder = 3)   #, norm = colors.Normalize(vmin=lower_frac_lim, vmax=100-lower_frac_lim))   
    ax.scatter(OM_CP[stable_constraint_1],DM_CP[stable_constraint_1],  
                s = scattersize_stab, color = consistency1color, label =  label_stable1, marker="s", zorder = 4)
    return sc

def set_x_y_lim_ticks(data, ax, ax_num, scale = 0.7):
    OM_CP = data[0]
    DM_CP = data[1]
    
    xmin = min(OM_CP)*scale
    xmax = max(OM_CP)/scale
    ymin = min(DM_CP)*scale
    ymax = max(DM_CP)/scale
    exp_x_min = math.ceil(math.log10(xmin))
    exp_x_max = math.floor(math.log10(xmax))
    exp_y_min = math.ceil(math.log10(ymin))
    exp_y_max = math.floor(math.log10(ymax))
    xticks = np.logspace(exp_x_min, exp_x_max, exp_x_max-exp_x_min+1)
    if ax_num == 3 or ax_num == 4:
        ax.set_xticks(xticks, [r"$10^{%i}$"%np.log10(x) for x in xticks])
    else:
        ax.set_xticks(xticks,["" for x in xticks])
    yticks = np.logspace(exp_y_min, exp_y_max, exp_y_max-exp_y_min+1)
    ax.set_yticks(yticks, [r"$10^{%i}$"%np.log10(x) for x in yticks])   
    # if ax_num == 1 or ax_num == 3:
    #     ax.set_yticks(yticks, [r"$10^{%i}$"%np.log10(x) for x in yticks])
    # else:
    #     ax.set_yticks(xticks,["" for x in xticks])
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])

def plot_stability_big(EoS, show = False):
    fig, [[ax1,ax2],[ax3,ax4]] = plt.subplots(2,2,figsize = (16,16))

    plt.subplots_adjust(wspace=0.12, hspace=0)   

    scale = 0.8

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax4.set_xscale("log")
    ax4.set_yscale("log")

    data = get_data("stability_results/%s_data_m_DM_%e.dat"%(EoS,500))
    sc = plot_stability(data, ax1)
    set_x_y_lim_ticks(data, ax1, 1, scale = scale)
    data = get_data("stability_results/%s_data_m_DM_%e.dat"%(EoS,1000))
    plot_stability(data, ax2)
    set_x_y_lim_ticks(data, ax2, 2, scale = scale)
    data = get_data("stability_results/%s_data_m_DM_%e.dat"%(EoS,2000))
    plot_stability(data, ax3)
    set_x_y_lim_ticks(data, ax3, 3, scale = scale)
    data = get_data("stability_results/%s_data_m_DM_%e.dat"%(EoS,4000))
    plot_stability(data, ax4)
    set_x_y_lim_ticks(data, ax4, 4, scale = scale)

    xtext = 0.02
    ytext = 0.98
    props = dict(facecolor='white', alpha=0.5)
    ax1.text(xtext, ytext, r'$m_{C} = %0.1f\,\text{GeV}$'%(0.5), ha='left', va='top', transform=ax1.transAxes, fontsize=24, bbox=props)
    ax2.text(xtext, ytext, r'$m_{C} = %i\,\text{GeV}$'%(1), ha='left', va='top', transform=ax2.transAxes, fontsize=24, bbox=props)
    ax3.text(xtext, ytext, r'$m_{C} = %i\,\text{GeV}$'%(2), ha='left', va='top', transform=ax3.transAxes, fontsize=24, bbox=props)
    ax4.text(xtext, ytext, r'$m_{C} = %i\,\text{GeV}$'%(4), ha='left', va='top', transform=ax4.transAxes, fontsize=24, bbox=props)

    ax1.grid(color = "lightgrey", alpha = 0.8, zorder = -7)
    ax2.grid(color = "lightgrey", alpha = 0.8, zorder = -7)
    ax3.grid(color = "lightgrey", alpha = 0.8, zorder = -7)
    ax4.grid(color = "lightgrey", alpha = 0.8, zorder = -7)

    xlabel = r"$p_{0,O}\,\text{[MeV/fm}^3\text{]}$"
    ylabel = r"$p_{0,D}\,\text{[MeV/fm}^3\text{]}$"

    ax3.set_xlabel(xlabel, fontsize = styles.fontsize)
    ax1.set_ylabel(ylabel, fontsize = styles.fontsize)
    ax3.set_ylabel(ylabel, fontsize = styles.fontsize)
    ax4.set_xlabel(xlabel, fontsize = styles.fontsize)

    tick_arr = [lower_frac_lim,2,3,4,6,8,max_dm_frac]
    
    # fig.subplots_adjust(right=0.94)
    cbar_ax = fig.add_axes([0.91, 0.2, 0.03, 0.6])
    cb = fig.colorbar(sc, cax=cbar_ax, fraction=0.5)
    # cb = plt.colorbar(sc)#, ticks=tick_arr)
    cb.Ticks = []
    cb.set_ticks(tick_arr)
    cb.set_label(r'$\text{Fraction of DM (%)}$', rotation=270, labelpad=15, fontsize = styles.fontsize)
    cb.ax.set_yticklabels([r"$\text{%i}$"%x for x in tick_arr])



    # ax1.text(0.01, 1.05, r"$\text{%s}$"%EoS_name(EoS), ha='left', va='top', transform=ax1.transAxes, fontsize=16)

    ax2.legend(loc="upper right",title = r"$\text{%s}$"%EoS_name(EoS), fontsize = styles.fontsize)
    fig.savefig('plots/stable_big_alt.pdf', bbox_inches = "tight")
    fig.savefig('plots/stable_big_alt.png', bbox_inches = "tight", dpi = 300)
    if show:
        plt.show()
    plt.close(fig)

# mass_list=[125, 190,250, 304.753414,  371.498572, 452.861832, 552.044757, 672.950096, 820.335356, 1000.000000, 1219.013654, 1485.994289, 1811.447329, 2208.179027, 2691.800385, 3281.341424, 4000.000000]

if __name__ == "__main__":
    plot_stability_big("light_EoS_III", show = False)




