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

scattersize = 11.5
scattersize_stab = 11.5

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
    OM_rad = np.genfromtxt(file, usecols = 2) # OM matter radius
    tot_mass = np.genfromtxt(file, usecols = 7) # total mass
    stability =  np.genfromtxt(file, usecols = 17) # compactness @R_D
    DM_mass = np.genfromtxt(file, usecols = 6) # OM mass
    dm_frac = np.array([100*DM_mass[x]/tot_mass[x] for x in range(len(tot_mass))])
    return OM_rad, tot_mass, stability, dm_frac

label_stable = r"$\text{stable}$"
label_stable1 = r"$\text{stable}<1\%$"

def plot_M_R(data, ax):
    OM_rad, tot_mass, stability, dm_frac = data
    stable = (stability == 1)
    stable_constraint = (stability == 1) & (dm_frac > lower_frac_lim) & (dm_frac <= max_dm_frac)
    stable_constraint_1 = (stability == 1) & (dm_frac <= lower_frac_lim)
    ax.scatter(OM_rad[stable],tot_mass[stable],  s = scattersize, label =  label_stable, color = stabilitycolor, marker="s",zorder=1)
    sc = ax.scatter(OM_rad[stable_constraint],tot_mass[stable_constraint],
                       c = dm_frac[stable_constraint],  s = scattersize, 
                       cmap=plt.get_cmap(colormp), norm = colors.LogNorm(vmin=lower_frac_lim, vmax=max_dm_frac), marker="s",zorder=2)   #, norm = colors.Normalize(vmin=lower_frac_lim, vmax=100-lower_frac_lim))   
    ax.scatter(OM_rad[stable_constraint_1],tot_mass[stable_constraint_1],  
                s = scattersize, color = consistency1color, label =  label_stable1, marker="s",zorder=3)
    return sc

def plot_M_R_big(EoS, show = False):
    fig, [[ax1,ax2],[ax3,ax4]] = plt.subplots(2,2,figsize = (16,12),sharex=False,sharey=False)
    # fig.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)   


    data = get_data("data/%s_data_m_DM_%e.dat"%(EoS,500))
    sc = plot_M_R(data, ax1)
    data = get_data("data/%s_data_m_DM_%e.dat"%(EoS,1000))
    plot_M_R(data, ax2)
    data = get_data("data/%s_data_m_DM_%e.dat"%(EoS,2000))
    plot_M_R(data, ax3)
    data = get_data("data/%s_data_m_DM_%e.dat"%(EoS,4000))
    plot_M_R(data, ax4)

    xtext = 0.98
    ytext = 0.02
    props = dict(facecolor='white', alpha=0.5)
    ax1.text(xtext, ytext, r'$m_{C} = %0.1f\,\text{GeV}$'%(0.5), ha='right', va='bottom', transform=ax1.transAxes, fontsize=24, bbox=props)
    ax2.text(xtext, ytext, r'$m_{C} = %i\,\text{GeV}$'%(1), ha='right', va='bottom', transform=ax2.transAxes, fontsize=24, bbox=props)
    ax3.text(xtext, ytext, r'$m_{C} = %i\,\text{GeV}$'%(2), ha='right', va='bottom', transform=ax3.transAxes, fontsize=24, bbox=props)
    ax4.text(xtext, ytext, r'$m_{C} = %i\,\text{GeV}$'%(4), ha='right', va='bottom', transform=ax4.transAxes, fontsize=24, bbox=props)

    draw_errors(ax1)
    draw_errors(ax2)
    draw_error_text(ax4)
    draw_errors(ax3)
    draw_errors(ax4)

    ax1.grid(color = "lightgrey", alpha = 0.8, zorder = -3)
    ax2.grid(color = "lightgrey", alpha = 0.8, zorder = -3)
    ax3.grid(color = "lightgrey", alpha = 0.8, zorder = -3)
    ax4.grid(color = "lightgrey", alpha = 0.8, zorder = -3)

    ax3.set_xlabel(r"$R_{O}\,\text{[km]}$", fontsize = styles.fontsize)
    ax4.set_xlabel(r"$R_{O}\,\text{[km]}$", fontsize = styles.fontsize)
    ax1.set_ylabel(r"$M_{tot}\,\text{[}M_\odot\text{]}$", fontsize = styles.fontsize)
    ax3.set_ylabel(r"$M_{tot}\,\text{[}M_\odot\text{]}$", fontsize = styles.fontsize)
    xticks = np.linspace(8,12,3)
    ax3.set_xticks(xticks, [r"$%i$"%x for x in xticks])
    xticks = np.linspace(8,14,4)
    ax4.set_xticks(xticks, [r"$%i$"%x for x in xticks])
    yticks = np.linspace(0.5,2.5,5)
    ax1.set_yticks(yticks, [r"$%1.1f$"%x for x in yticks])
    yticks = np.linspace(0.5,2,4)
    ax3.set_yticks(yticks, [r"$%1.1f$"%x for x in yticks])

    ax1.set_xticks(xticks,["" for x in xticks])
    ax2.set_xticks(xticks,["" for x in xticks])
    ax2.set_yticks(yticks,["" for x in yticks])
    ax4.set_yticks(yticks,["" for x in yticks])

    tick_arr = [lower_frac_lim,2,3,4,6,8,max_dm_frac]
    
    # fig.subplots_adjust(right=0.94)
    cbar_ax = fig.add_axes([0.91, 0.2, 0.03, 0.6])
    cb = fig.colorbar(sc, cax=cbar_ax, fraction=0.5)
    # cb = plt.colorbar(sc)#, ticks=tick_arr)
    cb.Ticks = []
    cb.set_ticks(tick_arr)
    cb.set_label(r'$\text{Fraction of DM (%)}$', rotation=270, labelpad=15, fontsize = styles.fontsize)
    cb.ax.set_yticklabels([r"$\text{%i}$"%x for x in tick_arr])

    ax1.set_xlim([8,14])
    ax2.set_xlim([8,14])
    ax3.set_xlim([8,14])
    ax4.set_xlim([8,14])
    ax1.set_ylim([0.5,2.5])
    ax2.set_ylim([0.5,2.5])
    ax3.set_ylim([0.5,2.5])
    ax4.set_ylim([0.5,2.5])

    # ax1.text(0.01, 1.05, r"$\text{%s}$"%EoS_name(EoS), ha='left', va='top', transform=ax1.transAxes, fontsize=16)

    ax1.legend(loc="upper left",title = r"$\text{%s}$"%EoS_name(EoS), fontsize = styles.fontsize)
    fig.savefig('plots/M_R_big.pdf', bbox_inches = "tight")
    if show:
        plt.show()
    plt.close(fig)

# mass_list=[125, 190,250, 304.753414,  371.498572, 452.861832, 552.044757, 672.950096, 820.335356, 1000.000000, 1219.013654, 1485.994289, 1811.447329, 2208.179027, 2691.800385, 3281.341424, 4000.000000]

if __name__ == "__main__":
    plot_M_R_big("light_EoS_II", show = False)




