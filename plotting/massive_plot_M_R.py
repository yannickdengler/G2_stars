#!/usr/bin/env python
# coding: utf-8

import os, sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
from ellipses import draw_errors
from ellipses import draw_error_text

# import styles

import math

plt.rcParams['figure.figsize'] = [10, 6] 
fontsize = 36
font = {'size'   : fontsize}
matplotlib.rc('font', **font)
plt.rcParams.update({
    "mathtext.fontset": "cm",   # Computer Modern
})

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
    ax.scatter(OM_rad[stable],tot_mass[stable],  s = scattersize, color = stabilitycolor, marker="s", zorder = 1) #, label =  label_stable)
    sc = ax.scatter(OM_rad[stable_constraint],tot_mass[stable_constraint],
                       c = dm_frac[stable_constraint],  s = scattersize, 
                       cmap=plt.get_cmap(colormp), norm = colors.LogNorm(vmin=lower_frac_lim, vmax=max_dm_frac), marker="s", zorder = 2)   #, norm = colors.Normalize(vmin=lower_frac_lim, vmax=100-lower_frac_lim))   
    ax.scatter(OM_rad[stable_constraint_1],tot_mass[stable_constraint_1],  
                s = scattersize, color = consistency1color, zorder = 3) #, label =  label_stable1, marker="s")
    return sc

def do_four_plots(axs, EoS):
    [ax1,ax2,ax3,ax4] = axs
    data = get_data("data/%s_data_m_DM_%e.dat"%(EoS,500))
    sc = plot_M_R(data, ax1)
    data = get_data("data/%s_data_m_DM_%e.dat"%(EoS,1000))
    plot_M_R(data, ax2)
    data = get_data("data/%s_data_m_DM_%e.dat"%(EoS,2000))
    plot_M_R(data, ax3)
    data = get_data("data/%s_data_m_DM_%e.dat"%(EoS,4000))
    plot_M_R(data, ax4)

    xtext = 0.97
    ytext = 0.03
    props = dict(facecolor='white', alpha=0.7)
    ax1.text(xtext, ytext, r'$m_{C} = %0.1f\,\text{GeV}$'%(0.5), ha='right', va='bottom', transform=ax1.transAxes, fontsize=fontsize, bbox=props)
    ax2.text(xtext, ytext, r'$m_{C} = %i\,\text{GeV}$'%(1), ha='right', va='bottom', transform=ax2.transAxes, fontsize=fontsize, bbox=props)
    ax3.text(xtext, ytext, r'$m_{C} = %i\,\text{GeV}$'%(2), ha='right', va='bottom', transform=ax3.transAxes, fontsize=fontsize, bbox=props)
    ax4.text(xtext, ytext, r'$m_{C} = %i\,\text{GeV}$'%(4), ha='right', va='bottom', transform=ax4.transAxes, fontsize=fontsize, bbox=props)

    ax1.grid(color = "lightgrey", alpha = 0.8, zorder = -7)
    ax2.grid(color = "lightgrey", alpha = 0.8, zorder = -7)
    ax3.grid(color = "lightgrey", alpha = 0.8, zorder = -7)
    ax4.grid(color = "lightgrey", alpha = 0.8, zorder = -7)

    draw_errors(ax1)
    draw_errors(ax2)
    # draw_error_text(ax4)
    draw_errors(ax3)
    draw_errors(ax4)

    return sc


def plot_M_R_massive(show = False):
    fig, axss = plt.subplots(5,4,figsize = (32,40),sharex=True,sharey=True)
    plt.subplots_adjust(wspace=0, hspace=0)   

    xlim = [8,15]
    ylim = [0.5,2.6]

    EoSs=["light_EoS_I","heavy_EoS_I","heavy_EoS_II","light_EoS_III","heavy_EoS_III"]

    for i in range(len(EoSs)):
        sc = do_four_plots(axss[i],EoSs[i])

    for tmp in axss:
        for tmptmp in tmp:
            ticks_tmp = [0,]
            tmptmp.set_xticks(ticks_tmp,["" for x in ticks_tmp])
            tmptmp.set_yticks(ticks_tmp,["" for x in ticks_tmp])

    xtext = 0.03
    ytext = 0.97
    props = dict(facecolor='white', alpha=0.7)

    for i in range(len(EoSs)):
        yticks = [1,1.5,2]
        axss[i,0].set_yticks(yticks, [r"$%1.1f$"%x for x in yticks])
        axss[i,0].set_ylabel(r"$M_{tot}\,\text{[}M_\odot\text{]}$", fontsize = fontsize)
        axss[i,0].text(xtext, ytext, r'%s'%EoS_name(EoSs[i]), ha='left', va='top', transform=axss[i,0].transAxes, fontsize=fontsize, bbox=props)

    for tmp in np.transpose(axss):
        xticks = [9,10,11,12,13,14]
        tmp[4].set_xticks(xticks, [r"$%i$"%x for x in xticks])
        tmp[4].set_xlabel(r"$R_{O}\,\text{[km]}$", fontsize = fontsize)

    tick_arr = [lower_frac_lim,2,3,4,6,8,max_dm_frac]
    
    cbar_ax = fig.add_axes([0.91, 0.2, 0.03, 0.6])
    cb = fig.colorbar(sc, cax=cbar_ax, fraction=0.5)
    cb.Ticks = []
    cb.set_ticks(tick_arr)
    cb.set_label(r'$\text{Fraction of DM (%)}$', rotation=270, labelpad=15, fontsize = fontsize)
    cb.ax.set_yticklabels([r"$\text{%i}$"%x for x in tick_arr])

    axss[1,1].set_xlim(xlim)
    axss[1,1].set_ylim(ylim)
    fig.savefig('plots/M_R_massive.pdf', bbox_inches = "tight", dpi = 100)
    fig.savefig('plots/M_R_massive.png', bbox_inches = "tight", dpi = 300)
    if show:
        plt.show()
    plt.close(fig)

if __name__ == "__main__":
    plot_M_R_massive(show = False)




