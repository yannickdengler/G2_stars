import numpy as np
import matplotlib.pyplot as plt
import matplotlib

import styles_PP

font = {'size'   : 12}
matplotlib.rc('font', **font)

fig, ax = plt.subplots()

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$n$")
# ax.set_ylabel("$p$")

ax.set_xlim([0.01,100])
ax.set_ylim([0.01,100])

# ax.axhline(1)
fontsize = 16

xpoints = [0.01,0.025,0.25,0.8,20,100]
ypoints = [0.02,0.3  ,0.5 ,7 ,15,100]

ax.plot(xpoints[:2],ypoints[:2],ls="dashed",color = "black")
ax.plot(xpoints[4:],ypoints[4:],ls="dashed",color = "black")
ax.plot(xpoints[1:5],ypoints[1:5],color = "black")

ax.scatter(xpoints[1:5],ypoints[1:5],marker="x",color="black")

ax.annotate(r"$\Gamma_{i-1},c_{i-1},K_{i-1}$",[0.1,0.422],[0.06,1],arrowprops=dict(facecolor='black',headwidth=6,width=0.4),color="darkgreen",ha="center", fontsize=fontsize)
ax.annotate(r"$\Gamma_{i},c_{i},K_{i}$",[0.419,1.75],[0.03,5],arrowprops=dict(facecolor='black',headwidth=6,width=0.4),color="darkgreen", fontsize=fontsize)
ax.annotate(r"$\Gamma_{i+1},c_{i+1},K_{i+1}$",[2.14,9.12],[0.04,30],arrowprops=dict(facecolor='black',headwidth=6,width=0.4),color="darkgreen", fontsize=fontsize)

ax.text(x=0.06,y=0.22,s=r"$i-1$",color="red", fontsize=fontsize)
ax.text(x=0.476,y=1.264,s=r"$i$",color="red", fontsize=fontsize)
ax.text(x=4.5,y=6.5,s=r"$i+1$",color="red", fontsize=fontsize)

ax.text(x=0.288,y=0.345,s=r"$\mu_{i-1},n_{i-1},P_{i-1}$",color="blue", fontsize=fontsize)
ax.text(x=0.76,y=4.2,s=r"$\mu_{i},n_{i},P_{i}$",color="blue", fontsize=fontsize)

ax.set_yticklabels([])
ax.set_xticklabels([])

plt.savefig("piecewise_poly_sketch.pdf", bbox_inches = "tight")
# plt.show()