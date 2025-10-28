import matplotlib.pyplot as plt
import numpy as np

NICER = np.transpose(np.genfromtxt("NICER.csv",delimiter=","))
GW170817 = np.transpose(np.genfromtxt("GW170817.csv",delimiter=","))

print(len(NICER))


plt.fill(NICER[0],NICER[1], alpha = 0.5)
plt.fill(GW170817[0],GW170817[1], alpha = 0.5)


PSR_low_up = [2.14-0.18,2.14+0.2]

plt.axhspan(PSR_low_up[0],PSR_low_up[1], color = "grey", alpha = 0.3)

plt.show()