import numpy as np
import matplotlib.pyplot as plt
import global_vars as g
# from scipy.interpolate import interp1d as inter
import scipy.interpolate


class interpolate1d(scipy.interpolate.interp1d):
    """Extend scipy interp1d to interpolate/extrapolate per axis in log space"""
    
    def __init__(self, x, y, *args, xspace='linear', yspace='linear', **kwargs):
        self.xspace = xspace
        self.yspace = yspace
        if self.xspace == 'log': x = np.log10(x)
        if self.yspace == 'log': y = np.log10(y)
        super().__init__(x, y, *args, **kwargs)
        
    def __call__(self, x, *args, **kwargs):
        if self.xspace == 'log': x = np.log10(x)
        if self.yspace == 'log':
            return 10**super().__call__(x, *args, **kwargs)
        else:
            return super().__call__(x, *args, **kwargs)

if __name__ == "__main__":
    neutron_mass = g.neutron_mass  # MeV
    hbarc = g.MeV_fm
    n0 = 0.16                                       # 0.16 fm^-3

    n_Heb, rho, P_soft, eps_soft, e_soft, R_soft, M_soft, P_immi, eps_immi, e_immi, R_immi, M_immi, P_hard, eps_hard, e_hard, R_hard, M_hard = np.transpose(np.genfromtxt("EoS/data/EoS_Hebeler"))
    n_I, P_I, eps_I, mu_I, c_s_I, R_I, M_I = np.transpose(np.genfromtxt("EoS/data/EoS_I"))
    n_II, P_II, eps_II, mu_II, c_s_II, R_II, M_II = np.transpose(np.genfromtxt("EoS/data/EoS_II"))
    n_III, P_III, eps_III, mu_III, c_s_III, R_III, M_III = np.transpose(np.genfromtxt("EoS/data/EoS_III"))

    n_Heb = n_Heb*n0
    n_I = n_I*n0
    n_II = n_II*n0
    n_III = n_III*n0
    
    mu_I = mu_I*1000
    mu_II = mu_II*1000
    mu_III = mu_III*1000
    
    for i in range(len(eps_I)):
        print(P_I[i],eps_I[i],n_I[i],mu_I[i],P_I[i]+eps_I[i],n_I[i]*mu_I[i])

    EoS_soft = interpolate1d(P_soft, eps_soft, xspace="log", yspace="log")
    EoS_immi = interpolate1d(P_immi, eps_immi, xspace="log", yspace="log")
    EoS_hard = interpolate1d(P_hard, eps_hard, xspace="log", yspace="log")
    n_soft = interpolate1d(P_soft, n_Heb, xspace="log", yspace="log")
    n_immi = interpolate1d(P_immi, n_Heb, xspace="log", yspace="log")
    n_hard = interpolate1d(P_hard, n_Heb, xspace="log", yspace="log")

    EoS_I = interpolate1d(P_I, eps_I, xspace="log", yspace="log")
    EoS_II = interpolate1d(P_II, eps_II, xspace="log", yspace="log")
    EoS_III = interpolate1d(P_III, eps_III, xspace="log", yspace="log")
    n_func_I = interpolate1d(P_I, n_I, xspace="log", yspace="log", fill_value="extrapolate")
    n_func_II = interpolate1d(P_II, n_II, xspace="log", yspace="log", fill_value="extrapolate")
    n_func_III = interpolate1d(P_III, n_III, xspace="log", yspace="log", fill_value="extrapolate")
    mu_func_I = interpolate1d(P_I, mu_I, xspace="log", yspace="log", fill_value="extrapolate")
    mu_func_II = interpolate1d(P_II, mu_II, xspace="log", yspace="log", fill_value="extrapolate")
    mu_func_III = interpolate1d(P_III, mu_III, xspace="log", yspace="log", fill_value="extrapolate")

    gamma_53 = 5./3

    K_soft = min(P_soft)/(min(eps_soft))**(gamma_53)              # nochmal überprüfen!!! Hararison Wheeler?
    K_immi = min(P_immi)/(min(eps_immi))**(gamma_53)
    K_hard = min(P_hard)/(min(eps_hard))**(gamma_53)
    
    def EoS_I_func(P):
        if P < 0.447:
            poly = (P/K_soft)**(1/gamma_53)
            if poly < P:
                return P, 0
            else: 
                return poly, 0
        elif P < 2.163:
            return EoS_soft(P), n_soft(P)
        elif P < max(P_I):
            return EoS_I(P), n_func_I(P)
        else:
            return max(eps_I)+(P-max(P_I)), -3

    def EoS_II_func(P):
        if P < 0.696:
            poly = (P/K_hard)**(1/gamma_53)
            if poly < P:
                return P, 0
            else: 
                return poly, 0
        elif P < 3.542:
            return EoS_hard(P), n_hard(P)
        elif P < max(P_II):
            return EoS_II(P), n_func_II(P)
        else:
            return max(eps_II)+(P-max(P_II)), -3
    
    def EoS_III_func(P):
        if P < 0.696:
            poly = (P/K_hard)**(1/gamma_53)
            if poly < P:
                return P, 0
            else: 
                return poly, 0
        elif P < 3.542:
            return EoS_hard(P), n_hard(P)
        elif P < max(P_III):
            return EoS_III(P), n_func_III(P)
        else:
            return max(eps_III)+(P-max(P_III)), -3
            
    P_arr = []

    # numP = [100,500,10]
    numP = [20,50,5]

    for P in np.logspace(-10,-2,numP[0], base = 10):
        P_arr.append(P)
    for P in np.logspace(-2,3,numP[1], base = 10):
        P_arr.append(P)
    for P in np.logspace(3,6,numP[2], base = 10):
        P_arr.append(P)

    with open("EoS/data/EoS_I_inter","w") as f:
        f.write("%e\t%e\t%e\t%e\n"%(0,0,0,neutron_mass))
        for P in P_arr:
            eps, n = EoS_I_func(P)
            if n == 0:
                mu = neutron_mass
                n = (eps+P)/mu
            elif n == -3:
                mu = mu_func_I(P)
                n = (eps+P)/mu
            else:
                mu = (eps+P)/n
            f.write("%e\t%e\t%e\t%e\n"%(P, eps, n, mu))
        f.write("%e\t%e\t%e\t%e\n"%(1e10,1e10,1e10,1e10))

    with open("EoS/data/EoS_II_inter","w") as f:
        f.write("%e\t%e\t%e\t%e\n"%(0,0,0,neutron_mass))
        for P in P_arr:
            eps, n = EoS_II_func(P)
            if n == 0:
                mu = neutron_mass
                n = (eps+P)/mu
            elif n == -3:
                mu = mu_func_II(P)
                n = (eps+P)/mu
            else:
                mu = (eps+P)/n
            f.write("%e\t%e\t%e\t%e\n"%(P, eps, n, mu))
        f.write("%e\t%e\t%e\t%e\n"%(1e10,1e10,1e10,1e10))

    with open("EoS/data/EoS_III_inter","w") as f:
        f.write("%e\t%e\t%e\t%e\n"%(0,0,0,neutron_mass))
        for P in P_arr:
            eps, n = EoS_III_func(P)
            if n == 0:
                mu = neutron_mass
                n = (eps+P)/mu
            elif n == -3:
                mu = mu_func_III(P)
                n = (eps+P)/mu
            else:
                mu = (eps+P)/n
            f.write("%e\t%e\t%e\t%e\n"%(P, eps, n, mu))
        f.write("%e\t%e\t%e\t%e\n"%(1e10,1e10,1e10,1e10))


    with open("EoS/data/EoS_I_inter_P_eps","w") as f:
        f.write("%e\t%e\n"%(1e-50,1e-50))
        for P in P_arr:
            eps, n = EoS_I_func(P)
            f.write("%e\t%e\n"%(P, eps))
        f.write("%e\t%e\n"%(1e30,1e30))

    with open("EoS/data/EoS_II_inter_P_eps","w") as f:
        f.write("%e\t%e\n"%(1e-50,1e-50))
        for P in P_arr:
            eps, n = EoS_II_func(P)
            f.write("%e\t%e\n"%(P, eps))
        f.write("%e\t%e\n"%(1e10,1e10))

    with open("EoS/data/EoS_III_inter_P_eps","w") as f:
        f.write("%e\t%e\n"%(1e-50,1e-50))
        for P in P_arr:
            eps, n = EoS_III_func(P)
            f.write("%e\t%e\n"%(P, eps))
        f.write("%e\t%e\n"%(1e10,1e10))

    plt.xscale("log")
    plt.yscale("log")

    P_I_inter, eps_I_inter, n_I_inter, mu_I_inter = np.transpose(np.genfromtxt("EoS/data/EoS_I_inter"))
    P_II_inter, eps_II_inter, n_II_inter, mu_II_inter = np.transpose(np.genfromtxt("EoS/data/EoS_II_inter"))
    P_III_inter, eps_III_inter, n_III_inter, mu_III_inter = np.transpose(np.genfromtxt("EoS/data/EoS_III_inter"))

    plt.scatter(P_soft, eps_soft, label = "soft", marker = "^", s = 1)
    plt.scatter(P_immi, eps_immi, label = "immi", marker = "<", s = 1)
    plt.scatter(P_hard, eps_hard, label = "hard", marker = "v", s = 1)

    plt.scatter(P_I, eps_I, label = "EoSI", color = "green")
    plt.scatter(P_II, eps_II, label = "EoSII", color = "red")#, linewidth = 3)
    plt.scatter(P_III, eps_III, label = "EoSIII", color = "blue")#, linewidth = 3)


    plt.plot(P_I_inter, eps_I_inter, label = "inter I", color = "lime")
    plt.plot(P_II_inter, eps_II_inter, label = "inter II", color = "orange")
    plt.plot(P_III_inter, eps_III_inter, label = "inter III", color = "cyan")
    plt.plot(P_I_inter, n_I_inter, color = "lime", ls = "dashed")
    plt.plot(P_II_inter, n_II_inter, color = "orange", ls = "dashed")
    plt.plot(P_III_inter, n_III_inter, color = "cyan", ls = "dashed")
    plt.plot(P_I_inter, mu_I_inter, color = "lime", ls = "dotted")
    plt.plot(P_II_inter, mu_II_inter, color = "orange", ls = "dotted")
    plt.plot(P_III_inter, mu_III_inter, color = "cyan", ls = "dotted")

    # plt.xlim([1e-5,1e-1])
    # plt.ylim([1e-3,1e-1])

    plt.legend()

    # plt.show()