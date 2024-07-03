'''
This file takes the G2 EoS (n(mu)) and
1. calculates the physical parameters (in terms of DM mass)
2. avarages it
3. performs the piecewise polytropes (PP)
and writes the interpolated result in a "look up table".

Every quantity is scaled with the DM mass mu(P=0, eps=0) = 1

There are 3 regimes: 
1. Free fermigas n = c*(mu-m_DM)^(3/2)     (Automatically done by PP)
2. PP
3. Fermi-Dirac: n = n_s/(exp(a-b*mu)+1)    (Calculated and then done by PP)
Inputs: low mu cut-off, high mu cut-off

Linear interpolation is faster later than calculating the PP
'''

import numpy as np
import global_vars as g
import matplotlib.pyplot as plt
from scipy.optimize import bisect

def delete_doubles():
    '''
    Deletes the doubles and avarages over the entries of same mu
    Also makes the EoS unitless by scaling it with a DM mass of (1 MeV)
    '''
    EoS_OG = np.transpose(np.genfromtxt(g.PATH+"EoS/data/EoS_G2_lat"))                # In the format [mu,n,nerr]

    EoS_res = []
    for i in range(len(EoS_OG)):
        EoS_res.append([])
    
    for i in range(len(EoS_OG[0])-1):
        if EoS_OG[0][i] == EoS_OG[0][i+1]:
            EoS_res[0].append(EoS_OG[0][i])                         
            EoS_res[1].append(0.5*(EoS_OG[1][i]+EoS_OG[1][i+1]))
            EoS_res[2].append(0.5*(EoS_OG[1][i]+EoS_OG[1][i+1]))
        elif EoS_OG[0][i] != EoS_OG[0][i-1]:
            EoS_res[0].append(EoS_OG[0][i])
            EoS_res[1].append(EoS_OG[1][i])
            EoS_res[2].append(EoS_OG[2][i])

    EoS_res[0].append(EoS_OG[0][len(EoS_OG[0])-1])                                  # add the last line
    EoS_res[1].append(EoS_OG[1][len(EoS_OG[0])-1])
    EoS_res[2].append(EoS_OG[2][len(EoS_OG[0])-1])

    for i in range(len(EoS_res[0])):
        EoS_res[0][i] = EoS_res[0][i]*g.num_const/g.mn_light_lat                     # mu_prime
        EoS_res[1][i] = EoS_res[1][i]/g.num_const/g.mn_light_lat**3                  # n_prime
        EoS_res[2][i] = EoS_res[2][i]/g.num_const/g.mn_light_lat**3                  # n_prime

    with open(g.PATH+"EoS/data/EoS_G2_prime","w") as f:
        for i in range(len(EoS_res[0])):
            f.write("%e\t%e\t%e\n"%(EoS_res[0][i],EoS_res[1][i],EoS_res[2][i]))

################## From this point everything is unitless and scaled by the DM mass

def fermi_dirac(mu,a,b):
    return g.lat_sat_den_prime/(np.exp(abs(a)-abs(b)*mu)+1)
def get_ab(mu1,mu2,n1,n2):
    a = np.log((n2*(g.lat_sat_den_prime-n1))/(n1*(g.lat_sat_den_prime-n2)))*mu1/(mu2-mu1)+np.log(g.lat_sat_den_prime/n1-1)
    b = np.log((n2*(g.lat_sat_den_prime-n1))/(n1*(g.lat_sat_den_prime-n2)))/(mu2-mu1)
    return a, b

def add_fermi_dirac(num_FD = 20):               # FD stands for Fermi-Dirac
    '''
    Extent the data from "delete_doubles" with an interpolation using fermi-dirac statistics
    '''
    EoS_in = np.transpose(np.genfromtxt(g.PATH+"EoS/data/EoS_G2_prime"))                # In the format [mu,n,nerr]

    leng = len(EoS_in[0])
    a,b=get_ab(EoS_in[0][leng-4],EoS_in[0][leng-2],EoS_in[1][leng-4],EoS_in[1][leng-2])                 # used "leng-4" because it is in decent agreement with 1312.5579

    xarr = np.linspace(EoS_in[0][leng-2]*1.1,EoS_in[0][leng-1], num_FD)
    yarr = fermi_dirac(xarr,a,b)

    with open("EoS/data/EoS_G2_prime_FD","w") as f: 
        for i in range(len(EoS_in[0])-1):
            if EoS_in[0][i] > 1:
                f.write("%e\t%e\n"%(EoS_in[0][i],EoS_in[1][i]))
        for i in range(len(xarr)):
            f.write("%e\t%e\n"%(xarr[i], yarr[i]))

def next_gamma(mu0, mu1, K0, n0, n1, gam0):
    print(mu0, mu1, K0, n0, n1, gam0)
    def zero_func(gamma):
        return mu1-mu0+K0*gamma/(1-gamma)*n0**(gam0-gamma)*(n1**(gamma-1)-n0**(gamma-1))              # should be >0 for gamma = 0
    for b in np.logspace(-3,10,200):
        # print(b)
        if zero_func(b) < 0:
            return bisect(f=zero_func, a=0, b=b)
    print("Error in next_gamma")
    return None

def next_K(n_0, gamma_1, P_0):
    return P_0/n_0**gamma_1

def next_c(mu_1, n_1, K_1, gamma_1):
    return mu_1 - K_1*gamma_1*n_1**(gamma_1-1)/(gamma_1-1)

def next_P(n_1, gamma_1, K_1):
    return K_1*n_1**gamma_1

def PP(gamma0 = 5./3, c0 = 1):
    EoS_in = np.transpose(np.genfromtxt(g.PATH+"EoS/data/EoS_G2_prime_FD"))

    mu_arr = EoS_in[0]
    n_arr = EoS_in[1]
    gamma_arr = np.zeros(len(mu_arr))
    c_arr = np.zeros(len(mu_arr))
    K_arr = np.zeros(len(mu_arr))
    P_arr = np.zeros(len(mu_arr))
    eps_arr = np.zeros(len(mu_arr))

    gamma_arr[0] = gamma0
    c_arr[0] = c0
    K_arr[0] = (mu_arr[0]-c_arr[0])*(gamma_arr[0]-1)*n_arr[0]**(1-gamma_arr[0])/gamma_arr[0]
    P_arr[0] = K_arr[0]*n_arr[0]**gamma_arr[0]
    for i in range(1, len(mu_arr)):
        print(i)
        gamma_arr[i] = next_gamma(mu_arr[i-1], mu_arr[i], K_arr[i-1], n_arr[i-1], n_arr[i], gamma_arr[i-1])

        K_arr[i] = next_K(n_arr[i-1], gamma_arr[i], P_arr[i-1])
        c_arr[i] = next_c(mu_arr[i], n_arr[i], K_arr[i], gamma_arr[i])
        P_arr[i] = next_P(n_arr[i], gamma_arr[i], K_arr[i])
        eps_arr[i] = mu_arr[i]*n_arr[i]-P_arr[i]
    
    with open("EoS/PP","w") as f:
        for i in range(len(mu_arr)):
            f.write("%e\t%e\t%e\t%e\t%e\t%e\t%e\n"%(mu_arr[i],n_arr[i],P_arr[i],eps_arr[i],gamma_arr[i],K_arr[i],c_arr[i]))

def get_index(n, n_arr):
    index = 0
    if n > n_arr[len(n_arr)-1]:
        index = 999
    else:
        for i in range(len(n_arr)):
            if n > n_arr[i] and n <= n_arr[i+1]:
                index = i+1
    return index

def mu_n(n, c, gamma, K):
    return c+gamma*K*n**(gamma-1)/(gamma-1)

# def n_mu(mu,c,gamma,K)

def P_n(n, K, gamma):
    return K*n**gamma

def eps_n(n, c, K, gamma):
        return c*n+K*n**gamma/(gamma-1)

def c_s_n(n, c, K, gamma):
    c_s = gamma*P_n(n, K, gamma)/(eps_n(n, c, K, gamma)+P_n(n, K, gamma))
    return c_s
    if c_s > 1:
        # print(n, c_s, "c_s > 1")
        return 1
    else:
        return c_s

def lin_inter(num = 10000):
    n_crit = 0.98561                                                                                       # n > n_crit -> c_s > 1!!        
    ind_crit = 0 
    mu_arr, n_arr, P_arr, eps_arr, gamma_arr, K_arr, c_arr = np.transpose(np.genfromtxt("EoS/PP"))
    
    # n_w = np.linspace(0, max(n_arr), num)
    n_w = np.linspace(0, 1.5, num)
    mu_w = np.zeros(num)
    P_w = np.zeros(num)
    eps_w = np.zeros(num)
    cs_w = np.zeros(num)
    for i in range(num):
        n_tmp = n_w[i]
        if n_tmp < n_crit:  
            index = get_index(n_tmp, n_arr)
            mu_w[i] = mu_n(n_tmp, c_arr[index], gamma_arr[index], K_arr[index])
            P_w[i] = P_n(n_tmp, K_arr[index], gamma_arr[index])
            eps_w[i] = eps_n(n_tmp,c_arr[index], K_arr[index], gamma_arr[index])
            cs_w[i] = c_s_n(n_tmp,c_arr[index], K_arr[index], gamma_arr[index])
            ind_crit = i
        else:
            index = get_index(n_w[ind_crit], n_arr)
            P_w[i] = K_arr[index] * n_w[i]**((eps_w[i-1]+P_w[i-1])/P_w[i-1])
            print(P_w[i], P_w[ind_crit])
            eps_w[i] = eps_w[ind_crit] + P_w[i] - P_w[ind_crit]
            cs_w[i] = 0.99999999999999999999999999
            mu_w[i] = (P_w[i]+eps_w[i])/n_w[i]

    with open("EoS/data/EoS_G2_PP_inter", "w") as f:
        for i in range(num):
            f.write("%f\t%f\t%f\t%f\t%f\n"%(n_w[i],mu_w[i],P_w[i],eps_w[i],cs_w[i]))


    with open("EoS/data/EoS_G2_PP_inter_P_eps", "w") as f:              # writes in unitless. P = P/m_F^4
        for i in range(num):
            f.write("%f\t%f\n"%(P_w[i],eps_w[i]))
        f.write("%f\t%f\n"%(1e30,1e30))

    # plt.xlabel("n")
    # plt.scatter(n_arr, mu_arr, label = "mu_OG")
    # plt.scatter(n_arr, P_arr, label = "P_PP")
    # plt.scatter(n_arr, eps_arr, label = "eps_PP")
    # plt.plot(n_w, mu_w, label = "mu inter")
    # plt.plot(n_w, P_w, label = "P inter")
    # plt.plot(n_w, eps_w, label = "eps inter")
    # plt.plot(n_w, cs_w, label = "cs inter")
    # plt.legend()
    # plt.show()












if __name__ == "__main__":
    delete_doubles()
    add_fermi_dirac()
    PP()
    lin_inter()
