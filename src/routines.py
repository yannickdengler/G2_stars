import integration as integrate
import functions as func
import global_vars as v
import numpy as np
import TOV
import sys
sys.path.insert(0, '/home/dengler_yannick/Documents/G2_stars/EoS/scripts')
import EoS

# almost_zero = 1e-70
almost_zero = 1e-20

################################### print EoS 1 fluid

def print_EoS_mu_1_fluid(mu_start, mu_end, steps):
    f = open("results/EoS_mu_1_fluid.out","w")
    f.close()
    for i in range(steps):
        mu = mu_start*(mu_end/mu_start)**(i/(steps-1))
        f = open("results/EoS_mu_1_fluid.out","a")
        f.write("%e\t%e\t%e\t%e\n"%(mu,EoS.P_mu(mu),EoS.eps_mu(mu),EoS.c_s_mu(mu)))
        f.close()

def print_EoS_P_1_fluid(P_start, P_end, steps):
    f = open("results/EoS_P_1_fluid.out","w")
    f.close()
    for i in range(steps):
        P = P_start*(P_end/P_start)**(i/(steps-1))
        # print(P)
        f = open("results/EoS_P_1_fluid.out","a")
        f.write("%e\t%e\t%e\t%e\n"%(P,EoS.mu_P(P),EoS.eps_P(P),EoS.c_s_P(P)))
        f.close()

################################### print EoS 2 fluid

def print_EoS_mu_2_fluid(mu_start, mu_end, steps):
    f = open("results/EoS_mu_2_fluid.out","w")
    f.close()
    for i in range(steps):
        mu = mu_start*(mu_end/mu_start)**(i/(steps-1))
        f = open("results/EoS_mu_2_fluid.out","a")
        f.write("%e\t%e\t%e\t%e\t%e\t%e\t%e\n"%(mu,EoS.P_OM_mu(mu), EoS.eps_OM_mu(mu), EoS.c_s_OM_mu(mu),EoS.P_DM_mu(mu), EoS.eps_DM_mu(mu), EoS.c_s_DM_mu(mu)))
        f.close()

def print_EoS_P_2_fluid(P_start, P_end, steps, pref = ""):
    f = open("results/EoS_P_2_fluid"+pref+".out","w")
    f.close()
    for i in range(steps):
        P = P_start*(P_end/P_start)**(i/(steps-1))
        f = open("results/EoS_P_2_fluid"+pref+".out","a")
        # f.write("%e\t%e\t%e\t%e\t%e\t%e\t%e\n"%(P,EoS.mu_OM_P(P), EoS.eps_OM_P(P), EoS.c_s_OM_P(P),EoS.mu_DM_P(P), EoS.eps_DM_P(P), EoS.c_s_DM_P(P)))
        f.write("%e\t%e\t%e\t%e\t%e\n"%(P, EoS.eps_OM_P(P), EoS.c_s_OM_P(P), EoS.eps_DM_P(P), EoS.c_s_DM_P(P)))
        f.close()

################################### single 1 fluid

def single_TIDAL_mu(init_mu, stepsize):
    f = open("results/single_mu_1.out","w")
    f.close()
    mu = init_mu
    y = 2
    r = stepsize
    M = 4*v.pi/3*EoS.eps_mu(mu)*r**3
    while( EoS.P_mu(mu) > 0):
        r_safe = r
        r, (mu,M,y), stepsize = integrate.adapt_stepsize_step(r=r, stepsize=stepsize, funcs=(func.d_mu_r, func.d_M_r_mu,func.d_y_r_mu), vals=(mu,M,y), limit=1e-3)
        # r, (mu,M,y) = integrate.Euler_step(r=r, stepsize=stepsize, funcs=(func.d_mu_r, func.d_M_r_mu,func.d_y_r_mu), vals=(mu,M,y))
        # r, (mu,M,y) = integrate.RK_4_step(r=r, stepsize=stepsize, funcs=(func.d_mu_r, func.d_M_r_mu,func.d_y_r_mu), vals=(mu,M,y))
        if r != r_safe:
            f = open("results/single_mu_1.out","a")
            P = EoS.P_mu(mu)
            eps = EoS.eps_mu(mu)
            c_s = EoS.c_s_mu(mu)
            f.write("%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n"%(r, stepsize, M, mu, P, eps, c_s, y, func.Q_mu(r, mu, M), func.F_mu(r, mu, M),func.d_y_r_mu(r, (mu,M,y))))
            f.close()
        if EoS.P_mu(mu) < almost_zero:
            mu = 0

def single_TIDAL_P(init_P, stepsize):
    f = open("results/single_P_1.out","w")
    f.close()
    P = init_P
    y = 2
    r = stepsize
    M = 4*v.pi/3*EoS.eps_P(P)*r**3
    while( P > 0):
        r_safe = r
        r, (P,M,y), stepsize = integrate.adapt_stepsize_step(r=r, stepsize=stepsize, funcs=(func.d_P_r, func.d_M_r_P,func.d_y_r_P), vals=(P,M,y), limit=1e-3)
        # r, (P,M,y) = integrate.Euler_step(r=r, stepsize=stepsize, funcs=(func.d_P_r, func.d_M_r_P,func.d_y_r_P), vals=(P,M,y))
        # r, (P,M,y) = integrate.RK_4_step(r=r, stepsize=stepsize, funcs=(func.d_P_r, func.d_M_r_P,func.d_y_r_P), vals=(P,M,y))
        if r != r_safe:
            f = open("results/single_P_1.out","a")
            eps = EoS.eps_P(P)
            c_s = EoS.c_s_P(P)
            f.write("%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n"%(r, stepsize, M, EoS.mu_P(P), P, eps, c_s, y, func.Q_P(r, P, M), func.F_P(r, P, M),func.d_y_r_P(r, (P,M,y))))
            f.close()
        if P < almost_zero:
            P = 0

################################### TOV_TIDAL 2 fluid asdasdasd

def single_TIDAL_mu_2_fluid(init_mu, stepsize, dm_om_ratio):
    f = open("results/single_mu_2.out","w")
    f.close()
    mu_OM = init_mu
    if dm_om_ratio < 0:
        mu_DM = -dm_om_ratio*init_mu
        mu_OM = 0
    else:
        mu_DM = dm_om_ratio*init_mu
    y = 2
    r = stepsize
    M_OM = 4*v.pi/3*EoS.eps_OM_mu(mu_OM)*r**3
    M_DM = 4*v.pi/3*EoS.eps_DM_mu(mu_DM)*r**3
    OM_fin = 1
    DM_fin = 1
    R_OM = 0
    R_DM = 0
    r_safe = 0
    while( EoS.P_OM_mu(mu_OM) > 0 or EoS.P_DM_mu(mu_DM) > 0 ):
        r_safe = r
        r, (mu_OM, mu_DM,M_OM, M_DM, y), stepsize = integrate.adapt_stepsize_step(r=r, stepsize=stepsize, funcs=(func.d_mu_r_OM, func.d_mu_r_DM, func.d_M_r_mu_OM, func.d_M_r_mu_DM, func.d_y_r_mu_2), vals=(mu_OM,mu_DM,M_OM,M_DM,y), limit=1e-3)
        # r, (mu_OM, mu_DM,M_OM, M_DM, y) = integrate.Euler_step(r=r, stepsize=stepsize, funcs=(func.d_mu_r_OM, func.d_mu_r_DM, func.d_M_r_mu_OM, func.d_M_r_mu_DM, func.d_y_r_mu_2), vals=(mu_OM,mu_DM,M_OM,M_DM,y))
        # r, (mu_OM, mu_DM,M_OM, M_DM, y) = integrate.RK_4_step(r=r, stepsize=stepsize, funcs=(func.d_mu_r_OM, func.d_mu_r_DM, func.d_M_r_mu_OM, func.d_M_r_mu_DM, func.d_y_r_mu_2), vals=(mu_OM,mu_DM,M_OM,M_DM,y))
        if r != r_safe:
            f = open("results/single_mu_2.out","a")
            P_OM = EoS.P_OM_mu(mu_OM)
            eps_OM = EoS.eps_OM_mu(mu_OM)
            c_s_OM = EoS.c_s_OM_mu(mu_OM)
            P_DM = EoS.P_DM_mu(mu_DM)
            eps_DM = EoS.eps_DM_mu(mu_DM)
            c_s_DM = EoS.c_s_DM_mu(mu_DM)
            M = M_OM + M_DM
            f.write("%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n"%(r, stepsize, M_OM, M_DM, mu_OM, mu_DM, P_OM, P_DM, eps_OM, eps_DM, c_s_OM, c_s_DM, y, func.Q_mu_2_fluid(r, mu_OM, mu_DM, M), func.F_mu_2_fluid(r, mu_OM, mu_DM, M),func.d_y_r_mu_2(r, (mu_OM, mu_DM,M_OM, M_DM, y))))
            f.close()
        if EoS.P_OM_mu(mu_OM) < almost_zero:
            mu_OM = 0
        if EoS.P_DM_mu(mu_DM) < almost_zero:
            mu_DM = 0
        print("P:",EoS.P_OM_mu(mu_OM),EoS.P_DM_mu(mu_DM))

def single_TIDAL_P_2_fluid(init_P, dm_om_ratio, init_stepsize=1e-20):
    f = open("results/single_P_2.out","w")
    f.close()
    P_OM = init_P
    if dm_om_ratio < 0:
        P_DM = -dm_om_ratio*init_P
        P_OM = 0
    else:
        P_DM = dm_om_ratio*init_P
    y = 2
    r = init_stepsize
    stepsize = init_stepsize
    M_OM = 4*v.pi/3*EoS.eps_OM_P(P_OM)*r**3
    M_DM = 4*v.pi/3*EoS.eps_DM_P(P_DM)*r**3
    OM_fin = 1
    DM_fin = 1
    R_OM = 0
    R_DM = 0
    r_safe = 0
    while( P_OM > 0 or P_DM > 0 ):
        r_safe = r
        r, (P_OM, P_DM,M_OM, M_DM, y), stepsize = integrate.adapt_stepsize_step(r=r, stepsize=stepsize, funcs=(func.d_P_r_OM, func.d_P_r_DM, func.d_M_r_P_OM, func.d_M_r_P_DM, func.d_y_r_P_2), vals=(P_OM,P_DM,M_OM,M_DM,y), limit=1e-4)
        # r, (P_OM, P_DM,M_OM, M_DM, y) = integrate.Euler_step(r=r, stepsize=stepsize, funcs=(func.d_P_r_OM, func.d_P_r_DM, func.d_M_r_P_OM, func.d_M_r_P_DM, func.d_y_r_P_2), vals=(P_OM,P_DM,M_OM,M_DM,y))
        # r, (P_OM, P_DM,M_OM, M_DM, y) = integrate.RK_4_step(r=r, stepsize=stepsize, funcs=(func.d_P_r_OM, func.d_P_r_DM, func.d_M_r_P_OM, func.d_M_r_P_DM, func.d_y_r_P_2), vals=(P_OM,P_DM,M_OM,M_DM,y))
        if P_OM < almost_zero:
            P_OM = 0
        if P_DM < almost_zero:
            P_DM = 0
        if stepsize < almost_zero or stepsize > 1/almost_zero:
            P_OM = 0
            P_DM = 0        
        # if r != r_safe:
        if True:
            f = open("results/single_TIDAL_P_2.out","a")
            mu_OM = EoS.mu_OM_P(P_OM)
            eps_OM = EoS.eps_OM_P(P_OM)
            c_s_OM = EoS.c_s_OM_P(P_OM)
            mu_DM = EoS.mu_DM_P(P_DM)
            eps_DM = EoS.eps_DM_P(P_DM)
            c_s_DM = EoS.c_s_DM_P(P_DM)
            M = M_OM + M_DM
            f.write("%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n"%(r, stepsize, M_OM, M_DM, mu_OM, mu_DM, P_OM, P_DM, eps_OM, eps_DM, c_s_OM, c_s_DM, y, func.Q_P_2_fluid(r, P_OM, P_DM, M), func.F_P_2_fluid(r, P_OM, P_DM, M),func.d_y_r_P_2(r, (P_OM, P_DM,M_OM, M_DM, y))))
            f.close()

def single_P_2_fluid(init_P, dm_om_ratio, init_stepsize=1e-20):
    f = open("results/single_P_2.out","w")
    f.close()
    P_OM = init_P
    if dm_om_ratio < 0:
        P_DM = -dm_om_ratio*init_P
        P_OM = 0
    else:
        P_DM = dm_om_ratio*init_P
    r = init_stepsize
    stepsize = init_stepsize
    M_OM = 4*v.pi/3*EoS.eps_OM_P(P_OM)*r**3
    M_DM = 4*v.pi/3*EoS.eps_DM_P(P_DM)*r**3
    OM_fin = 1
    DM_fin = 1
    R_OM = 0
    R_DM = 0
    r_safe = 0
    while( P_OM > 0 or P_DM > 0 ):
        r_safe = r
        r, (P_OM, P_DM,M_OM, M_DM), stepsize = integrate.adapt_stepsize_step(r=r, stepsize=stepsize, funcs=(func.d_P_r_OM, func.d_P_r_DM, func.d_M_r_P_OM, func.d_M_r_P_DM), vals=(P_OM,P_DM,M_OM,M_DM), limit=1e-3)
        if P_OM < almost_zero:
            P_OM = 0
        if P_DM < almost_zero:
            P_DM = 0
        if stepsize < almost_zero or stepsize > 1/almost_zero:
            P_OM = 0
            P_DM = 0        
        # if r != r_safe:
        if True:
            f = open("results/single_P_2.out","a")
            mu_OM = EoS.mu_OM_P(P_OM)
            eps_OM = EoS.eps_OM_P(P_OM)
            c_s_OM = EoS.c_s_OM_P(P_OM)
            mu_DM = EoS.mu_DM_P(P_DM)
            eps_DM = EoS.eps_DM_P(P_DM)
            c_s_DM = EoS.c_s_DM_P(P_DM)
            M = M_OM + M_DM
            f.write("%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n"%(r, stepsize, M_OM, M_DM, mu_OM, mu_DM, P_OM, P_DM, eps_OM, eps_DM, c_s_OM, c_s_DM))
            f.close()


def M_R_1_fluid(var_start, var_end, steps, var, init_stepsize=1e-10):
    f = open("results/M_R_%s_1_fluid.out"%var,"w")
    f.close()
    for i in range(steps):
        init_var = var_start*(var_end/var_start)**(i/(steps-1))
        print(init_var)
        if var == "P":
            (init_var,R_fin,M_fin) = TOV.TOV_P(init_var, init_stepsize)
        if var == "mu":
            (init_var,R_fin,M_fin) = TOV.TOV_mu(init_var, init_stepsize)
        f = open("results/M_R_%s_1_fluid.out"%var,"a")
        f.write("%e\t%e\t%e\n"%(init_var,R_fin,M_fin))
        f.close()

def M_R_TIDAL_1_fluid(var_start, var_end, steps, var, init_stepsize=1e-10):
    f = open("results/M_R_%s_1_fluid.out"%var,"w")
    f.close()
    for i in range(steps):
        init_var = var_start*(var_end/var_start)**(i/(steps-1))
        print(init_var)
        if var == "P":
            (init_var,R_fin,M_fin,y) = TOV.TOV_TIDAL_P(init_var, init_stepsize)
        if var == "mu":
            (init_var,R_fin,M_fin,y) = TOV.TOV_TIDAL_mu(init_var, init_stepsize)
        f = open("results/M_R_%s_1_fluid.out"%var,"a")
        C = M_fin/R_fin
        k_2 = func.k_2_y(y, C)
        f.write("%e\t%e\t%e\t%e\t%e\t%e\n"%(init_var,R_fin,M_fin,y,C, k_2))
        f.close()

def M_R_2_fluid(var_start, var_end, steps, dm_om_ratio, var, init_stepsize=1e-10):
    f = open("results/M_R_%s_dm_om_%1.1e_2_fluid.out"%(var,dm_om_ratio),"w")
    f.close()
    if steps == 1:
        init_var = var_start
        print(init_var)
        if var == "P":    
            R_OM, R_DM , M_OM, M_DM = TOV.TOV_2_fluid_P(init_var, init_stepsize, dm_om_ratio)
        if var == "mu":
            R_OM, R_DM , M_OM, M_DM = TOV.TOV_2_fluid_mu(init_var, init_stepsize, dm_om_ratio)
        f = open("results/M_R_%s_dm_om_%1.1e_2_fluid.out"%(var,dm_om_ratio),"a")
        f.write("%e\t%e\t%e\t%e\t%e\n"%(init_var,R_OM,R_DM,M_OM,M_DM))
        f.close()
    else:
        for i in range(steps):
            init_var = var_start*(var_end/var_start)**(i/(steps-1))
            print(init_var)
            if var == "P":    
                R_OM, R_DM , M_OM, M_DM = TOV.TOV_2_fluid_P(init_var, init_stepsize, dm_om_ratio)
            if var == "mu":
                R_OM, R_DM , M_OM, M_DM = TOV.TOV_2_fluid_mu(init_var, init_stepsize, dm_om_ratio)
            f = open("results/M_R_%s_dm_om_%1.1e_2_fluid.out"%(var,dm_om_ratio),"a")
            f.write("%e\t%e\t%e\t%e\t%e\n"%(init_var,R_OM,R_DM,M_OM,M_DM))
            f.close()

def M_R_TIDAL_2_fluid(var_start, var_end, steps, dm_om_ratio, var, init_stepsize=1e-10):
    f = open("results/M_R_TIDAL_%s_dm_om_%1.3e_M_%1.3e_2_fluid.out"%(var,dm_om_ratio, v.mass_DM),"w")
    f.close()
    for i in range(steps):
        init_var = var_start*(var_end/var_start)**(i/(steps-1))
        print("%i/%i steps! init var: %e, init_var with units: %e"%(i+1, steps, init_var, init_var*v.conv_P))
        if var == "P":    
            R_OM, R_DM , M_OM, M_DM,y = TOV.TOV_TIDAL_2_fluid_P(init_var, init_stepsize, dm_om_ratio)
        if var == "mu":
            R_OM, R_DM , M_OM, M_DM,y = TOV.TOV_TIDAL_2_fluid_mu(init_var, init_stepsize, dm_om_ratio)
        f = open("results/M_R_TIDAL_%s_dm_om_%1.3e_M_%1.3e_2_fluid.out"%(var,dm_om_ratio, v.mass_DM),"a")
        C = (M_OM+M_DM)/max(R_OM,R_DM, 1e-100)
        k_2 = func.k_2_y(y, C)
        f.write("%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n"%(init_var,R_OM,R_DM,M_OM,M_DM,y,C,k_2))
        f.close()

def write_stability_file(filename_template, pref = ""):
    print(filename_template)
    data = np.transpose(np.genfromtxt("results"+"%s"%pref+"/M_R_TIDAL_"+filename_template+".out"))
    P_0 = []
    P_1 = []
    P_0.append(0)
    P_1.append(0)
    stable = False
    for i in range(1, len(data[0])-1):
        if (data[3][i] <= data[3][i+1]) and (data[4][i] <= data[4][i+1]):
            if not stable:
                stable = True
                P_0.append(data[0][i])
        elif (data[3][i] >= data[3][i+1]) or (data[4][i] >= data[4][i+1]):
            if stable:
                stable = False
                P_1.append(data[0][i])
        if stable:
            if i == len(data[0])-2:
                P_1.append(data[0][i])
    with open("results"+"%s"%pref+"/stability_"+filename_template, "w") as f:
        for i in range(len(P_0)):
            f.write("%e\t%e\n"%(P_0[i],P_1[i]))
