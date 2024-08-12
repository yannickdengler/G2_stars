import integration as integrate
import functions as func
import EoS
import global_vars as v

# almost_zero = 1e-70
almost_zero = 1e-20

################################### print EoS 2 fluid

def print_EoS_mu_2_fluid(mu_start, mu_end, steps):
    f = open("results/EoS_mu_2_fluid.out","w")
    f.close()
    for i in range(steps):
        mu = mu_start*(mu_end/mu_start)**(i/(steps-1))
        f = open("results/EoS_mu_2_fluid.out","a")
        f.write("%e\t%e\t%e\t%e\t%e\t%e\t%e\n"%(mu,EoS.P_OM_mu(mu), EoS.eps_OM_mu(mu), EoS.c_s_OM_mu(mu),EoS.P_DM_mu(mu), EoS.eps_DM_mu(mu), EoS.c_s_DM_mu(mu)))
        f.close()

def print_EoS_P_2_fluid(P_start, P_end, steps):
    f = open("results/EoS_P_2_fluid.out","w")
    f.close()
    for i in range(steps):
        P = P_start*(P_end/P_start)**(i/(steps-1))
        f = open("results/EoS_P_2_fluid.out","a")
        f.write("%e\t%e\t%e\t%e\t%e\t%e\t%e\n"%(mu,EoS.mu_OM_P(P), EoS.eps_OM_P(P), EoS.c_s_OM_P(P),EoS.mu_DM_P(P), EoS.eps_DM_P(P), EoS.c_s_DM_P(P)))
        f.close()

################################### TOV 1 fluid

def TOV_mu(init_mu, stepsize):
    mu = init_mu
    r = stepsize
    M = 4*v.pi/3*EoS.eps_mu(mu)*r**3
    P = EoS.P_mu(mu)
    while( EoS.P_mu(mu) > 0 and mu != 0):
        r, (mu,M), stepsize = integrate.adapt_stepsize_step(r=r, stepsize=stepsize, funcs=(func.d_mu_r, func.d_M_r_mu), vals=(mu,M), limit=1e-3)
        # r, (mu,M) = integrate.Euler_step(r=r, stepsize=stepsize, funcs=(func.d_mu_r, func.d_M_r_mu), vals=(mu,M))
        # r, (mu,M) = integrate.RK_4_step(r=r, stepsize=stepsize, funcs=(func.d_mu_r, func.d_M_r_mu), vals=(mu,M))
        if EoS.P_mu(mu) < almost_zero or mu < almost_zero:
            mu = 0
    return (init_mu,r,M)

def TOV_P(init_P, stepsize):
    P = init_P
    r = stepsize
    M = 4*v.pi/3*EoS.eps_P(P)*r**3
    while( P > 0):
        r, (P,M), stepsize = integrate.adapt_stepsize_step(r=r, stepsize=stepsize, funcs=(func.d_P_r, func.d_M_r_P), vals=(P,M), limit=1e-3)
        # r, (mu,M) = integrate.Euler_step(r=r, stepsize=stepsize, funcs=(func.d_P_r, func.d_M_r_P), vals=(P,M))
        # r, (mu,M) = integrate.RK_4_step(r=r, stepsize=stepsize, funcs=(func.d_P_r, func.d_M_r_P), vals=(P,M))
        if P < almost_zero:
            P = 0
    return (init_P,r,M)

################################### TOV 2 fluid

def TOV_2_fluid_mu(init_mu, stepsize, dm_om_ratio):
    mu_OM = init_mu
    if dm_om_ratio < 0:
        mu_DM = -dm_om_ratio*init_mu
        mu_OM = 0
    else:
        mu_DM = dm_om_ratio*init_mu
    r = stepsize
    M_OM = 4*v.pi/3*EoS.eps_OM_mu(mu_OM)*r**3
    M_DM = 4*v.pi/3*EoS.eps_DM_mu(mu_DM)*r**3
    OM_fin = 0
    DM_fin = 0
    R_OM = 0
    R_DM = 0
    while( EoS.P_OM_mu(mu_OM) > 0 or EoS.P_DM_mu(mu_DM) > 0 ):
        r, (mu_OM, mu_DM,M_OM, M_DM), stepsize = integrate.adapt_stepsize_step(r=r, stepsize=stepsize, funcs=(func.d_mu_r_OM, func.d_mu_r_DM, func.d_M_r_mu_OM, func.d_M_r_mu_DM), vals=(mu_OM,mu_DM,M_OM,M_DM), limit=1e-3)
        # r, (mu_OM, mu_DM,M_OM, M_DM) = integrate.Euler_step(r=r, stepsize=stepsize, funcs=(func.d_mu_r_OM, func.d_mu_r_DM, func.d_M_r_mu_OM, func.d_M_r_mu_DM), vals=(mu_OM,mu_DM,M_OM,M_DM))
        # r, (mu_OM, mu_DM,M_OM, M_DM) = integrate.RK_4_step(r=r, stepsize=stepsize, funcs=(func.d_mu_r_OM, func.d_mu_r_DM, func.d_M_r_mu_OM, func.d_M_r_mu_DM), vals=(mu_OM,mu_DM,M_OM,M_DM))
        
        if EoS.P_OM_mu(mu_OM) < almost_zero and OM_fin == 0:
            mu_OM = 0
            R_OM = r
            OM_fin = 1
        if EoS.P_DM_mu(mu_DM) < almost_zero and DM_fin == 0:
            mu_DM = 0
            R_DM = r
            DM_fin = 1
    return (R_OM, R_DM , M_OM, M_DM)

def TOV_2_fluid_P(init_P, stepsize, dm_om_ratio):
    P_OM = init_P
    if dm_om_ratio < 0:
        P_DM = -dm_om_ratio*init_P
        P_OM = 0
    else:
        P_DM = dm_om_ratio*init_P
    r = stepsize
    M_OM = 4*v.pi/3*EoS.eps_OM_P(P_OM)*r**3
    M_DM = 4*v.pi/3*EoS.eps_DM_P(P_DM)*r**3
    OM_fin = 0
    DM_fin = 0
    R_OM = 0
    R_DM = 0
    while( P_OM > 0 or P_DM > 0 ):
        # print(stepsize*v.con_radius, P_OM*v.con_pressure, P_DM*v.con_pressure, M_OM*v.con_mass, M_DM*v.con_mass, OM_fin, DM_fin)
        r, (P_OM, P_DM,M_OM, M_DM), stepsize = integrate.adapt_stepsize_step(r=r, stepsize=stepsize, funcs=(func.d_P_r_OM, func.d_P_r_DM, func.d_M_r_P_OM, func.d_M_r_P_DM), vals=(P_OM,P_DM,M_OM,M_DM), limit=1e-2)
        # r, (P_OM, P_DM,M_OM, M_DM) = integrate.Euler_step(r=r, stepsize=stepsize, funcs=(func.d_P_r_OM, func.d_P_r_DM, func.d_M_r_P_OM, func.d_M_r_P_DM), vals=(P_OM,P_DM,M_OM,M_DM))
        # r, (P_OM, P_DM,M_OM, M_DM) = integrate.RK_4_step(r=r, stepsize=stepsize, funcs=(func.d_P_r_OM, func.d_P_r_DM, func.d_M_r_P_OM, func.d_M_r_P_DM), vals=(P_OM,P_DM,M_OM,M_DM))
        
        if stepsize < almost_zero or stepsize > 1/almost_zero:
            P_OM = 0
            P_DM = 0    
        if P_OM < almost_zero and OM_fin == 0:
            P_OM = 0
            R_OM = r
            OM_fin = 1
        if P_DM < almost_zero and DM_fin == 0:
            P_DM = 0
            R_DM = r
            DM_fin = 1 
    return (R_OM, R_DM , M_OM, M_DM)

################################### TOV_TIDAL 1 fluid

def TOV_TIDAL_mu(init_mu, stepsize):
    mu = init_mu
    y = 2
    r = stepsize
    M = 4*v.pi/3*EoS.eps_mu(mu)*r**3
    P = EoS.P_mu(mu)
    while( EoS.P_mu(mu) > 0):
        r, (mu,M,y), stepsize = integrate.adapt_stepsize_step(r=r, stepsize=stepsize, funcs=(func.d_mu_r, func.d_M_r_mu,func.d_y_r_mu), vals=(mu,M,y), limit=1e-3)
        # r, (mu,M,y) = integrate.Euler_step(r=r, stepsize=stepsize, funcs=(func.d_mu_r, func.d_M_r_mu,func.d_y_r_mu), vals=(mu,M,y))
        # r, (mu,M,y) = integrate.RK_4_step(r=r, stepsize=stepsize, funcs=(func.d_mu_r, func.d_M_r_mu,func.d_y_r_mu), vals=(mu,M,y))
        if EoS.P_mu(mu) < almost_zero:
            mu = 0
    return (init_mu,r,M,y)

def TOV_TIDAL_P(init_P, stepsize):
    P = init_P
    y = 2
    r = stepsize
    M = 4*v.pi/3*EoS.eps_P(P)*r**3
    while( P > 0):
        r, (P,M,y), stepsize = integrate.adapt_stepsize_step(r=r, stepsize=stepsize, funcs=(func.d_P_r, func.d_M_r_P,func.d_y_r_P), vals=(P,M,y), limit=1e-3)
        # r, (P,M,y) = integrate.Euler_step(r=r, stepsize=stepsize, funcs=(func.d_P_r, func.d_M_r_P,func.d_y_r_P), vals=(P,M,y))
        # r, (P,M,y) = integrate.RK_4_step(r=r, stepsize=stepsize, funcs=(func.d_P_r, func.d_M_r_P,func.d_y_r_P), vals=(P,M,y))
        if P < almost_zero:
            P = 0
    return (init_P,r,M,y)

################################### TOV_TIDAL 2 fluid asdasdasd

def TOV_TIDAL_2_fluid_mu(init_mu, stepsize, dm_om_ratio):
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
    OM_fin = 0
    DM_fin = 0
    R_OM = 0
    R_DM = 0
    # print("mu, P",mu_OM, EoS.P_OM_mu(mu_OM))
    while( EoS.P_OM_mu(mu_OM) > 0 or EoS.P_DM_mu(mu_DM) > 0 ):
        r, (mu_OM, mu_DM,M_OM, M_DM, y), stepsize = integrate.adapt_stepsize_step(r=r, stepsize=stepsize, funcs=(func.d_mu_r_OM, func.d_mu_r_DM, func.d_M_r_mu_OM, func.d_M_r_mu_DM, func.d_y_r_mu_2), vals=(mu_OM,mu_DM,M_OM,M_DM,y), limit=1e-3)
        # r, (mu_OM, mu_DM,M_OM, M_DM, y) = integrate.Euler_step(r=r, stepsize=stepsize, funcs=(func.d_mu_r_OM, func.d_mu_r_DM, func.d_M_r_mu_OM, func.d_M_r_mu_DM, func.d_y_r_mu_2), vals=(mu_OM,mu_DM,M_OM,M_DM,y))
        # r, (mu_OM, mu_DM,M_OM, M_DM, y) = integrate.RK_4_step(r=r, stepsize=stepsize, funcs=(func.d_mu_r_OM, func.d_mu_r_DM, func.d_M_r_mu_OM, func.d_M_r_mu_DM, func.d_y_r_mu_2), vals=(mu_OM,mu_DM,M_OM,M_DM,y))
        if EoS.P_OM_mu(mu_OM) < almost_zero and OM_fin == 0:
            mu_OM = 0
            R_OM = r
            OM_fin = 1
        if EoS.P_DM_mu(mu_DM) < almost_zero and DM_fin == 0:
            mu_DM = 0
            R_DM = r
            DM_fin = 1
    return (R_OM, R_DM, M_OM, M_DM, y)

def TOV_TIDAL_2_fluid_P(init_P, stepsize, dm_om_ratio):
    P_OM = init_P
    if dm_om_ratio < 0:
        P_DM = -dm_om_ratio*init_P
        P_OM = 0
    else:
        P_DM = dm_om_ratio*init_P
    y = 2
    r = stepsize
    M_OM = 4*v.pi/3*EoS.eps_OM_P(P_OM)*r**3
    M_DM = 4*v.pi/3*EoS.eps_DM_P(P_DM)*r**3
    OM_fin = 0
    DM_fin = 0
    R_OM = 0
    R_DM = 0
    while( P_OM > 0 or P_DM > 0 ):
        r, (P_OM, P_DM,M_OM, M_DM, y), stepsize = integrate.adapt_stepsize_step(r=r, stepsize=stepsize, funcs=(func.d_P_r_OM, func.d_P_r_DM, func.d_M_r_P_OM, func.d_M_r_P_DM, func.d_y_r_P_2), vals=(P_OM,P_DM,M_OM,M_DM,y), limit=1e-3)
        # r, (P_OM, P_DM,M_OM, M_DM, y) = integrate.Euler_step(r=r, stepsize=stepsize, funcs=(func.d_P_r_OM, func.d_P_r_DM, func.d_M_r_P_OM, func.d_M_r_P_DM, func.d_y_r_P_2), vals=(P_OM,P_DM,M_OM,M_DM,y))
        # r, (P_OM, P_DM,M_OM, M_DM, y) = integrate.RK_4_step(r=r, stepsize=stepsize, funcs=(func.d_P_r_OM, func.d_P_r_DM, func.d_M_r_P_OM, func.d_M_r_P_DM, func.d_y_r_P_2), vals=(P_OM,P_DM,M_OM,M_DM,y))
        
        if P_OM < almost_zero and OM_fin == 0:
            P_OM = 0
            R_OM = r
            OM_fin = 1
        if P_DM < almost_zero and DM_fin == 0:
            P_DM = 0
            R_DM = r
            DM_fin = 1
    return (R_OM, R_DM , M_OM, M_DM, y)

