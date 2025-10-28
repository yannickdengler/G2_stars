import numpy as np
import EoS
import integrate
import func
import warnings
import multiprocessing
import sys

n_cores = 8

almost_zero = 1e-50

from EoS import eps_DM as eps_DM
from EoS import n_of_P_DM as n_DM
from EoS import c_s2_DM as c_s2_DM

from EoS import eps_OM as eps_OM
from EoS import n_OM as n_OM
from EoS import c_s2_OM as c_s2_OM
from EoS import rescale_OM as rescale_OM

def print_EoS():
    f = open("results/EoS.out","w")
    Parr_OM = np.logspace(-6.5, np.log10(0.005746510987081789), 200)
    Parr_DM = np.logspace(-6.5, np.log10(0.929846126816694), 200)

    for P_OM, P_DM in zip(Parr_OM,Parr_DM):
        f.write("%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n"%(P_OM,n_OM(P_OM),eps_OM(P_OM),c_s2_OM(P_OM),P_DM,n_DM(P_DM),eps_DM(P_DM),c_s2_DM(P_DM)))
    f.close()

def single(p0_OM, p0_DM):
    f = open("results/single_TIDAL_P_2.out","w")
    stepsize = 0.03
    P_OM = p0_OM
    P_DM = p0_DM
    r = 0.0001
    y = 2
    M_OM = 4*np.pi/3*eps_OM(P_OM)*r**3
    M_DM = 4*np.pi/3*eps_DM(P_DM)*r**3
    N_OM = 4*np.pi/3*n_OM(P_OM)*r**3
    N_DM = 4*np.pi/3*n_DM(P_OM)*r**3
    OM_fin = 0
    DM_fin = 0
    R_OM = 0
    R_DM = 0
    r_safe = 0
    while(P_OM > 0 or P_DM > 0):
        if r != r_safe:
            f = open("results/single_TIDAL_P_2.out","a")
            e_OM = eps_OM(P_OM)
            c_s_OM = c_s2_OM(P_OM)
            e_DM = eps_DM(P_DM)
            c_s_DM = c_s2_DM(P_DM)
            M = M_OM + M_DM
            f.write("%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n"%(r, stepsize, M_OM, M_DM, N_OM, N_DM, P_OM, P_DM, e_OM, e_DM, c_s_OM, c_s_DM, y, func.Q_P_2_fluid(r, P_OM, P_DM, M), func.F_P_2_fluid(r, P_OM, P_DM, M),func.d_y_r_P_2(r,(P_OM,P_DM,M_OM,M_DM,N_OM,N_DM,y))))
            f.close()
        r_safe = r
        # r, (P_OM,P_DM,M_OM,M_DM,N_OM,N_DM,y), stepsize = integrate.adapt_stepsize_step(r=r, stepsize=stepsize, funcs=(func.d_P_r_OM, func.d_P_r_DM, func.d_M_r_P_OM, func.d_M_r_P_DM, func.d_N_r_OM, func.d_N_r_DM, func.d_y_r_P_2), vals=(P_OM,P_DM,M_OM,M_DM,N_OM,N_DM,y),limit=0.1)
        r, (P_OM,P_DM,M_OM,M_DM,N_OM,N_DM,y) = integrate.RK_4_step(r=r, stepsize=stepsize, funcs=(func.d_P_r_OM, func.d_P_r_DM, func.d_M_r_P_OM, func.d_M_r_P_DM, func.d_N_r_OM, func.d_N_r_DM, func.d_y_r_P_2), vals=(P_OM,P_DM,M_OM,M_DM,N_OM,N_DM,y))
        y = min([max([y,-3]),2])
        if P_OM < almost_zero and OM_fin == 0:
            P_OM = 0
            R_OM = r_safe
            OM_fin = 1
        if P_DM < almost_zero and DM_fin == 0:
            P_DM = 0
            R_DM = r_safe
            DM_fin = 1  
    return (R_OM, R_DM , M_OM, M_DM, N_OM, N_DM, y)

def TOV(p0_OM, p0_DM, stepsize=0.01):
    # stepsize = 0.005                # test 1 - 284 s
    # # stepsize = 0.03                # test 2 - 37 s
    # # stepsize = 0.01                # test 3 - 47 s
    P_OM = p0_OM
    P_DM = p0_DM
    r = 0.0001
    y = 2
    M_OM = 4*np.pi/3*eps_OM(P_OM)*r**3
    M_DM = 4*np.pi/3*eps_DM(P_DM)*r**3
    N_OM = 4*np.pi/3*n_OM(P_OM)*r**3
    N_DM = 4*np.pi/3*n_DM(P_DM)*r**3
    OM_fin = 0
    DM_fin = 0
    R_OM = 0
    R_DM = 0
    C_OM = 0
    C_DM = 0
    r_safe = 0
    while( P_OM > 0 or P_DM > 0 ):
        if y < -1 or y > 2:
            print(y)
        r_safe = r
        # r, (P_OM,P_DM,M_OM,M_DM,N_OM,N_DM,y), stepsize = integrate.adapt_stepsize_step(r=r, stepsize=stepsize, funcs=(func.d_P_r_OM, func.d_P_r_DM, func.d_M_r_P_OM, func.d_M_r_P_DM, func.d_N_r_OM, func.d_N_r_DM, func.d_y_r_P_2), vals=(P_OM,P_DM,M_OM,M_DM,N_OM,N_DM,y), limit = 0.2)
        r, (P_OM,P_DM,M_OM,M_DM,N_OM,N_DM,y) = integrate.RK_4_step(r=r, stepsize=stepsize, funcs=(func.d_P_r_OM, func.d_P_r_DM, func.d_M_r_P_OM, func.d_M_r_P_DM, func.d_N_r_OM, func.d_N_r_DM, func.d_y_r_P_2), vals=(P_OM,P_DM,M_OM,M_DM,N_OM,N_DM,y))
        if y > 2:
            y=2
        elif y < -1:
            y=-1
        if P_OM < almost_zero and OM_fin == 0:
            P_OM = 0
            R_OM = r_safe
            C_OM = (M_OM+M_DM)/R_OM
            OM_fin = 1
        if P_DM < almost_zero and DM_fin == 0:
            P_DM = 0
            R_DM = r_safe
            C_DM = (M_OM+M_DM)/R_DM
            DM_fin = 1
    return (R_OM, R_DM, M_OM, M_DM, N_OM, N_DM, y, C_OM, C_DM)

m_neutron = 939.5654205

def compute_one_pair(p0_OM, p0_DM):
    if p0_OM == 0 and p0_DM == 0:
        return "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n" % (
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    else:
        R_OM, R_DM, M_OM, M_DM, N_OM, N_DM, y, C_OM, C_DM = TOV(p0_OM, p0_DM, float(sys.argv[2]))
        C = (M_OM+M_DM)/max(R_OM,R_DM, 1e-100)
        k_2 = func.k_2_y(y, C)
        return "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n" % (
            p0_OM,p0_DM,R_OM,R_DM, max(R_OM,R_DM),M_OM,M_DM,
            M_OM+M_DM,N_OM,N_DM,y,C,k_2, C_OM, C_DM,
            eps_OM(p0_OM), eps_DM(p0_DM))

def worker(args):
    p0_OM, p0_DM = args
    return compute_one_pair(p0_OM, p0_DM,)

def M_R_parallel(p0_OM_min, p0_OM_max, p0_DM_min, p0_DM_max, steps_OM=40, steps_DM=40, m_DM = m_neutron, pref=""):
    rescale_OM(m_DM=m_DM)
    p_0_OM_arr = [0,]
    for P_0 in np.logspace(np.log10(p0_OM_min),np.log10(p0_OM_max),steps_OM):
        p_0_OM_arr.append(P_0)
    p_0_DM_arr = [0,]
    for P_0 in np.logspace(np.log10(p0_DM_min),np.log10(p0_DM_max),steps_DM):
        p_0_DM_arr.append(P_0)

    tasks = [(p0_OM, p0_DM) for p0_OM in p_0_OM_arr for p0_DM in p_0_DM_arr]

    with multiprocessing.Pool(n_cores) as pool:
        results = pool.map(worker, tasks)

    # Write everything at once
    with open("results/"+pref+"_M_R_mDM_%e.out"%m_DM, "w") as f:
        f.writelines(results)

def M_R(p0_OM_min, p0_OM_max, p0_DM_min, p0_DM_max, steps_OM=40, steps_DM=40, m_DM = m_neutron, pref="",stepsize=0.01):
    rescale_OM(m_DM=m_DM)
    p_0_OM_arr = [0,]
    for P_0 in np.logspace(np.log10(p0_OM_min),np.log10(p0_OM_max),steps_OM):
        p_0_OM_arr.append(P_0)
    p_0_DM_arr = [0,]
    for P_0 in np.logspace(np.log10(p0_DM_min),np.log10(p0_DM_max),steps_DM):
        p_0_DM_arr.append(P_0)
    # with open("results/M_R_"+pref+".out","w") as f:
    with open("results/"+pref+"_M_R_mDM_%e.out"%m_DM,"w") as f:
        for p0_OM in p_0_OM_arr:
            for p0_DM in p_0_DM_arr:    
                if p0_OM == 0 and p0_DM == 0:
                    f.write("%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n"%(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
                else:
                    R_OM, R_DM, M_OM, M_DM, N_OM, N_DM, y, C_OM, C_DM = TOV(p0_OM, p0_DM, stepsize)
                    C = (M_OM+M_DM)/max(R_OM,R_DM, 1e-100)
                    k_2 = func.k_2_y(y, C)
                    f.write("%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n"%(p0_OM,p0_DM,R_OM,R_DM, max(R_OM,R_DM),M_OM,M_DM, M_OM+M_DM,N_OM,N_DM,y,C,k_2, C_OM, C_DM, eps_OM(p0_OM), eps_DM(p0_DM)))

# if __name__ =="__main__":
    # res = single(1e-3,2e-3)<
    # print(res)
    # print(TOV(1e-3,2e-3))
    # print(func.k_2_y(res[6], (res[2]+res[3])/max([res[0],res[1]]) ))

    # print_EoS()

    # M_R(5e-6, 0.018, 0.00004, 6.757564e-03, pref="_big",steps_OM=110,steps_DM=100)    
    # M_R(5e-6, 1e-2, 0.00004, 0.007, pref="runtime_test_1",steps_OM=25,steps_DM=20)    
    # M_R(5e-6, 1e-2, 0.00004, 0.007, pref="runtime_test_1",steps_OM=25,steps_DM=20)