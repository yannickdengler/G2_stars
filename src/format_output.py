import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, ticker

m_neutron = 939.5654205                                                              # in MeV
MeV_fm = 197.32698  
MeV_solar = 8.97*1e-61                                                      # 1 MeV in solar masses
Planck_Mass = 1.095*1e-38                                                   # in solar masses
Planck_Length = 1.6163e-38                                                  # in km

# con_mass = pow(Planck_Mass, 3)*pow(mass_DM*MeV_solar, -2)
# con_radius = pow(Planck_Mass, 2)*pow(mass_DM*MeV_solar, -2)*Planck_Length  
# con_pressure = mass_DM**4/MeV_fm**3
# con_num_den = mass_DM**3/MeV_fm**3
def con_mass(mass_DM):
    return pow(Planck_Mass, 3)*pow(mass_DM*MeV_solar, -2)
def con_radius(mass_DM):
    return pow(Planck_Mass, 2)*pow(mass_DM*MeV_solar, -2)*Planck_Length  
def con_pressure(mass_DM):
    return mass_DM**4/MeV_fm**3
def con_num_den(mass_DM):
    return mass_DM**3/MeV_fm**3

def get_data(pref, num_OM, num_DM,m_DM=m_neutron):                 # TODO: REMEMBER C_OM, C_DM
    data = np.transpose(np.genfromtxt("results/"+pref+"M_R_mDM_%e.out"%m_DM))
    print(pref, m_DM)

    data[:2] = data[:2]*con_pressure(m_DM)
    data[2:5] = data[2:5]*con_radius(m_DM)
    data[5:8] = data[5:8]*con_mass(m_DM)
    data[8:10] = data[8:10]*con_num_den(m_DM)
    data[15:17] = data[15:17]*con_pressure(m_DM)

    data = data.reshape((len(data),num_OM+1,num_DM+1))

    data_pure_OM = data[:,1:,0]
    data_pure_DM = data[:,0,1:]
    data = data[:,1:,1:]

    return data_pure_OM, data_pure_DM, data


def three_point_derivative(x0,y0,x1,y1,x2,y2):
    return (y1-y0)/(x1-x0)+(x1-x0)*((y2-y1)/(x2-x1)-(y1-y0)/(x1-x0))/(x2-x0)

def calc_stability_kappa(data):
    ((p0_OM_OM,p0_DM_OM,R_OM_OM,R_DM_OM,R_max_OM,M_OM_OM,M_DM_OM,M_tot_OM,N_OM_OM,N_DM_OM,y_OM,C_tot_OM,k_2_OM,C_OM_OM,C_DM_OM,eps0_OM_OM,eps0_DM_OM),(p0_OM_DM,p0_DM_DM,R_OM_DM,R_DM_DM,R_max_DM,M_OM_DM,M_DM_DM,M_tot_DM,N_OM_DM,N_DM_DM,y_DM,C_tot_DM,k_2_DM,C_OM_DM,C_DM_DM,eps0_OM_DM,eps0_DM_DM),(p0_OM,p0_DM,R_OM,R_DM,R_max,M_OM,M_DM,M_tot,N_OM,N_DM,y,C_tot,k_2,C_OM,C_DM,eps0_OM,eps0_DM)) = data
    num_OM = len(p0_OM_OM)
    num_DM = len(p0_DM_DM)

    matrix = np.zeros((num_OM-2,num_DM-2, 2,2))
    k1 = np.zeros((num_OM-2,num_DM-2))
    k2 = np.zeros((num_OM-2,num_DM-2))
    stable = np.zeros((num_OM-2,num_DM-2))

    stable_OM = np.zeros(num_OM-2)
    stable_DM = np.zeros(num_DM-2)

    for j in range(1,num_OM-1):
        if M_OM_OM[j+1] > M_OM_OM[j-1]:
            stable_OM[j-1] = int(1)
        else:
            stable_OM[j-1] = int(0)
    for i in range(1,num_DM-1):
        if M_DM_DM[i+1] > M_DM_DM[i-1]:
            stable_DM[i-1] = int(1)
        else:
            stable_DM[j-1] = int(0)

    for i in range(1,num_OM-1):
        for j in range(1,num_DM-1):
            matrix[i-1,j-1,0,0] = three_point_derivative(eps0_OM[i-1][j],N_OM[i-1][j],eps0_OM[i][j],N_OM[i][j],eps0_OM[i+1][j],N_OM[i+1][j])
            matrix[i-1,j-1,0,1] = three_point_derivative(eps0_DM[i][j-1],N_OM[i][j-1],eps0_DM[i][j],N_OM[i][j],eps0_DM[i][j+1],N_OM[i][j+1])
            matrix[i-1,j-1,1,0] = three_point_derivative(eps0_OM[i-1][j],N_DM[i-1][j],eps0_OM[i][j],N_DM[i][j],eps0_OM[i+1][j],N_DM[i+1][j])
            matrix[i-1,j-1,1,1] = three_point_derivative(eps0_DM[i][j-1],N_DM[i][j-1],eps0_DM[i][j],N_DM[i][j],eps0_DM[i][j+1],N_DM[i][j+1])
            eigval, eigvec = np.linalg.eig(matrix[i-1,j-1])
            k1[i-1,j-1], k2[i-1,j-1] = eigval
            if k1[i-1,j-1] > 0 and k2[i-1,j-1] > 0:
                stable[i-1,j-1]=int(1)
            else:
                stable[i-1,j-1]=int(0)
    return stable, stable_OM, stable_DM

def created_outfile_Axel(pref="", num_OM=10, num_DM=10,m_DM=m_neutron):
    data = get_data(pref, num_OM, num_DM,m_DM=m_DM)
    stable, stable_OM, stable_DM = calc_stability_kappa(data)
    data_OM, data_DM, data_full = data


    with open("stability_results/%sdata_m_DM_%e.dat"%(pref,m_DM),"w") as file:
        for j in range(num_OM-2):
            for k in range(num_DM-2):
                for i in range(20):
                    if i == 17:
                        file.write("%i\t"%stable[j,k])
                    elif i == 18:
                        file.write("%e\t"%m_DM)
                    elif i == 19:
                        file.write("%e\t"%(2*data_full[12,j+1,k+1]/(3*data_full[11,j+1,k+1]**5)))             # unitless Lambda
                    else:
                        file.write("%e\t"%data_full[i,j+1,k+1])
                file.write("\n")
    with open("stability_results/%spure_OM_data_m_DM_%e.dat"%(pref,m_DM),"w") as file:
        for j in range(num_OM-2):
            for i in range(20):
                if i == 17:
                    file.write("%i\t"%stable_OM[j])
                elif i == 18:
                    file.write("%e\t"%m_DM)
                elif i == 19:
                    file.write("%e\t"%(2*data_OM[12,j+1]/(3*data_OM[11,j+1]**5)))             # unitless Lambda
                else:
                    file.write("%e\t"%data_OM[i,j+1])
            file.write("\n")
    with open("stability_results/%spure_DM_data_m_DM_%e.dat"%(pref,m_DM),"w") as file:
        for k in range(num_DM-2):
            for i in range(20):
                if i == 17:
                    file.write("%i\t"%stable_DM[k])
                elif i == 18:
                    file.write("%e\t"%m_DM)
                elif i == 19:
                    file.write("%e\t"%(2*data_DM[12,k+1]/(3*data_DM[11,k+1]**5)))             # unitless Lambda
                else:
                    file.write("%e\t"%data_DM[i,k+1])
            file.write("\n")

if __name__ == "__main__":
    for m_DM in (500, 1000.000000, 2000.000000,4000):
        for pref in ["light_EoS_I_","light_EoS_II_","light_EoS_III_", "heavy_EoS_I_", "heavy_EoS_II_", "heavy_EoS_III_"]:
            created_outfile_Axel(pref=pref, num_OM=100, num_DM=100, m_DM=m_DM)