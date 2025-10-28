import TOV
import sys

MeV_fm = 197.32698
def con_pressure(mass_DM):
    return mass_DM**4/MeV_fm**3

m_neutron = 939.5654205
m_DM = float(sys.argv[1])

if sys.argv[3][10:] == "I":
    print(sys.argv[3], m_DM)
    p_low_OM = 5.030643560643215e-6*m_neutron**4/m_DM**4
    p_high_OM = 0.01859823271207*m_neutron**4/m_DM**4
elif sys.argv[3][10:] == "II":
    print(sys.argv[3], m_DM)
    p_low_OM = 6.367580732002786e-6*m_neutron**4/m_DM**4
    p_high_OM = 0.0171003357007*m_neutron**4/m_DM**4
elif sys.argv[3][10:] == "III":
    print(sys.argv[3], m_DM)
    p_low_OM = 5.030643560643215e-6*m_neutron**4/m_DM**4
    p_high_OM = 0.0218441625408*m_neutron**4/m_DM**4
else:
    print("wrong EoS!!")
    exit()

p_low_DM = 0.000004
p_high_DM = 0.03

# TOV.M_R(p_low_OM, p_high_OM, p_low_DM, p_high_DM, pref=sys.argv[3],steps_OM=80,steps_DM=80,m_DM=m_DM, stepsize = float(sys.argv[2]))
TOV.M_R_parallel(p_low_OM, p_high_OM, p_low_DM, p_high_DM, pref=sys.argv[3],steps_OM=100,steps_DM=100,m_DM=m_DM)