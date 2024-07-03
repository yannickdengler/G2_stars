import sys

PATH = "/home/dengler_yannick/Documents/G2_stars/"
num_const = 3
lat_const = 0.343               # in fm
MeV_fm = 197.32698  
mn_light_lat = 1.63              # m_N*a (m_DM in our case) = m_lat        
neutron_mass = 939.5654205                                                  # in MeV
mass_DM = float(sys.argv[1])   # in MeV
conv_P = mass_DM**4/MeV_fm**3
pi = 3.1415926535897932384626433832795028841971693993
MeV_solar = 8.97*1e-61                                                      # 1 MeV in solar masses
Planck_Mass = 1.095*1e-38                                                   # in solar masses
Planck_Length = 1.6163e-38                                                  # in km
con_mass = pow(Planck_Mass, 3)*pow(mass_DM*MeV_solar, -2)
con_radius = pow(Planck_Mass, 2)*pow(mass_DM*MeV_solar, -2)*Planck_Length               #1.769 * 1e-76# *neutron_mass_fm**2

# rat_ferm_neutron = fermion_mass/neutron_mass
lat_sat_den_prime = 1.077566