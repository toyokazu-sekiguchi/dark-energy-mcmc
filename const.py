import numpy as np
from scipy import special

# constants in SI units

# astronomical units
km = 1e3 # km [m]
pc = 3.085675e16 # parsec [m]
Mpc = 1e6*pc
yr = 3.1557600e7 # year [s]
Gyr = 1e9*yr

# physical constants
c = 2.99792458e8 # speed of light [m/s]
m_H = 1.6735575e-27 # hydrogen atom mass [kg]
m_He = 6.6464764e-27 # helium4 atom mass [kg]
sigmaT = 6.6524587158e-29 # Thomson cross section [m^2]
G = 6.67430e-11 # Newton constant [m^3/kg/s^2]
kB = 1.38064852e-23 # Boltzmann constant [m^2*kg/s^2/K] 
eV = 1.602176634e-19 # eV [m^2*kg/s^2]
hbar = 1.054571817e-34 # reduced Planck constant [m^2*kg/s]

# cosmological units
BigH = 100*km/Mpc # H0/h = 100km/s/Mpc in [/s]
rhoch2 = 3*(c*BigH)**2/(8*np.pi*G) # critical density devided by h^2 [kg/m^3]
TCMB = 2.7255*kB # CMB temperature [m^2*kg/s^2]
TCnuB = TCMB*(4/11)**(1/3) # CnuB temperature [m^2*kg/s^2]

# numerical values useful for massive neutrino
nu_energy = 120/7/np.pi**4
nu_number = 180*special.zeta(3)/7/np.pi**4
nnu_standard = 3.046
m2nu21 = 7.39e-5 # [eV^2]
m2nu32 = 2.45e-3 # [eV^2]

