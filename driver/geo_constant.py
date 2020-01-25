### some constants used in atmospheric science ###
from numpy import *
# earth radius
ra=6378137.0  # meter
rb=6356752.3141 # meter
rab = ra/rb
c=2.99792E8 # m/s
# gravity constant
g0 = 9.80665  # in m s^{-2}

# gas constant for dry air
Rd=287.04 # m^2 s-2 K^-1

# ratio between  R and Cd, the heating rate with pressure fixed.
kesi=0.286 # R/Cd
# plancl constant
h=6.6260755E-34     # J s
# Boltzmann's constant
kb=1.3807E-19 # cm^3 mb K^-1 molecular-1

# Avogadro's number 
An=6.02213E23 # molecules mole^-1
# dry air mass
mg=28.966 # g mole^-1
# h2o mass
mh2o=18.02 # g mole^-1
# o3 mass
mo3=48.0   # g mole^-1
# earth angular velocity
omega=7.292115E-5 # s-1
# parameters for the geoconstant
mco2=44.0 #g mole^-1
mc=12 # g mole^-1
mco=28.0 # g mole^-1



gm=3.986005e14
j2=0.0010826256
j4=-0.0000023709122
#
kb_j=1.3807e-23 # j k-1
h_o_k=h/kb_j  # Hz-1 k
h_ok_mhz=1.0E6*h_o_k

h_o_kc=h_o_k*3.0e8 # h*c/k
E0=(1.0E-3/An)*c*c # mc^2 of H atom  in J 
# some derivated constants
hk=1.0E6*h_o_k # MHz-1 K used as h*v/kT v--> freq in MHz T in K 
hkc=100.0*h_o_kc # cm K used as h*c*E(in cm-1)/kT 
kE0=kb_j/E0
wd0=sqrt(2.0*log(2.0)*kE0) # in sqrt(K-1)
# 1/kb in nm-2 km-1 hPa-1 K

rkb=1.0/kb*(1.0/(1.E14*1.0E-5*1.0E-5))
tspace=2.735 # K

kgCO_s_2_GtC_y=mc*3600.0*365*24.0/(1.0e12*mco)
kgCO2_s_2_GtC_y=mc*3600.0*365*24.0/(1.0e12*mco2)
kg_s_2_Gt_y=3600.0*365*24.0/1.0e12

EMFACTCO2CO = 12.068


