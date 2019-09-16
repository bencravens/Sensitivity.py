import numpy as np

# Defining constants 
#gravity at the poles (m s^(-2))
g = 9.832
#Ratio of keel depth and sail height (Worby (2008))
R_h = 4.4
#Ratio of average distance between keels and average distance between sails
R_d = 1.0
#Weight variables of sails and keels (how many sails vs how many keels)
#alpha in Tsamados Here naming W_s for weighting of sails
W_s = 0.5 
#beta in Tsamados. Here naming W_k for weighting of keels
W_k = 0.5 
#Slope of sails(rad) (Worby (2008))
alpha_r = 0.45
#Slope of keels(rad) (Worby (2008)) 
alpha_k = 0.45 
#Attenuation parameter in sheltering function (given by Tsamados)
s_l = 0.18 
#Ice concentration (fraction of 1.0)
#Get input from user
A = float(input("Input ice concentration: "))
#Sail height(m) (Worby (2008)) 
H_s = 0.57 
#height of keel
H_k = 4.4*H_s
#Ratio of aice to ardg, Ridged ice area fraction
R_f = 1.0/0.11
#dimensionless scaling coefficent for sails
c_ra = 0.1
#dimensionless scaling coefficent for keels
c_kw = 0.3
#dimensionless scaling coefficent for skin drag(atmosphere)
c_sa = 0.0005
#dimensionless scaling coefficent for skin drag(ocean)
c_sw = 0.002
#atmospheric skin drag tunable parameter (from sail height) (Tsamados)
m_a = 20.0
#oceanic skin drag tunable parameter (from keel depth) (Tsamados)
m_w = 10.0
#density of air (kg/m^3) (Wikipedia) (STP)
rho_a = 1.225
#density of ocean water (kg/m^3) (Wikipedia) (STP)
rho_w = 1025
#von karman constant
kappa = 0.40
#specific heat capacity of salt water (J/ kg K)
cp = 4000
#speed (m/s) of the atmosphere flowing above the ice at 10m (reference height used in Tsamados)
#source: Rodrigo et al 
u_a = 6.0
#speed (m/s) of the ocean flowing below the ice at 10m (reference height used in Tsamados)
#source: NEMO output
u_w = 0.03 
#delta_T is the temperature difference between the sea ice bottom and the surface of the sea. Source: CICE output
delta_T = -0.75
#Enthalpy of vapourization for water. J/kg. Source: https://en.wikipedia.org/wiki/Enthalpy_of_vaporization#Other_common_substances
L_v = 2257e3
#specific gas constant for dry air (J/(kg K)). Source: https://en.wikipedia.org/wiki/Gas_constant
gas_d = 287.058 
#specific gas constant for water vapour (J/(kg L)) Source: https://www.engineeringtoolbox.com/individual-universal-gas-constant-d_588.html
gas_v = 461.52
#Reference temperature (K)
T_0 = 273.15
#mean air temperature at ross sea (K) in February. Source: https://web.archive.org/web/20101009114203/http://coolantarctica.com/Antarctica fact file/antarctica environment/climate_graph/vostok_south_pole_mcmurdo.htm
T_z0 = 263.65
#assume wet lapse rate of 5 deg
T_10 = T_z0 - 0.05
#relative humidity (percentage, dimensionless.) Source: https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2000JC000411
p_r = 1.0 
#latent heat of vapourization of water (J/kg) source: https://link.springer.com/referenceworkentry/10.1007%2F978-90-481-2642-2_327
yotta = 2260e3
#roughness lengths
z_0a = 0.16e-3
z_0w = 0.0061
my_lambda = np.log(z_0a/10)
