import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from decimal import Decimal 
import csv
import os

"""
python script I am using to calculate the momentum flux tau
of an idealized ice floe with properties taken from the literature.
I will calculate the momentum flux as a function of the form drag for
sails and keels, floe edges, and the atmospheric and oceanic skin drag.
"""

#want to save images into diagrams directory
os.chdir('./diagrams')

#Defining constants 
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
#roughness length of atmosphere (given by TSAMADOS paper)
z_0a = 0.16e-3 
#roughness length of ocean (given by TSAMADOS paper)
z_0w = 0.0061
#Ice concentration (fraction of 1.0)
A = 0.5 
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
#lambda for atmosphere
my_lambda = np.log(z_0a/10)
#latent heat of vapourization of water (J/kg) source: https://link.springer.com/referenceworkentry/10.1007%2F978-90-481-2642-2_327
yotta = 2260e3

#calculating pressure at roughness length and 10m reference height...
pressure = lambda h: 101325*(1 - 2.25577e-5*h)**(5.2558)
P_z0 = pressure(z_0a)
P_10 = pressure(10)
print("pressure at z_0 is {} Pa".format(P_z0))
print("pressure at 10m is {} Pa".format(P_10))

#calculating potential temperature at roughness length and 10m reference height...
potential_temp = lambda P,T : T*(1e5/P)**(0.286)
theta_z0 = potential_temp(P_z0,T_z0)
theta_10 = potential_temp(P_10,T_10)
print("potential temperature at z0 is {} K".format(theta_z0))
print("potential temperature at 10m is {} K".format(theta_10))
print("potential temperature ratio at 10m is {}".format(theta_z0/theta_10))

#calculating specific humidity at roughness length and 10m reference height... 
#need to calculate saturated vapor pressure
e_s = lambda T: 611*np.exp((17.67*(T-T_0))/(T-29.65))
e_s_z0 = e_s(T_z0)
e_s_10 = e_s(T_10)

#calculating w. Then h_s = w/(1+w)
mass_mixing = lambda P, e_s: (e_s*gas_d)/(gas_v*(P - e_s))
w_z0 = mass_mixing(P_z0,e_s_z0)
w_10 = mass_mixing(P_10,e_s_10)
print("mass mixing at z0 is {}".format(w_z0))
print("mass mixing at 10m is {}".format(w_10))
h_s_z0 = w_z0/(1 + w_z0)
h_s_10 = w_10/(1 + w_10)
print("specific humidity at z0 is {}".format(h_s_z0))
print("specific humidity at 10m is {}".format(h_s_10))
print("specific humidity ratio at 10m is {}".format(h_s_z0/h_s_10))

#now calculating stability
#gamma = (gamma_c/(C_d**2))*(s1*C_s + s2*C_l)
#constant out the front
gamma_c = (kappa*g*z_0a)/(u_a**2)
s_1 = (theta_10 - theta_z0)/(theta_z0*(1 + 0.606*h_s_z0))
s_2 = (h_s_10 - h_s_z0)/(1/0.606 + h_s_z0)
print("s_1 is {}, s_2 is {}".format(s_1,s_2))
print("gamma_c is {}".format(gamma_c))

def stability_fun(c_d, c_s, c_l):
    #given transfer coefficients calculate stability
    return (gamma_c/(c_d**2))*(c_s*s_1 + c_l*s_2)

#get sail height from sail spacing
def H_s_fun(D_s):
    return D_s/(2.0*R_f*(W_s/math.tan(alpha_r) + (W_k/math.tan(alpha_k))*(R_h/R_d)))

#Keel height(m)
def H_k_fun(D_k):
    return R_h*H_s

#get sail spacing from sail height
def D_s_fun(H_s):
    return 2.0*H_s*R_f*(W_s/math.tan(alpha_r) + (W_k/math.tan(alpha_k))*(R_h/R_d))

#Distance between keels(m) from keel height. Use D_k = (R_d/R_h)*alpha* H_s = R_d * alpha * H_k, where alpha is given by H_s = alpha *  D_s
def D_k_fun(H_k):
    return (R_d)*2.0*H_k*R_f*(W_s/math.tan(alpha_r) + (W_k/math.tan(alpha_k))*(R_h/R_d))

#Keel height(m)
def H_k_fun(H_s):
    return R_h*H_s

#Distance between keels(m)
def D_k_fun(D_s):
    return R_d*D_s

#Now defining individual drag components
#first sheltering function
def S_c_fun(D,H):
    return abs(1.0 - np.exp(-s_l*D/H))

#Now the form drag coefficient from sails
def C_dar_fun(H_s,D_s,heat=False):
    return 0.5*c_ra*S_c_fun(D_s,H_s)*(H_s/D_s)*A*pow((np.log(H_s/z_0a)/np.log(10/z_0a)),2.0)

#Form drag coefficient from keels
def C_dwr_fun(H_k,D_k,heat=False):
    return 0.5*c_kw*S_c_fun(D_k,H_k)*(H_k/D_k)*A*pow((np.log(H_k/z_0w)/np.log(10/z_0w)),2.0)

#Skin drag coefficient for atmosphere
def C_das_fun(H_s,D_s):
    return A*(1-m_a*(H_s/D_s))*(c_sa)

#Skin drag coefficient for ocean
def C_dws_fun(H_k,D_k):
    return A*(1-m_w*(H_k/D_k))*(c_sw)

#momentum flux from atmosphere to ice
def tau_atmosphere_fun(C_d):
    #tau = C_d * rho * u^2
    return C_d*rho_a*(u_a**2)

#momentum flux from ocean to ice
def tau_ocean_fun(C_d):
    #tau = C_d * rho * u^2
    return 0.5*C_d*rho_w*(u_w**2)

#heat flux from ocean to ice 
def F_ocean(C_d):
    # 1 J/s = 1 W so no unit conversion needed...
    #need momentum
    mag_tau = np.abs(tau_ocean_fun(C_d))
    return -C_d*cp*rho_w*delta_T*math.sqrt(mag_tau/rho_w)

#sensible heat transfer(atmospheric)
def F_atmosphere_sensible(c):
    return rho_w*cp*c*(theta_10-theta_z0)*u_w

def F_atmosphere_latent(c):
    return rho_w*yotta*c*(h_s_10 - h_s_z0)*u_w

def F_atmosphere_total(c):
    return F_atmosphere_sensible(c) + F_atmosphere_latent(c)

#get compact form of stability
def chi_fun(stability):
    return (1 - 16*stability)**(0.25)

#get psi's (integrated flux profiles) from stability 
def psi_fun(stability, chi):
    if stability<0:
        psi_m = 2*np.log(0.5*(1 + chi)) + np.log(0.5*(1 + chi**2)) - 2*np.arctan(chi) + (math.pi/2)
        psi_s = 2*np.log(0.5*(1 + chi**2)) 
    else:
        psi_m = -(0.7*stability + 0.75*(stability - 14.3)*np.exp(-0.35*stability) + 10.7)
        psi_s = psi_m        
    return psi_m, psi_s 

def iterdrag(c_d):
    c_d = c_d
    c_s = c_d
    c_l = c_d
    for i in range(5):
        stability = stability_fun(c_d,c_s,c_l)
        chi = chi_fun(stability)
        psi_m, psi_s = psi_fun(stability,chi)
        #now updating coefficients based on psi
        print("iter {}, cd is {}".format(i,c_d))
        c_d = c_d/(1 + (c_d/kappa)*(my_lambda - psi_m))
        c_s = c_s/(1 + (c_s/kappa)*(my_lambda - psi_s))
        c_l = c_s
    return c_d, c_s, c_l

#DEFAULT EXCHANGE COEFFICIENTS
C_d_ocean_momentum = 0.006 
C_d_ocean_heat = 0.006
C_d_a_init = 1.31e-3
print("DEFAULT C_D_A: {}".format(C_d_a_init))
C_d_atmosphere,C_s_atmosphere,C_l_atmosphere = iterdrag(C_d_a_init)
print("Iterated C_D_A: {}".format(C_d_atmosphere))
#for non form drag scheme heat transfer, roughness length is 20% of momentum case
z_0a = 0.2*z_0a
foo,C_s_atmosphere,baz = iterdrag(C_d_a_init)
z_0a = 5*z_0a

#DEFAULT MOMENTUM AND HEAT TRANSFER
tau_ocean_default = tau_ocean_fun(C_d_ocean_momentum)
tau_atmosphere_default = tau_atmosphere_fun(C_d_atmosphere)
heat_ocean_default = F_ocean(C_d_ocean_heat)
heat_atmosphere_default = F_atmosphere_total(C_d_atmosphere)

#convert to nice string to put in plot title
C_d_atmosphere_default = '{:.2e}'.format(C_d_atmosphere)
C_s_atmosphere_default = '{:.2e}'.format(C_s_atmosphere)
tau_ocean_default = '{:.2e}'.format(tau_ocean_default)
tau_atmosphere_default = '{:.2e}'.format(tau_atmosphere_default)
heat_ocean_default = '{:.2e}'.format(heat_ocean_default)
heat_atmosphere_default = '{:.2e}'.format(heat_atmosphere_default)

#vectorizing all functions
tau_atmosphere_fun = np.vectorize(tau_atmosphere_fun)
tau_ocean_fun = np.vectorize(tau_ocean_fun)
F_ocean = np.vectorize(F_ocean)
psi_fun = np.vectorize(psi_fun)
iterdrag = np.vectorize(iterdrag)

def plot_2d(H_s,D_s):
    #function to plot dependence of drag coefficents on ridge height and distance between ridges. Because the distance between the ridges depends on the ridge height, when I change the distance between ridges and hold the ridge height constant I am effectively changing the geometry of the ridges, making them flatter/sharper. 
    
    #need to hold these constant whilst calculating skin drag
    #Sail height(m) (Worby (2008)) 
    H_s_constant = 0.57 
    #height of keel
    H_k_constant = 4.4*H_s_constant
    
    #FIRST DOING ATMOSPHERIC CALCULATIONS
    #initializing coefficient
    formsail = np.zeros(len(H_s))
    skinsail = np.zeros(len(H_s))
    #making form drag
    for (index,h) in enumerate(H_s):
        formsail[index] = C_dar_fun(h,D_s_fun(h)) #we want the only part of D_s changing to be H_s, thus put H_s into D_s_fun
    #making skin drag
    for (index,d) in enumerate(D_s):
        skinsail[index] = C_das_fun(H_s_constant,d)
    #total neutral drag coefficient
    totalsail = formsail + skinsail
    
    #now ITERATE
    c_d_a, c_s_a, c_l_a = iterdrag(totalsail)
    tau_sail = [tau_atmosphere_fun(i) for i in c_d_a]
    f_atmosphere_sensible = [F_atmosphere_sensible(i) for i in c_s_a]
    f_atmosphere_latent = [F_atmosphere_latent(i) for i in c_l_a]
    f_atmosphere = np.add(f_atmosphere_sensible,f_atmosphere_latent)

    plt.title("How atmospheric drag coeff changes w sail height. {} def".format(C_d_atmosphere_default))
    plt.plot(H_s,c_d_a,label="Sails (Atmospheric)")
    plt.xlabel("Sail height (m)")
    plt.ylabel("Force (N)")
    plt.savefig('atmosphere_drag_2d.png')
    plt.close('all')

    plt.title("How atmospheric heat transfer changes with sail height. {} def".format(C_s_atmosphere_default))
    plt.plot(H_s,c_s_a,label="Sails (Atmospheric)")
    plt.xlabel("Sail height (m)")
    plt.ylabel("Heat transfer (W/m^3)")
    plt.savefig('atmosphere_heat_2d.png')
    plt.close('all')

    #NOW DOING OCEANIC CALCULATIONS
    #MOMENTUM
    H_k = [H_k_fun(h) for h in H_s]
    D_k = [D_k_fun(d) for d in D_s]
    formkeel = np.zeros(len(H_s))
    skinkeel = np.zeros(len(H_s))
    #making form drag
    for (index,h) in enumerate(H_k):
        formkeel[index] = C_dwr_fun(h,D_k_fun(h)) #we want the only part of D_k changing to be H_k, thus put H_k into D_k_fun
    #making skin drag
    for (index,d) in enumerate(D_k):
        skinkeel[index] = C_dws_fun(H_k_constant,d)
    totalkeel = formkeel+skinkeel
    tau_keel = [tau_ocean_fun(i) for i in totalkeel]
    oceanheatflux = np.zeros(len(tau_keel))
  
    #plotting
    plt.title("How ocean momentum transfer changes with keel depth. {} def".format(tau_ocean_default))
    plt.plot(H_k,tau_keel,label="Keels (Oceanic)")
    plt.xlabel("Keel depth (m)")
    plt.ylabel("Momentum transfer (N)")
    plt.savefig('ocean_drag_2d.png')
    plt.close('all')

    #HEAT
    #calculating heat flux
    heat_flux_ocean = [F_ocean(i) for i in totalkeel]

    #plotting heat flux
    plt.title("How ocean heat flux changes with keel depth. {} def".format(heat_ocean_default))
    plt.plot(H_k, heat_flux_ocean, label="Keels (Oceanic)")
    plt.xlabel("Keel depth(m)")
    plt.ylabel("Heat transfer (W/m^3)")
    plt.savefig('ocean_heat_2d.png')
    plt.close('all')

    for (index,c_d) in enumerate(np.linspace(0,0.006,len(oceanheatflux))):
        oceanheatflux[index] = F_ocean(c_d)
    
    #plotting heat flux
    plt.plot(np.linspace(0,0.006,len(oceanheatflux)), oceanheatflux)
    plt.title("Heat flux from ocean to ice as a function of C_d_ocean.")
    plt.xlabel("C_d")
    plt.ylabel("F_ocn (W/m^3)")
    plt.savefig('ocean_heat_2d_cd.png')
    plt.close('all')

    #Store momentum and heat transfer in CSV file
    with open('heat_momentum.csv', mode='w') as file:
        file_writer = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        file_writer.writerow(["variable", "ocean or atmosphere?", "min", "max", "units"])
        file_writer.writerow(['tau', 'atmosphere', np.min(tau_sail), np.max(tau_sail),'N/m^2'])
        file_writer.writerow(['tau', 'ocean', np.min(tau_keel), np.max(tau_keel),'N/m^2'])
        file_writer.writerow(['f_ocn', 'ocean', np.min(oceanheatflux), np.max(oceanheatflux), 'W/m^2'])
        file_writer.writerow(['c_d', 'atmosphere', np.min(c_d_a), np.max(c_d_a), '-'])
        file_writer.writerow(['c_s', 'atmosphere', np.min(c_s_a), np.max(c_s_a), '-'])
    #stuff to save plots if needed later
    """
    plt.legend(loc="center right")
    if skin:
        fig.savefig("/home/bcravens/Desktop/honours2019/thesis/2D_plot_skin.png")
    elif form:
        fig.savefig("/home/bcravens/Desktop/honours2019/thesis/2D_plot_form.png")
    """

def plot_3d(H_s,D_s):
    #function to plot 3d surface of dependence of drag coefficient on height of sails and distance between sails
    H_k_temp = [H_k_fun(i) for i in H_s]
    D_k_temp = [D_k_fun(j) for j in D_s]
    #now making 3d mesh
    x_s,y_s = np.meshgrid(H_s,D_s)
    x_k,y_k = np.meshgrid(H_k_temp,D_k_temp)
    
    #calculating coefficients (form and skin)
    totalsailform = C_dar_fun(x_s,y_s)
    totalkeelform = C_dwr_fun(x_k,y_k)
    totalsailskin = C_das_fun(x_s,y_s)
    totalkeelskin = C_dws_fun(x_k,y_k)
    totalkeel = totalkeelform + totalkeelskin
    totalsail = totalsailform + totalsailskin

    #iterating
    c_d_a, c_s_a, c_l_a = iterdrag(totalsail) 

    #calculating momentum (just scalar multiple of total neutral coefficients)
    tau_atmosphere = [tau_atmosphere_fun(i) for i in c_d_a]
    tau_ocean = [tau_ocean_fun(i) for i in totalkeel]
    f_ocean = [F_ocean(i) for i in totalkeel]
    f_atmosphere_sensible = [F_atmosphere_sensible(i) for i in c_s_a]
    f_atmosphere_latent = [F_atmosphere_latent(i) for i in c_l_a]
    f_atmosphere_total = np.add(f_atmosphere_sensible,f_atmosphere_latent)
    #max newtons
    tau_max = np.max(tau_atmosphere)
    tau_max1 = np.max(tau_ocean)


    ###PLOTTING ATMOSPHERIC MOMENTUM CO-EFF### 
    cp = plt.contourf(x_s,y_s,tau_atmosphere,levels=np.linspace(0,np.max(tau_atmosphere),100))
    plt.colorbar(cp)
    cp.set_clim(0,np.max(tau_atmosphere))
    plt.xlabel("Height of sails(m)")
    plt.ylabel("Distance between sails(m)")
    plt.title("Atmospheric drag coefficient. {} default".format(tau_atmosphere_default))
    plt.savefig('atmosphere_momentum_3d.png')
    plt.close('all')
    
    ###PLOTTING OCEANIC MOMENTUM###
    momentum_arg = totalkeel
    cp = plt.contourf(x_k,y_k,momentum_arg,levels=np.linspace(0,np.max(momentum_arg),100))
    plt.colorbar(cp)
    cp.set_clim(0,np.max(momentum_arg))
    plt.xlabel("Height of keels(m)")
    plt.ylabel("Distance between keels(m)")
    plt.title("Oceanic drag coefficient. 0.006 default")
    plt.savefig('ocean_momentum_3d.png')
    plt.close('all')

    ###OCEANIC HEAT TRANSFER
    heat_arg = totalkeel
    cp = plt.contourf(x_k,y_k,totalkeel,levels=np.linspace(0,np.max(heat_arg),100),cmap='inferno')
    plt.colorbar(cp)
    cp.set_clim(0,np.max(heat_arg))
    plt.xlabel("Height of keels(m)")
    plt.ylabel("Distance between keels(m)")
    plt.title("Ocean to ice heat transfer coefficient. 0.006 default.")
    plt.savefig('ocean_heat_3d.png')
    plt.close('all')
    
    ###PLOTTING ATMOSPHERIC HEAT TRANSFER###
    cp = plt.contourf(x_s,y_s,c_s_a,levels=np.linspace(0,np.max(np.abs(c_s_a))),cmap='inferno')
    plt.colorbar(cp)
    cp.set_clim(0,np.max(c_s_a))
    plt.xlabel("Height of sails(m)")
    plt.ylabel("Distance between sails(m)")
    plt.title("Atmospheric heat transfer coefficient. {} default".format(C_s_atmosphere_default))
    plt.savefig('atmosphere_heat_3d.png')
    plt.close('all')

#a normal range for sail height is something like 0.3-4.0 m, do a bit more than this (source: Tsamados 3 a (1)) below eq 25
diff_h = 4.0-0.3 
diff_d = 500-30

plot_2d(np.linspace(0.3*0.9,4.0*1.1,50),np.linspace(30.0*0.9,500.0*1.1,50))
plot_3d(np.linspace(0.3*0.9,4.0*1.1,50),np.linspace(30.0*0.9,500.0*1.1,50))
