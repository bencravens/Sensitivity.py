import math
from math import isclose
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from decimal import Decimal 
import csv
import os
from constants import *

"""
calculate the momentum and heat flux for some sea ice gridcell with properties taken from the literature. Plot as a function of the geometric properties of the sea ice in the gridcell (scalar field of ridge height and spacing)
"""
class ice:
    
    def pressure(self,h):
        return 101325*(1 - 2.25577e-5*h)**(5.2558)       

    def potential_temp(self,P,T):
        return T*(1e5/P)**(0.286)

    def specific_humidity(self,P,T):
        #need to first calculate saturated vapour pressure
        e_s = 611*np.exp((17.67*(T-T_0))/(T-29.65))
        #now calculate mass_mixing
        w = (e_s*gas_d)/(gas_v*(P - e_s))
        #now calculate specific humidity
        return w/(1 + w)    

    #get sail height from sail spacing
    def H_s_fun(self,D_s):
        return D_s/(2.0*R_f*(W_s/math.tan(alpha_r) + (W_k/math.tan(alpha_k))*(R_h/R_d)))

    #Keel height(m)
    def H_k_fun(self,D_k):
        return R_h*H_s

    #get sail spacing from sail height
    def D_s_fun(self,H_s):
        return 2.0*H_s*R_f*(W_s/math.tan(alpha_r) + (W_k/math.tan(alpha_k))*(R_h/R_d))

    #Distance between keels(m) from keel height. Use D_k = (R_d/R_h)*alpha* H_s = R_d * alpha * H_k, where alpha is given by H_s = alpha *  D_s
    def D_k_fun(self,H_k):
        return (R_d)*2.0*H_k*R_f*(W_s/math.tan(alpha_r) + (W_k/math.tan(alpha_k))*(R_h/R_d))

    #Keel height(m)
    def H_k_fun(self,H_s):
        return R_h*H_s

    #Distance between keels(m)
    def D_k_fun(self,D_s):
        return R_d*D_s

    #Now defining individual drag components
    #first sheltering function
    def S_c_fun(self,D,H):
        return abs(1.0 - np.exp(-s_l*D/H))

    #Now the form drag coefficient from sails
    def C_dar_fun(self,H_s,D_s,z_0):
        return 0.5*c_ra*self.S_c_fun(D_s,H_s)*(H_s/D_s)*A*pow((np.log(H_s/z_0)/np.log(10/z_0)),2.0)

    #Form drag coefficient from keels
    def C_dwr_fun(self,H_k,D_k,z_0):
        return 0.5*c_kw*self.S_c_fun(D_k,H_k)*(H_k/D_k)*A*pow((np.log(H_k/z_0)/np.log(10/z_0)),2.0)

    #Skin drag coefficient for atmosphere
    def C_das_fun(self,H_s,D_s):
        return A*(1-m_a*(H_s/D_s))*(c_sa)

    #Skin drag coefficient for ocean
    def C_dws_fun(self,H_k,D_k):
        return A*(1-m_w*(H_k/D_k))*(c_sw)

    #momentum flux from atmosphere to ice
    def tau_atmosphere_fun(self,C_d):
        #tau = C_d * rho * u^2
        return C_d*rho_a*(u_a**2)

    #momentum flux from ocean to ice
    def tau_ocean_fun(self,C_d):
        #tau = C_d * rho * u^2
        return 0.5*C_d*rho_w*(u_w**2)

    #heat flux from ocean to ice 
    def F_ocean(self,C_d):
        # 1 J/s = 1 W so no unit conversion needed...
        #need momentum
        mag_tau = np.abs(self.tau_ocean_fun(C_d))
        return -C_d*cp*rho_w*delta_T*math.sqrt(mag_tau/rho_w)

    #sensible heat transfer(atmospheric)
    def F_atmosphere_sensible(self,c):
        return rho_w*cp*c*(theta_10-theta_z0)*u_w

    def F_atmosphere_latent(self,c):
        return rho_w*yotta*c*(h_s_10 - h_s_z0)*u_w

    def F_atmosphere_total(self,c):
        return F_atmosphere_sensible(c) + F_atmosphere_latent(c)
  
    def stability_fun(self,c_d, c_s, c_l):
        #given transfer coefficients calculate stability
        return (self.gamma_c/(c_d**2))*(c_s*self.s_1 + c_l*self.s_2)

    #get compact form of stability
    def chi_fun(self, stability):
        #STABILITY MAGNITUDE NEEDS TO BE LESS THAN TEN
        print(stability)
        assert abs(stability)<10, "error: STABILITY IS TOO BIG!"
        chi = (1 - 16*stability)**(0.25) 
        #chi needs to be >= 1.
        print(chi)
        assert chi>=1, "error: chi is too small"
        return chi

    #get psi's (integrated flux profiles) from stability 
    def psi_fun(self, stability):
        assert abs(stability)<10, "error: STABILITY TOO BIG"
        if stability<0:
            chi = self.chi_fun(stability)
            psi_m = 2*np.log(0.5*(1 + chi)) + np.log(0.5*(1 + chi**2)) - 2*np.arctan(chi) + (math.pi/2)
            psi_s = 2*np.log(0.5*(1 + chi**2)) 
        else:
            psi_m = -(0.7*stability + 0.75*(stability - 14.3)*np.exp(-0.35*stability) + 10.7)
            psi_s = psi_m        
        return psi_m, psi_s 

    #iterate the drag coefficients 5 times. This is the overall function called 
    #when wanting to calculate drag coefficients.  
    def iterdrag(self, c_d):
        c_d = c_d
        c_s = c_d
        c_l = c_d
        for i in range(5):
            stability = self.stability_fun(c_d,c_s,c_l)
            psi_m, psi_s = self.psi_fun(stability)
            #now updating coefficients based on psi
            c_d = c_d/(1 + (c_d/kappa)*(my_lambda - psi_m))
            c_s = c_s/(1 + (c_s/kappa)*(my_lambda - psi_s))
            c_l = c_s
        return c_d, c_s, c_l
        
    def __init__(self,H_s,D_s,z_0a,z_0w,T_z0):
        #sail height
        self.H_s = H_s
        #distance between sails
        self.D_s = D_s
        #atmospheric roughness length
        self.z_0a = z_0a
        #oceanic roughness length
        self.z_0w = z_0w        
        #temperature at roughness length
        self.T_z0 = 263.65
        #temperature at 10m
        self.T_10 = T_z0 - 0.05 #assume wet lapse rate of 5 deg per km
        #pressure at roughness length
        self.P_z0 = self.pressure(z_0a) 
        #pressure at 10m
        self.P_10 = self.pressure(10)
        #potential temperature at roughness length
        self.theta_z0 = self.potential_temp(self.P_z0,self.T_z0)
        #potential temperature at 10m
        self.theta_10 = self.potential_temp(self.P_10,self.T_10)
        #specific humidity
        self.h_s_z0 = self.specific_humidity(self.P_z0,self.T_z0)
        self.h_s_10 = self.specific_humidity(self.P_10,self.T_10)
        #Height of keels
        self.H_k = self.H_k_fun(self.H_s)
        #Distance between keels
        self.D_k = self.D_k_fun(self.D_s)
        #Constants that go into stability calculation
        self.gamma_c = (kappa*g*self.z_0a)/(u_a**2)
        self.s_1 = (self.theta_10 - self.theta_z0)/(self.theta_z0*(1 + 0.606*self.h_s_z0))
        self.s_2 = (self.h_s_10 - self.h_s_z0)/(1/0.606 + self.h_s_z0) 
 
        #VECTORIZE FUNCTIONS SO THAT WE CAN PASS MATRICES TO THEM
        self.stability_fun = np.vectorize(self.stability_fun)
        self.chi_fun = np.vectorize(self.chi_fun)
        self.psi_fun = np.vectorize(self.psi_fun)
        self.iterdrag = np.vectorize(self.iterdrag)
        
        #making default drag coefficients
        C_d_atmosphere,C_s_atmosphere,C_l_atmosphere = self.iterdrag(C_d_a_init)
        #now storing them in a nice string format to put in plot titles
        self.C_d_atmosphere_default = '{:.2e}'.format(C_d_atmosphere)
        self.C_s_atmosphere_default = '{:.2e}'.format(C_s_atmosphere)

        #making 3d mesh for plotting
        self.x_s,self.y_s = np.meshgrid(self.H_s,self.D_s)
        self.x_k,self.y_k = np.meshgrid(self.H_k,self.D_k)
        
        #making scalar fields for 3d plotting
        #first momentum transfer coefficient
        self.momentum_sail = self.C_dar_fun(self.x_s,self.y_s,self.z_0a) + self.C_das_fun(self.x_s,self.y_s)
        self.momentum_keel = self.C_dwr_fun(self.x_k,self.y_k,self.z_0w) + self.C_dws_fun(self.x_k,self.y_k)
        self.momentum_sail, foo, bar = self.iterdrag(self.momentum_sail) 

        #now heat transfer coefficient
        #roughness length for heat transfer is 0.2* (source: UM documentation)
        self.heat_sail = self.C_dar_fun(self.x_s,self.y_s,0.2*self.z_0a) + self.C_das_fun(self.x_s,self.y_s)
        self.heat_keel = self.C_dwr_fun(self.x_k,self.y_k,0.2*self.z_0w) + self.C_dws_fun(self.x_k,self.y_k)
        foo, self.heat_sail, bar = self.iterdrag(self.heat_sail) 
 
    def plot_3d(self):
        #function to plot 3d surface of dependence of drag coefficient on height of sails and distance between sails
        #printing min, max, default
        print("Atmospheric momentum:")
        print("min {}\n max {}\n default {}".format(np.min(self.momentum_sail),np.max(self.momentum_sail),self.C_d_atmosphere_default))
        cp = plt.contourf(self.x_s,self.y_s,self.momentum_sail,levels=np.linspace(0,np.max(self.momentum_sail),100))
        plt.colorbar(cp)
        cp.set_clim(0,np.max(self.momentum_sail))
        plt.xlabel("Height of sails(m)")
        plt.ylabel("Distance between sails(m)")
        plt.title("Atmospheric drag coeff. {} default. {} aice".format(self.C_d_atmosphere_default,A))
        plt.savefig('atmosphere_momentum_3d_{}.png'.format(A))
        plt.close('all')
        
        ###PLOTTING OCEANIC MOMENTUM###
        momentum_arg = self.momentum_keel 
        print("Oceanic momentum:")
        print("min {}\n max {}\n default {}".format(np.min(momentum_arg),np.max(momentum_arg),0.006))
        cp = plt.contourf(self.x_k,self.y_k,momentum_arg,levels=np.linspace(0,np.max(momentum_arg),100))
        plt.colorbar(cp)
        cp.set_clim(0,np.max(momentum_arg))
        plt.xlabel("Height of keels(m)")
        plt.ylabel("Distance between keels(m)")
        plt.title("Oceanic drag coeff. 0.006 default, {} aice".format(A))
        plt.savefig('ocean_momentum_3d_{}.png'.format(A))
        plt.close('all')

        ###OCEANIC HEAT TRANSFER
        heat_arg = self.heat_keel 
        print("Oceanic heat:")
        print("min {}\n max {}\n default {}".format(np.min(heat_arg),np.max(heat_arg),0.006))
        cp = plt.contourf(self.x_k,self.y_k,heat_arg,levels=np.linspace(0,np.max(heat_arg),100),cmap='inferno')
        plt.colorbar(cp)
        cp.set_clim(0,np.max(heat_arg))
        plt.xlabel("Height of keels(m)")
        plt.ylabel("Distance between keels(m)")
        plt.title("Ocean to ice heat transfer coeff. 0.006 default. {} aice".format(A))
        plt.savefig('ocean_heat_3d_{}.png'.format(A))
        plt.close('all')
        
        ###PLOTTING ATMOSPHERIC HEAT TRANSFER###
        print("Atmospheric heat:")
        print("min {}\n max {}\n default {}".format(np.min(self.heat_sail),np.max(self.heat_sail),self.C_s_atmosphere_default))
        cp = plt.contourf(self.x_s,self.y_s,self.heat_sail,levels=np.linspace(0,np.max(np.abs(self.heat_sail))),cmap='inferno')
        plt.colorbar(cp)
        cp.set_clim(0,np.max(self.heat_sail))
        plt.xlabel("Height of sails(m)")
        plt.ylabel("Distance between sails(m)")
        plt.title("Atmospheric heat transfer coeff. {} default. {} aice".format(self.C_s_atmosphere_default,A))
        plt.savefig('atmosphere_heat_3d_{}.png'.format(A))
        plt.close('all')
    
    def display(self):
        print("pressure at z_0 is {}".format(self.P_z0))
        print("pressure at 10 is {}".format(self.P_10))
        print("theta at z_0 is {}".format(self.theta_z0))
        print("theta at 10 is {}".format(self.theta_10))
        print("specific humidity at z_0 is {}".format(self.h_s_z0))
        print("specific humidity at 10m is {}".format(self.h_s_10))

    def testfun(self):
        #this function tests the calculated properties of the sea ice to check that they are realistic...
        #basing these tests on expected values for these parameters as laid out in Tsamados paper (cited in report...)
        print("Beginning tests.")
        assert np.min(self.H_s)>0,"ERROR. Sails must have a height..."
        assert np.max(self.H_s)<5,"ERROR. Sails are too tall!"        
        assert np.min(self.D_s)>20,"ERROR. Sails are too close together!"
        assert np.max(self.D_s)<550,"ERROR. Sails are too far apart!"
        assert 1e-5 < self.z_0a < 1e-2,"ERROR. Atmospheric roughness length is too large/small"
        assert 1e-5 < self.z_0w < 1e-2,"ERROR. Oceanic roughness length is too large/small"
        assert self.T_10<self.T_z0, "ERROR.. temperature should decrease with altitude. Chech T_10 and T_z0"
        assert self.P_10<self.P_z0, "ERROR.. pressure should decrease with altitude"
        assert self.P_z0>1.013e5, "ERROR. Pressure at rougnness length is too low. Should be >1e5 (ambient pressure is 101325"
        #potential temperature literature vals taken from https://tinyurl.com/thetavals
        assert 260<self.theta_z0<265, "ERROR. potential temperature at z_0 is not valid."
        assert 260<self.theta_10<265, "ERROR. potential temperature at 10m is not valid."
        #humidity calculator results given by http://www.humcal.com/index.php
        assert isclose(1.83e-3,self.h_s_z0,abs_tol=0.05e-3), "ERROR. specific humidity at z0 is incorrect."
        assert isclose(1.83e-3,self.h_s_10,abs_tol=0.05e-3), "ERROR. specific humidity at z0 is incorrect."
        #atmospheric humidity decreases with altitude...
        assert self.h_s_z0>self.h_s_10, "ERROR: humidity should decrease with altitude!"
        assert np.mean(self.momentum_keel) > np.mean(self.momentum_sail), "Oceanic drag coefficient should be large than atmospheric (source: Tsamados)"
        assert np.mean(self.heat_keel) > np.mean(self.heat_sail), "Oceanic heat transfer should be larger than atmospheric (source: Tsamados)"
        print("All tests passed!")

my_H_s = np.linspace(0.3,4.0,50)  
my_D_s = np.linspace(30.0,500.0,50)
testice = ice(my_H_s,my_D_s,0.16e-3,0.0061,263.65)
testice.display()
testice.testfun()
testice.plot_3d()
