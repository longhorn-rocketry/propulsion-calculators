from rocketcea.cea_obj import CEA_Obj
from math import exp, pi, sqrt, tan, sin, radians
from numpy import log
import RP1_prop as rp_prop
import matplotlib.pyplot as plt

#Created by Akshay Kulkarni - 02/17/2019
#This script is intended to perform combustion chamber heat transfer analysis.
#Please beware that the parameters are based on ideal assumptions and equilibrium calculations based on NASACEA outputs.

#Constants
R = 1545.32 #ft-lbf/lbmol-R
g = 32.2 #ft/s^2

#Optimization Variables
pcc = 500 #Combustion Chamber Pressure (psia)
pamb = 14.7 #Ambient Atmospheric Pressure (psia)
mr = 6.3 #oxidizer to fuel mass ratio
mass_flow = 3.63 #mass flow rate (lbm/s)

#CEA Variables
ispObj = CEA_Obj( oxName='N2O', fuelName='HTPB')
eps = ispObj.get_eps_at_PcOvPe(pcc, mr, pcc/pamb); #Optimim Expansion Ratio

#Chamber Properties
c_properties = ispObj.get_Chamber_MolWt_gamma(pcc, mr, eps) #tuple of chamber properties
c_molweight = c_properties[0] #Molecular Weight - in lb/lbmol
c_gamma = c_properties[1] #Ratio of Specific Heats
c_temp = ispObj.get_Tcomb(pcc, mr) #Chamber Combustion Temperature (Rankine)
c_cp = ispObj.get_Chamber_Cp(pcc, mr)
c_star = ispObj.get_Cstar(pcc, mr)

#Throat Properties
t_properties = ispObj.get_Throat_MolWt_gamma(pcc, mr, eps)
t_molweight = t_properties[0] #Molecular Weight - in lb/lbmol
t_gamma = t_properties[1] #Ratio of Specific Heats
t_temp = 0.91*c_temp #RPE - general correlation for Combustion Temp to Throat Temp

#Augmented Constants
R_specific_comb = (R/c_molweight)
R_specific_comb_in = (R/(12*c_molweight)) #Universal Gas Constant for Chamber

#Throat
t_area = (mass_flow/pcc)*sqrt((c_temp*R_specific_comb_in)/c_gamma)*(1+(((c_gamma-1)/2)**((c_gamma+1)/(2*(c_gamma-1))))); #Throat Area (in^2) - refer: https://www.grc.nasa.gov/www/k-12/airplane/rktthsum.html
#Assuming Chamber Temperature and Chamber Pressure are Stagnation Pressures - velocity in CC assumed to be negligible
t_radius = sqrt(t_area/pi) #Throat Radius (in)
t_diameter = 2*t_radius
t_area_cm = t_area*6.4516 #Throat Area (cm^2)

#Combustion Chamber Parameters
c_area = t_area_cm*((8.0*(t_diameter)**(-0.6)) + 1.25)*0.155 #Chamber Area (in^2) Obtained from: https://www.youtube.com/watch?v=b9NcdDhhZio
c_diameter = 2*sqrt(c_area/pi)
c_area_ratio = c_area/t_area #Chamber to Throat Area Ratio - verify with RPE Table 3-2

lc_eqn = exp((0.029*log(t_radius*2)**2) + 0.47*log(t_radius*2) + 1.94) #Chamber Cylindrical LC - Referenced from Fig 1.7: http://www.braeunig.us/space/propuls.htm
#Reference Fig 4.9 H&H
theta = radians(45) #Combustion Chamber - Throat Angle (Ranges from 25-45 degrees)
c_volume_total = t_area*((lc_eqn*c_area_ratio) + 0.333*sqrt(t_area/pi)*(1/tan(theta))*((c_area_ratio**(0.333))-1)) #H&H Eqn 4.5
c_surface_area = 2*lc_eqn*sqrt(pi*c_area_ratio*t_area) + (1/sin(theta))*(c_area_ratio-1)*t_area #H&H Eqn.4.6

#Ablative Wall Thermal Properties
w_k = 0.928 #Thermal Conductivity of ceramic wall (BTU/(hr-F-ft) - https://www.azom.com/properties.aspx?ArticleID=133
w_k_adjusted = (w_k)/(12*3600) #Thermal Conductivity of ceramic wall (BTU/(s-F-in)
w_thickness = 0.004 #Wall Thickness (in)
fg_k = 0.696 #Thermal Conductivity of fiberglass wall (BTU/(hr-F-ft) - https://www.engineeringtoolbox.com/thermal-conductivity-d_429.html
fg_k_adjusted = (fg_k)/(12*3600) #Thermal Conductivity of ceramic wall (BTU/(s-F-in)
fg_thickness = 0.35 #Wall Thickness (in)

#Aluminum Heat Thermal Properties Sink
c_amb_temp = 530.67 #Case Ambient Temperature (R)
c_thickness = 0.375 #Case Thickness (in)
c_k = 137 #Thermal Conductivity of Aluminum Case (BTU/(hr-F-ft)
c_k_adjusted = (c_k)/(12*3600) #Thermal Conductivity of Aluminum Case (BTU/(s-F-in) #

#Case Dimensions and Properties
c_case_density = 0.098 #Density of 6061 Aluminum (lb/in3)
c_postcomb_len = 4 #Post-combustion Chamber Length (in)
c_case_innerdia = 4.625 #Combustion Chamber Diameter (in)
c_case_surfarea = (pi*c_case_innerdia)*(c_postcomb_len) #Heated Surface Area (in2)
c_case_outerdia = 5 #Case Outer Diameter (in)
c_case_length = 6 #Case Length (in)
c_case_vol = pi*(((c_case_outerdia**2)-(c_case_innerdia**2))/4)*c_case_length #Case Volume (in3)
c_case_spht = 0.215 #Specific Heat Capacity of Aluminum (BTU/lb-F)

c_case_temp = c_amb_temp
flight_time = 0
time = []
temp = []
heatflux = []

for x in range(1,400000):
    dt = 0.00001
    flight_time = flight_time + dt

    #Combustion Chamber Heat Transfer
    c_local_mach = 0.2  # Chamber Mach Number

    # Combustion Chamber Heat Transfer Calcs
    c_Pr = (4 * c_gamma) / ((9 * c_gamma) - 5)  # Chamber Prandtl Number
    c_mu = (46.6 * 10 ** (-10)) * (c_molweight ** 0.5) * (c_temp ** 0.6)  # Chamber mu
    c_sigma = 1 / (((0.5 * (c_case_temp / c_temp) * (1 + ((c_gamma - 1) / (2)) * (c_local_mach ** 2)) + 0.5) ** 0.68) * ((1 + ((c_gamma - 1) / (2)) * (c_local_mach ** 2)) ** (0.12)))  # Chamber Sigma
    c_hg = (((0.026) / (t_diameter ** 0.2)) * (((c_mu ** 0.2) * c_cp) / (c_Pr ** 0.6)) * (((pcc * g) / (c_star)) ** 0.8) * ((t_diameter / (1.5 * t_radius)) ** 0.1)) * ((t_area / c_area) ** 0.9) * c_sigma  # Chamber wall heat transfer coefficient

    c_U = 1/((1/c_hg) + (w_thickness/w_k_adjusted) + (fg_thickness/fg_k_adjusted))
    c_q = c_U*c_case_surfarea*(c_temp - c_case_temp)

    c_case_temp = ((c_q*dt)/(c_case_density*c_case_vol*c_case_spht)) + c_case_temp
    time.append(flight_time)
    temp.append(c_case_temp)
    heatflux.append(c_q)

plt.plot (time, temp)
plt.grid()
plt.xlabel('Time [s]');
plt.ylabel('Temp [R]');
plt.title('Case Temperature vs Time');
plt.show()
print("")
plt.plot (time, heatflux)
plt.grid()
plt.xlabel('Time [s]');
plt.ylabel('Heat Flow [BTU/s]');
plt.title('Case Heat Flow vs Time');
plt.show()

print ("Input Parameters")
print ("Ceramic Thermal Conductivity (BTU/(hr-F-ft) = ", w_k)
print ("Ceramic Thickness (in) = ", w_thickness)
print ("Ablative Thermal Conductivity (BTU/(hr-F-ft) = ", fg_k)
print ("Ablative Thickness (in) = ", fg_thickness)
print ("Case Initial Temperature (R) = ", c_amb_temp)
print ("Case Specific Heat (BTU/(hr-F-ft) = ", c_case_spht)
print ("")
print ("Combustion Chamber Geometry")
print ("Combustion Chamber Length (in) = ", c_postcomb_len)
print ("Combustion Chamber Surface Area (in^2) = ", c_case_surfarea)
print ("Combustion Chamber Volume (in^3) = ", c_case_vol)

print("")
print ("Chamber Heat Transfer Coefficient (Btu/in^2-s-F) = ",c_hg)
print ("Combustion Chamber Heat Flux (Btu/in^2-s) = ", c_q)
print ("Combustion Chamber Temperature (R)", c_temp)
print ("Chamber Final Temperature (R) = ", c_case_temp)
print ("Chamber Diameter (in) = ", c_diameter)
print ("Chamber Length (in) = ", lc_eqn)

