from rocketcea.cea_obj import CEA_Obj
from math import exp, pi, sqrt, tan, radians
from numpy import log
import RP1_prop as rp_prop

#Created by Akshay Kulkarni - 02/09/2019
#This script is intended to determine nozzle and combustion chamber dimensions based on specific inputs.
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
isp = ispObj.estimate_Ambient_Isp(pcc, mr, eps, pamb)

#Chamber Properties
c_properties = ispObj.get_Chamber_MolWt_gamma(pcc, mr, eps) #tuple of chamber properties
c_molweight = c_properties[0] #Molecular Weight - in lb/lbmol
c_gamma = c_properties[1] #Ratio of Specific Heats
c_temp = ispObj.get_Tcomb(pcc, mr) #Chamber Combustion Temperature (Rankine)
c_star = ispObj.get_Cstar(pcc, mr)
c_cp = ispObj.get_Chamber_Cp(pcc, mr, eps)
c_density = ispObj.get_Chamber_Density(pcc, mr, eps)

#Throat Properties
t_prop = ispObj.get_IvacCstrTc_ThtMwGam(pcc, mr, eps) #Gas Properties in the Throat
t_gamma = t_prop[4]
t_temp = (2*c_temp)/(t_gamma+1) #Eqn 3-22 from RPE, Using Throat Gamma (Rankine)
t_pressureratio = ispObj.get_Throat_PcOvPe(pcc, mr) #Chamber to Throat Pressure Ratio
t_pressure = pcc/t_pressureratio #Pressure at the throat (psi) - verified with RPE pg.57 - Pc/Pt = 0.56

#Augmented Constants
R_specific_comb = (R/c_molweight)
R_specific_comb_in = (R/(12*c_molweight)) #Universal Gas Constant for Chamber

#Nozzle Parameters
#Throat
t_area = (mass_flow/pcc)*sqrt((c_temp*R_specific_comb_in)/c_gamma)*(1+(((c_gamma-1)/2)**((c_gamma+1)/(2*(c_gamma-1))))); #Throat Area (in^2) - refer: https://www.grc.nasa.gov/www/k-12/airplane/rktthsum.html
#Assuming Chamber Temperature and Chamber Pressure are Stagnation Pressures - velocity in CC assumed to be negligible
t_radius = sqrt(t_area/pi) #Throat Radius (in)
t_diameter = 2*2.54*t_radius #Throat Diameter (cm)
t_area_cm = t_area*6.4516 #Throat Area (cm^2)

#Exit
e_pressureratio = ispObj.get_PcOvPe(pcc, mr, eps) #Ratio of Pressure from Combustion Chamber to Exit
e_pressure_actual = pcc/e_pressureratio #Actual Exit Pressure (psi)
e_area = eps*t_area #Nozzle Exit Area (in^2)
e_prop = ispObj.get_IvacCstrTc_exitMwGam(pcc, mr, eps) #Nozzle Exit Properties
e_molweight = e_prop[3] #Molecular Weight - in lb/lbmol
e_gamma = e_prop[4] #Nozzle Exit Gamma - this parameter changes significantly from the chamber gamma
e_mach = ispObj.get_MachNumber(pcc, mr, eps) #Exit Mach Number

#Augmented Constants
R_specific_exit = (R/e_molweight) #Universal Gas Constant for Nozzle Exit

#Thrust Calcs
e_vel_HH = sqrt(((2*g*e_gamma*R_specific_exit*c_temp)/(e_gamma-1))*(1- ((pamb)/(pcc))**((e_gamma-1)/(e_gamma)))) #H&H Eqn 1.18
thrust_HH = (mass_flow*(e_vel_HH))/(32.2) + (e_pressure_actual - pamb)*e_area


print ("Input Parameters")
print ("Chamber Pressure (psi) = ", pcc)
print ("Expansion Ratio = ", eps)
print ("Mass Flow Rate (lb/s) = ", mass_flow)
print ("Mass Ratio = ", mr)
print ("Ambient Pressure (psi) = ", pamb)
print ("")
print ("Calculated Parameters")
print ("Chamber Temperature (R) = ", c_temp)
print ("R for Exit Gases (ft-lbf/lbm R) = ", R_specific_exit)
print ( "Throat Area (in^2) = ", t_area )
print ("Throat Radius (in) = ", t_radius)
print ('C* (ft/s) = ', c_star)
print ("Exit Mach Number = ", e_mach)
print ("Isp = ", isp)
print ("Chamber Cp = ", c_cp)
print ("Chamber Molecular Weight = ", c_molweight)
print ("Chamber Density (lbm/ft^3) = ", c_density)
print ("")
print ("Outputs")
print ("Exit Velocity ft(s) = ", e_vel_HH)
print ("Thrust (lbf)", thrust_HH)
