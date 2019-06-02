#include <math.h>
#include "nitrous_oxide.h" /* header file */
#include <iostream>
using namespace std;
const float pCrit = 72.51f; /* critical pressure, Bar Abs */
const float rhoCrit = 452.0f; /* critical density, kg/m3 */
const float tCrit = 309.57f; /* critical temperature, Kelvin (36.42 Centigrade) */
const float ZCrit = 0.28f; /* critical compressibility factor */
const float gamma = 1.3; /* average over subcritical range */
/* ranges of function validity */
const float CENTIGRADE_TO_KELVIN = 273.5;
const float lower_temp_limit = -90.0 + CENTIGRADE_TO_KELVIN;
const float upper_temp_limit = 36.4 + CENTIGRADE_TO_KELVIN;
static int dd;
static double bob, rab, shona, Tr;
/* signum of a number, used below */
short int SGN(double bob)
{
 short int signum;
 if (bob >= 0.0)
 signum = 1;
 else
 signum = -1;
 return (signum);
}
/* Nitrous oxide vapour pressure, Bar */
double nox_vp(double T_Kelvin)
{
 const float p[4] = {1.0f, 1.5f, 2.5f, 5.0f};
 const float b[4] = {-6.71893f, 1.35966f, -1.3779f, -4.051f};
 Tr = T_Kelvin / tCrit;
 rab = 1.0 - Tr;
 shona = 0.0;
 for (dd = 0; dd < 4; dd++)
 shona += b[dd] * pow(rab,p[dd]);
 bob = pCrit * exp((shona / Tr));
 return(bob);
}
/* Nitrous oxide saturated liquid density, kg/m3 */
double nox_Lrho(double T_Kelvin)
{
 const float b[4] = {1.72328f, -0.8395f, 0.5106f, -0.10412f};
 Tr = T_Kelvin / tCrit;
 rab = 1.0 - Tr;
 shona = 0.0;
 for (dd = 0; dd < 4; dd++)
 shona += b[dd] * pow(rab,((dd+1) / 3.0));
 bob = rhoCrit * exp(shona);
 return(bob);
}
/* Nitrous oxide saturated vapour density, kg/m3 */
double nox_Vrho(double T_Kelvin)
{
 const float b[5] = {-1.009f, -6.28792f, 7.50332f, -7.90463f, 0.629427f};
 Tr = T_Kelvin / tCrit;
 rab = (1.0 / Tr) - 1.0;
 shona = 0.0;
 for (dd = 0; dd < 5; dd++)
 shona += b[dd] * pow(rab,((dd+1) / 3.0));
 bob = rhoCrit * exp(shona);
 return(bob);
}
/* Nitrous liquid Enthalpy (Latent heat) of vaporisation, J/kg */
double nox_enthV(double T_Kelvin)
{
 const float bL[5] = {-200.0f, 116.043f, -917.225f, 794.779f, -589.587f};
 const float bV[5] = {-200.0f, 440.055f, -459.701f, 434.081f, -485.338f};
 double shonaL, shonaV;
 Tr = T_Kelvin / tCrit;
 rab = 1.0 - Tr;
 shonaL = bL[0];
 shonaV = bV[0];
 for (dd = 1; dd < 5; dd++)
 {
 shonaL += bL[dd] * pow(rab,(dd / 3.0)); /* saturated liquid enthalpy */
 shonaV += bV[dd] * pow(rab,(dd / 3.0)); /* saturated vapour enthalpy */
 }

 bob = (shonaV - shonaL) * 1000.0; /* net during change from liquid to vapour */
 return(bob);
}
/* Nitrous saturated liquid isobaric heat capacity, J/kg K */
double nox_CpL(double T_Kelvin)
{
 const float b[5] = {2.49973f, 0.023454f, -3.80136f, 13.0945f, -14.518f};
 Tr = T_Kelvin / tCrit;
 rab = 1.0 - Tr;
 shona = 1.0 + b[1] / rab;
 for (dd = 1; dd < 4; dd++)
 shona += b[(dd+1)] * pow(rab,dd);
 bob = b[0] * shona * 1000.0; /* convert from KJ to J */
 return(bob);
}
/* liquid nitrous thermal conductivity, W/m K */
double nox_KL(double T_Kelvin)
{
 const float b[4] = {72.35f, 1.5f, -3.5f, 4.5f};
 /* max. 10 deg C */
 if (T_Kelvin > 283.15)
 Tr = 283.15 / tCrit;
 else
 Tr = T_Kelvin / tCrit;
 rab = 1.0 - Tr;
 shona = 1.0 + b[3] * rab;
 for (dd = 1; dd < 3; dd++)
 shona += b[dd] * pow(rab,(dd / 3.0));
 bob = b[0] * shona / 1000; /* convert from mW to W */
 return(bob);
}
/* nitrous temperature based on pressure (bar) */
double nox_on_press(double P_Bar_abs)
{
 const float p[4] = {1.0f, 1.5f, 2.5f, 5.0f};
 const float b[4] = {-6.71893f, 1.35966f, -1.3779f, -4.051f};
 double pp_guess, step, tempK;
 step = -1.0;
 tempK = (tCrit - 0.1) - step;
 do /* iterative loop */
 {
 do
 {
 tempK += step;
 Tr = tempK / tCrit;
 rab = 1.0 - Tr;
 shona = 0.0;
 for (dd = 0; dd < 4; dd++)
 shona += b[dd] * pow(rab,p[dd]);
 pp_guess = pCrit * exp((shona / Tr));
 }
 while ( ((pp_guess - P_Bar_abs) * SGN(step)) < 0.0);
 step = step / (-2.0); /* reduce step size */
 }
 while (fabs((pp_guess - P_Bar_abs)) > 0.01);
 bob = tempK;
 return(bob); /* return temperature */
}
double bar_to_kPa(double press_bar)
{
    return press_bar * 100;
}

double kg_to_lbs(double mass_kg)
{
    return 2.204622622 * mass_kg;
}

double bar_to_psi(double press_bar)
{
    return press_bar * 14.503773773;
}
