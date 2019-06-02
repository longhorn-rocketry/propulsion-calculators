#ifndef NITROUS_OXIDE_H_INCLUDED
#define NITROUS_OXIDE_H_INCLUDED
/* header file nitrous oxide.h */
/* Nitrous oxide vapour pressure, Bar */
double nox_vp(double T_Kelvin);
/* Nitrous oxide saturated liquid density */
double nox_Lrho(double T_Kelvin);
/* Nitrous oxide saturated vapour density */
double nox_Vrho(double T_Kelvin);
/* Nitrous oxide Enthalpy (Latent heat) of vapourisation */
double nox_enthV(double T_Kelvin);
/* Nitrous oxide saturated liquid isobaric heat capacity */
double nox_CpL(double T_Kelvin);
/* Nitrous oxide liquid thermal conductivity, W/m K */
double nox_KL(double T_Kelvin);
/* mean temperature K based on pressure */
double nox_on_press(double P_Bar_abs);
/* bar to kPa conversion */
double bar_to_kPa(double press_bar);
/* bar to psi conversion */
double bar_to_psi(double press_bar);
/* kg to pounds conversion */
double kg_to_lbs(double mass_kg);


#endif // NITROUS_OXIDE_H_INCLUDED
