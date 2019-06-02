#include "nitrous_oxide.h" /* header file for nitrous oxide property calcs subroutines */
#include "hybrid_tank.h" /* header file for hybrid tank properties */
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
using namespace std;
hybrid_tank hybrid;
const float nox_pCrit = 72.51f; /* critical pressure, Bar Abs */
const float nox_rhoCrit = 452.0f; /* critical density, kg/m3 */
const float nox_tCrit = 309.57f; /* critical temperature, Kelvin (36.42 Centigrade) */
const float nox_ZCrit = 0.28f; /* critical compressibility factor */
const float nox_gamma = 1.3; /* average over subcritical range */
const double delta_time = .001;
/* prototypes */
static double injector_model(double upstream_pressure, double downstream_pressure);
//static bool first = true;
/* square */
double SQR(double bob)
{
 return((bob * bob));
}
#define CENTIGRADE_TO_KELVIN 273.15 // to Kelvin
#define BAR_TO_PASCALS 100000.0
#define PASCALS_TO_BAR (1.0 / BAR_TO_PASCALS)
/* calculate injector pressure drop (Bar) and mass flowrate (kg/sec) */
static double injector_model(double upstream_pressure, double downstream_pressure)
{
 double mass_flowrate;
 double pressure_drop;
 pressure_drop = upstream_pressure - downstream_pressure; /* Bar */
 /* reality check */
 if (pressure_drop < 0.00001)
 pressure_drop = 0.00001;
 /* is injector pressure drop lower than 20 percent of chamber pressure? */
 if ((pressure_drop / hybrid.chamber_pressure_bar) < 0.2)
 hybrid.hybrid_fault = 3; // too low for safety
 /* Calculate fluid flowrate through the injector, based on the */
 /* total-pressure loss factor between the tank and combustion chamber */
 /* (injector_loss_coefficient includes K coefficient and orifice cross-sectional areas) */
 mass_flowrate =
 sqrt((2.0 * hybrid.tank_liquid_density * pressure_drop / hybrid.injector_loss_coefficient));
 return(mass_flowrate); /* kg/sec */
}
/* Equilibrium (instantaneous boiling) tank blowdown model */
/* Empty tank of liquid nitrous */
void Nitrous_tank_liquid(void)
{
 double bob;
 double Chamber_press_bar_abs;
 double delta_outflow_mass, deltaQ, deltaTemp;
 double Enth_of_vap;
 double Spec_heat_cap;
 double tc;
 static double lagged_bob = 0.0;
 /* blowdown simulation using nitrous oxide property calcs subroutines */
 /* update last-times values, O = 'old' */
 hybrid.Omdot_tank_outflow = hybrid.mdot_tank_outflow;
 Enth_of_vap = nox_enthV(hybrid.tank_fluid_temperature_K); /* Get enthalpy (latent heat) of
vaporisation */
 Spec_heat_cap = nox_CpL(hybrid.tank_fluid_temperature_K); /* Get specific heat capacity
of the liquid nitrous */
 /* Calculate the heat removed from the liquid nitrous during its vaporisation */
 deltaQ = hybrid.vaporised_mass_old * Enth_of_vap;//PICK ARBITRARY NON ZERO VALUE FOR Mv
 /* temperature drop of the remaining liquid nitrous due to losing this heat */
 deltaTemp = -(deltaQ / (hybrid.tank_liquid_mass * Spec_heat_cap));
 hybrid.tank_fluid_temperature_K += deltaTemp; /* update fluid temperature */
 /* reality checks */
 if (hybrid.tank_fluid_temperature_K < (-90.0 + CENTIGRADE_TO_KELVIN))
 {
 hybrid.tank_fluid_temperature_K = (-90.0 + CENTIGRADE_TO_KELVIN); /* lower limit */
 hybrid.hybrid_fault = 1;
 }
 else if (hybrid.tank_fluid_temperature_K > (36.0 + CENTIGRADE_TO_KELVIN))
 {
 hybrid.tank_fluid_temperature_K = (36.0 + CENTIGRADE_TO_KELVIN); /* upper limit */
 hybrid.hybrid_fault = 2;
 }
 /* get current nitrous properties */
 hybrid.tank_liquid_density = nox_Lrho(hybrid.tank_fluid_temperature_K);
 hybrid.tank_vapour_density = nox_Vrho(hybrid.tank_fluid_temperature_K);
 hybrid.tank_pressure_bar = nox_vp(hybrid.tank_fluid_temperature_K); /* vapour pressure,
Bar abs */
 Chamber_press_bar_abs = hybrid.chamber_pressure_bar; /* Bar Abs */
 /* calculate injector pressure drop and mass flowrate */
 hybrid.mdot_tank_outflow = injector_model(hybrid.tank_pressure_bar, Chamber_press_bar_abs);
  /* integrate mass flowrate using Addams second order integration formula */
 /* (my preferred integration formulae, feel free to choose your own.) */
 /* Xn=X(n-1) + DT/2 * ((3 * Xdot(n-1) - Xdot(n-2)) */
 /* O infront of a variable name means value from previous timestep (Old) */
 delta_outflow_mass = 0.5 * delta_time * (3.0 * hybrid.mdot_tank_outflow - hybrid.Omdot_tank_outflow);
 /* drain the tank based on flowrates only */
     hybrid.tank_propellant_contents_mass -= delta_outflow_mass; /* update mass within tank
    for next iteration */
    hybrid.old_liquid_nox_mass -= delta_outflow_mass; /* update liquid mass within tank for next
    iteration */

     /* now the additional effects of phase changes */
     /* The following equation is applicable to the nitrous tank, containing saturated nitrous: */
     /* tank_volume = liquid_nox_mass / liquid_nox_density + gaseous_nox_mass /
    gaseous_nox_density */
     /* Rearrage this equation to calculate current liquid_nox_mass */
     bob = (1.0 / hybrid.tank_liquid_density) - (1.0 / hybrid.tank_vapour_density);
     hybrid.tank_liquid_mass = (hybrid.tank_volume - (hybrid.tank_propellant_contents_mass /
    hybrid.tank_vapour_density)) / bob;
    //cout<<"bob "<<bob<<endl<<"liquid density "<<hybrid.tank_liquid_density<<endl<<"vol "<<hybrid.tank_volume<<endl<<"nox mass "<<hybrid.tank_propellant_contents_mass<<endl<<"vapor density "<<hybrid.tank_vapour_density<<endl;
     hybrid.tank_vapour_mass = hybrid.tank_propellant_contents_mass -
    hybrid.tank_liquid_mass;
     /* update for next iteration */
     bob = hybrid.old_liquid_nox_mass - hybrid.tank_liquid_mass;
     /* Add a 1st-order time lag (of 0.15 seconds) to aid numerical */
     /* stability (this models the finite time required for boiling) */
     tc = delta_time / 0.15;
    lagged_bob = tc * (bob - lagged_bob) + lagged_bob; // 1st-order lag
     hybrid.vaporised_mass_old = lagged_bob;
     // Check for model fault at nearly zero liquid oxidiser mass
     // If this occurs, use the fault flag to trigger burnout
     if (hybrid.tank_liquid_mass > hybrid.old_liquid_nox_mass)
     hybrid.hybrid_fault = 9;
     /* update tank contents for next iteration */
     hybrid.old_liquid_nox_mass = hybrid.tank_liquid_mass;
}
/* Linear interpolation routine, with limiters added */
/* incase x isn't within the range range x1 to x2 */
/* Limits updated to allow descending x values
 x = input value
 x1 = MINIMUM BOUNDS ( minimum x )
 y1 = MINIMUM VALUE ( output value at MINIMUM BOUNDS )
 x2 = MAXIMUM BOUNDS ( maximum x )
 y2 = MAXIMUM VALUE ( output value at MAXIMUM BOUNDS )
*/
double LinearInterpolate(double x, double x1, double y1, double x2, double y2)
{
 // This procedure extrapolates the y value for the x position
 // on a line defined by x1,y1; x2,y2
 double c, m, y; // the constants to find in y=mx+b
 if ( (x1 < x2) && ((x <= x1) || (x >= x2)) ) // ascending x values
 {
 if (x <= x1) return y1; else return y2;
 }
 else if ( (x1 > x2) && ((x >= x1) || (x <= x2)) ) // descending x values
 {
 if (x >= x1) return y1; else return y2;
 }
 else
 {
 m = (y2 - y1) / (x2 - x1); // calculate the gradient m
 c = y1 - m * x1; // calculate the y-intercept c
 y = m * x + c; // the final calculation
 return(y);
 }
}

/* Compressibility factor of vapour (subcritical) on the saturation line */
double compressibility_factor(double P_Bar_abs, double pCrit, double ZCrit)
{
 double Z;
 Z = LinearInterpolate(P_Bar_abs, 0.0, 1.0, pCrit, ZCrit);
 return(Z);
}
/* subroutine to model the tank emptying of vapour only */
/* Isentropic vapour-only blowdown model */
void subcritical_tank_no_liquid()
{
 static int Aim, OldAim;
 static bool first = true;
 double bob;
 double current_Z, current_Z_guess;
 double delta_outflow_mass;
 double step;
 static double initial_vapour_density, initial_vapour_mass, initial_vapour_pressure_bar;
 static double initial_vapour_temp_K, initial_Z;
 // capture initial conditions
 if (first == true)
 {
 initial_vapour_temp_K = hybrid.tank_fluid_temperature_K;
 initial_vapour_mass = hybrid.tank_vapour_mass;
 initial_vapour_pressure_bar = hybrid.tank_pressure_bar;
 initial_vapour_density = hybrid.tank_vapour_density;
 initial_Z = compressibility_factor(hybrid.tank_pressure_bar, nox_pCrit, nox_ZCrit);
 hybrid.Omdot_tank_outflow = 0.0; // reset
 first = false;
 }
 // calculate injector pressure drop and mass flowrate
 hybrid.mdot_tank_outflow = injector_model(hybrid.tank_pressure_bar, hybrid.chamber_pressure_bar);
 /* integrate mass flowrate using Addams second order integration formula */
 /* Xn=X(n-1) + DT/2 * ((3 * Xdot(n-1) - Xdot(n-2)) */
delta_outflow_mass
 = 0.5 * delta_time * (3.0 * hybrid.mdot_tank_outflow - hybrid.Omdot_tank_outflow);
 // drain the tank based on flowrates only
 hybrid.tank_propellant_contents_mass -= delta_outflow_mass; // update mass within tank for next iteration
 // drain off vapour
hybrid.tank_vapour_mass -= delta_outflow_mass; // update vapour mass within tank for next iteration
 // initial guess
 current_Z_guess = compressibility_factor(hybrid.tank_pressure_bar, nox_pCrit, nox_ZCrit);
 step = 1.0 / 0.9; // initial step size
 OldAim = 2; Aim = 0; // flags used below to home-in
 // recursive loop to get correct compressibility factor
 do
 {
 // use isentropic relationships
 bob = nox_gamma - 1.0;
 hybrid.vapour_temperature_K = initial_vapour_temp_K
 * pow(((hybrid.tank_vapour_mass * current_Z_guess) / (initial_vapour_mass * initial_Z)), bob);
 bob = nox_gamma / (nox_gamma - 1.0);
 hybrid.tank_pressure_bar = initial_vapour_pressure_bar
 * pow((hybrid.vapour_temperature_K / initial_vapour_temp_K), bob);
 current_Z = compressibility_factor(hybrid.tank_pressure_bar, nox_pCrit, nox_ZCrit);
 OldAim = Aim;
 if (current_Z_guess < current_Z)
 {
 current_Z_guess *= step;
 Aim = 1;
 }
 else
 {
 current_Z_guess /= step;
 Aim = -1;
 }
 /* check for overshoot of target, and if so, reduce step nearer to 1.0 */
 if (Aim == -OldAim)
 step = sqrt(step);
 /* leave loop upon convergence to required accuracy */
 }
 while ( ((current_Z_guess / current_Z) > 1.000001) ||
 ((current_Z_guess / current_Z) < (1.0 / 1.000001)) );
 bob = 1.0 / (nox_gamma - 1.0);
 hybrid.tank_vapour_density = initial_vapour_density
 * pow((hybrid.vapour_temperature_K / initial_vapour_temp_K), bob);
}
/* Subroutine to initialise main program variables */
/* Gets called only once, prior to firing */
void initialise_hybrid_engine(void)
{
 hybrid.hybrid_fault = 0;
 hybrid.tank_vapour_mass = 0.0;
 hybrid.mdot_tank_outflow = 0.0;
 /* set either initial nitrous vapour (tank) pressure */
 /* or initial nitrous temperature (deg Kelvin) */
 if (hybrid.press_or_temp == true)
 hybrid.tank_fluid_temperature_K = nox_on_press(hybrid.initial_tank_pressure); /* set tank
pressure */
 else
 { /* set nitrous temperature */
 hybrid.tank_fluid_temperature_K = hybrid.initial_fluid_propellant_temp
 + CENTIGRADE_TO_KELVIN;
 hybrid.initial_tank_pressure = nox_vp(hybrid.tank_fluid_temperature_K);
 }
 /* reality check */
 if (hybrid.tank_fluid_temperature_K > (36.0 + CENTIGRADE_TO_KELVIN))
 {
 hybrid.tank_fluid_temperature_K = 36.0 + CENTIGRADE_TO_KELVIN;
 hybrid.hybrid_fault = 2;
 }
 /* get initial nitrous properties */
 hybrid.tank_liquid_density = nox_Lrho(hybrid.tank_fluid_temperature_K);
 hybrid.tank_vapour_density = nox_Vrho(hybrid.tank_fluid_temperature_K);
 /* base the nitrous vapour volume on the tank percentage ullage (gas head-space) */
 hybrid.tank_vapour_volume = (hybrid.initial_ullage / 100.0) * hybrid.tank_volume;
 hybrid.tank_liquid_volume = hybrid.tank_volume - hybrid.tank_vapour_volume;
 hybrid.tank_liquid_mass = hybrid.tank_liquid_density * hybrid.tank_liquid_volume;
 hybrid.tank_vapour_mass = hybrid.tank_vapour_density * hybrid.tank_vapour_volume;
 hybrid.tank_propellant_contents_mass = hybrid.tank_liquid_mass
 + hybrid.tank_vapour_mass; /* total mass within tank
*/
 /* initialise values needed later */
 hybrid.old_liquid_nox_mass = hybrid.tank_liquid_mass;
 //hybrid.old_vapour_nox_mass = hybrid.tank_vapour_mass;
 hybrid.initial_liquid_propellant_mass = hybrid.tank_liquid_mass;
 hybrid.initial_vapour_propellant_mass = hybrid.tank_vapour_mass;
 /* guessed initial value of amount of nitrous vaporised per iteration */
 /* in the nitrous tank blowdown model (actual value is not important) */
 hybrid.vaporised_mass_old = 0.001;

 /* individual injector orifice total loss coefficent K2 */
 double Ainj = M_PI * SQR((hybrid.orifice_diameter / 2.0)); /* orifice cross sectional area */
  hybrid.injector_loss_coefficient
 = (hybrid.orifice_k2_coefficient / (SQR((hybrid.orifice_number * Ainj)) ) )
 * PASCALS_TO_BAR;
}
int main()
{
//    cout<<nox_CpL(275)<<endl;
//    cout<<nox_KL(275)<<endl;
//    cout<<nox_Lrho(275)<<endl;
//    cout<<nox_vp(275)<<endl;
//    cout<<nox_Vrho(275)<<endl;
//    cout<<nox_on_press(33.0948)<<endl;
//    setprecision(6);
//    initialise_hybrid_engine(); //Initialize engine properties
//    cout<<hybrid.initial_tank_pressure<<endl;
//    cout<<"liq: "<<hybrid.tank_liquid_mass<<endl<<"fault; "<<hybrid.hybrid_fault<<endl;
//    Nitrous_tank_liquid();
//    cout<<hybrid.tank_pressure_bar<<endl;
//    cout<<"liq: "<<hybrid.tank_liquid_mass<<endl<<"fault; "<<hybrid.hybrid_fault<<endl;
//    Nitrous_tank_liquid();
//    cout<<hybrid.tank_pressure_bar<<endl;
//    cout<<"liq: "<<hybrid.tank_liquid_mass<<endl<<"fault; "<<hybrid.hybrid_fault<<endl;
//    Nitrous_tank_liquid();
//    cout<<hybrid.tank_pressure_bar<<endl;
//    cout<<"liq: "<<hybrid.tank_liquid_mass<<endl<<"fault; "<<hybrid.hybrid_fault<<endl;
//    Nitrous_tank_liquid();
//    cout<<hybrid.tank_pressure_bar<<endl;
//    cout<<"liq: "<<hybrid.tank_liquid_mass<<endl<<"fault; "<<hybrid.hybrid_fault<<endl;
//    Nitrous_tank_liquid();
//    cout<<hybrid.tank_pressure_bar<<endl;
//    cout<<"liq: "<<hybrid.tank_liquid_mass<<endl<<"fault; "<<hybrid.hybrid_fault<<endl;
//    Nitrous_tank_liquid();
//    cout<<hybrid.tank_pressure_bar<<endl;
//    cout<<"liq: "<<hybrid.tank_liquid_mass<<endl<<"fault; "<<hybrid.hybrid_fault<<endl;
//    Nitrous_tank_liquid();
//    cout<<hybrid.tank_pressure_bar<<endl;
//    cout<<"liq: "<<hybrid.tank_liquid_mass<<endl<<"fault; "<<hybrid.hybrid_fault<<endl;
//    Nitrous_tank_liquid();
//    cout<<hybrid.tank_pressure_bar<<endl;
//    cout<<"liq: "<<hybrid.tank_liquid_mass<<endl<<"fault; "<<hybrid.hybrid_fault<<endl;
//    Nitrous_tank_liquid();
//    cout<<hybrid.tank_pressure_bar<<endl;
//    cout<<"liq: "<<hybrid.tank_liquid_mass<<endl<<"fault; "<<hybrid.hybrid_fault<<endl;
    int i = 0; //index

    /* Initialize data vectors */
    vector<double> pressure_kPa(0);
    vector<double>pressure_psi(0);
    vector<double> mass_flow(0);
    vector<double>temp(0);
    vector<double>liquid_mass(0);
    vector<double>vapor_mass(0);

    initialise_hybrid_engine();

    /* Tank liquid phase */
    while(hybrid.tank_liquid_mass>.01)
    {
        Nitrous_tank_liquid();
        //cout<<"Time: "<<i*.001<<" Tl: "<<hybrid.tank_fluid_temperature_K<<" P: "<<hybrid.tank_pressure_bar<<" ml: "<<hybrid.tank_liquid_mass<<endl;
        pressure_kPa.push_back(bar_to_kPa(hybrid.tank_pressure_bar));
        pressure_psi.push_back(bar_to_psi(hybrid.tank_pressure_bar));
        mass_flow.push_back(hybrid.mdot_tank_outflow);
        temp.push_back(hybrid.tank_fluid_temperature_K);
        liquid_mass.push_back(hybrid.tank_liquid_mass);
        vapor_mass.push_back(hybrid.tank_vapour_mass);
        i++;
    }

    /* Tank vapor only phase */
    for(int x = 0; hybrid.mdot_tank_outflow > .002; x++)
    {
        subcritical_tank_no_liquid();
        //cout<<"Time: "<<i*.001<<" Tv: "<<hybrid.vapour_temperature_K<<" P: "<<hybrid.tank_pressure_bar<<" mv: "<<hybrid.tank_vapour_mass<<endl;
        pressure_kPa.push_back(bar_to_kPa(hybrid.tank_pressure_bar));
        pressure_psi.push_back(bar_to_psi(hybrid.tank_pressure_bar));
        mass_flow.push_back(hybrid.mdot_tank_outflow);
        temp.push_back(hybrid.vapour_temperature_K);
        liquid_mass.push_back(hybrid.tank_liquid_mass);
        vapor_mass.push_back(hybrid.tank_vapour_mass);
        i++;
    }

    ofstream myout;
    /* data for excel */
    myout.open("data.csv");
    myout << "time (s), Pressure (kPa), Mass Flow Rate (kg/s), Temperature (K), Liquid Mass (kg), Vapor Mass (kg)" << endl; //x-axis
    for(int index = 0; index < i; index++)
    {
        myout << index * .001 << ", " << pressure_kPa[index] << ", " << mass_flow[index] << ", " << temp[index] << ", " << liquid_mass[index] << ", " << vapor_mass[index] << endl;
    }
    myout.close();
    /* data for gnuplot */
    myout.open("PvT.dat"); //PvT
    for(int index = 0; index < i; index++)
    {
        myout << index * .001 << " " << pressure_psi[index] << endl;
    }
    myout.close();
    myout.open("MvT.dat");
    for(int index = 0; index < i; index++)
    {
        myout << index * .001 <<  " " << kg_to_lbs(mass_flow[index]) << endl;
    }
    myout.close();
}
