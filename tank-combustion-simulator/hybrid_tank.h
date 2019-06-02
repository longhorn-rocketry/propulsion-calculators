#ifndef HYBRID_TANK_H_INCLUDED
#define HYBRID_TANK_H_INCLUDED
#include <iostream>

class hybrid_tank
{
private:

public:
    double chamber_pressure_bar; //update
    int hybrid_fault;
    double tank_liquid_density;
    double injector_loss_coefficient; //update (K/(N*Ainj)^2)
    double tank_fluid_temperature_K;
    double tank_liquid_mass;
    double tank_pressure_bar;
    double tank_vapour_density;
    double tank_volume;
    double tank_propellant_contents_mass;
    double tank_vapour_mass;
    bool press_or_temp;
    double initial_tank_pressure;
    double initial_fluid_propellant_temp;
    double tank_vapour_volume;
    double initial_ullage;
    double tank_liquid_volume;
    double initial_liquid_propellant_mass;
    double initial_vapour_propellant_mass;
    double orifice_diameter;
    double orifice_k2_coefficient;
    double orifice_number;
    double mdot_tank_outflow;
    double Omdot_tank_outflow;
    double vaporised_mass_old;
    double old_liquid_nox_mass;
    double vapour_temperature_K;

    hybrid_tank();
};


#endif // HYBRID_TANK_H_INCLUDED
