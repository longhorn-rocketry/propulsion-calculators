#include <iostream>
#include "hybrid_tank.h"
using namespace std;

hybrid_tank :: hybrid_tank()
{
    press_or_temp = true; //
    initial_tank_pressure = 57.2265; //830 psi to bar
    tank_volume = 0.0084995122110964; //518.6720581 in3 to m3
    orifice_diameter = 0.0015113; //0.0595 in to m
    orifice_number = 22;
    chamber_pressure_bar = 34.4738; //500 psi to bar
    initial_ullage = 85;
    orifice_k2_coefficient = 30.0;
}
