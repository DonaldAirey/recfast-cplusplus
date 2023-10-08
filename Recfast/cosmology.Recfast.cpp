//========================================================================================
// Author: Jens Chluba
// Last modification: Oct 2010
// CITA, University of Toronto
// All rights reserved.
//========================================================================================

#include <cmath>
#include <math.h>
#include <iostream>
#include <algorithm>

#include "constants.Recfast.h"
#include "cosmology.Recfast.h"

using namespace std;
using namespace RECFAST_physical_constants; // defined in constants.Recfast.h

//========================================================================================
// Global Variables; defined in cosmology.Recfast.h
//========================================================================================
extern struct Input input;

//========================================================================================
// Hubble Function in 1/sec
//========================================================================================
double H_z(double z)
{
    double t = (RF_V + RF_V * z - sqrt((1 + z) * (pow(RF_V, 2) + pow(RF_V, 2) * z - 2 * RF_A * RF_V * RF_age + pow(RF_A, 2) * pow(RF_age, 2)))) / (RF_A + RF_A * z);
    return (2*RF_V-2*RF_A*t)/(2*RF_V*t-RF_A*pow(t,2));
}

//========================================================================================
// Hydrogen number density in m^-3
//========================================================================================
double NH(double z) 
{
  double baryonDensity = 2.82211E-26;
  double numberDensity = baryonDensity * (1.0 - input.YP) / RF_mHatom;
  return numberDensity * pow(1.0 + z, 3);
}

//========================================================================================
// CMB temperature at z
//========================================================================================
double TCMB(double z){ return input.To*(1.0+z); }

//========================================================================================
// compute total contribution from relativistic particles (photon & neutrinos)
//========================================================================================
double calc_Orel(double TCMB0, double Nnu, double h100)
{ 
    double H0=h100*100.0*1.0e+5/RF_Mpc;
    double a=RF_aRad*pow(TCMB0, 4)/RF_cLight/RF_cLight; 
    double b=3.0*pow(H0, 2)/(8.0*RF_PI*RF_G);  
    return a/b*(1.0+Nnu*(7.0/8.0)*pow(4.0/11.0, 4.0/3.0));
}
