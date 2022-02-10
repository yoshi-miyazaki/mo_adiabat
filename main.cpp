/*
  adiabat.cpp
  
  created by Yoshi Miyazaki on Feb 9, 2022
 */

#include "CGgibbsmin.h"
#include "gcell.h"


int main(){
    /* set the mantle potential tempeature ->
       and  solve for the temperature structure of magma ocean

       This code can only be used for Earth (for now). */


    // input potential temperature
    double Tpotl;
    cin >> Tpotl;
    
    
    /* --------------------------------------------------------
       model set up */
    int    n     = 200;
    double Psurf = 0;
    double rmin  = 3400e3; /* radius: bottom of magma-ocean */
    double rmax  = 6300e3; /* radius: surface               */
    
    /* --------------------------------------------------------
       composition */
    string           file_species     = "molecule_simple.txt";
    string           file_composition = "n_pyrolite.txt";    
    MoleculeData_G   list_species(file_species, 273.15, 1e5);
    initcomp         n_init(file_composition, list_species);
    
    /* --------------------------------------------------------
       initial & boundary condition */
    // double Tpotl = 2460;   /* potential temperature of magma-ocean */

    magma_ocean mo_system(n, rmin, rmax, Tpotl, n_init);


}
