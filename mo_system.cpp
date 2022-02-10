/*
  mo_system.cpp
  - collection of grid for magma ocean solidification
  
  created by Yoshi Miyazaki on Feb 9, 2022
*/

#include "gcell.h"

double density(double P){
    double rho0, K0d, K0;
    if (P < 24e9){ // enstatite
        rho0 = 3300;  K0d = 5;    K0 = 125e9;
    }else{ // bridgmanite
        rho0 = 4100;  K0d = 3.97; K0 = 247e9;
    } // add ppv for super-Earths

    double rho = rho0*pow(1.+K0d*P/K0, 1/K0d);

    return rho;
}

magma_ocean::magma_ocean(const int& n0, const double& rmin, const double& rmax, 
                         const double& Tpotl, const initcomp& obj_init){
    
    /* construct collection of grid-cell "gcell"
       ... with the minimum set of information
       1. n0   : num. of "gcell" 
       2. rmin : bottom  radius of magma ocean
       3. rmax : surface radius of magma ocean 
       4. Tsurf: surface (and initial) temperature
       5. obj_init: object "initcomp" containing initial composition & list of species */
    
    n = n0;   magma.resize(n0); 
    
    double     dr = (rmax - rmin)/(n-1);
    double     T  = Tpotl;
    double     P  = 0;               /* assume that MO reaches the surface */
    string     list_species = obj_init.getlist();
    tensor1d<double> n_init = obj_init.getmole();
    
    for (int i=(n-1); i>=0; i--){
        double r  = rmin + dr*i;
        double rhom = density(P);
        double dP   = rhom*g*dr;
        tensor1d<double>  n_here = n_init*square(r/rmax);
        gcell  thisg(r, T, P, dP, n_here, list_species);
        magma[i] = thisg;
        magma[i].setrho(rhom);
        magma[i].setrhos(0.0);
        magma[i].setrhol(0.0);
	
	/* dP = int rho g dz */
        int div = 100;
        for (int j=div; j>=0; j--){
            r -= dr/div;
            rhom = density(P);
            dP = rhom*g*(dr/div);
            P += dP;
        }
        
    }
    therm_structure(Tpotl);  /* solve for whole magma ocean tempeature
                                gradient assmuing adiabat gradient */
}
void magma_ocean::therm_structure(double Tpotl){
    /* given the mantle potential tempeature, solve for the temperature structure
       of magma ocean up to the depth of r_rigid */
    double T  = Tpotl, dT = 0.0;
    
    for (int i=(n-1); i>=0; i--){

	/* output results */
        cout << i << "\t" << setw(7) << T << "\t" << magma[i].getP()/1e9 << endl;
        
        magma[i].setT(T);     /* set T calced from adiabat gradient */
        magma[i].Gibbsmin_CG();
        
        /* T & P for the next grid */
        double dP = magma[i].getdP();
        dT = magma[i].adiabat()*dP;   /* adiabat */
        T += dT;
        
        /* b/c dT is the critical part of this model
           and yet is based on a very unstable calculation, 
           the results are double checked: */
        double mf   = magma[i].melt_fraction();
        double mf1  = 0.0;
        double plus = 0.;
        if (i != (n-1)){ mf1  = magma[i+1].melt_fraction(); }
	
	while ((i<30 and (mf1-mf) > 0.025) or (dT<5 or dT>150)){
	    plus++;
	    dT= magma[i].recalc_adiabat(plus)*dP;
	    //cout<< "Rev: dT = " << dT << endl;
        }
        
    }
}
