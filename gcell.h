/*
  gcell.h
  - header file for class "mantle" to solve magma ocean solidication
  
  Created by Yoshi Miyazaki on Apr 17, 2017
 */

#ifndef MO_gcell
#define MO_gcell

#include "CGgibbsmin.h"
#include "const.h"
#include <omp.h>


const double CpSi = 1000;   /* heat capacity of Fe metal */
const double CpFe = 800;    /* heat capacity of Fe metal */
const double g    = 9.8;    /* gravitational constant */
const double Mcore  = (1.835+0.09675)*1e24;


/*--------------------------------------------------------------------------
 // Def of class: gcell
 ---------------------------------------------------------------------------*/
class gcell{
 public:
    gcell();
    gcell(const double r0, const double T0, const double P0, const double dP0,
          const tensor1d<double> n_init, const string s_list);
    gcell operator=(const gcell&);  /* assginment */
    
    void             setT(double T0){ T = T0;};
    void             setP(double P0){ P = P0;};
    void             setrho(double rho0) { rho  = rho0; };
    void             setrhos(double rho0){ rhos = rho0; };
    void             setrhol(double rho0){ rhol = rho0; };
    void             set_list(string _s){ list_txt = _s; };
    void             set_meltC (tensor1d<double> _v){ n_melt = _v; };
    void             set_solidC(tensor1d<double> _v){ n_solid= _v; };
    
    double           getr(){return r;};
    double           getT(){return T;};
    double           getP(){return P;};
    double           getdP(){return dP;};
    double           getrho()    {return rho;   };
    double           getrhos()   {return rhos;  };
    double           getrhol()   {return rhol;  };
    double           getalpha()  {return alpha; };
    double           getCp()     {return Cp;    };
    Matrix<double>   get_massm() {return massm; };
    tensor1d<double> get_meltC() {return n_melt; };
    tensor1d<double> get_solidC(){return n_solid;};
    tensor1d<double> get_totalC(){return n_melt+n_solid;};
    tensor1d<double> get_masvec(){return mass_spec;};
    
    void             Gibbsmin_CG();
    double           melt_fraction();
    double           adiabat();
    
    double           recalc_adiabat(double plus = 10);

 private:
    double           r;
    double           T;
    double           P;
    double           dP;
    double           rho;
    double           rhos;
    double           rhol;
    double           alpha;
    double           Cp;        /* in J/(K*kg) */
    
    string           list_txt;
    Matrix<double>   massm;
    tensor1d<double> n_melt;
    tensor1d<double> n_solid;
    tensor1d<double> mass_spec;
    int              m_comp;    /* size of n_melt & n_solid */
};

/*--------------------------------------------------------------------------
 // Def of class: magma_ocean
 ---------------------------------------------------------------------------*/
class magma_ocean{
 public:
    magma_ocean(const int&, const double&, const double&, const double&, const initcomp&);
    void  next(double);
    void  output(string);
    
    void  setTsurf(double Ts){ Tsurf = Ts; cout << "Tsurf is " << Tsurf << endl;}
    void  setTcore(double Tc){ Tcore = Tc; }
    
    double getTcore(){ return Tcore; }
    
 private:
    int               n;
    tensor1d<gcell>   magma;
    tensor1d<double>  v_perc;
    
    int              isep;
    int              isol;
    double           rrigid;
    double           rsolid;
    double           Tsurf;
    double           Tcore;
    
    void    therm_structure(double);
    void    grain_settle();
    void    porous(double);
    void    compact(double);
    void    homogenize();
    double  heat_capacity(double, double);
    double  rigid_boundary();
    double  solid_boundary();
};

#endif
