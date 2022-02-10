//
//  gibbsminCG.h
//  Gibbs
//
//  Created by Yoshi Miyazaki on 2015/04/25.
//  Modified on 2017/4/8

#ifndef Gibbs_ConjugateG_h
#define Gibbs_ConjugateG_h

#include <cstdlib>
#include <cfloat>
#include <time.h>
#include "const.h"
#include "gibbse.h"
#include "initcomp.h"
#include "massb_oxide.h"

const double LIM  = 2.22044604925031e-16*2;
const int print = 0;

/*--------------------------------------------------------------------------
  define `melt_mixing` class
 ---------------------------------------------------------------------------*/
class melt_mixing{
 public:
    melt_mixing(string s1 = "", string s2 = "", string s3 = "");
    melt_mixing(string s1, string s2, string s3, tensor1d<double>);
    
    void  add_endmember(string, string, string);
    void  add_Margules(double, double, tensor1d<double>&);
    
    bool    is_insystem(string);
    void    reset(){ nmole = 0; };
    void    add_mole(string, double);   /* store mole amount of end member */
    double  tot_mole();                 /* return total amount in this solid-solution */
    double  activity(string);
 private:
    tensor1d<string>    end_members;
    tensor1d<double>    Wmix;
    tensor1d<double>    nmole;
};

/*--------------------------------------------------------------------------
  define `solution` class
 ---------------------------------------------------------------------------*/
class solution{
 public:
    solution(string s1 = "", string s2 = "", int no_sites=0);
    solution(string s1, string s2, string s3, int no_sites);
    solution(string s1, string s2, string s3, string s4, int no_sites);
    solution(string, double, int, string, double, int, string, double, int, int);
    solution(string, double, int, string, double, int, string, double, int, string, double, int, int);
    
    int     get_nosites(){ return no_sites; }
    bool    is_insystem(string);
    void    reset(){ nmole = 0; };
    void    add_mole(string, double);   /* store mole amount of end member */
    double  tot_mole();                 /* return total amount in this solid-solution */
    double  config(string, double);     /* return configuration entropy               */
 private:
    tensor1d<string>    end_members;
    int                 no_sites;
     
    tensor1d<int>       group;
    tensor1d<double>    contribution;
    tensor1d<double>    nmole;
};

/*--------------------------------------------------------------------------
  define `gibbsminCG` class
 ---------------------------------------------------------------------------*/
class gibbsminCG{
 public:
    gibbsminCG(MoleculeData_G&, tensor1d<double>&);
    gibbsminCG(MoleculeData_G&, tensor1d<double>&, tensor1d<double>);
    
    /* output result */
    tensor1d<double> getnbest(){return nbest;};
    tensor1d<int>    getphase(){return phase;};
    double           getGbest(){return gibbsE(nbest);};
    double           getG_in() {return gibbsE(n_in);};
    void             result(double, double, string&);
    double           melt_mfrac();
    double           melt_vfrac();
    double           meltfrac(double, double);
    double           meltfrac(double, double, string&);

 private:
    /* pressure and temperature */
    double           T;
    double           P;

    /* calculate Gibbs energy and gradient */
    void             conjugateG();      /* minimize function gibbsE from n_in */
    double           gibbsE(tensor1d<double>&);
    tensor1d<double> grad(tensor1d<double>, int);
    
    /* solid solution and melt mixing */
    melt_mixing        meltSi;
    tensor1d<solution> ss_list;   /* store the list of solid solution */
    void               create_sslist();
    
    /* initial n & massbalance */
    tensor1d<double> n_in;
    tensor1d<double> atomic;
    Matrix<double>   massm;
    
    /* gibbs energy data & phase data */
    MoleculeData_G   list;
    tensor1d<double> gibbse;
    tensor1d<int>    phase;
    tensor1d<double> getGibbslist(MoleculeData_G&);
    tensor1d<int>    getphaselist(MoleculeData_G&);
    
    /* conjugate gradient */
    double           G;
    tensor1d<double> nbest;
    int              dore;
    
    void             mass_balance(tensor1d<double>&, tensor1d<int>&);
    Matrix<double>   compose_projection(Matrix<double>&);
    void             zerogas(tensor1d<double>&, tensor1d<double>&, tensor1d<int>&);
    void             avoidneg(tensor1d<double>&, tensor1d<double>&, tensor1d<int>&);
    tensor1d<double> shrink  (tensor1d<double>&, tensor1d<double>&, double&);
    tensor1d<double> perturb_zero(tensor1d<double>&, tensor1d<double>&);
    void             bisection(tensor1d<double>, tensor1d<double>&, tensor1d<double>&, double, tensor1d<int>&, int&, int&);
    bool             ck_smaller(tensor1d<double>&, tensor1d<double>&, tensor1d<double>&, tensor1d<int>&);
    double           dG_direction(tensor1d<double>&, tensor1d<int>&, tensor1d<double>);

    /* projection matrix P = I - B (BT B)-1 BT */
    Matrix<double>   projv;
    Matrix<double>   projv_saved;    /* reduce calculation */
    tensor1d<int>    exist_saved;    /* if diminish is same as previous */
};

#endif
