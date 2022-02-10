/*
  gcell.cpp
  - grid for magma ocean solidification
  
  created by Yoshi on 2017, Apr 20
*/

#include "gcell.h"

/* gcell class ---------------------------------------------------*/
gcell::gcell(){
    T = 0; P = 0; dP = 0; n_melt.resize(0); n_solid.resize(0); list_txt = "";
}
gcell::gcell(const double r0, const double T0, const double P0, const double dP0,
             const tensor1d<double> n_init, const string s_list){
    r = r0;  T = T0;  P = P0;  dP = dP0;
    
    list_txt = s_list;
    /* for initialization, put all components in solid */
    n_solid  = n_init;   m_comp = n_solid.size();  
    n_melt.resize(0.0,m_comp);  /* and put 0 for all melt components */
    
    /* store molar mass of each species */
    MoleculeData_G _md(list_txt,T0,P0);
    mass_spec.resize(0.0,m_comp);
    for (int i=0; i<m_comp; i++){ mass_spec[i] = _md[i].getWeight(); }
    
    /* store mass balance matrix */
    massBalance    _mb(_md);
    massm = _mb.get_massBalance();
    
}
gcell gcell::operator=(const gcell &copy){
    /* copy and assign */
    r        = copy.r;
    T        = copy.T;      P  = copy.P;    dP       = copy.dP;   
    alpha    = copy.alpha;  Cp = copy.Cp;   rho      = copy.rho;
    list_txt = copy.list_txt;
    n_melt   = copy.n_melt;
    n_solid  = copy.n_solid;
    mass_spec= copy.mass_spec;
    massm    = copy.massm;
    m_comp   = copy.m_comp;
    
    return *this;
}
double gcell::melt_fraction(){
    /* re-run Gibbs free energy min*/
    /* calc volume fraction of melt
       ... applicable ONLY for the MgO-FeO-SiO2 ternary system */
    MoleculeData_G   slist(list_txt, T, P);
    massBalance      massbo(slist);
    Matrix<double>   massm  = massbo.get_massBalance().transpose();
    
    /* reassure n_melt & n_solid are stored correctly */
    tensor1d<double> ntotl = get_totalC();
    tensor1d<int>    phase = slist.getphase();
    n_melt = 0; n_solid = 0;
    for (int j=0; j<m_comp; j++){
        if (phase[j] == 1){ n_melt[j]  = ntotl[j];}
        else              { n_solid[j] = ntotl[j];}
    }
    
    tensor1d<double> e_mel = massm*n_melt, e_sol = massm*n_solid;
    double v_mel =(e_mel[0]*11.3 + e_mel[1]*22.67 + e_mel[2]*12.5)*1.03;
    double v_sol = e_sol[0]*11.3 + e_sol[1]*22.67 + e_sol[2]*12.5;
    //cout << "mfrac: " << v_mel << " " << v_sol << " " << e_mel[0] << " " << e_mel[1];
    //cout << " " << e_mel[2] << "/" << e_sol[0] << " " << e_sol[1] << " " << e_sol[2] << endl;
    
    return v_mel/(v_mel+v_sol);
}
void gcell::Gibbsmin_CG(){
    /* calculate equilibrium composition in the given T&P composition */
    MoleculeData_G    slist(list_txt, T, P);
    tensor1d<double>  n_init = n_melt + n_solid;
    gibbsminCG        res(slist, n_init);
    double            G0 = res.getGbest()*(R*T);
    //cout << "C: "<<n_melt[7]<<"  "<<n_melt[8]<<"  "<<n_melt[9]<<endl;
    
    /* store result */
    tensor1d<double> nbest = res.getnbest();
    tensor1d<int>    phase = res.getphase();
    int m_comp = nbest.size();
    if (m_comp != n_melt.size()){
        cout << "ERROR: no. of species not expected! " << m_comp;
        cout << " but, " << n_melt.size() << endl;
        for (int j=0; j<m_comp; j++){
            cout << slist[j].getMoleculeName()<<"\t"<<n_melt[j]<<" -> "<<nbest[j]<<endl;
        }
    }
    n_melt = 0; n_solid = 0;  /* reset composition */
    for (int j=0; j<m_comp; j++){
        if (phase[j] == 1){ n_melt[j]  = nbest[j];}
        else              { n_solid[j] = nbest[j];}
    }
    //cout << "D: "<< n_melt[7]<<"  "<<n_melt[8]<<"  "<<n_melt[9]<<endl;
    
    /* perturb to calc alpha, rho, and Cp */
    double eT = 0.1,  eP = 1e5;
    MoleculeData_G    slistP(list_txt, T, P+eP);
    gibbsminCG        resP(slistP, nbest);          double GP = resP.getGbest()*(R*T);
    
    MoleculeData_G    slistT(list_txt, T+eT, P);
    gibbsminCG        resT(slistT, nbest);          double GT = resT.getGbest()*(R*(T+eT));
    
    MoleculeData_G    slistTP(list_txt, T+eT, P+eP);
    gibbsminCG        resTP(slistTP, nbest);        double GTP = resTP.getGbest()*(R*(T+eT));

    MoleculeData_G    slistTT(list_txt, T+2*eT, P);
    gibbsminCG        resTT(slistTT, nbest);        double GTT = resTT.getGbest()*(R*(T+2*eT));
    
    double  V0 =  (GP-G0)/eP,  VT =  (GTP-GT)/eP,  umass = n_init*mass_spec;
    double  S0 = -(GT-G0)/eT,  ST = -(GTT-GT)/eT;
    alpha = (VT-V0)/(eT*(V0+VT)/2.0);   rho = umass/V0;   Cp = T*(ST-S0)/(eT*umass);
    
}
double gcell::adiabat(){
    /* return adiabat gradient (al*T)/(rho*Cp) */
    if (alpha < 0 or Cp < 0){
        double dT = recalc_adiabat(1.0);
    }
    
    return (alpha*T)/(rho*Cp);
}
double gcell::recalc_adiabat(double plus){
    double pert  = plus*10*((double)rand()/RAND_MAX*2 - 1.0);
    double pertP = plus*((double)rand()/RAND_MAX*2 - 1.0)*1e8;
    double eT = 0.1 + 0.1*(double)rand()/RAND_MAX,  eP = 1e4 + 1e5*(double)rand()/RAND_MAX;
    tensor1d<double> n_init = n_melt + n_solid;
    cout << "eT = " << eT<< " @" << T+pert << "K, P=" << (P+pertP)/1e9 << endl;
    
    MoleculeData_G    slist(list_txt, T+pert, P+pertP);
    gibbsminCG        res(slist, n_init);
    tensor1d<double>  nbest = res.getnbest();
    tensor1d<int>     phase = res.getphase();
    double            G0 = res.getGbest()*(R*(T+pert));
    
    int m_comp = nbest.size();
    n_melt = 0; n_solid = 0;  /* reset composition */
    for (int j=0; j<m_comp; j++){
        if (phase[j] == 1){ n_melt[j]  = nbest[j];}
        else              { n_solid[j] = nbest[j];}
    }
    
    MoleculeData_G    slistP(list_txt, T+pert, P+pertP+eP);
    gibbsminCG        resP(slistP, nbest);      double GP = resP.getGbest()*(R*(T+pert));
    tensor1d<double>  nP = resP.getnbest();
    
    MoleculeData_G    slistT(list_txt, T+pert+eT, P+pertP);
    gibbsminCG        resT(slistT, nbest);      double GT = resT.getGbest()*(R*(T+pert+eT));
    tensor1d<double>  nT = resT.getnbest();
    
    MoleculeData_G    slistTP(list_txt, T+pert+eT, P+pertP+eP);
    gibbsminCG        resTP(slistTP, nbest);    double GTP = resTP.getGbest()*(R*(T+pert+eT));
    tensor1d<double>  nTP = resTP.getnbest();
    
    MoleculeData_G    slistTT(list_txt, T+pert+2*eT, P+pertP);
    gibbsminCG        resTT(slistTT, nbest);    double GTT = resTT.getGbest()*(R*(T+pert+2*eT));
    tensor1d<double>  nTT = resTT.getnbest();
    
    double  V0 =  (GP-G0)/eP,  VT =  (GTP-GT)/eP,  umass = n_init*mass_spec;
    double  S0 = -(GT-G0)/eT,  ST = -(GTT-GT)/eT;
    alpha = (VT-V0)/(eT*V0);   Cp = T*(ST-S0)/(eT*umass);   rho = umass/V0;
    
    cout << "\t S0: " << S0 << "\t ST: " << ST << endl;
    cout << "++: " << plus << " , alpha = " << alpha << " , Cp = " << Cp << " ,rho = " << rho << endl;
    
    cout << "grid com: " << endl;
    for (int j=0; j<(int)nbest.size(); j++){ cout << nbest[j]<<"\t"<<nP[j]<<"\t"<<nT[j]<<"\t"<<nTP[j]<<"\t"<<nTT[j]<<endl; }
    
    if (alpha < 0 or Cp < 0){ plus++; return recalc_adiabat(plus); }
    else  { cout << "res:" << (alpha*T)/(rho*Cp) << endl; return (alpha*T)/(rho*Cp); }
}
