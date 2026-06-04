#pragma once
#include <vector>
#include <complex>

using namespace std;

//value at grid center
//--------------------------------
struct GridCenter{ 

    //coordinate (special definition)
    //--------------------------------
    vector<double> x;
    vector<double> r;
    //--------------------------------

    //number densitty
    //--------------------------------
    vector<vector<double> > rhoe;
    vector<vector<double> > rhoe_old; 
   
    vector<vector<double> > rhoi;
    vector<vector<double> > rhoi_old;

    vector<vector<double> > rhon;
    vector<vector<double> > rhon_old;

    vector<vector<double> > rhom;
    vector<vector<double> > rhom_old;
    //--------------------------------
    
    //bulk velocity
    //--------------------------------
    vector<vector<double> > Uip;
    vector<vector<double> > Uip_old;
    //--------------------------------

    //energy
    //--------------------------------
    vector<vector<double> > rhoeps;
    vector<vector<double> > rhoeps_old;
    //--------------------------------

    //collision frequency
    //--------------------------------
    vector<vector<double> > nu_m; //電子運動量移動周波数
    vector<vector<double> > nu_m1; //電子運動量移動周波数 (有効衝突追加)
    vector<vector<double> > nui_m; //イオン衝突周波数
    vector<vector<double> > nu_elas; //弾性周波数
    vector<vector<double> > nu_super; //超弾性周波数
    vector<vector<double> > nu_ionz; //電離衝突周波数
    vector<vector<double> > nu_ionzStep; //段階電離衝突周波数
    vector<vector<double> > nu_exc; //励起周波数
    vector<vector<double> > nu_excReso; //励起(共鳴)周波数
    vector<vector<double> > nu_excMeta; //励起(準安定)周波数
    vector<vector<double> > nu_excStep; //段階励起周波数
    vector<vector<double> > nu_ano; //異常衝突周波数
    vector<vector<double> > nu_ano_IAT; //異常衝突周波数
    //--------------------------------

    //production rate
    //--------------------------------
    vector<vector<double> > rate_ionize;
    //--------------------------------

    //energy exchange
    //--------------------------------
    vector<vector<double> > rate_eloss; //電子エネルギー獲得量
    vector<vector<double> > delossdeps; //電子エネルギーロスの微係数
    //--------------------------------

    //electron mobility
    //--------------------------------
    vector<vector<double> > mu_para;
    vector<vector<double> > mu_perp;
    vector<vector<double> > mu_min;
    vector<vector<double> > mu_xx;
    vector<vector<double> > mu_rr;
    vector<vector<double> > mu_xr;
    vector<vector<double> > mu_inv_xx;
    vector<vector<double> > mu_inv_rr;
    vector<vector<double> > mu_inv_xr;
    //--------------------------------
    
    //temperature
    //--------------------------------
    vector<vector<double> > Te;
    vector<vector<double> > Te_old;
    //--------------------------------

    //field
    //--------------------------------
    vector<vector<double> > phi;
    vector<vector<double> > phi_old;

    vector<vector<double> > Ap;
    vector<vector<double> > Bx;
    vector<vector<double> > Br;

    vector<vector<complex<double> > > E1p;
    //--------------------------------

    //current
    //--------------------------------
    vector<vector<complex<double> > > J1p;
    //--------------------------------

    //permittivity
    //--------------------------------
    vector<vector<double> > epsr;
    //--------------------------------

    //electron Mach number
    //--------------------------------
    vector<vector<double> > Mache;
    //--------------------------------

    //power absorption density
    //--------------------------------
    vector<vector<double> > Pabs;
    //--------------------------------

    //judge boundary
    //--------------------------------
    vector<vector<int> > jdgBnd_Ep;
    vector<vector<int> > jdgBnd_flc;
    //--------------------------------

    //region index
    //--------------------------------
    vector<vector<int> > i_flc_bl;
    vector<vector<int> > j_flc_bl;
    //--------------------------------

    //index-convergion
    //--------------------------------
    vector<vector<int> > kp;   //convergion of (i,j) → k (E-field)
    //--------------------------------

    //constructor
    //--------------------------------
    explicit GridCenter(int ni,int nj)
        :
        //coordinate
        //--------------------------------
        x(ni+2,0.0),
        r(nj+2,0.0),
        //--------------------------------

        //number densitty
        //--------------------------------
        rhoe(ni+2,vector<double>(nj+2,0.0)),
        rhoe_old(ni+2,vector<double>(nj+2,0.0)),

        rhoi(ni+2,vector<double>(nj+2,0.0)),
        rhoi_old(ni+2,vector<double>(nj+2,0.0)),

        rhom(ni+2,vector<double>(nj+2,0.0)),
        rhom_old(ni+2,vector<double>(nj+2,0.0)),

        rhon(ni+2,vector<double>(nj+2,0.0)),
        rhon_old(ni+2,vector<double>(nj+2,0.0)),
        //--------------------------------

        //bulk velocity
        //--------------------------------
        Uip(ni+2,vector<double>(nj+2,0.0)),
        Uip_old(ni+2,vector<double>(nj+2,0.0)),
        //--------------------------------

        //energy
        //--------------------------------
        rhoeps(ni+2,vector<double>(nj+2,0.0)),
        rhoeps_old(ni+2,vector<double>(nj+2,0.0)),
        //--------------------------------

        //collision frequency
        //--------------------------------
        nu_m(ni+2,vector<double>(nj+2,0.0)), //電子運動量移動周波数
        nu_m1(ni+2,vector<double>(nj+2,0.0)), //電子運動量移動周波数 (有効衝突追加)
        nui_m(ni+2,vector<double>(nj+2,0.0)), //イオン衝突周波数
        nu_elas(ni+2,vector<double>(nj+2,0.0)), //弾性周波数
        nu_super(ni+2,vector<double>(nj+2,0.0)), //超弾性周波数
        nu_ionz(ni+2,vector<double>(nj+2,0.0)), //電離衝突周波数
        nu_ionzStep(ni+2,vector<double>(nj+2,0.0)), //段階電離衝突周波数
        nu_exc(ni+2,vector<double>(nj+2,0.0)), //励起周波数
        nu_excReso(ni+2,vector<double>(nj+2,0.0)), //励起(共鳴)周波数
        nu_excMeta(ni+2,vector<double>(nj+2,0.0)), //励起(準安定)周波数
        nu_excStep(ni+2,vector<double>(nj+2,0.0)), //段階励起周波数
        nu_ano(ni+2,vector<double>(nj+2,0.0)), //異常衝突周波数
        nu_ano_IAT(ni+2,vector<double>(nj+2,0.0)), //異常衝突周波数
        //--------------------------------

        //production rate
        //--------------------------------
        rate_ionize(ni+2,vector<double>(nj+2,0.0)),
        //--------------------------------

        //energy exchange
        //--------------------------------
        rate_eloss(ni+2,vector<double>(nj+2,0.0)), //電子エネルギー獲得量
        delossdeps(ni+2,vector<double>(nj+2,0.0)), //電子エネルギーロスの微係数
        //--------------------------------

        //electron mobility
        //--------------------------------
        mu_para(ni+2,vector<double>(nj+2,0.0)),
        mu_perp(ni+2,vector<double>(nj+2,0.0)),
        mu_min(ni+2,vector<double>(nj+2,0.0)),
        mu_xx(ni+2,vector<double>(nj+2,0.0)),
        mu_rr(ni+2,vector<double>(nj+2,0.0)),
        mu_xr(ni+2,vector<double>(nj+2,0.0)),
        mu_inv_xx(ni+2,vector<double>(nj+2,0.0)),
        mu_inv_rr(ni+2,vector<double>(nj+2,0.0)),
        mu_inv_xr(ni+2,vector<double>(nj+2,0.0)),
        //--------------------------------

        //temperature
        //--------------------------------
        Te(ni+2,vector<double>(nj+2,0.0)),
        Te_old(ni+2,vector<double>(nj+2,0.0)),
        //--------------------------------

        //field
        //--------------------------------
        phi(ni+2,vector<double>(nj+2,0.0)),
        phi_old(ni+2,vector<double>(nj+2,0.0)),

        Ap(ni+2,vector<double>(nj+2,0.0)),
        Bx(ni+2,vector<double>(nj+2,0.0)),
        Br(ni+2,vector<double>(nj+2,0.0)),

        E1p(ni+2,vector<complex<double> >(nj+2,complex<double>(0.0,0.0))),
        //--------------------------------

        //current
        //--------------------------------
        J1p(ni+2,vector<complex<double> >(nj+2,complex<double>(0.0,0.0))),
        //--------------------------------

        //permittivity
        //--------------------------------
        epsr(ni+2,vector<double>(nj+2,0.0)),
        //--------------------------------

        //electron Mach number
        //--------------------------------
        Mache(ni+2,vector<double>(nj+2,0.0)),
        //--------------------------------
        
        //power absorption density
        //--------------------------------
        Pabs(ni+2,vector<double>(nj+2,0.0)),
        //--------------------------------

        //judge boundary
        //--------------------------------
        jdgBnd_Ep(ni+2,vector<int>(nj+2,0.0)),
        jdgBnd_flc(ni+2,vector<int>(nj+2,0.0)),
        //--------------------------------

        //region index
        //--------------------------------
        i_flc_bl(ni+2,vector<int>(nj+2,0.0)),
        j_flc_bl(ni+2,vector<int>(nj+2,0.0)),
        //--------------------------------

        //index-convergion
        //--------------------------------
        kp(ni+2,vector<int>(nj+2,0.0))   //convergion of (i,j) → k (E-field)
        //--------------------------------
    {}
    //--------------------------------
};
//--------------------------------

//value at grid interface
//--------------------------------
struct GridInterfaceX{ 

    //flux
    //--------------------------------
    vector<vector<double> > rhoUex;
    vector<vector<double> > rhoUex_old;
    
    vector<vector<double> > rhoUix;
    vector<vector<double> > rhoUix_old;
    
    vector<vector<double> > rhoUnx;
    vector<vector<double> > rhoUnx_old;
    
    vector<vector<double> > rhoUmx;
    vector<vector<double> > rhoUmx_old;
    
    vector<vector<double> > nUex;
    vector<vector<double> > nUex_old;
    //--------------------------------

    //energy flux
    //--------------------------------
    vector<vector<double> > Gx;
    vector<vector<double> > Gx_old;
    //--------------------------------

    //wall flux
    //--------------------------------
    vector<vector<double> > rhoUex_wall;

    vector<vector<double> > rhoUix_wall;
    
    vector<vector<double> > rhoUnx_wall;
    
    vector<vector<double> > rhoUmx_wall;
    
    vector<vector<double> > nUex_wall;
    //--------------------------------

    //energy flux
    //--------------------------------
    vector<vector<double> > Gx_wall;
    //--------------------------------

    //bulk velocity
    //--------------------------------
    vector<vector<double> > Uix;
    vector<vector<double> > Uix_old;
    //--------------------------------

    //field
    //--------------------------------
    vector<vector<double> > Ex;

    vector<vector<complex<double> > > E1x;
    //--------------------------------

    //current
    //--------------------------------
    vector<vector<complex<double> > > J1x;
    //--------------------------------

    //judge boundary
    //--------------------------------
    vector<vector<int> > jdgBnd_Ex;

    vector<vector<int> > jdgBnd_flx;
    //--------------------------------

    //region index
    //--------------------------------
    vector<vector<int> > i_flx_bl;
    vector<vector<int> > j_flx_bl;
    //--------------------------------

    //index-convergion
    //--------------------------------
    vector<vector<int> > kx;   //convergion of (i,j) → k (E-field)
    //--------------------------------

    //constructor
    //--------------------------------
    explicit GridInterfaceX(int ni,int nj)
        :

        //flux
        //--------------------------------
        rhoUex(ni+2,vector<double>(nj+2,0.0)),
        rhoUex_old(ni+2,vector<double>(nj+2,0.0)),
        
        rhoUix(ni+2,vector<double>(nj+2,0.0)),
        rhoUix_old(ni+2,vector<double>(nj+2,0.0)),
        
        rhoUnx(ni+2,vector<double>(nj+2,0.0)),
        rhoUnx_old(ni+2,vector<double>(nj+2,0.0)),
        
        rhoUmx(ni+2,vector<double>(nj+2,0.0)),
        rhoUmx_old(ni+2,vector<double>(nj+2,0.0)),
        
        nUex(ni+2,vector<double>(nj+2,0.0)),
        nUex_old(ni+2,vector<double>(nj+2,0.0)),
        //--------------------------------
    
        //energy flux
        //--------------------------------
        Gx(ni+2,vector<double>(nj+2,0.0)),
        Gx_old(ni+2,vector<double>(nj+2,0.0)),
        //--------------------------------
    
        //wall flux
        //--------------------------------
        rhoUex_wall(ni+2,vector<double>(nj+2,0.0)),
    
        rhoUix_wall(ni+2,vector<double>(nj+2,0.0)),
       
        rhoUnx_wall(ni+2,vector<double>(nj+2,0.0)),
        
        rhoUmx_wall(ni+2,vector<double>(nj+2,0.0)),
        
        nUex_wall(ni+2,vector<double>(nj+2,0.0)),
        //--------------------------------
    
        //energy flux
        //--------------------------------
        Gx_wall(ni+2,vector<double>(nj+2,0.0)),
        //--------------------------------
    
        //bulk velocity
        //--------------------------------
        Uix(ni+2,vector<double>(nj+2,0.0)),
        Uix_old(ni+2,vector<double>(nj+2,0.0)),
        //--------------------------------
    
        //field
        //--------------------------------
        Ex(ni+2,vector<double>(nj+2,0.0)),

        E1x(ni+2,vector<complex<double> >(nj+2,complex<double>(0.0,0.0))),
        //--------------------------------
    
        //current
        //--------------------------------
        J1x(ni+2,vector<complex<double> >(nj+2,complex<double>(0.0,0.0))),
        //--------------------------------

        //judge boundary
        //--------------------------------
        jdgBnd_Ex(ni+2,vector<int>(nj+2,0.0)),
        jdgBnd_flx(ni+2,vector<int>(nj+2,0.0)),
        //--------------------------------

        //region index
        //--------------------------------
        i_flx_bl(ni+2,vector<int>(nj+2,0.0)),
        j_flx_bl(ni+2,vector<int>(nj+2,0.0)),
        //--------------------------------

        //index-convergion
        //--------------------------------
        kx(ni+2,vector<int>(nj+2,0.0))   //convergion of (i,j) → k (E-field)
        //--------------------------------
        
    {}
    //--------------------------------
};
//--------------------------------


//value at grid interface
//--------------------------------
struct GridInterfaceR{ 

    //flux
    //--------------------------------
    vector<vector<double> > rhoUer;
    vector<vector<double> > rhoUer_old;

    vector<vector<double> > rhoUir;
    vector<vector<double> > rhoUir_old;

    vector<vector<double> > rhoUnr;
    vector<vector<double> > rhoUnr_old;

    vector<vector<double> > rhoUmr;
    vector<vector<double> > rhoUmr_old;

    vector<vector<double> > nUer;
    vector<vector<double> > nUer_old;
    //--------------------------------

    //energy flux
    //--------------------------------
    vector<vector<double> > Gr;
    vector<vector<double> > Gr_old;
    //--------------------------------

    //wall flux
    //--------------------------------
    vector<vector<double> > rhoUer_wall;

    vector<vector<double> > rhoUir_wall;
 
    vector<vector<double> > rhoUnr_wall;
    
    vector<vector<double> > rhoUmr_wall;
    
    vector<vector<double> > nUer_wall;
    //--------------------------------

    //energy flux
    //--------------------------------
    vector<vector<double> > Gr_wall;
    //--------------------------------

    //bulk velocity
    //--------------------------------
    vector<vector<double> > Uir;
    vector<vector<double> > Uir_old;
    //--------------------------------

    //field
    //--------------------------------
    vector<vector<double> > Er;

    vector<vector<complex<double> > > E1r;
    //--------------------------------

    //current
    //--------------------------------
    vector<vector<complex<double> > > J1r;
    //--------------------------------

    //judge boundary
    //--------------------------------
    vector<vector<int> > jdgBnd_Er;

    vector<vector<int> > jdgBnd_flr;
    //--------------------------------

    //region index
    //--------------------------------
    vector<vector<int> > i_flr_bl;
    vector<vector<int> > j_flr_bl;
    //--------------------------------

    //index-convergion
    //--------------------------------
    vector<vector<int> > kr;   //convergion of (i,j) → k (E-field)
    //--------------------------------

    //constructor
    //--------------------------------
    explicit GridInterfaceR(int ni,int nj)
        :

        //flux
        //--------------------------------
        rhoUer(ni+2,vector<double>(nj+2,0.0)),
        rhoUer_old(ni+2,vector<double>(nj+2,0.0)),
    
        rhoUir(ni+2,vector<double>(nj+2,0.0)),
        rhoUir_old(ni+2,vector<double>(nj+2,0.0)),
    
        rhoUnr(ni+2,vector<double>(nj+2,0.0)),
        rhoUnr_old(ni+2,vector<double>(nj+2,0.0)),
    
        rhoUmr(ni+2,vector<double>(nj+2,0.0)),
        rhoUmr_old(ni+2,vector<double>(nj+2,0.0)),
    
        nUer(ni+2,vector<double>(nj+2,0.0)),
        nUer_old(ni+2,vector<double>(nj+2,0.0)),
        //--------------------------------
    
        //energy flux
        //--------------------------------
        Gr(ni+2,vector<double>(nj+2,0.0)),
        Gr_old(ni+2,vector<double>(nj+2,0.0)),
        //--------------------------------
    
        //wall flux
        //--------------------------------
        rhoUer_wall(ni+2,vector<double>(nj+2,0.0)),
    
        rhoUir_wall(ni+2,vector<double>(nj+2,0.0)),
    
        rhoUnr_wall(ni+2,vector<double>(nj+2,0.0)),
        
        rhoUmr_wall(ni+2,vector<double>(nj+2,0.0)),
    
        nUer_wall(ni+2,vector<double>(nj+2,0.0)),
        //--------------------------------
    
        //energy flux
        //--------------------------------
        Gr_wall(ni+2,vector<double>(nj+2,0.0)),
        //--------------------------------
    
        //bulk velocity
        //--------------------------------
        Uir(ni+2,vector<double>(nj+2,0.0)),
        Uir_old(ni+2,vector<double>(nj+2,0.0)),
        //--------------------------------
    
        //field
        //--------------------------------
        Er(ni+2,vector<double>(nj+2,0.0)),
    
        E1r(ni+2,vector<complex<double> >(nj+2,complex<double>(0.0,0.0))),
        //--------------------------------
    
        //current
        //--------------------------------
        J1r(ni+2,vector<complex<double> >(nj+2,complex<double>(0.0,0.0))),
        //--------------------------------

        //judge boundary
        //--------------------------------
        jdgBnd_Er(ni+2,vector<int>(nj+2,0.0)),
        jdgBnd_flr(ni+2,vector<int>(nj+2,0.0)),
        //--------------------------------

        //region index
        //--------------------------------
        i_flr_bl(ni+2,vector<int>(nj+2,0.0)),
        j_flr_bl(ni+2,vector<int>(nj+2,0.0)),
        //--------------------------------

        //index-convergion
        //--------------------------------
        kr(ni+2,vector<int>(nj+2,0.0))    //convergion of (i,j) → k (E-field)
        //--------------------------------
        
    {}
    //--------------------------------
};
//--------------------------------


//value at k-space
//--------------------------------
struct GridK{ 
    vector<int> ikx; //convergion of k→ i (Ex)
    vector<int> jkx; //convergion of k→ j (Ex)
    
    vector<int> ikr; //convergion of k→ i (Er)
    vector<int> jkr; //convergion of k→ j (Er)

    vector<int> ikp; //convergion of k→ i (Ep)
    vector<int> jkp; //convergion of k→ j (Ep)
};
//--------------------------------

//boundary setting for microwave
//--------------------------------
struct MicrowaveBC{ 

    vector<int> iBndWx; 
    vector<int> jBndWx;
    vector<int> sBndWx;

    vector<int> iBndEx;
    vector<int> jBndEx;
    vector<int> sBndEx;

    vector<int> iBndSx;
    vector<int> jBndSx;
    vector<int> sBndSx;

    vector<int> iBndNx;
    vector<int> jBndNx;
    vector<int> sBndNx;

    vector<int> iBndWr;
    vector<int> jBndWr;
    vector<int> sBndWr;

    vector<int> iBndEr;
    vector<int> jBndEr;
    vector<int> sBndEr;

    vector<int> iBndSr;
    vector<int> jBndSr;
    vector<int> sBndSr;

    vector<int> iBndNr;
    vector<int> jBndNr;
    vector<int> sBndNr;

    vector<int> iBndWp;
    vector<int> jBndWp;
    vector<int> sBndWp;

    vector<int> iBndEp;
    vector<int> jBndEp;
    vector<int> sBndEp;

    vector<int> iBndSp;
    vector<int> jBndSp;
    vector<int> sBndSp;

    vector<int> iBndNp;
    vector<int> jBndNp;
    vector<int> sBndNp;
};
//--------------------------------

//for input value of BOLSIG+
//--------------------------------
struct BolsigVec{
    
    vector<double> nu_elas_N; //弾性周波数/nn <σelas*vth>
    vector<double> nu_super_N; //超弾性周波数/nn <σsuper*vth>
    vector<double> nu_ionz_N; //イオン化周波数/nn <σionz*vth>
    vector<double> nu_ionzStep_N; //段階イオン化周波数/nn <σionz*vth>
    vector<double> nu_exc_N; //励起周波数/nn <σexc*vth>
    vector<double> nu_excReso_N; //共鳴励起周波数/nn <σexc*vth>
    vector<double> nu_excMeta_N; //準安定励起周波数/nn <σexc*vth>
    vector<double> nu_excStep_N; //段階励起周波数/nn <σexc*vth>
    vector<double> rate_eloss_n_N; //エネルギーロスレート/nn
    vector<double> rate_eloss_m_N; //エネルギーロスレート/nn
    vector<double> nui_m_N; //イオン運動量移動周波数/nn <σionz*vth>

};
//--------------------------------

//for SEE table
//--------------------------------
struct SeeVec{
    
    vector<double> coefEISEE_eb_table; //EI elastic-backscattering
    vector<double> coefEISEE_rd_table; //EI rediffusion
    vector<double> coefEISEE_ts_table; //EI true-secondary
};
//--------------------------------
