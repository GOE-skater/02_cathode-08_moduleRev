#pragma once
#include <string>
#include <iostream>
#include <vector>
#include <sys/stat.h>

#include "constants.hpp"
#include "params.hpp"
#include "miniFuncs.hpp"
#include "arrays.hpp"
#include "solverModule.hpp"

using namespace std;

//*****************************************************************
//**                                                             **
//**           class FluidModule                                 **
//**                                                             **
//*****************************************************************
class FluidModule{
    private:

    public:
        void update_transport_coef(Params &pm, GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr, BolsigVec &bo);
};

//*****************************************************************
//**                                                             **
//**           void update_transport_coef                        **
//**                                                             **
//*****************************************************************
void FluidModule::update_transport_coef(Params &pm,GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr, BolsigVec &bo){
    
    //calculate maximum electron density to avoid zero-division in Te calc.
    //------------------------------------
    double rhoe_max = 0.0;
    double ratio_rhoe = 0.001;
    for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
        for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
            for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                if(rhoe_max < gc.rhoe[i][j]){
                    rhoe_max = gc.rhoe[i][j];
                }
            }
        }
    }
    //------------------------------------

    //IAT model
    //------------------------------------
    for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
        for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
            for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                
                //BC - set adjacent flag (adjacent = 0, others = 1)
                //adjacent "cell" to the wall for setting zero-flux on the wall
                //------------------------------------
                //left
                double bLx_wall  = double(i!=gc.i_flc_bl[0][0] || iblock != 0); //z0
                bLx_wall = bLx_wall*double(i!=gc.i_flc_bl[2][0] || j>gc.j_flc_bl[0][0]-1); //z1
                bLx_wall = bLx_wall*double(i!=gc.i_flc_bl[3][0] || j>gc.j_flc_bl[2][0]-1); //z2
                bLx_wall = bLx_wall*double(i!=gc.i_flc_bl[4][0] || j<gc.j_flc_bl[1][0]); //z4
                double bLx  = bLx_wall;
                
                //right
                double bRx_wall  = double(i!=gc.i_flc_bl[1][1] || j<gc.j_flc_bl[1][0]); //z3
                bRx_wall = bRx_wall*double(i!=gc.i_flc_bl[4][1]); //z5
                double bRx  = bRx_wall;
                
                //lower
                double bLr_wall  = double(j!=gc.j_flc_bl[0][0] || i>gc.i_flc_bl[0][1]); //x0
                bLr_wall = bLr_wall*double(j!=gc.j_flc_bl[2][0] || i>gc.i_flc_bl[2][1]); //x2
                double bLr_cen = double(j!=gc.j_flc_bl[3][0]); //x6
                double bLr  = bLr_wall*bLr_cen;
                
                //upper
                double bRr_wall  = double(j!=gc.j_flc_bl[0][1] || iblock > 2); //x1
                bRr_wall = bRr_wall*double(j!=gc.j_flc_bl[3][1] || i<gc.i_flc_bl[1][1]+1 || i>gc.i_flc_bl[3][1]); //x4
                double bRr_open  = (j!=gc.j_flc_bl[4][1] || iblock != 4); //x5
                double bRr = bRr_wall*bRr_open;
                //------------------------------------

                double ut = sqrt(2.0*ph::Boltz*gc.Te[i][j]/pm.masse);

                //radial position at upper and lower cell interface
                //------------------------------------
                double rL = (gc.r[j]+gc.r[j-1])/2.0;
                double rR = (gc.r[j]+gc.r[j+1])/2.0;
                //------------------------------------

                //flux at cell center
                //------------------------------------
                double rhoUex_tmp = (gx.rhoUex[i+1][j]*bRx_wall + gx.rhoUex[i][j]*bLx_wall)/2.0;
                double rhoUer_tmp = (rR*gr.rhoUer[i][j+1]*bRr_wall + rL*gr.rhoUer[i][j]*bLr_wall)/(2.0*gc.r[j]);
                //double rhoUex_tmp = (rhoUex[i][j]   + rhoUex[i+1][j]  )/2.0;
                //double rhoUer_tmp = (rL*rhoUer[i][j]   + rR*rhoUer[i][j+1]  )/(2.0*r[j]);
                //------------------------------------
                //velocity at cell center
                //------------------------------------
                double Uex_tmp = rhoUex_tmp/(gc.rhoe[i][j] + rhoe_max*ratio_rhoe);
                double Uer_tmp = rhoUer_tmp/(gc.rhoe[i][j] + rhoe_max*ratio_rhoe);

                double Uix_tmp = (gx.Uix[i][j] + gx.Uix[i+1][j])/2.0;
                double Uir_tmp = (rL*gr.Uir[i][j] + rR*gr.Uir[i][j+1])/(2.0*gc.r[j]);    
                //------------------------------------

                gc.Mache[i][j] = sqrt(pow(Uex_tmp,2) + pow(Uer_tmp,2))/(ut + 100.0);


                /*
                //Model of Mikellides 2007
                //------------------------------------
                //Mache[i][j] = sqrt(pow(Uex_tmp-Uix_tmp,2) + pow(Uer_tmp-Uir_tmp,2))/(ut + 100.0);
                double fM = 1.0;
                double alpha = 0.171;
                double b = 1.25;

                if(Mache[i][j]<1.3){
                    fM = alpha*pow(Mache[i][j],b);
                }else{
                    fM = 1.0;
                }
                double omega_pe =sqrt(rhoe[i][j]*e0*e0/(masse*eps0));
                nu_ano_IAT[i][j] = pow(masse/massi,1.0/3.0)*fM*sqrt(2.0)*omega_pe;
                //------------------------------------
                */
    
                
                //Model of Sagdeev only
                //------------------------------------
                double cs = sqrt(2.0*ph::Boltz*(gc.Te[i][j] + pm.Ti)/pm.massi);
                double omega_pi =sqrt(gc.rhoi[i][j]*ph::e0*ph::e0/(pm.massi*ph::eps0));
                double ud = sqrt(pow(Uex_tmp-Uix_tmp,2) + pow(Uer_tmp-Uir_tmp,2));
                
                gc.nu_ano_IAT[i][j] = 0.01*omega_pi*ud/cs*gc.Te[i][j]/pm.Ti;
                //------------------------------------

                //std::cout << Mache[i][j] << "," <<nu_ano_IAT[i][j] << ","<<0.01*omega_pi*ud/cs*Te[i][j]/Ti << std::endl;
            }
        }
    }

    //collision frequencies    
    for (int i=0;i<pm.ni+2;i++){
        for (int j=0;j<pm.nj+2;j++){

            //****************** 電子衝突計算 ******************
            gc.nu_elas[i][j]     = gc.rhon[i][j]*func_Te(bo.nu_elas_N,     gc.Te[i][j],pm.dTe); //弾性衝突
            gc.nu_exc[i][j]      = gc.rhon[i][j]*func_Te(bo.nu_exc_N,      gc.Te[i][j],pm.dTe); //励起衝突
            gc.nu_excReso[i][j]  = gc.rhon[i][j]*func_Te(bo.nu_excReso_N,  gc.Te[i][j],pm.dTe); //励起(共鳴)衝突
            gc.nu_excMeta[i][j]  = gc.rhon[i][j]*func_Te(bo.nu_excMeta_N,  gc.Te[i][j],pm.dTe); //励起(準安定)衝突
            gc.nu_ionz[i][j]     = gc.rhon[i][j]*func_Te(bo.nu_ionz_N,     gc.Te[i][j],pm.dTe); //イオン化衝突
            gc.nu_ionzStep[i][j] = gc.rhom[i][j]*func_Te(bo.nu_ionzStep_N, gc.Te[i][j],pm.dTe); //段階イオン化衝突
            gc.nu_excStep[i][j]  = gc.rhom[i][j]*func_Te(bo.nu_excStep_N,  gc.Te[i][j],pm.dTe); //段階励起衝突
            gc.nu_super[i][j]    = gc.rhom[i][j]*func_Te(bo.nu_super_N,    gc.Te[i][j],pm.dTe); //超弾性衝突
            
            //運動量移動周波数
            double nu_m_real = gc.rhon[i][j]*func_Te(bo.nu_elas_N,gc.Te[i][j],pm.dTe);
            gc.nu_m[i][j] = nu_m_real
                + double(pm.icon_Bohm)*gc.nu_ano[i][j] 
                + double(pm.icon_Sagdeev)*gc.nu_ano_IAT[i][j];
            gc.nu_m1[i][j] = (nu_m_real + pm.nu_eff); //マイクロ波用

            //イオン生成レート
            gc.rate_ionize[i][j] = gc.rhoe[i][j]*(gc.nu_ionz[i][j] + gc.nu_ionzStep[i][j]);
            //****************** イオン衝突計算 ******************
            double vv = 0.0;
            if(i!=0 && i!=pm.ni+1 && j!=0 && j!=pm.nj+1){
                double rL = (gc.r[j]+gc.r[j-1])/2.0;
                double rR = (gc.r[j]+gc.r[j+1])/2.0;

                double Uix_tmp = (gx.Uix[i][j] + gx.Uix[i+1][j])/2.0;
                double Uir_tmp = (rL*gr.Uir[i][j] + rR*gr.Uir[i][j+1])/(2.0*gc.r[j]);
                double Uip_tmp = gc.Uip[i][j];
                vv = Uix_tmp*Uix_tmp + Uir_tmp*Uir_tmp + Uip_tmp*Uip_tmp;
            }
            double Ei_tmp = pm.massi*vv;
            
            //イオン衝突周波数
            gc.nui_m[i][j] = gc.rhon[i][j]*func_Ei(bo.nui_m_N, Ei_tmp, pm.dEi, pm.massi); 

            //****************** 電子エネルギー交換 ******************
            gc.rate_eloss[i][j] = gc.rhon[i][j]*func_Te(bo.rate_eloss_n_N,gc.Te[i][j],pm.dTe) 
                + gc.rhom[i][j]*func_Te(bo.rate_eloss_m_N,gc.Te[i][j],pm.dTe);
            
            //移動度の微係数の計算
            gc.delossdeps[i][j] = 2.0/(3.0*ph::Boltz)*gc.rhon[i][j]*func_dTe(bo.rate_eloss_n_N,gc.Te[i][j],pm.dTe)
                + 2.0/(3.0*ph::Boltz)*gc.rhom[i][j]*func_dTe(bo.rate_eloss_m_N,gc.Te[i][j],pm.dTe);
            //****************** 電子移動度計算 ******************
            gc.mu_para[i][j] = ph::e0/(pm.masse*gc.nu_m[i][j]+1e-100);
            double Bx_tmp = gc.Bx[i][j];
            double Br_tmp = gc.Br[i][j];
            double nu_m_tmp = gc.nu_m[i][j];

            double Bmag = sqrt(Bx_tmp*Bx_tmp + Br_tmp*Br_tmp);
            double Omega = ph::e0*Bmag/(pm.masse*nu_m_tmp+1e-100);
            double Omega_x = ph::e0*Bx_tmp/(pm.masse*nu_m_tmp+1e-100);
            double Omega_r = ph::e0*Br_tmp/(pm.masse*nu_m_tmp+1e-100);
            double deno = 1.0 + Omega*Omega;
            double deno1 = 2.0 + Omega*Omega;
    
            gc.mu_xx[i][j] = gc.mu_para[i][j]*(1.0 + Omega_x*Omega_x)/deno;
            gc.mu_rr[i][j] = gc.mu_para[i][j]*(1.0 + Omega_r*Omega_r)/deno;
            gc.mu_xr[i][j] = gc.mu_para[i][j]*(Omega_x*Omega_r      )/deno;
            gc.mu_perp[i][j] = gc.mu_para[i][j]/deno;
            gc.mu_min[i][j] = gc.mu_para[i][j]*2.0/deno;
        }
    }
}
