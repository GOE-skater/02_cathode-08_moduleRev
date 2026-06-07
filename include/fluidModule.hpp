#pragma once
#include <string>
#include <iostream>
#include <vector>
#include <sys/stat.h>

#include "constants.hpp"
#include "params.hpp"
#include "miniFuncs.hpp"
#include "arrays.hpp"
#include "solverFuncs.hpp"

using namespace std;

//*****************************************************************
//**                                                             **
//**           class FluidModule                                 **
//**                                                             **
//*****************************************************************
class FluidModule{
    private:

    public:
        void check_CFL(Params &pm,GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr);
        
        void solve_Uix_tmp(Params &pm,GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr);
        void solve_Uir_tmp(Params &pm,GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr);
        void solve_Uip_tmp(Params &pm,GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr);
        void solve_rhoi_constTe(Params &pm,GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr);
        void correct_Ui_constTe(Params &pm,GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr);

        void solve_phi_couple_wdTe_wSEE_PC(Params &pm,GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr, SeeVec &se);
        void solve_rhoe_wdTe_wSEE_PC(Params &pm,GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr, SeeVec &se);
        void solve_Te_wdTe_wSEE_PC(Params &pm,GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr, SeeVec &se);
        void update_rhon(Params &pm,GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr); 
        void update_rhom(Params &pm,GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr);

        void update_transport_coef(Params &pm, GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr, BolsigVec &bo);
};

//*****************************************************************
//**                                                             **
//**           void check_CFL                                    **
//**                                                             **
//*****************************************************************
void FluidModule::check_CFL(Params &pm,GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr){
    

    double rhoe_max = 0.0;
    for (int i=0;i<pm.ni+2;i++){
        for (int j=0;j<pm.nj+2;j++){
            if(rhoe_max < gc.rhoe[i][j]){
                rhoe_max = gc.rhoe[i][j];
            }
        }
    }

    double mu_max = 0.0;
    for (int i=0;i<pm.ni+2;i++){
        for (int j=0;j<pm.nj+2;j++){
            if(mu_max < gc.mu_para[i][j]){
                mu_max = gc.mu_para[i][j];
            }
        }
    }

    double Te_max = 0.0;
    for (int i=0;i<pm.ni+2;i++){
        for (int j=0;j<pm.nj+2;j++){
            if(Te_max < gc.Te[i][j]){
                Te_max = gc.Te[i][j];
            }
        }
    }

    double ue_max = 0.0;
    for (int i=1;i<pm.ni+1;i++){
        for (int j=1;j<pm.nj+1;j++){
            double rL = (gc.r[j]+gc.r[j-1])/2.0;
            double rR = (gc.r[j]+gc.r[j+1])/2.0;

            double uex_tmp = (gx.rhoUex[i+1][j] + gx.rhoUex[i][j])/2.0/(gc.rhoe[i][j]+1e-100);
            double uer_tmp = (rR*gr.rhoUer[i][j+1] + rL*gr.rhoUer[i][j])/(2.0*gc.r[j])/(gc.rhoe[i][j]+1e-100);
            
            double ue_tmp = sqrt(uex_tmp*uex_tmp + uer_tmp*uer_tmp);
            if(ue_max < ue_tmp){
                ue_max = ue_tmp;
            }
        }
    }

    double vth = sqrt(8.0*ph::Boltz*Te_max/(M_PI*pm.masse));
    ue_max = ue_max + vth;

    double dt_min = ph::eps0/(ph::e0*rhoe_max*mu_max);

    if(pm.flag_adp_dt == 1){ //アダプティブdt (dtを計算)
        pm.dt = dt_min*pm.CFL;
    }else{ //固定dt (CFLを計算)
        pm.CFL = pm.dt/dt_min;
    }


} // end check_CFL


//*****************************************************************
//**                                                             **
//**           void solve_Uix_tmp                                **
//**                                                             **
//*****************************************************************
void FluidModule::solve_Uix_tmp(Params &pm,GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr){

    vector<vector<double> > aW(pm.ni+2,vector<double>(pm.nj+2,0.0));
    vector<vector<double> > aE(pm.ni+2,vector<double>(pm.nj+2,0.0));
    vector<vector<double> > aS(pm.ni+2,vector<double>(pm.nj+2,0.0));
    vector<vector<double> > aN(pm.ni+2,vector<double>(pm.nj+2,0.0));
    vector<vector<double> > aP(pm.ni+2,vector<double>(pm.nj+2,0.0));
    vector<vector<double> > b (pm.ni+2,vector<double>(pm.nj+2,0.0));

    //time step for ion
    //------------------------------------
    double dt_i = pm.dt*double(pm.ndt_i);
    //------------------------------------
    
    //Reserve the past values
    //------------------------------------
    for (int i=0;i<=pm.ni+1;i++){
        for (int j=0;j<=pm.nj+1;j++){
            gx.Uix_old[i][j]  = gx.Uix[i][j];
            gr.Uir_old[i][j]  = gr.Uir[i][j];
            gc.Uip_old[i][j]  = gc.Uip[i][j];
            gc.rhoi_old[i][j] = gc.rhoi[i][j];
        }
    }
    //------------------------------------

    //Make coefficient for implicit x-directional ion momentum conservation for 5 point stencil
    //=================================================================================
    for (int iblock=0;iblock<pm.n_bl-1;iblock++){
        for (int i=gx.i_flx_bl[iblock][0];i<=gx.i_flx_bl[iblock][1];i++){
            for (int j=gx.j_flx_bl[iblock][0];j<=gx.j_flx_bl[iblock][1];j++){
                //radial position at upper and lower cell interface
                //------------------------------------
                double rL = (gc.r[j]+gc.r[j-1])/2.0;
                double rR = (gc.r[j]+gc.r[j+1])/2.0;
                //------------------------------------

                //normalized radial position at upper and lower cell interface
                //------------------------------------
                double qL = rL/gc.r[j];
                double qR = rR/gc.r[j];
                //------------------------------------
                
                //values for upwind discretization
                //*****************************************************
                //velocity at cell (Uix control volume) interface
                //------------------------------------
                double Uix_Lx = (gx.Uix_old[i][j] + gx.Uix_old[i-1][j])/2.0;
                double Uix_Rx = (gx.Uix_old[i][j] + gx.Uix_old[i+1][j])/2.0;
                double Uir_Lr = (gr.Uir_old[i][j] + gr.Uir_old[i-1][j])/2.0;
                double Uir_Rr = (gr.Uir_old[i][j+1] + gr.Uir_old[i-1][j+1])/2.0;
                //------------------------------------

                //values at cell (Uix control volume) center
                //------------------------------------
                double rhoi_tmp = (gc.rhoi[i][j] + gc.rhoi[i-1][j])/2.0;
                double rate_ionize_tmp = (gc.rate_ionize[i][j] + gc.rate_ionize[i-1][j])/2.0;
                double nui_m_tmp = (gc.nui_m[i][j] + gc.nui_m[i-1][j])/2.0;
                //------------------------------------
                
                //density gradient at cell (Uix control volume) center
                //------------------------------------
                double dndx = (gc.rhoi[i][j]-gc.rhoi[i-1][j])/pm.dx;
                //------------------------------------

                //Lorentz-force at cell (Uix control volume) center
                //------------------------------------
                double Uip_tmp = (gc.Uip_old[i][j] + gc.Uip_old[i-1][j])/2.0;
                double Br_tmp = (gc.Br[i][j]+gc.Br[i-1][j])/2.0;
                double lorentz = - ph::e0*Uip_tmp*Br_tmp;
                //------------------------------------
                //*****************************************************

                //coefficient for ion momentum conservation for 5 point stencil
                //*****************************************************
                
                //vertical and horizental 5 stencil
                //------------------------------------
                aE[i][j] = -(Uix_Rx- fabs(Uix_Rx))/(2.0*pm.dx);
                aW[i][j] =  (Uix_Lx + fabs(Uix_Lx))/(2.0*pm.dx);
                aN[i][j] = -qR*(Uir_Rr - fabs(Uir_Rr))/(2.0*pm.dr);
                aS[i][j] =  qL*(Uir_Lr + fabs(Uir_Lr))/(2.0*pm.dr);
                aP[i][j] = 1.0/dt_i  + rate_ionize_tmp/(rhoi_tmp + 1e-100) + nui_m_tmp
                    + aE[i][j] + aW[i][j] + aN[i][j] + aS[i][j];
                //------------------------------------

                //RHS
                //------------------------------------
                b[i][j] = gx.Uix_old[i][j]/dt_i + (ph::e0*gx.Ex[i][j] + lorentz)/pm.massi
                    - ph::Boltz*pm.Ti/(pm.massi*rhoi_tmp+1e-100)*dndx;
                //------------------------------------
                //*****************************************************
            }
        }
    }

    //Set BCs for coefficients
    //=====================================================================
    //left-wall-normal-BC (z0) energy-conservation
    //------------------------------------
    for (int j=gx.j_flx_bl[0][0];j<=gx.j_flx_bl[0][1];j++){
        int i=gx.i_flx_bl[0][0];
        //calculate velocity at wall
        //------------------------------------
        double delta_phi = (gc.phi[i][j] + gc.phi[i-1][j])/2.0;
        double Ui_wall = -sqrt(pow(gx.Uix_old[i][j],2) + fmax(2.0*ph::e0*delta_phi/pm.massi,0.0));
        //------------------------------------
        
        //modify coefficient
        //------------------------------------
        b[i][j] = b[i][j] + aW[i][j]*Ui_wall;
        aW[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------
    
    //left-wall-normal-BC (z1) energy-conservation
    //------------------------------------
    for (int j=gx.j_flx_bl[2][0];j<=gx.j_flx_bl[0][0]-1;j++){
        int i=gx.i_flx_bl[2][0];
        //calculate velocity at wall
        //------------------------------------
        double delta_phi = (gc.phi[i][j] + gc.phi[i-1][j])/2.0;
        double Ui_wall = -sqrt(pow(gx.Uix_old[i][j],2) + fmax(2.0*ph::e0*delta_phi/pm.massi,0.0));
        //------------------------------------
        
        //modify coefficient
        //------------------------------------
        b[i][j] = b[i][j] + aW[i][j]*Ui_wall;
        aW[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //left-wall-normal-BC (z2) energy-conservation
    //------------------------------------
    for (int j=gx.j_flx_bl[3][0];j<=gx.j_flx_bl[2][0]-1;j++){
        int i=gx.i_flx_bl[3][0];
        //calculate velocity at wall
        //------------------------------------
        double delta_phi = (gc.phi[i][j] + gc.phi[i-1][j])/2.0;
        double Ui_wall = -sqrt(pow(gx.Uix_old[i][j],2) + fmax(2.0*ph::e0*delta_phi/pm.massi,0.0));
        //------------------------------------
        
        //modify coefficient
        //------------------------------------
        b[i][j] = b[i][j] + aW[i][j]*Ui_wall;
        aW[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //left-wall-normal-BC (z4) energy-conservation
    //------------------------------------
    for (int j=gx.j_flx_bl[1][0];j<=gx.j_flx_bl[4][1];j++){
        int i=gx.i_flx_bl[4][0];
        //calculate velocity at wall
        //------------------------------------
        double delta_phi = (gc.phi[i][j] + gc.phi[i-1][j])/2.0;
        double Ui_wall = -sqrt(pow(gx.Uix_old[i][j],2) + fmax(2.0*ph::e0*delta_phi/pm.massi,0.0));
        //------------------------------------
        
        //modify coefficient
        //------------------------------------
        b[i][j] = b[i][j] + aW[i][j]*Ui_wall;
        aW[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //right-wall-normal-BC (z3) energy-conservation
    //------------------------------------
    for (int j=gx.j_flx_bl[1][0];j<=gx.j_flx_bl[1][1];j++){
        int i=gx.i_flx_bl[1][1];
        //calculate velocity at wall
        //------------------------------------
        double delta_phi = (gc.phi[i][j] + gc.phi[i-1][j])/2.0;
        double Ui_wall = sqrt(pow(gx.Uix_old[i][j],2) + fmax(2.0*ph::e0*delta_phi/pm.massi,0.0));
        //------------------------------------

        //modify coefficient
        //------------------------------------
        b[i][j] = b[i][j] + aE[i][j]*Ui_wall;
        aE[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //right-wall-normal-BC (z5) energy-conservation
    //------------------------------------
    for (int j=gx.j_flx_bl[4][0];j<=gx.j_flx_bl[4][1];j++){
        int i=gx.i_flx_bl[4][1];
        //calculate velocity at wall
        //------------------------------------
        double delta_phi = (gc.phi[i][j] + gc.phi[i-1][j])/2.0 - pm.V_bias;
        double Ui_wall = sqrt(pow(gx.Uix_old[i][j],2) + fmax(2.0*ph::e0*delta_phi/pm.massi,0.0));
        //------------------------------------

        //modify coefficient
        //------------------------------------
        b[i][j] = b[i][j] + aE[i][j]*Ui_wall;
        aE[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //lower-wall-tangential-BC (x0) zero-flux
    //------------------------------------
    for (int i=gx.i_flx_bl[0][0];i<=gx.i_flx_bl[0][1];i++){
        int j=gx.j_flx_bl[0][0];
        //modify coefficient
        //------------------------------------
        aS[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //lower-wall-tangential-BC (x2) zero-flux
    //------------------------------------
    for (int i=gx.i_flx_bl[2][0];i<=gx.i_flx_bl[2][1];i++){
        int j=gx.j_flx_bl[2][0];
        //modify coefficient
        //------------------------------------
        aS[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //upper-wall-tangential-BC (x1) zero-flux
    //------------------------------------
    for (int i=gx.i_flx_bl[0][0];i<=gx.i_flx_bl[1][1];i++){
        int j=gx.j_flx_bl[0][1];
        //modify coefficient
        //------------------------------------
        aN[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------
    
    //upper-wall-tangential-BC (x4) zero-flux
    //------------------------------------
    for (int i=gx.i_flx_bl[1][1]+1;i<=gx.i_flx_bl[3][1];i++){
        int j=gx.j_flx_bl[3][1];
        //modify coefficient
        //------------------------------------
        aN[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //upper-open-tangential-BC (x5) zero-flux
    //------------------------------------
    for (int i=gx.i_flx_bl[4][0];i<=gx.i_flx_bl[4][1];i++){
        int j=gx.j_flx_bl[4][1];
        //modify coefficient
        //------------------------------------
        aN[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------
    //=====================================================================

    //Solver
    //=====================================================================
    if(pm.flag_iter_Ui == 0){
        //SOR
        //------------------------------------
        double alpha = 1.9;
        solver_SOR(aP,aE,aW,aN,aS,b,gx.i_flx_bl,gx.j_flx_bl,pm.n_bl-1,alpha,pm.maxITR_SOR_Ui,pm.error_cnv_SOR_Ui,0,gx.Uix);
        //------------------------------------
    }else{
        //Semi-coarse Multi-Grid by Hypre
        //------------------------------------
        solver_SMG(aP,aE,aW,aN,aS,b,gx.i_flx_bl,gx.j_flx_bl,pm.n_bl-1,pm.maxITR_SOR_Ui,pm.error_cnv_SOR_Ui,0,gx.Uix);
        //------------------------------------

        //for (int iblock=0;iblock<pm.n_bl-1;iblock++){
        //    std::cout << iblock << "("<<gx.i_flx_bl[iblock][0] << ","<<gx.j_flx_bl[iblock][0] <<") -> ("<<gx.i_flx_bl[iblock][1]<<","<<gx.j_flx_bl[iblock][1]<<")"<<std::endl;
        //}
    }
    //=====================================================================

    //Reflect BCs for ghost cell
    //=====================================================================
    //left-wall-normal-BC (chamber) energy-conservation
    //****************************************
    //left-wall-normal-BC (z0) energy-conservation
    //------------------------------------
    for (int j=gx.j_flx_bl[0][0];j<=gx.j_flx_bl[0][1];j++){
        int i=gx.i_flx_bl[0][0];
        //calculate velocity at wall
        //------------------------------------
        double delta_phi = (gc.phi[i][j] + gc.phi[i-1][j])/2.0;
        double Ui_wall = -sqrt(pow(gx.Uix_old[i][j],2) + fmax(2.0*ph::e0*delta_phi/pm.massi,0.0));
        //------------------------------------

        //reflect BC
        //------------------------------------
        gx.Uix[i-1][j] = Ui_wall;
        //------------------------------------
    }
    //------------------------------------
    
    //left-wall-normal-BC (z1) energy-conservation
    //------------------------------------
    for (int j=gx.j_flx_bl[2][0];j<=gx.j_flx_bl[0][0]-1;j++){
        int i=gx.i_flx_bl[2][0];
        //calculate velocity at wall
        //------------------------------------
        double delta_phi = (gc.phi[i][j] + gc.phi[i-1][j])/2.0;
        double Ui_wall = -sqrt(pow(gx.Uix_old[i][j],2) + fmax(2.0*ph::e0*delta_phi/pm.massi,0.0));
        //------------------------------------

        //reflect BC
        //------------------------------------
        gx.Uix[i-1][j] = Ui_wall;
        //------------------------------------
    }
    //------------------------------------

    //left-wall-normal-BC (z2) energy-conservation
    //------------------------------------
    for (int j=gx.j_flx_bl[3][0];j<=gx.j_flx_bl[2][0]-1;j++){
        int i=gx.i_flx_bl[3][0];
        //calculate velocity at wall
        //------------------------------------
        double delta_phi = (gc.phi[i][j] + gc.phi[i-1][j])/2.0;
        double Ui_wall = -sqrt(pow(gx.Uix_old[i][j],2) + fmax(2.0*ph::e0*delta_phi/pm.massi,0.0));
        //------------------------------------

        //reflect BC
        //------------------------------------
        gx.Uix[i-1][j] = Ui_wall;
        //------------------------------------
    }
    //------------------------------------

    //left-wall-normal-BC (z4) energy-conservation
    //------------------------------------
    for (int j=gx.j_flx_bl[1][0];j<=gx.j_flx_bl[4][1];j++){
        int i=gx.i_flx_bl[4][0];
        //calculate velocity at wall
        //------------------------------------
        double delta_phi = (gc.phi[i][j] + gc.phi[i-1][j])/2.0;
        double Ui_wall = -sqrt(pow(gx.Uix_old[i][j],2) + fmax(2.0*ph::e0*delta_phi/pm.massi,0.0));
        //------------------------------------

        //reflect BC
        //------------------------------------
        gx.Uix[i-1][j] = Ui_wall;
        //------------------------------------
    }
    //------------------------------------

    //right-wall-normal-BC (z3) energy-conservation
    //------------------------------------
    for (int j=gx.j_flx_bl[1][0];j<=gx.j_flx_bl[1][1];j++){
        int i=gx.i_flx_bl[1][1];
        //calculate velocity at wall
        //------------------------------------
        double delta_phi = (gc.phi[i][j] + gc.phi[i-1][j])/2.0;
        double Ui_wall = sqrt(pow(gx.Uix_old[i][j],2) + fmax(2.0*ph::e0*delta_phi/pm.massi,0.0));
        //------------------------------------
        
        //reflect BC
        //------------------------------------
        gx.Uix[i+1][j] = Ui_wall;
        //------------------------------------
    }
    //------------------------------------

    //right-wall-normal-BC (z5) energy-conservation
    //------------------------------------
    for (int j=gx.j_flx_bl[4][0];j<=gx.j_flx_bl[4][1];j++){
        int i=gx.i_flx_bl[4][1];
        //calculate velocity at wall
        //------------------------------------
        double delta_phi = (gc.phi[i][j] + gc.phi[i-1][j])/2.0 - pm.V_bias;
        double Ui_wall = sqrt(pow(gx.Uix_old[i][j],2) + fmax(2.0*ph::e0*delta_phi/pm.massi,0.0));
        //------------------------------------
        
        //reflect BC
        //------------------------------------
        gx.Uix[i+1][j] = Ui_wall;
        //------------------------------------
    }
    //------------------------------------

    //lower-wall-tangential-BC (x0) zero-flux
    //  [Attention] convex boundary is specified once by normal-BC
    //------------------------------------
    for (int i=gx.i_flx_bl[0][0];i<=gx.i_flx_bl[0][1]-1;i++){
        int j=gx.j_flx_bl[0][0];
        //reflect BC
        //------------------------------------
        gx.Uix[i][j-1] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //lower-wall-tangential-BC (x2) zero-flux
    //  [Attention] convex boundary is specified once by normal-BC
    //------------------------------------
    for (int i=gx.i_flx_bl[2][0];i<=gx.i_flx_bl[2][1]-1;i++){
        int j=gx.j_flx_bl[2][0];
        //reflect BC
        //------------------------------------
        gx.Uix[i][j-1] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //upper-wall-tangential-BC (x1) zero-flux
    //------------------------------------
    for (int i=gx.i_flx_bl[0][0];i<=gx.i_flx_bl[1][1];i++){
        int j=gx.j_flx_bl[0][1];
        //reflect BC
        //------------------------------------
        gx.Uix[i][j-1] = 0.0;
        //------------------------------------
    }
    //------------------------------------
    
    //upper-wall-tangential-BC (x4) zero-flux
    //  [Attention] convex boundary is specified once by normal-BC
    //------------------------------------
    for (int i=gx.i_flx_bl[1][1]+2;i<=gx.i_flx_bl[3][1]-1;i++){
        int j=gx.j_flx_bl[3][1];
        //reflect BC
        //------------------------------------
        gx.Uix[i][j+1] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //upper-open-tangential-BC (x5) zero-flux
    //------------------------------------
    for (int i=gx.i_flx_bl[4][0];i<=gx.i_flx_bl[4][1];i++){
        int j=gx.j_flx_bl[4][1];
        //reflect BC
        //------------------------------------
        gx.Uix[i][j+1] = 0.0;
        //------------------------------------
    }
    //------------------------------------
    //=====================================================================

    //check coefficients
    if(0==1){
        std::ofstream outputfile1("results/Uix_coef0.csv");
        //std::ofstream outputfile1_1(char1+"_tmp"+char2+char_csv);
        //outputfile1<<"x,rho,Ui,rhoUix,Ue,rhoUex,rhoUex_E,rhoUex_D,E,phi,psi,rho_th,U_th,E_th,phi_th,psi_th,rate_ionize,divi,dive1,dive2,dive3, div_poisson, LHS, RHS,zero" << std::endl;
        //outputfile1<<"x,rhon,rho,rhoe(test),Uix,rhoUix,Uex,Uey,rhoUex,rhoUey,Te,heat_flux,Bz,E,phi,nu_m,nu_en,nu_ei,nu_ionz,nu_exc,nu_wall,nu_ano,mue_parae,mue_perp,Halle,jd,Id,rate_eloss,rate_ionize,divi,dive,divn,engy_LHS,engy_LHS1,engy_LHS2,engy_LHS3,engy_LHS4,engy_RHS,engy_RHS1,engy_RHS2,engy_RHS3,zero" << std::endl;
        outputfile1<<"i,j,x,r,aP,aW,aE,aN,aS,b,zero" << std::endl;
        //outputfile1_1 << "i,j,x,r,rhom,rhoUmx,rhoUmr,nabla_rhoUm,rhon,rhoUnx,rhoUnr,nabla_rhoUn,zero" << std::endl;
        
        for(int i=1;i<=pm.ni;i++){
            for(int j=1;j<=pm.nj;j++){
                outputfile1 << i << ","<< j << "," << gc.x[i]<< ","<< gc.r[j]
                    << "," << aP[i][j] 
                    << "," << aW[i][j] 
                    << "," << aE[i][j] 
                    << "," << aN[i][j] 
                    << "," << aS[i][j] 
                    << "," << b[i][j]  
                    << "," << 0.0
                    <<std::endl;
            }
        }
    }

}

//*****************************************************************
//**                                                             **
//**           void solve_Uir_tmp                                **
//**                                                             **
//*****************************************************************
void FluidModule::solve_Uir_tmp(Params &pm,GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr){

    vector<vector<double> > aW(pm.ni+2,vector<double>(pm.nj+2,0.0)); //係数
    vector<vector<double> > aE(pm.ni+2,vector<double>(pm.nj+2,0.0)); //係数
    vector<vector<double> > aS(pm.ni+2,vector<double>(pm.nj+2,0.0)); //係数
    vector<vector<double> > aN(pm.ni+2,vector<double>(pm.nj+2,0.0)); //係数
    vector<vector<double> > aP(pm.ni+2,vector<double>(pm.nj+2,0.0)); //係数
    vector<vector<double> > b (pm.ni+2,vector<double>(pm.nj+2,0.0)); //ソース

    //time step for ion
    //------------------------------------
    double dt_i = pm.dt*double(pm.ndt_i);
    //------------------------------------

    //Make coefficient for implicit r-directional ion momentum conservation for 5 point stencil
    //------------------------------------
    for (int iblock=0;iblock<pm.n_bl-1;iblock++){
        for (int i=gr.i_flr_bl[iblock][0];i<=gr.i_flr_bl[iblock][1];i++){
            for (int j=gr.j_flr_bl[iblock][0];j<=gr.j_flr_bl[iblock][1];j++){
                
                //radial position at upper and lower cell interface
                //------------------------------------
                double r_tmp = (gc.r[j]+gc.r[j-1])/2.0;
                double rL = (gc.r[j-1]+gc.r[j-2])/2.0;
                double rR = (gc.r[j+1]+gc.r[j])/2.0;
                //------------------------------------

                //values for upwind discretization
                //*****************************************************
                //velocity at cell (Uir control volume) interface
                //------------------------------------
                double Uix_Lx = (gc.r[j]*gx.Uix_old[i][j]   + gc.r[j-1]* gx.Uix_old[i][j-1]  )/(2.0*r_tmp);
                double Uix_Rx = (gc.r[j]*gx.Uix_old[i+1][j] + gc.r[j-1]* gx.Uix_old[i+1][j-1])/(2.0*r_tmp);
                double Uir_Lr = (r_tmp*gr.Uir_old[i][j] + rL*gr.Uir_old[i][j-1])/(2.0*gc.r[j-1]);
                double Uir_Rr = (r_tmp*gr.Uir_old[i][j] + rR*gr.Uir_old[i][j+1])/(2.0*gc.r[j]);
                //------------------------------------

                //values at cell (Uir control volume) center
                //------------------------------------
                double rhoi_tmp = (gc.rhoi[i][j] + gc.rhoi[i][j-1])/2.0;
                double rate_ionize_tmp = (gc.rate_ionize[i][j] + gc.rate_ionize[i][j-1])/2.0;
                double nui_m_tmp = (gc.nui_m[i][j] + gc.nui_m[i][j-1])/2.0;
                //------------------------------------

                //density gradient at cell (Uir control volume) center
                //------------------------------------
                double dndr = (gc.rhoi[i][j]-gc.rhoi[i][j-1])/pm.dr;
                //------------------------------------

                //Lorentz & centrifugal-force at cell (Uir control volume) center
                //------------------------------------
                double Bx_tmp = (gc.r[j]*gc.Bx[i][j] + gc.r[j-1]*gc.Bx[i][j-1])/(2.0*r_tmp);
                double Uip_tmp = (gc.r[j]*gc.Uip_old[i][j] + gc.r[j-1]*gc.Uip_old[i][j-1])/(2.0*r_tmp);
                double lorentz = ph::e0*Uip_tmp*Bx_tmp;
                double centrifugal = -Uip_tmp*Uip_tmp/r_tmp;
                //------------------------------------
                //*****************************************************
    
                //coefficient for ion momentum conservation for 5 point stencil
                //*****************************************************

                //vertical and horizental 5 stencil
                //------------------------------------
                aE[i][j] = -(Uix_Rx - fabs(Uix_Rx))/(2.0*pm.dx);
                aW[i][j] =  (Uix_Lx + fabs(Uix_Lx))/(2.0*pm.dx);
                aN[i][j] = -gc.r[j  ]/r_tmp*(Uir_Rr - fabs(Uir_Rr))/(2.0*pm.dr);
                aS[i][j] =  gc.r[j-1]/r_tmp*(Uir_Lr + fabs(Uir_Lr))/(2.0*pm.dr);
                aP[i][j] = 1.0/dt_i + rate_ionize_tmp/(rhoi_tmp + 1e-100) + nui_m_tmp
                    + aE[i][j] + aW[i][j] + aN[i][j] + aS[i][j];
                //------------------------------------
    
                //RHS
                //------------------------------------
                b[i][j] = gr.Uir_old[i][j]/dt_i + (ph::e0*gr.Er[i][j] + lorentz)/pm.massi 
                    - ph::Boltz*pm.Ti/(pm.massi*rhoi_tmp+1e-100)*dndr + centrifugal;
                //------------------------------------
                //*****************************************************
            }
        }
    }

    //Set BCs for coefficients
    //=====================================================================
    //left-wall-tangential-BC (z0) zero-flux
    //------------------------------------
    for (int j=gr.j_flr_bl[0][0];j<=gr.j_flr_bl[0][1];j++){
        int i=gr.i_flr_bl[0][0];
        //modify coefficient
        //------------------------------------
        aW[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //left-wall-tangential-BC (z1) zero-flux
    //------------------------------------
    for (int j=gr.j_flr_bl[2][0];j<=gr.j_flr_bl[0][0]-1;j++){
        int i=gr.i_flr_bl[2][0];
        //modify coefficient
        //------------------------------------
        aW[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //left-wall-tangential-BC (z2) zero-flux
    //------------------------------------
    for (int j=gr.j_flr_bl[3][0];j<=gr.j_flr_bl[2][0]-1;j++){
        int i=gr.i_flr_bl[3][0];
        //modify coefficient
        //------------------------------------
        aW[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //left-wall-tangential-BC (z4) zero-flux
    //------------------------------------
    for (int j=gr.j_flr_bl[1][0];j<=gr.j_flr_bl[4][1];j++){
        int i=gr.i_flr_bl[4][0];
        //modify coefficient
        //------------------------------------
        aW[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //right-wall-tangential-BC (z3) zero-flux
    //------------------------------------
    for (int j=gr.j_flr_bl[1][0];j<=gr.j_flr_bl[1][1];j++){
        int i=gr.i_flr_bl[1][1];
        //modify coefficient
        //------------------------------------
        aE[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------
    
    //right-wall-tangential-BC (z5) zero-flux
    //------------------------------------
    for (int j=gr.j_flr_bl[4][0];j<=gr.j_flr_bl[4][1];j++){
        int i=gr.i_flr_bl[4][1];
        //modify coefficient
        //------------------------------------
        aE[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //lower-wall-normal-BC (x0) energy-conservation
    //------------------------------------
    for (int i=gr.i_flr_bl[0][0];i<=gr.i_flr_bl[0][1];i++){
        int j=gr.j_flr_bl[0][0];
        //calculate velocity at wall
        //------------------------------------
        double delta_phi = (gc.phi[i][j] + gc.phi[i][j-1])/2.0;
        double Ui_wall = -sqrt(pow(gr.Uir_old[i][j],2) + fmax(2.0*ph::e0*delta_phi/pm.massi,0.0));
        //------------------------------------

        //modify coefficient
        //------------------------------------
        b[i][j] = b[i][j] + aS[i][j]*Ui_wall;
        aS[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //lower-wall-normal-BC (x2) energy-conservation
    //------------------------------------
    for (int i=gr.i_flr_bl[2][0];i<=gr.i_flr_bl[2][1];i++){
        int j=gr.j_flr_bl[2][0];
        //calculate velocity at wall
        //------------------------------------
        double delta_phi = (gc.phi[i][j] + gc.phi[i][j-1])/2.0;
        double Ui_wall = -sqrt(pow(gr.Uir_old[i][j],2) + fmax(2.0*ph::e0*delta_phi/pm.massi,0.0));
        //------------------------------------

        //modify coefficient
        //------------------------------------
        b[i][j] = b[i][j] + aS[i][j]*Ui_wall;
        aS[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------
    
    //centerline-normal-BC (x6) zero-flux
    //------------------------------------
    for (int i=gr.i_flr_bl[3][0];i<=gr.i_flr_bl[4][1];i++){ 
        int j = gr.j_flr_bl[3][0];
        //modify coefficient
        //------------------------------------
        aS[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //upper-wall-nornmal-BC (x1) energy-conservation
    //------------------------------------
    for (int i=gr.i_flr_bl[0][0];i<=gr.i_flr_bl[1][1];i++){
        int j=gr.j_flr_bl[0][1];
        //calculate velocity at wall
        //------------------------------------
        double delta_phi = (gc.phi[i][j] + gc.phi[i][j-1])/2.0;
        double Ui_wall = sqrt(pow(gr.Uir_old[i][j],2) + fmax(2.0*ph::e0*delta_phi/pm.massi,0.0));
        //------------------------------------
        
        //modify coefficient
        //------------------------------------
        b[i][j] = b[i][j] + aN[i][j]*Ui_wall;
        aN[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //upper-wall-nornmal-BC (x4) energy-conservation
    //------------------------------------
    for (int i=gr.i_flr_bl[1][1]+1;i<=gr.i_flr_bl[3][1];i++){
        int j=gr.j_flr_bl[3][1];
        //calculate velocity at wall
        //------------------------------------
        double delta_phi = (gc.phi[i][j] + gc.phi[i][j-1])/2.0;
        double Ui_wall = sqrt(pow(gr.Uir_old[i][j],2) + fmax(2.0*ph::e0*delta_phi/pm.massi,0.0));
        //------------------------------------
        
        //modify coefficient
        //------------------------------------
        b[i][j] = b[i][j] + aN[i][j]*Ui_wall;
        aN[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //upper-open-nornmal-BC (x5) continuity
    //------------------------------------
    for (int i=gr.i_flr_bl[4][0];i<=gr.i_flr_bl[4][1];i++){
        int j=gr.j_flr_bl[4][1];
        //radial position at upper and lower cell interface
        //------------------------------------
        double r_tmp = (gc.r[j]+gc.r[j-1])/2.0;
        double rR = (gc.r[j+1]+gc.r[j])/2.0;
        //------------------------------------

        /*
        //modify coefficient (allow backflow)
        //------------------------------------
        aP[i][j] = aP[i][j] - aN[i][j]*r_tmp/rR;
        aN[i][j] = 0.0;
        //------------------------------------
        */

        //modify coefficient (block backflow)
        //------------------------------------
        b[i][j] = b[i][j] + aN[i][j]*r_tmp/rR*fmax(gr.Uir_old[i][j],0);
        aN[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------
    //=====================================================================

    //Solver
    //=====================================================================
    if(pm.flag_iter_Ui == 0){
        //SOR
        //------------------------------------
        double alpha = 1.9; //acceleration factor
        solver_SOR(aP,aE,aW,aN,aS,b,gr.i_flr_bl,gr.j_flr_bl,pm.n_bl-1,alpha,pm.maxITR_SOR_Ui,pm.error_cnv_SOR_Ui,0,gr.Uir);
        //------------------------------------
    }else{
        //Semi-coarse Multi-Grid by Hypre
        //------------------------------------
        solver_SMG(aP,aE,aW,aN,aS,b,gr.i_flr_bl,gr.j_flr_bl,pm.n_bl-1,pm.maxITR_SOR_Ui,pm.error_cnv_SOR_Ui,0,gr.Uir);
        //------------------------------------
    }
    //=====================================================================

    //Reflect BCs for ghost cell
    //=====================================================================
    //left-wall-tangential-BC (z0) zero-flux
    //------------------------------------
    for (int j=gr.j_flr_bl[0][0];j<=gr.j_flr_bl[0][1];j++){
        int i=gr.i_flr_bl[0][0];
        //reflect BC
        //------------------------------------
        gr.Uir[i-1][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //left-wall-tangential-BC (z1) zero-flux
    //------------------------------------
    for (int j=gr.j_flr_bl[2][0];j<=gr.j_flr_bl[0][0]-2;j++){
        int i=gr.i_flr_bl[2][0];
        //reflect BC
        //------------------------------------
        gr.Uir[i-1][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //left-wall-tangential-BC (z2) zero-flux
    //------------------------------------
    for (int j=gr.j_flr_bl[3][0];j<=gr.j_flr_bl[2][0]-2;j++){
        int i=gr.i_flr_bl[3][0];
        //reflect BC
        //------------------------------------
        gr.Uir[i-1][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //left-wall-tangential-BC (z4) zero-flux
    //------------------------------------
    for (int j=gr.j_flr_bl[1][0]+1;j<=gr.j_flr_bl[4][1];j++){
        int i=gr.i_flr_bl[4][0];
        //reflect BC
        //------------------------------------
        gr.Uir[i-1][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //right-wall-tangential-BC (z3) zero-flux
    //------------------------------------
    for (int j=gr.j_flr_bl[1][0]+1;j<=gr.j_flr_bl[1][1];j++){
        int i=gr.i_flr_bl[1][1];
        //reflect BC
        //------------------------------------
        gr.Uir[i+1][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------
    
    //right-wall-tangential-BC (z5) zero-flux
    //------------------------------------
    for (int j=gr.j_flr_bl[4][0];j<=gr.j_flr_bl[4][1];j++){
        int i=gr.i_flr_bl[4][1];
        //reflect BC
        //------------------------------------
        gr.Uir[i+1][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //lower-wall-normal-BC (x0) energy-conservation
    //------------------------------------
    for (int i=gr.i_flr_bl[0][0];i<=gr.i_flr_bl[0][1];i++){
        int j=gr.j_flr_bl[0][0];
        //calculate velocity at wall
        //------------------------------------
        double delta_phi = (gc.phi[i][j] + gc.phi[i][j-1])/2.0;
        //------------------------------------

        //modify coefficient
        //------------------------------------
        double Ui_wall = -sqrt(pow(gr.Uir_old[i][j],2) + fmax(2.0*ph::e0*delta_phi/pm.massi,0.0));
        gr.Uir[i][j-1] = Ui_wall;
        //------------------------------------
    }
    //------------------------------------

    //lower-wall-normal-BC (x2) energy-conservation
    //------------------------------------
    for (int i=gr.i_flr_bl[2][0];i<=gr.i_flr_bl[2][1];i++){
        int j=gr.j_flr_bl[2][0];
        //calculate velocity at wall
        //------------------------------------
        double delta_phi = (gc.phi[i][j] + gc.phi[i][j-1])/2.0;
        //------------------------------------

        //modify coefficient
        //------------------------------------
        double Ui_wall = -sqrt(pow(gr.Uir_old[i][j],2) + fmax(2.0*ph::e0*delta_phi/pm.massi,0.0));
        gr.Uir[i][j-1] = Ui_wall;
        //------------------------------------
    }
    //------------------------------------
    
    //centerline-normal-BC (x6) zero-flux
    //------------------------------------
    for (int i=gr.i_flr_bl[3][0];i<=gr.i_flr_bl[4][1];i++){ 
        int j = gr.j_flr_bl[3][0];
        //modify coefficient
        //------------------------------------
        gr.Uir[i][j-1] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //upper-wall-nornmal-BC (x1) energy-conservation
    //------------------------------------
    for (int i=gr.i_flr_bl[0][0];i<=gr.i_flr_bl[1][1];i++){
        int j=gr.j_flr_bl[0][1];
        //calculate velocity at wall
        //------------------------------------
        double delta_phi = (gc.phi[i][j] + gc.phi[i][j-1])/2.0;
        double Ui_wall = sqrt(pow(gr.Uir_old[i][j],2) + fmax(2.0*ph::e0*delta_phi/pm.massi,0.0));
        //------------------------------------
        
        //modify coefficient
        //------------------------------------
        gr.Uir[i][j+1] = Ui_wall;
        //------------------------------------
    }
    //------------------------------------

    //upper-wall-nornmal-BC (x4) energy-conservation
    //------------------------------------
    for (int i=gr.i_flr_bl[1][1]+1;i<=gr.i_flr_bl[3][1];i++){
        int j=gr.j_flr_bl[3][1];
        //calculate velocity at wall
        //------------------------------------
        double delta_phi = (gc.phi[i][j] + gc.phi[i][j-1])/2.0;
        double Ui_wall = sqrt(pow(gr.Uir_old[i][j],2) + fmax(2.0*ph::e0*delta_phi/pm.massi,0.0));
        //------------------------------------
        
        //modify coefficient
        //------------------------------------
        gr.Uir[i][j+1] = Ui_wall;
        //------------------------------------
    }
    //------------------------------------

    //upper-open-nornmal-BC (x5) continuity
    //------------------------------------
    for (int i=gr.i_flr_bl[4][0];i<=gr.i_flr_bl[4][1];i++){
        int j=gr.j_flr_bl[4][1];
        //radial position at upper and lower cell interface
        //------------------------------------
        double r_tmp = (gc.r[j]+gc.r[j-1])/2.0;
        double rR = (gc.r[j+1]+gc.r[j])/2.0;
        //------------------------------------

        /*
        //modify coefficient (allow backflow)
        //------------------------------------
        gr.Uir[i][j+1] = r_tmp/rR*gr.Uir[i][j];
        //------------------------------------
        */
        
        //modify coefficient (block backflow)
        //------------------------------------
        gr.Uir[i][j+1] = r_tmp/rR*fmax(gr.Uir_old[i][j],0);
        //------------------------------------
    }
    //------------------------------------
    //=====================================================================
}

//*****************************************************************
//**                                                             **
//**           void solve_Uip_tmp                                **
//**                                                             **
//*****************************************************************
void FluidModule::solve_Uip_tmp(Params &pm,GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr){

    vector<vector<double> > aW(pm.ni+2,vector<double>(pm.nj+2,0.0)); //係数
    vector<vector<double> > aE(pm.ni+2,vector<double>(pm.nj+2,0.0)); //係数
    vector<vector<double> > aS(pm.ni+2,vector<double>(pm.nj+2,0.0)); //係数
    vector<vector<double> > aN(pm.ni+2,vector<double>(pm.nj+2,0.0)); //係数
    vector<vector<double> > aP(pm.ni+2,vector<double>(pm.nj+2,0.0)); //係数
    vector<vector<double> > b (pm.ni+2,vector<double>(pm.nj+2,0.0)); //ソース

    //time step for ion
    //------------------------------------
    double dt_i = pm.dt*double(pm.ndt_i);
    //------------------------------------

    //Make coefficient for implicit p-directional ion momentum conservation for 5 point stencil
    //=================================================================================
    for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
        for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
            for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                
                //radial position at upper and lower cell interface
                //------------------------------------
                double rL = (gc.r[j]+gc.r[j-1])/2.0;
                double rR = (gc.r[j]+gc.r[j+1])/2.0;
                //------------------------------------
                
                //normalized radial position at upper and lower cell interface
                //------------------------------------
                double qL = rL/gc.r[j];
                double qR = rR/gc.r[j];
                //------------------------------------
                
                //values for upwind discretization
                //*****************************************************
                //velocity at cell (Uip control volume) interface
                //------------------------------------
                double Uix_Lx = gx.Uix_old[i][j];
                double Uix_Rx = gx.Uix_old[i+1][j];
                double Uir_Lr = gr.Uir_old[i][j];
                double Uir_Rr = gr.Uir_old[i][j+1];
                //------------------------------------

                //Lorentz-force at cell (Uip control volume) center
                //------------------------------------
                double Uix_tmp = (gx.Uix_old[i+1][j] + gx.Uix_old[i][j])/2.0;
                double Uir_tmp = (rL*gr.Uir_old[i][j] + rR*gr.Uir_old[i][j+1])/(2.0*gc.r[j]);
                double lorentz = ph::e0*(Uix_tmp*gc.Br[i][j] - Uir_tmp*gc.Bx[i][j]);
                double centrifugal = gc.Uip[i][j]*Uir_tmp/gc.r[j];
                //------------------------------------
                //*****************************************************

                //coefficient for ion momentum conservation for 5 point stencil
                //*****************************************************

                //vertical and horizental 5 stencil
                //------------------------------------
                aE[i][j] = -(Uix_Rx - fabs(Uix_Rx))/(2.0*pm.dx);
                aW[i][j] =  (Uix_Lx + fabs(Uix_Lx))/(2.0*pm.dx);
                aN[i][j] = -qR*(Uir_Rr - fabs(Uir_Rr))/(2.0*pm.dr);
                aS[i][j] =  qL*(Uir_Lr + fabs(Uir_Lr))/(2.0*pm.dr);
                aP[i][j] = 1.0/dt_i  + gc.rate_ionize[i][j]/(gc.rhoi[i][j] + 1e-100) + gc.nui_m[i][j]
                    + aE[i][j] + aW[i][j] + aN[i][j] + aS[i][j];
                //------------------------------------

                //RHS
                //------------------------------------
                b[i][j] = gc.Uip_old[i][j]/dt_i + lorentz/pm.massi + centrifugal;
                //------------------------------------
                //*****************************************************
            }
        }
    }

    //Set BCs for coefficients
    //=====================================================================
    //left-wall-tangential-BC (z0) continuity
    //------------------------------------
    for (int j=gc.j_flc_bl[0][0];j<=gc.j_flc_bl[0][1];j++){
        int i=gc.i_flc_bl[0][0];
        //modify coefficient
        //------------------------------------
        aP[i][j] = aP[i][j] - aW[i][j];
        aW[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //left-wall-tangential-BC (z1) continuity
    //------------------------------------
    for (int j=gc.j_flc_bl[2][0];j<=gc.j_flc_bl[0][0]-1;j++){
        int i=gc.i_flc_bl[2][0];
        //modify coefficient
        //------------------------------------
        aP[i][j] = aP[i][j] - aW[i][j];
        aW[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //left-wall-tangential-BC (z2) continuity
    //------------------------------------
    for (int j=gc.j_flc_bl[3][0];j<=gc.j_flc_bl[2][0]-1;j++){
        int i=gc.i_flc_bl[3][0];
        //modify coefficient
        //------------------------------------
        aP[i][j] = aP[i][j] - aW[i][j];
        aW[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //left-wall-tangential-BC (z4) continuity
    //------------------------------------
    for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[4][1];j++){
        int i=gc.i_flc_bl[4][0];
        //modify coefficient
        //------------------------------------
        aP[i][j] = aP[i][j] - aW[i][j];
        aW[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //right-wall-tangential-BC (z3) continuity
    //------------------------------------
    for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[1][1];j++){
        int i=gc.i_flc_bl[1][1];
        //modify coefficient
        //------------------------------------
        aP[i][j] = aP[i][j] - aE[i][j];
        aE[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //right-wall-tangential-BC (z5) continuity
    //------------------------------------
    for (int j=gc.j_flc_bl[4][0];j<=gc.j_flc_bl[4][1];j++){
        int i=gc.i_flc_bl[4][1];
        //modify coefficient
        //------------------------------------
        aP[i][j] = aP[i][j] - aE[i][j];
        aE[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //lower-wall-tangential-BC (x0) continuity
    //------------------------------------
    for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[0][1];i++){
        int j=gc.j_flc_bl[0][0];
        //modify coefficient
        //------------------------------------
        aP[i][j] = aP[i][j] - aS[i][j]*gc.r[j]/gc.r[j-1];
        aS[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //lower-wall-tangential-BC (x2) continuity
    //------------------------------------
    for (int i=gc.i_flc_bl[2][0];i<=gc.i_flc_bl[2][1];i++){
        int j=gc.j_flc_bl[2][0];
        //modify coefficient
        //------------------------------------
        aP[i][j] = aP[i][j] - aS[i][j]*gc.r[j]/gc.r[j-1];
        aS[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //upper-wall-tangential-BC (x1) continuity
    //------------------------------------
    for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[1][1];i++){
        int j=gc.j_flc_bl[0][1];
        //modify coefficient
        //------------------------------------
        aP[i][j] = aP[i][j] - aN[i][j]*gc.r[j]/gc.r[j+1];
        aN[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //upper-wall-tangential-BC (x4) continuity
    //------------------------------------
    for (int i=gc.i_flc_bl[1][1]+1;i<=gc.i_flc_bl[3][1];i++){
        int j=gc.j_flc_bl[3][1];
        //modify coefficient
        //------------------------------------
        aP[i][j] = aP[i][j] - aN[i][j]*gc.r[j]/gc.r[j+1];
        aN[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //upper-open-tangential-BC (x5) continuity
    //------------------------------------
    for (int i=gc.i_flc_bl[4][0];i<=gc.i_flc_bl[4][1];i++){
        int j=gc.j_flc_bl[4][1];
        //modify coefficient
        //------------------------------------
        aP[i][j] = aP[i][j] - aN[i][j]*gc.r[j]/gc.r[j+1];
        aN[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------

    //Solver
    //=====================================================================
    if(pm.flag_iter_Ui == 0){
        //SOR
        //------------------------------------
        double alpha = 1.9; //acceleration factor
        solver_SOR(aP,aE,aW,aN,aS,b,gc.i_flc_bl,gc.j_flc_bl,pm.n_bl-1,alpha,pm.maxITR_SOR_Ui,pm.error_cnv_SOR_Ui,0,gc.Uip);
    }else{
        //Semi-coarse Multi-Grid by Hypre
        //------------------------------------
        solver_SMG(aP,aE,aW,aN,aS,b,gc.i_flc_bl,gc.j_flc_bl,pm.n_bl-1,pm.maxITR_SOR_Ui,pm.error_cnv_SOR_Ui,0,gc.Uip);
        //------------------------------------
    }
    //=====================================================================

    //Set BCs for coefficients
    //=====================================================================
    //left-wall-tangential-BC (z0) continuity
    //------------------------------------
    for (int j=gc.j_flc_bl[0][0];j<=gc.j_flc_bl[0][1];j++){
        int i=gc.i_flc_bl[0][0];
        //reflect BC
        //------------------------------------
        gc.Uip[i-1][j] = gc.Uip[i][j];
        //------------------------------------
    }
    //------------------------------------

    //left-wall-tangential-BC (z1) continuity
    //------------------------------------
    for (int j=gc.j_flc_bl[2][0];j<=gc.j_flc_bl[0][0]-1;j++){
        int i=gc.i_flc_bl[2][0];
        //reflect BC
        //------------------------------------
        gc.Uip[i-1][j] = gc.Uip[i][j];
        //------------------------------------
    }
    //------------------------------------

    //left-wall-tangential-BC (z2) continuity
    //------------------------------------
    for (int j=gc.j_flc_bl[3][0];j<=gc.j_flc_bl[2][0]-1;j++){
        int i=gc.i_flc_bl[3][0];
        //reflect BC
        //------------------------------------
        gc.Uip[i-1][j] = gc.Uip[i][j];
        //------------------------------------
    }
    //------------------------------------

    //left-wall-tangential-BC (z4) continuity
    //------------------------------------
    for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[4][1];j++){
        int i=gc.i_flc_bl[4][0];
        //reflect BC
        //------------------------------------
        gc.Uip[i-1][j] = gc.Uip[i][j];
        //------------------------------------
    }
    //------------------------------------

    //right-wall-tangential-BC (z3) continuity
    //------------------------------------
    for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[1][1];j++){
        int i=gc.i_flc_bl[1][1];
        //reflect BC
        //------------------------------------
        gc.Uip[i+1][j] = gc.Uip[i][j];
        //------------------------------------
    }
    //------------------------------------

    //right-wall-tangential-BC (z5) continuity
    //------------------------------------
    for (int j=gc.j_flc_bl[4][0];j<=gc.j_flc_bl[4][1];j++){
        int i=gc.i_flc_bl[4][1];
        //reflect BC
        //------------------------------------
        gc.Uip[i+1][j] = gc.Uip[i][j];
        //------------------------------------
    }
    //------------------------------------

    //lower-wall-tangential-BC (x0) continuity
    //------------------------------------
    for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[0][1];i++){
        int j=gc.j_flc_bl[0][0];
        //reflect BC
        //------------------------------------
        gc.Uip[i][j-1] = gc.r[j]/gc.r[j-1]*gc.Uip[i][j];
        //------------------------------------
    }
    //------------------------------------

    //lower-wall-tangential-BC (x2) continuity
    //------------------------------------
    for (int i=gc.i_flc_bl[2][0];i<=gc.i_flc_bl[2][1];i++){
        int j=gc.j_flc_bl[2][0];
        //reflect BC
        //------------------------------------
        gc.Uip[i][j-1] = gc.r[j]/gc.r[j-1]*gc.Uip[i][j];
        //------------------------------------
    }
    //------------------------------------

    //upper-wall-tangential-BC (x1) continuity
    //------------------------------------
    for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[1][1];i++){
        int j=gc.j_flc_bl[0][1];
        //reflect BC
        //------------------------------------
        gc.Uip[i][j+1] = gc.r[j]/gc.r[j+1]*gc.Uip[i][j];
        //------------------------------------
    }
    //------------------------------------

    //upper-wall-tangential-BC (x4) continuity
    //------------------------------------
    for (int i=gc.i_flc_bl[1][1]+1;i<=gc.i_flc_bl[3][1];i++){
        int j=gc.j_flc_bl[3][1];
        //reflect BC
        //------------------------------------
        gc.Uip[i][j+1] = gc.r[j]/gc.r[j+1]*gc.Uip[i][j];
        //------------------------------------
    }
    //------------------------------------

    //upper-open-tangential-BC (x5) continuity
    //------------------------------------
    for (int i=gc.i_flc_bl[4][0];i<=gc.i_flc_bl[4][1];i++){
        int j=gc.j_flc_bl[4][1];
        //reflect BC
        //------------------------------------
        gc.Uip[i][j+1] = gc.r[j]/gc.r[j+1]*gc.Uip[i][j];
        //------------------------------------
    }
    //------------------------------------
}

//*****************************************************************
//**                                                             **
//**           void solve_rhoi                                   **
//**                                                             **
//*****************************************************************
void FluidModule::solve_rhoi_constTe(Params &pm,GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr){

    vector<vector<double> > aW(pm.ni+2,vector<double>(pm.nj+2,0.0)); //係数
    vector<vector<double> > aE(pm.ni+2,vector<double>(pm.nj+2,0.0)); //係数
    vector<vector<double> > aS(pm.ni+2,vector<double>(pm.nj+2,0.0)); //係数
    vector<vector<double> > aN(pm.ni+2,vector<double>(pm.nj+2,0.0)); //係数
    vector<vector<double> > aP(pm.ni+2,vector<double>(pm.nj+2,0.0)); //係数
    vector<vector<double> > b (pm.ni+2,vector<double>(pm.nj+2,0.0)); //ソース

    //time step for ion
    //------------------------------------
    double dt_i = pm.dt*double(pm.ndt_i);
    //------------------------------------

    //Make coefficient for implicit r-directional ion continuity equation for 5 point stencil
    //=================================================================================
    for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
        for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
            for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                
                //radial position at upper and lower cell interface
                //------------------------------------
                double rL = (gc.r[j]+gc.r[j-1])/2.0;
                double rR = (gc.r[j]+gc.r[j+1])/2.0;
                //------------------------------------

                //normalized radial position at upper and lower cell interface
                //------------------------------------
                double qL = rL/gc.r[j];
                double qR = rR/gc.r[j];
                //------------------------------------

                //ion bulk velocity at cell interface
                //------------------------------------
                double Uix_Lx = gx.Uix[i][j];
                double Uix_Rx = gx.Uix[i+1][j];
                double Uir_Lr = gr.Uir[i][j];
                double Uir_Rr = gr.Uir[i][j+1];
                //------------------------------------

                //ion momentum frequency at cell interface
                //------------------------------------
                double nui_m_Lx = (gc.nui_m[i][j] + gc.nui_m[i-1][j])/2.0;
                double nui_m_Rx = (gc.nui_m[i][j] + gc.nui_m[i+1][j])/2.0;
                double nui_m_Lr = (gc.nui_m[i][j] + gc.nui_m[i][j-1])/2.0;
                double nui_m_Rr = (gc.nui_m[i][j] + gc.nui_m[i][j+1])/2.0;
                //------------------------------------

                //diffusion coefficient (numerical param.) at cell interface
                //------------------------------------
                double Te_rep = pm.Te_rep_eV*ph::e0/ph::Boltz;
                double nu_star = 1.0/pm.rmax*sqrt(ph::Boltz*Te_rep/pm.massi);
                double D_Lx = ph::Boltz*(pm.Ti + Te_rep)/(pm.massi*(nui_m_Lx+nu_star)+1e-100);
                double D_Rx = ph::Boltz*(pm.Ti + Te_rep)/(pm.massi*(nui_m_Rx+nu_star)+1e-100);
                double D_Lr = ph::Boltz*(pm.Ti + Te_rep)/(pm.massi*(nui_m_Lr+nu_star)+1e-100);
                double D_Rr = ph::Boltz*(pm.Ti + Te_rep)/(pm.massi*(nui_m_Rr+nu_star)+1e-100);
                //------------------------------------

                //set adjacent flag (adjacent = 0, others = 1)
                //*****************************************************
                //adjacent "cell" to the wall, open, or central axis BC
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
                //w/o x5
                double bRr = bRr_wall;
                //------------------------------------
                //*****************************************************


                //BC - normal flux component to the wall boundaries
                //*****************************************************
                //kinetic vetor flux divide by density (for implicit method)
                //------------------------------------
                double GammaPerN_wall_Lx = 0.0;
                double GammaPerN_wall_Rx = 0.0;
                double GammaPerN_wall_Lr = 0.0;
                double GammaPerN_wall_Rr = 0.0;
                if(bLx_wall == 0){
                    double vth = sqrt(8.0*ph::Boltz*pm.Ti/(M_PI*pm.massi));
                    double um = sqrt(2.0*ph::Boltz*pm.Ti/pm.massi);
                    double utx = gx.Uix[i][j]/(um + 1e-100);
                    GammaPerN_wall_Lx = -0.25*vth*(exp(-utx*utx) + utx*sqrt(M_PI)*(erf(utx) - 1.0));
                    //if(i == gc.i_flc_bl[2][0]) std::cout << "before,"<< j << ","<<GammaPerN_wall_Lx<<endl;
                }
                if(bRx_wall == 0){
                    double vth = sqrt(8.0*ph::Boltz*pm.Ti/(M_PI*pm.massi));
                    double um = sqrt(2.0*ph::Boltz*pm.Ti/pm.massi);
                    double utx = gx.Uix[i+1][j]/(um + 1e-100);
                    GammaPerN_wall_Rx = 0.25*vth*(exp(-utx*utx) + utx*sqrt(M_PI)*(erf(utx) + 1.0));
                }
                if(bLr_wall == 0){
                    double vth = sqrt(8.0*ph::Boltz*pm.Ti/(M_PI*pm.massi));
                    double um = sqrt(2.0*ph::Boltz*pm.Ti/pm.massi);
                    double utr = gr.Uir[i][j]/(um + 1e-100);
                    GammaPerN_wall_Lr = -0.25*vth*(exp(-utr*utr) + utr*sqrt(M_PI)*(erf(utr) - 1.0));
                }
                if(bRr_wall == 0){
                    double vth = sqrt(8.0*ph::Boltz*pm.Ti/(M_PI*pm.massi));
                    double um = sqrt(2.0*ph::Boltz*pm.Ti/pm.massi);
                    double utr = gr.Uir[i][j+1]/(um + 1e-100);
                    GammaPerN_wall_Rr = 0.25*vth*(exp(-utr*utr) + utr*sqrt(M_PI)*(erf(utr) + 1.0));
                }
                //------------------------------------
                //*****************************************************

                //coefficient for ion momentum conservation for 5 point stencil
                //*****************************************************
                
                //vertical and horizental 5 stencil
                //------------------------------------
                aP[i][j] = (Uix_Rx + fabs(Uix_Rx))/(2.0*pm.dx)*bRx - (Uix_Lx - fabs(Uix_Lx))/(2.0*pm.dx)*bLx
                    + qR*(Uir_Rr + fabs(Uir_Rr))/(2.0*pm.dr)*bRr - qL*(Uir_Lr - fabs(Uir_Lr))/(2.0*pm.dr)*bLr
                    + ((D_Rx*bRx + D_Lx*bLx)/(pm.dx*pm.dx) + (qR*D_Rr*bRr + qL*D_Lr*bLr)/(pm.dr*pm.dr))
                    -    GammaPerN_wall_Lx/pm.dx //wall-flux boundary at left  wall (implicit)
                    +    GammaPerN_wall_Rx/pm.dx //wall-flux boundary at right wall (implicit)
                    - qL*GammaPerN_wall_Lr/pm.dr //wall-flux boundary at lower wall (implicit)
                    + qR*GammaPerN_wall_Rr/pm.dr //wall-flux boundary at upper wall (implicit)
                    + 1.0/dt_i; //time term
                aE[i][j] = (-(Uix_Rx - fabs(Uix_Rx))/(2.0*pm.dx) + D_Rx/(pm.dx*pm.dx))*bRx;
                aW[i][j] = ( (Uix_Lx + fabs(Uix_Lx))/(2.0*pm.dx) + D_Lx/(pm.dx*pm.dx))*bLx;
                aN[i][j] = (-qR*(Uir_Rr - fabs(Uir_Rr))/(2.0*pm.dr) + qR*D_Rr/(pm.dr*pm.dr))*bRr;
                aS[i][j] = ( qL*(Uir_Lr + fabs(Uir_Lr))/(2.0*pm.dr) + qL*D_Lr/(pm.dr*pm.dr))*bLr;
                //------------------------------------

                //RHS
                //------------------------------------
                double divDflux = (D_Rx*gc.rhoi_old[i+1][j]*bRx - (D_Rx*bRx + D_Lx*bLx)*gc.rhoi_old[i][j] + D_Lx*gc.rhoi_old[i-1][j]*bLx)/(pm.dx*pm.dx)
                    + (qR*D_Rr*gc.rhoi_old[i][j+1]*bRr - (qR*D_Rr*bRr + qL*D_Lr*bLr)*gc.rhoi_old[i][j] + qL*D_Lr*gc.rhoi_old[i][j-1]*bLr)/(pm.dr*pm.dr);
                
                //double divDflux = (D_Rx*gc.rhoi_old[i+1][j] - (D_Rx + D_Lx)*gc.rhoi_old[i][j] + D_Lx*gc.rhoi_old[i-1][j])/(pm.dx*pm.dx)
                //    + (qR*D_Rr*gc.rhoi_old[i][j+1] - (qR*D_Rr + qL*D_Lr)*gc.rhoi_old[i][j] + qL*D_Lr*gc.rhoi_old[i][j-1])/(pm.dr*pm.dr);

                b[i][j] = gc.rhoi_old[i][j]/dt_i + gc.rate_ionize[i][j] - divDflux;
                //------------------------------------
            }
        }
    }

    //upper-open-tangential-BC (x5) zero-gradient
    //------------------------------------
    for (int i=gc.i_flc_bl[4][0];i<=gc.i_flc_bl[4][1];i++){
        int j=gc.j_flc_bl[4][1];
        //modify coefficient
        //------------------------------------
        aP[i][j] = aP[i][j] - aN[i][j];
        aN[i][j] = 0.0;
        //------------------------------------
    }
    //------------------------------------
    
    //Solver
    //=====================================================================
    if(pm.flag_iter_rhoi==0){
        //SOR
        //------------------------------------
        double alpha = 1.9;
        solver_SOR(aP,aE,aW,aN,aS,b,gc.i_flc_bl,gc.j_flc_bl,pm.n_bl-1,alpha,pm.maxITR_SOR_rhoi,pm.error_cnv_SOR_rhoi,0,gc.rhoi);
        //------------------------------------
    }else{
        //SOR
        //------------------------------------
        solver_SMG(aP,aE,aW,aN,aS,b,gc.i_flc_bl,gc.j_flc_bl,pm.n_bl-1,pm.maxITR_SOR_rhoi,pm.error_cnv_SOR_rhoi,0,gc.rhoi);
        //------------------------------------
    }
    //=====================================================================

    //left-wall-BC (z0)
    //------------------------------------
    for (int j=gc.j_flc_bl[0][0];j<=gc.j_flc_bl[0][1];j++){
        int i=gc.i_flc_bl[0][0];

        double vth = sqrt(8.0*ph::Boltz*pm.Ti/(M_PI*pm.massi));
        double um = sqrt(2.0*ph::Boltz*pm.Ti/pm.massi);
        double utx = gx.Uix[i][j]/(um + 1e-100);
        double GammaPerN_wall_Lx = -0.25*vth*(exp(-utx*utx) + utx*sqrt(M_PI)*(erf(utx) - 1.0));
        gx.rhoUix_wall[i][j] = gc.rhoi[i][j]*GammaPerN_wall_Lx;
    }

    //left-wall-BC (z1)
    //------------------------------------
    for (int j=gc.j_flc_bl[2][0];j<=gc.j_flc_bl[0][0]-1;j++){
        int i=gc.i_flc_bl[2][0];
        
        double vth = sqrt(8.0*ph::Boltz*pm.Ti/(M_PI*pm.massi));
        double um = sqrt(2.0*ph::Boltz*pm.Ti/pm.massi);
        double utx = gx.Uix[i][j]/(um + 1e-100);
        double GammaPerN_wall_Lx = -0.25*vth*(exp(-utx*utx) + utx*sqrt(M_PI)*(erf(utx) - 1.0));
        //std::cout << "after,"<< j << ","<< gc.rhoi[i][j] << ","<<GammaPerN_wall_Lx << ",";
        gx.rhoUix_wall[i][j] = gc.rhoi[i][j]*GammaPerN_wall_Lx;
        //std::cout <<gx.rhoUix_wall[i][j] <<endl;
    }
    //------------------------------------

    //left-wall-BC (z2)
    //------------------------------------
    for (int j=gc.j_flc_bl[3][0];j<=gc.j_flc_bl[2][0]-1;j++){
        int i=gc.i_flc_bl[3][0];
        
        double vth = sqrt(8.0*ph::Boltz*pm.Ti/(M_PI*pm.massi));
        double um = sqrt(2.0*ph::Boltz*pm.Ti/pm.massi);
        double utx = gx.Uix[i][j]/(um + 1e-100);
        double GammaPerN_wall_Lx = -0.25*vth*(exp(-utx*utx) + utx*sqrt(M_PI)*(erf(utx) - 1.0));
        gx.rhoUix_wall[i][j] = gc.rhoi[i][j]*GammaPerN_wall_Lx;
    }
    //------------------------------------

    //left-wall-BC (z4)
    //------------------------------------
    for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[4][1];j++){
        int i=gc.i_flc_bl[4][0];
        double vth = sqrt(8.0*ph::Boltz*pm.Ti/(M_PI*pm.massi));
        double um = sqrt(2.0*ph::Boltz*pm.Ti/pm.massi);
        double utx = gx.Uix[i][j]/(um + 1e-100);
        double GammaPerN_wall_Lx = -0.25*vth*(exp(-utx*utx) + utx*sqrt(M_PI)*(erf(utx) - 1.0));
        gx.rhoUix_wall[i][j] = gc.rhoi[i][j]*GammaPerN_wall_Lx;
    }
    //------------------------------------

    //right-wall-BC (z3)
    //------------------------------------
    for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[1][1];j++){
        int i=gc.i_flc_bl[1][1];
        double vth = sqrt(8.0*ph::Boltz*pm.Ti/(M_PI*pm.massi));
        double um = sqrt(2.0*ph::Boltz*pm.Ti/pm.massi);
        double utx = gx.Uix[i+1][j]/(um + 1e-100);
        double GammaPerN_wall_Rx = 0.25*vth*(exp(-utx*utx) + utx*sqrt(M_PI)*(erf(utx) + 1.0));
        gx.rhoUix_wall[i+1][j] = gc.rhoi[i][j]*GammaPerN_wall_Rx;
    }
    //------------------------------------

    //right-wall-BC (z5)
    //------------------------------------
    for (int j=gc.j_flc_bl[4][0];j<=gc.j_flc_bl[4][1];j++){
        int i=gc.i_flc_bl[4][1];
        double vth = sqrt(8.0*ph::Boltz*pm.Ti/(M_PI*pm.massi));
        double um = sqrt(2.0*ph::Boltz*pm.Ti/pm.massi);
        double utx = gx.Uix[i+1][j]/(um + 1e-100);
        double GammaPerN_wall_Rx = 0.25*vth*(exp(-utx*utx) + utx*sqrt(M_PI)*(erf(utx) + 1.0));
        gx.rhoUix_wall[i+1][j] = gc.rhoi[i][j]*GammaPerN_wall_Rx;
    }
    //=====================================================================================

    //lower-wall-BC (x0)
    //------------------------------------
    for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[0][1];i++){
        int j=gc.j_flc_bl[0][0];
        double vth = sqrt(8.0*ph::Boltz*pm.Ti/(M_PI*pm.massi));
        double um = sqrt(2.0*ph::Boltz*pm.Ti/pm.massi);
        double utr = gr.Uir[i][j]/(um + 1e-100);
        double GammaPerN_wall_Lr = -0.25*vth*(exp(-utr*utr) + utr*sqrt(M_PI)*(erf(utr) - 1.0));
        gr.rhoUir_wall[i][j] = gc.rhoi[i][j]*GammaPerN_wall_Lr;
    }
    //------------------------------------

    //lower-wall-BC (x2)
    //------------------------------------
    for (int i=gc.i_flc_bl[2][0];i<=gc.i_flc_bl[2][1];i++){
        int j=gc.j_flc_bl[2][0];
        double vth = sqrt(8.0*ph::Boltz*pm.Ti/(M_PI*pm.massi));
        double um = sqrt(2.0*ph::Boltz*pm.Ti/pm.massi);
        double utr = gr.Uir[i][j]/(um + 1e-100);
        double GammaPerN_wall_Lr = -0.25*vth*(exp(-utr*utr) + utr*sqrt(M_PI)*(erf(utr) - 1.0));
        gr.rhoUir_wall[i][j] = gc.rhoi[i][j]*GammaPerN_wall_Lr;
    }
    //------------------------------------

    //upper-wall-BC (x1)
    //------------------------------------
    for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[1][1];i++){
        int j=gc.j_flc_bl[0][1];
        double vth = sqrt(8.0*ph::Boltz*pm.Ti/(M_PI*pm.massi));
        double um = sqrt(2.0*ph::Boltz*pm.Ti/pm.massi);
        double utr = gr.Uir[i][j+1]/(um + 1e-100);
        double GammaPerN_wall_Rr = 0.25*vth*(exp(-utr*utr) + utr*sqrt(M_PI)*(erf(utr) + 1.0));
        gr.rhoUir_wall[i][j+1] = gc.rhoi[i][j]*GammaPerN_wall_Rr;
    }
    //------------------------------------

    //upper-wall-BC (x4)
    //------------------------------------
    for (int i=gc.i_flc_bl[1][1]+1;i<=gc.i_flc_bl[3][1];i++){
        int j=gc.j_flc_bl[3][1];
        double vth = sqrt(8.0*ph::Boltz*pm.Ti/(M_PI*pm.massi));
        double um = sqrt(2.0*ph::Boltz*pm.Ti/pm.massi);
        double utr = gr.Uir[i][j+1]/(um + 1e-100);
        double GammaPerN_wall_Rr = 0.25*vth*(exp(-utr*utr) + utr*sqrt(M_PI)*(erf(utr) + 1.0));
        gr.rhoUir_wall[i][j+1] = gc.rhoi[i][j]*GammaPerN_wall_Rr;
    }

    //------------------------------------

    //upper-open-BC (x5)
    //------------------------------------
    for (int i=gc.i_flc_bl[4][0];i<=gc.i_flc_bl[4][1];i++){
        int j=gc.j_flc_bl[4][1];
        double rhoUir_Rr = gc.rhoi[i][j]*gr.Uir[i][j+1];
        double Gamma_open_Rr = rhoUir_Rr;
        gr.rhoUir[i][j+1] =  Gamma_open_Rr;
        gr.rhoUir_wall[i][j+1] = 0.0;
    }
    //------------------------------------

    //Don’t update unless the value is positive.
    //------------------------------------
    int ncount = 0;
    for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
        for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
            for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                if(gc.rhoi[i][j]<0.0){
                    //cout << "Warning negative density at (" << i <<","<<j<<") ne = " << gc.rhoi[i][j]<< endl; 
                    ncount ++;
                    gc.rhoi[i][j] = fmax(gc.rhoi_old[i][j],0.0);
                }
            }
        }
    }
    //if(ncount!=0){
    //    cout << "negative rhoi count = " << ncount << endl;
    //}

    //upper-open-tangential-BC (x5) zero-gradient
    //------------------------------------
    for (int i=gc.i_flc_bl[4][0];i<=gc.i_flc_bl[4][1];i++){
        int j=gc.j_flc_bl[4][1];
        //reflect BC
        //------------------------------------
        gc.rhoi[i][j+1] =  gc.rhoi[i][j];
        //------------------------------------
    }
    //------------------------------------
}


//*****************************************************************
//**                                                             **
//**           void correct_Ui_constTe                           **
//**                                                             **
//*****************************************************************
void FluidModule::correct_Ui_constTe(Params &pm,GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr){

    double ratio = 1.0;
    double Te_rep = 0.1*ph::e0/ph::Boltz; //代表温度は3eVに固定
    double nu_star = 0.1/pm.rmax*sqrt(ph::Boltz*Te_rep/pm.massi); //代表衝突周波数

    //time step for ion
    //------------------------------------
    double dt_i = pm.dt*double(pm.ndt_i);
    //------------------------------------
    
    //update Uix
    //=================================================================================
    for (int iblock=0;iblock<pm.n_bl-1;iblock++){
        for (int i=gx.i_flx_bl[iblock][0];i<=gx.i_flx_bl[iblock][1];i++){
            for (int j=gx.j_flx_bl[iblock][0];j<=gx.j_flx_bl[iblock][1];j++){
                //diffusion coefficient (numerical param.)
                //------------------------------------
                double Te_rep = pm.Te_rep_eV*ph::e0/ph::Boltz;
                double nu_star = 1.0/pm.rmax*sqrt(ph::Boltz*Te_rep/pm.massi);
                double nui_m_tmp = (gc.nui_m[i][j]+gc.nui_m[i-1][j])/2.0;
                double Te_tmp = Te_rep;
                double D_tmp = ph::Boltz*(pm.Ti + Te_tmp)/(pm.massi*(nui_m_tmp+nu_star)+1e-100);
                //----------------------------------------------------
                
                //density at cell (Uix control volume) center
                //----------------------------------------------------
                double rhoi_tmp = (gc.rhoi[i][j]+gc.rhoi[i-1][j])/2.0;
                //----------------------------------------------------

                //density gradient at cell (Uix control volume) center
                //------------------------------------
                double dndx_old = (gc.rhoi_old[i][j] - gc.rhoi_old[i-1][j])/pm.dx; //density gradient calculation
                double dndx = (gc.rhoi[i][j] - gc.rhoi[i-1][j])/pm.dx;
                //------------------------------------

                //temporary flux
                //------------------------------------
                double rhoUix_tmp =  gc.rhoi[i-1][j]*(gx.Uix[i][j] + fabs(gx.Uix[i][j]))/2.0
                    + gc.rhoi[i][j]*(gx.Uix[i][j] - fabs(gx.Uix[i][j]))/2.0;
                //------------------------------------

                //Update Uix
                //----------------------------------------------------
                gx.rhoUix[i][j] = rhoUix_tmp - D_tmp*(dndx - dndx_old);
                gx.Uix[i][j] = gx.Uix[i][j] - D_tmp/(rhoi_tmp+1e-100)*(dndx - dndx_old);
                //----------------------------------------------------
            }
        }
    }
    //=================================================================================

    //update Uir
    //=================================================================================
    for (int iblock=0;iblock<pm.n_bl-1;iblock++){
        for (int i=gr.i_flr_bl[iblock][0];i<=gr.i_flr_bl[iblock][1];i++){
            for (int j=gr.j_flr_bl[iblock][0];j<=gr.j_flr_bl[iblock][1] + int(iblock == 4);j++){
                //diffusion coefficient (numerical param.)
                //------------------------------------
                double Te_rep = pm.Te_rep_eV*ph::e0/ph::Boltz;
                double nu_star = 1.0/pm.rmax*sqrt(ph::Boltz*Te_rep/pm.massi);
                double r_tmp = (gc.r[j]+gc.r[j-1])/2.0;
                double nui_m_tmp = (gc.nui_m[i][j] + gc.nui_m[i][j-1])/2.0;
                double Te_tmp = Te_rep;
                double D_tmp = ph::Boltz*(pm.Ti + Te_tmp)/(pm.massi*(nui_m_tmp+nu_star)+1e-100);
                //----------------------------------------------------
                
                //density at cell (Uir control volume) center
                //----------------------------------------------------
                double rhoi_tmp = (gc.rhoi[i][j] + gc.rhoi[i][j-1])/2.0;
                //----------------------------------------------------

                //density gradient at cell (Uir control volume) center
                //------------------------------------
                double dndr_old = (gc.rhoi_old[i][j] - gc.rhoi_old[i][j-1])/pm.dr;
                double dndr     = (gc.rhoi[i][j]     - gc.rhoi[i][j-1]    )/pm.dr;
                //------------------------------------

                //temporary flux
                //------------------------------------
                double rhoUir_tmp =  gc.rhoi[i][j-1]*(gr.Uir[i][j] + fabs(gr.Uir[i][j]))/2.0
                    + gc.rhoi[i][j]*(gr.Uir[i][j] - fabs(gr.Uir[i][j]))/2.0;
                //------------------------------------

                //Update Uir
                //----------------------------------------------------
                gr.rhoUir[i][j] = rhoUir_tmp - D_tmp*(dndr - dndr_old);
                gr.Uir[i][j] = gr.Uir[i][j] - D_tmp/(rhoi_tmp+1e-100)*(dndr - dndr_old);
                //----------------------------------------------------
            }
        }
    }
    //=================================================================================

    //update Uip
    //=================================================================================
    // No process due to absense of density gradient in azimuthal direction
    //=================================================================================

    //calculation of particle balance
    //------------------------------------
    if(pm.itime%pm.ndiv_fout == 0 || pm.itime == pm.ntime){
        
        //volume integrated values
        //------------------------------------
        double dNidt_sum         = 0.0;
        double nabla_rhoUi_sum   = 0.0;
        double rate_d_ionz_sum   = 0.0;
        double rate_s_ionz_sum   = 0.0;
        //------------------------------------
    
        //calculate volume integrated of each term
        //------------------------------------
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){
            for (int i=gx.i_flx_bl[iblock][0];i<=gx.i_flx_bl[iblock][1];i++){
                for (int j=gx.j_flx_bl[iblock][0];j<=gx.j_flx_bl[iblock][1];j++){
                    
                    //radial position at upper and lower cell interface
                    //------------------------------------
                    double rR = (gc.r[j]+gc.r[j+1])/2.0;
                    double rL = (gc.r[j]+gc.r[j-1])/2.0;
                    //------------------------------------

                    //normalized radial position at upper and lower cell interface
                    //------------------------------------
                    double qR = (gc.r[j]+gc.r[j+1])/2.0/gc.r[j];
                    double qL = (gc.r[j]+gc.r[j-1])/2.0/gc.r[j];
                    //------------------------------------

                    //BC - set flag for adjacent to the wall, open, or central axis (adjacent = 0, others = 1)
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
                    //w/o x5
                    double bRr = bRr_wall;
                    //------------------------------------

                    double nabla_rhoUi = (gx.rhoUix[i+1][j]*bRx_wall - gx.rhoUix[i][j]*bLx_wall)/pm.dx 
                                       + (qR*gr.rhoUir[i][j+1]*bRr_wall -qL*gr.rhoUir[i][j]*bLr_wall)/pm.dr
                                       + (gx.rhoUix_wall[i+1][j] - gx.rhoUix_wall[i][j])/pm.dx 
                                       + (qR*gr.rhoUir_wall[i][j+1] - qL*gr.rhoUir_wall[i][j])/pm.dr;
                    
                    double volume = 2.0*M_PI*gc.r[j]*pm.dr*pm.dx;
                    
                    dNidt_sum         = dNidt_sum         + volume*(gc.rhoi[i][j] - gc.rhoi_old[i][j])/pm.dt;
                    nabla_rhoUi_sum   = nabla_rhoUi_sum   + volume*nabla_rhoUi;
                    rate_d_ionz_sum   = rate_d_ionz_sum   + volume*gc.rhoe[i][j]*gc.nu_ionz[i][j];
                    rate_s_ionz_sum   = rate_s_ionz_sum   + volume*gc.rhoe[i][j]*gc.nu_ionzStep[i][j];
                }
            }
        }

        //output particle balance for ion
        //------------------------------------
        string char1="results/particle_balance_ion";
        string char2=to_string(pm.nOut);
        string char_csv=".csv";
        ofstream outputfile9(char1+char2+char_csv);
        outputfile9 << setprecision(numeric_limits<double>::max_digits10) << scientific;
        outputfile9<< "dN/dt (pcl/s), loss_to_bnd (pcl/s), direct-iz (pcl/s), sw-iz (pcl/s), sum (pcl/s)"<< endl;
        double sum = -dNidt_sum + (-nabla_rhoUi_sum + rate_d_ionz_sum + rate_s_ionz_sum);
        outputfile9<< dNidt_sum << "," << -nabla_rhoUi_sum << "," << rate_d_ionz_sum << ","<< rate_s_ionz_sum 
            << "," << -dNidt_sum + (-nabla_rhoUi_sum + rate_d_ionz_sum + rate_s_ionz_sum)<< endl;
        outputfile9.close();
        //------------------------------------
    }
    //------------------------------------

    //Caluculate residuals
    //------------------------------------
    calcRes(pm.error_rhoi, gc.rhoi,gc.rhoi_old,gc.i_flc_bl, gc.j_flc_bl, pm.n_bl-1, 1e-100, pm.flag_error);
    calcRes(pm.error_Uix, gx.Uix, gx.Uix_old, gx.i_flx_bl, gx.j_flx_bl, pm.n_bl-1, 1e-100, pm.flag_error);
    calcRes(pm.error_Uir, gr.Uir, gr.Uir_old, gr.i_flr_bl, gr.j_flr_bl, pm.n_bl-1, 1e-100, pm.flag_error);
    calcRes(pm.error_Uip, gc.Uip, gc.Uip_old, gc.i_flc_bl, gc.j_flc_bl, pm.n_bl-1, 1e-100, pm.flag_error);
    //------------------------------------

    //Caluculate residuals
    //=====================================================================================
    //ion density
    //------------------------------------
    pm.error_rhoi = 0.0;
    if(pm.flag_error == 0){ //normalized maximum
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
            for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
                for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                    double error_tmp = fabs((gc.rhoi[i][j] - gc.rhoi_old[i][j])/(gc.rhoi[i][j] + gc.rhoi_old[i][j]+1e-100)*2.0/pm.CFL);
                    if(error_tmp > pm.error_rhoi){
                        pm.error_rhoi = error_tmp;
                    }
                }
            }
        }
    }else{ //normalized L2 norm
        double ncount = 0;
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
            for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
                for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                    pm.error_rhoi += pow(gc.rhoi[i][j]-gc.rhoi_old[i][j],2)/(pow(gc.rhoi_old[i][j],2)+1e10)/pm.CFL;
                    ncount++;
                }
            }
        }
        pm.error_rhoi = sqrt(pm.error_rhoi/double(ncount));
    }
    //------------------------------------
    
    //x-directional ion bulk velocity
    //------------------------------------
    pm.error_Uix = 0.0;
    if(pm.flag_error == 0){ //normalized maximum
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
            for (int i=gx.i_flx_bl[iblock][0];i<=gx.i_flx_bl[iblock][1] + 1;i++){
                for (int j=gx.j_flx_bl[iblock][0];j<=gx.j_flx_bl[iblock][1];j++){
                    double error_tmp = fabs((gx.Uix[i][j] - gx.Uix_old[i][j])/(gx.Uix[i][j] + gx.Uix_old[i][j]+1e-100)*2.0/pm.CFL);
                    if(error_tmp > pm.error_Uix){
                        pm.error_Uix = error_tmp;
                    }
                }
            }
        }
    }else{ //normalized L2 norm
        double ncount = 0;
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
            for (int i=gx.i_flx_bl[iblock][0];i<=gx.i_flx_bl[iblock][1] + 1;i++){
                for (int j=gx.j_flx_bl[iblock][0];j<=gx.j_flx_bl[iblock][1];j++){
                    pm.error_Uix += pow(gx.Uix[i][j]-gx.Uix_old[i][j],2)/(pow(gx.Uix_old[i][j],2)+1e-6)/pm.CFL;
                    ncount++;
                }
            }
        }
        pm.error_Uix = sqrt(pm.error_Uix/double(ncount));
    }
    //------------------------------------

    //r-directional ion bulk velocity
    //------------------------------------
    pm.error_Uir = 0.0;
    if(pm.flag_error == 0){ //normalized maximum
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){
            for (int i=gr.i_flr_bl[iblock][0];i<=gr.i_flr_bl[iblock][1];i++){ 
                for (int j=gr.j_flr_bl[iblock][0]+1;j<=gr.j_flr_bl[iblock][1];j++){
                    double error_tmp = fabs((gr.Uir[i][j] - gr.Uir_old[i][j])/(gr.Uir[i][j] + gr.Uir_old[i][j]+1e-100)*2.0/pm.CFL);
                    if(error_tmp > pm.error_Uir){
                        pm.error_Uir = error_tmp;
                    }
                }
            }
        }
    }else{ //normalized L2 norm
        double ncount = 0;
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){
            for (int i=gr.i_flr_bl[iblock][0];i<=gr.i_flr_bl[iblock][1];i++){ 
                for (int j=gr.j_flr_bl[iblock][0]+1;j<=gr.j_flr_bl[iblock][1];j++){
                    pm.error_Uir += pow(gr.Uir[i][j]-gr.Uir_old[i][j],2)/(pow(gr.Uir_old[i][j],2)+1e-6)/pm.CFL;
                    ncount++;
                }
            }
        }
        pm.error_Uir = sqrt(pm.error_Uir/double(ncount));
    }
    //------------------------------------

    //p-directional ion bulk velocity
    //------------------------------------
    pm.error_Uip = 0.0;
    if(pm.flag_error == 0){ //normalized maximum
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
            for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
                for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                    double error_tmp = fabs((gc.Uip[i][j] - gc.Uip_old[i][j])/(gc.Uip[i][j] + gc.Uip_old[i][j]+1e-100)*2.0/pm.CFL);
                    if(error_tmp > pm.error_Uip){
                        pm.error_Uip = error_tmp;
                    }
                }
            }
        }
    }else{ //normalized L2 norm
        double ncount = 0;
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){
            for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
                for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                    pm.error_Uip += pow(gc.Uip[i][j]-gc.Uip_old[i][j],2)/(pow(gc.Uip_old[i][j],2)+1e-6)/pm.CFL;
                    ncount++;
                }
            }
        }
        pm.error_Uip = sqrt(pm.error_Uip/double(ncount));
    }
    //------------------------------------
    //=====================================================================================
}


//*****************************************************************
//**                                                             **
//**           void solve_phi_couple_wdTe_wSEE_PC()              **
//**                                                             **
//*****************************************************************
void FluidModule::solve_phi_couple_wdTe_wSEE_PC(Params &pm,GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr, SeeVec &se){

    std::vector<std::vector<double> > aW(pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //係数
    std::vector<std::vector<double> > aE(pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //係数
    std::vector<std::vector<double> > aS(pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //係数
    std::vector<std::vector<double> > aN(pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //係数
    std::vector<std::vector<double> > aP(pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //係数
    std::vector<std::vector<double> > b (pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //ソース
    std::vector<std::vector<double> > phi_new (pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //ソース
    std::vector<std::vector<double> > phi_old_tmp(pm.ni+2,std::vector<double>(pm.nj+2,0.0));
    std::vector<std::vector<double> > nUex_old_tmp(pm.ni+2,std::vector<double>(pm.nj+2,0.0));
    std::vector<std::vector<double> > nUer_old_tmp(pm.ni+2,std::vector<double>(pm.nj+2,0.0));
    std::vector<std::vector<double> > dphi(pm.ni+2,std::vector<double>(pm.nj+2,0.0));
    std::vector<std::vector<double> > dnUex(pm.ni+2,std::vector<double>(pm.nj+2,0.0));
    std::vector<std::vector<double> > dnUer(pm.ni+2,std::vector<double>(pm.nj+2,0.0));

    std::vector<std::vector<double> > Ex_old(pm.ni+2,std::vector<double>(pm.nj+2,0.0));
    std::vector<std::vector<double> > Er_old(pm.ni+2,std::vector<double>(pm.nj+2,0.0));

    //Numerical treatment for exponential overflow and zero-division
    //------------------------------------
    double lim_exp = 1.0e2; //expの制限
    double EPS = 1.0e-3; //ゼロ割回避用の数字
    //------------------------------------

    //Acceleration factor for PC method
    //------------------------------------
    double c_mu = 1.0;
    //------------------------------------

    //convergence criteria
    //------------------------------------
    int flag_error2 = pm.flag_error;
    //------------------------------------

    //output duration
    //------------------------------------
    int ndiv_out = 100;
    //------------------------------------

    //scaling factor of wall flux
    //------------------------------------
    double alpha_wflux = 1.0;
    //------------------------------------

    //Reserve the past values
    //------------------------------------
    for (int i=0;i<=pm.ni+1;i++){ 
        for (int j=0;j<=pm.nj+1;j++){ 
            gc.phi_old[i][j] = gc.phi[i][j];
            Ex_old[i][j] = gx.Ex[i][j];
            Er_old[i][j] = gr.Er[i][j];
            gx.nUex_old[i][j] = gx.nUex[i][j];
            gr.nUer_old[i][j] = gr.nUer[i][j];

            //gx.nUex_old[i][j] = gx.rhoUex_old[i][j];
            //gr.nUer_old[i][j] = gr.rhoUer_old[i][j];
            //gx.nUex[i][j] = gx.rhoUex[i][j];
            //gr.nUer[i][j] = gr.rhoUer[i][j];
        }
    }
    //------------------------------------

    //Initial prediction of correction
    //------------------------------------
    for (int i=0;i<=pm.ni+1;i++){ 
        for (int j=0;j<=pm.nj+1;j++){ 
            dphi[i][j] = gc.phi[i][j]-gc.phi_old[i][j];
        }
    }
    //------------------------------------

    //set BCs for PC method
    //=====================================================================
    //left-wall-BC (z0) potential = 0
    //------------------------------------
    for (int j=gc.j_flc_bl[0][0];j<=gc.j_flc_bl[0][1];j++){
        int i=gc.i_flc_bl[0][0];
        gc.phi[i-1][j] = - gc.phi[i][j];
    }
    //------------------------------------

    //left-wall-BC (z1) potential = 0
    //------------------------------------
    for (int j=gc.j_flc_bl[2][0];j<=gc.j_flc_bl[0][0]-1;j++){
        int i=gc.i_flc_bl[2][0];
        gc.phi[i-1][j] = - gc.phi[i][j];
    }
    //------------------------------------

    //left-wall-BC (z2) potential = 0
    //------------------------------------
    for (int j=gc.j_flc_bl[3][0];j<=gc.j_flc_bl[2][0]-1;j++){
        int i=gc.i_flc_bl[3][0];
        gc.phi[i-1][j] = - gc.phi[i][j];
    }
    //------------------------------------

    //left-wall-BC (z4) potential = 0
    //------------------------------------
    for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[4][1];j++){
        int i=gc.i_flc_bl[4][0];
        gc.phi[i-1][j] = - gc.phi[i][j];
    }
    //------------------------------------

    //right-wall-BC (z3) potential = 0
    //------------------------------------
    for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[1][1];j++){
        int i=gc.i_flc_bl[1][1];
        gc.phi[i+1][j] = -gc.phi[i][j];
    }
    //------------------------------------
    
    //right-wall-BC (z5) potential = pm.V_bias
    //------------------------------------
    for (int j=gc.j_flc_bl[4][0];j<=gc.j_flc_bl[4][1];j++){
        int i=gc.i_flc_bl[4][1];
        gc.phi[i+1][j] = 2.0*pm.V_bias - gc.phi[i][j];
    }
    //------------------------------------


    //lower-wall-BC (x0) potential = 0
    //------------------------------------
    for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[0][1];i++){
        int j=gc.j_flc_bl[0][0];
        gc.phi[i][j-1] = - gc.phi[i][j];
    }
    //------------------------------------

    //lower-wall-BC (x2) potential = 0
    //------------------------------------
    for (int i=gc.i_flc_bl[2][0];i<=gc.i_flc_bl[2][1];i++){
        int j=gc.j_flc_bl[2][0];
        gc.phi[i][j-1] = - gc.phi[i][j];
    }
    //------------------------------------

    //upper-wall-BC (x1) potential = 0
    //------------------------------------
    for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[1][1];i++){
        int j=gc.j_flc_bl[0][1];
        gc.phi[i][j+1] = - gc.phi[i][j];
    }
    //------------------------------------

    //upper-wall-BC (x4) potential = 0
    //------------------------------------
    for (int i=gc.i_flc_bl[1][1]+1;i<=gc.i_flc_bl[3][1];i++){
        int j=gc.j_flc_bl[3][1];
        gc.phi[i][j+1] = - gc.phi[i][j];
    }
    //------------------------------------

    //upper-wall-BC (x5) grad(phi) = 0
    //------------------------------------
    for (int i=gc.i_flc_bl[4][0];i<=gc.i_flc_bl[4][1];i++){
        int j=gc.j_flc_bl[4][1];
        gc.phi[i][j+1] = gc.phi[i][j];
    }
    //------------------------------------
    //=====================================================================

    //PC method iteration
    //------------------------------------
    int ncount = 0; //number of iteration of PC method
    int flag_end = 0;
    do{
        
        ncount ++;
        pm.itime_PC_phi++;

        //Reserve previous-step values
        //------------------------------------
        for (int i=0;i<=pm.ni+1;i++){ 
            for (int j=0;j<=pm.nj+1;j++){ 
                phi_old_tmp[i][j] = gc.phi[i][j];
                nUex_old_tmp[i][j] = gx.nUex[i][j];
                nUer_old_tmp[i][j] = gr.nUer[i][j];
            }
        }
        //------------------------------------

        //Predict flux
        //==============================================
        //for rhoUex
        //------------------------------------
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
            for (int i=gx.i_flx_bl[iblock][0];i<=gx.i_flx_bl[iblock][1];i++){
                for (int j=gx.j_flx_bl[iblock][0];j<=gx.j_flx_bl[iblock][1];j++){
                    //radial position at upper and lower cell interface
                    //------------------------------------
                    double rL = (gc.r[j]+gc.r[j-1])/2.0;
                    double rR = (gc.r[j]+gc.r[j+1])/2.0;
                    //------------------------------------

                    //moblity at cell interface
                    //------------------------------------
                    double mu_para_tmp = (gc.mu_para[i][j] + gc.mu_para[i-1][j])/2.0;
                    double mu_perp_tmp = (gc.mu_perp[i][j] + gc.mu_perp[i-1][j])/2.0;
                    
                    double mu_star_tmp = mu_perp_tmp*c_mu;
                    //------------------------------------

                    //calculate inverse moblity tensor at cell interface
                    //------------------------------------
                    double mu_rr_tmp = (gc.mu_rr[i][j] + gc.mu_rr[i-1][j])/2.0;
                    double mu_xr_tmp = (gc.mu_xr[i][j] + gc.mu_xr[i-1][j])/2.0;

                    double mu_inv_xx_tmp = mu_rr_tmp/(mu_para_tmp*mu_perp_tmp);
                    double mu_inv_xr_tmp = -mu_xr_tmp/(mu_para_tmp*mu_perp_tmp);
                    //------------------------------------

                    //values at cell interface
                    //------------------------------------
                    double rhoe_tmp = (gc.rhoe[i][j] + gc.rhoe[i-1][j])/2.0;
                    double dphidx = (gc.phi[i][j] - gc.phi[i-1][j])/pm.dx;
                    double dpedx = ph::Boltz*(gc.rhoe[i][j]*gc.Te[i][j] - gc.rhoe[i-1][j]*gc.Te[i-1][j])/pm.dx;
                    //------------------------------------
                
                    //BC - set adjacent flag (adjacent = 0, others = 1)
                    //adjacent "cell" to the wall for setting zero-flux on the wall
                    //------------------------------------
                    //lower
                    double bLr_wall  = double(j!=gc.j_flc_bl[0][0] || i>gc.i_flc_bl[0][1]); //x0
                    bLr_wall = bLr_wall*double(j!=gc.j_flc_bl[2][0] || i>gc.i_flc_bl[2][1]); //x2
                    double bLr_cen = double(j!=gc.j_flc_bl[3][0]); //x6
                    double bLr  = bLr_wall*bLr_cen;

                    //lower-left
                    double bLrD_wall  = double(j!=gc.j_flc_bl[0][0] || i-1>gc.i_flc_bl[0][1]); //x0
                    bLrD_wall = bLrD_wall*double(j!=gc.j_flc_bl[2][0] || i-1>gc.i_flc_bl[2][1]); //x2
                    double bLrD_cen = double(j!=gc.j_flc_bl[3][0]); //x6
                    double bLrD  = bLrD_wall*bLrD_cen;

                    //upper
                    double bRr_wall  = double(j!=gc.j_flc_bl[0][1] || iblock > 2); //x1
                    bRr_wall = bRr_wall*double(j!=gc.j_flc_bl[3][1] || i<gc.i_flc_bl[1][1]+1 || i>gc.i_flc_bl[3][1]); //x4
                    //double bRr_open  = (j!=gc.j_flc_bl[4][1] || iblock != 4); //x5
                    double bRr = bRr_wall;

                    //upper-left
                    double bRrD_wall  = double(j!=gc.j_flc_bl[0][1] || iblock > 2); //x1
                    bRrD_wall = bRrD_wall*double(j!=gc.j_flc_bl[3][1] || i-1<gc.i_flc_bl[1][1]+1 || i-1>gc.i_flc_bl[3][1]); //x4
                    //double bRrD_open  = (j!=gc.j_flc_bl[4][1] || iblock != 4); //x5
                    double bRrD = bRrD_wall;
                    //------------------------------------

                    //interpolation of transverse flux component (should be old value)
                    //------------------------------------
                    double nUer_tmp = (rR*nUer_old_tmp[i][j+1]*bRr + rL*nUer_old_tmp[i][j]*bLr
                        + rR*nUer_old_tmp[i-1][j+1]*bRrD + rL*nUer_old_tmp[i-1][j]*bLrD)/(4.0*gc.r[j]);
                    //------------------------------------

                    //update flux
                    //------------------------------------
                    gx.nUex[i][j] = nUex_old_tmp[i][j] 
                        - (mu_inv_xx_tmp*nUex_old_tmp[i][j] + mu_inv_xr_tmp*nUer_tmp)*mu_star_tmp 
                        + mu_star_tmp*(rhoe_tmp*dphidx - dpedx/ph::e0);
                    //------------------------------------
                }
            }     
        }
        //------------------------------------

        //left-wall-BC (z0)
        //------------------------------------
        for (int j=gc.j_flc_bl[0][0];j<=gc.j_flc_bl[0][1];j++){
            int i=gc.i_flc_bl[0][0];
            double vth = sqrt(8.0*ph::Boltz*gc.Te[i][j]/(M_PI*pm.masse));
            double delta_phi = fmax(-Ex_old[i][j]*pm.dx/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerN_wall_Lx = -0.25*alpha_wflux*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Lx = -pm.coefIISEE_ts*gc.rhoi[i][j]*gx.Uix[i][j] //IISEE
                            -(coefEISEE_eb + coefEISEE_rd + coefEISEE_ts)*gc.rhoe_old[i][j]*GammaPerN_wall_Lx //EISEE
                            -pm.coefMISEE_ts*gx.rhoUmx[i][j]; //MISEE
            gx.nUex[i][j] = gc.rhoe[i][j]*GammaPerN_wall_Lx + Gamma_wallSEE_Lx;
        }
        //------------------------------------
    
        //left-wall-BC (z1)
        //------------------------------------
        for (int j=gc.j_flc_bl[2][0];j<=gc.j_flc_bl[0][0]-1;j++){
            int i=gc.i_flc_bl[2][0];
            double vth = sqrt(8.0*ph::Boltz*gc.Te[i][j]/(M_PI*pm.masse));
            double delta_phi = fmax(-Ex_old[i][j]*pm.dx/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerN_wall_Lx = -0.25*alpha_wflux*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Lx = -pm.coefIISEE_ts*gc.rhoi[i][j]*gx.Uix[i][j] //IISEE
                            -(coefEISEE_eb + coefEISEE_rd + coefEISEE_ts)*gc.rhoe_old[i][j]*GammaPerN_wall_Lx //EISEE
                            -pm.coefMISEE_ts*gx.rhoUmx[i][j]; //MISEE
            gx.nUex[i][j] = gc.rhoe[i][j]*GammaPerN_wall_Lx + Gamma_wallSEE_Lx;
        }
        //------------------------------------
    
        //left-wall-BC (z2)
        //------------------------------------
        for (int j=gc.j_flc_bl[3][0];j<=gc.j_flc_bl[2][0]-1;j++){
            int i=gc.i_flc_bl[3][0];
            double vth = sqrt(8.0*ph::Boltz*gc.Te[i][j]/(M_PI*pm.masse));
            double delta_phi = fmax(-Ex_old[i][j]*pm.dx/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerN_wall_Lx = -0.25*alpha_wflux*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Lx = -pm.coefIISEE_ts*gc.rhoi[i][j]*gx.Uix[i][j] //IISEE
                            -(coefEISEE_eb + coefEISEE_rd + coefEISEE_ts)*gc.rhoe_old[i][j]*GammaPerN_wall_Lx //EISEE
                            -pm.coefMISEE_ts*gx.rhoUmx[i][j]; //MISEE
            gx.nUex[i][j] = gc.rhoe[i][j]*GammaPerN_wall_Lx + Gamma_wallSEE_Lx;
        }
        //------------------------------------
    
        //left-wall-BC (z4)
        //------------------------------------
        for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[4][1];j++){
            int i=gc.i_flc_bl[4][0];
            double vth = sqrt(8.0*ph::Boltz*gc.Te[i][j]/(M_PI*pm.masse));
            double delta_phi = fmax(-Ex_old[i][j]*pm.dx/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerN_wall_Lx = -0.25*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Lx = -pm.coefIISEE_ts*gc.rhoi[i][j]*gx.Uix[i][j] //IISEE
                            -(coefEISEE_eb + coefEISEE_rd + coefEISEE_ts)*gc.rhoe_old[i][j]*GammaPerN_wall_Lx //EISEE
                            -pm.coefMISEE_ts*gx.rhoUmx[i][j]; //MISEE
            gx.nUex[i][j] = gc.rhoe[i][j]*GammaPerN_wall_Lx + Gamma_wallSEE_Lx;
        }
        //------------------------------------
    
        //right-wall-BC (z3)
        //------------------------------------
        for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[1][1];j++){
            int i=gc.i_flc_bl[1][1];
            double vth = sqrt(8.0*ph::Boltz*gc.Te[i][j]/(M_PI*pm.masse));
            double delta_phi =  fmax(Ex_old[i+1][j]*pm.dx/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerN_wall_Rx = 0.25*alpha_wflux*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Rx = -pm.coefIISEE_ts*gc.rhoi[i][j]*gx.Uix[i+1][j] //IISEE
                            -(coefEISEE_eb + coefEISEE_rd + coefEISEE_ts)*gc.rhoe_old[i][j]*GammaPerN_wall_Rx //EISEE
                            -pm.coefMISEE_ts*gx.rhoUmx[i+1][j]; //MISEE
            gx.nUex[i+1][j] =  gc.rhoe[i][j]*GammaPerN_wall_Rx + Gamma_wallSEE_Rx;
        }
        //------------------------------------
    
        //right-wall-BC (z5)
        //------------------------------------
        for (int j=gc.j_flc_bl[4][0];j<=gc.j_flc_bl[4][1];j++){
            int i=gc.i_flc_bl[4][1];
            double vth = sqrt(8.0*ph::Boltz*gc.Te[i][j]/(M_PI*pm.masse));
            double delta_phi =  fmax(Ex_old[i+1][j]*pm.dx/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerN_wall_Rx = 0.25*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Rx = -pm.coefIISEE_ts*gc.rhoi[i][j]*gx.Uix[i+1][j] //IISEE
                            -(coefEISEE_eb + coefEISEE_rd + coefEISEE_ts)*gc.rhoe_old[i][j]*GammaPerN_wall_Rx //EISEE
                            -pm.coefMISEE_ts*gx.rhoUmx[i+1][j]; //MISEE
            gx.nUex[i+1][j] =  gc.rhoe[i][j]*GammaPerN_wall_Rx + Gamma_wallSEE_Rx;
        }
        //=====================================================================================

        //for rhoUer
        //------------------------------------
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){
            for (int i=gr.i_flr_bl[iblock][0];i<=gr.i_flr_bl[iblock][1];i++){ 
                for (int j=gr.j_flr_bl[iblock][0];j<=gr.j_flr_bl[iblock][1];j++){
                    
                    //radial position at cell interface
                    //------------------------------------
                    double r_tmp = (gc.r[j] + gc.r[j-1])/2.0;
                    //------------------------------------

                    //moblity at cell interface
                    //------------------------------------
                    double mu_para_tmp = (gc.mu_para[i][j] + gc.mu_para[i][j-1])/2.0;
                    double mu_perp_tmp = (gc.mu_perp[i][j] + gc.mu_perp[i][j-1])/2.0;

                    double mu_star_tmp = mu_perp_tmp*c_mu;
                    //------------------------------------

                    //calculate inverse moblity tensor at cell interface
                    //------------------------------------
                    double mu_xx_tmp = (gc.mu_xx[i][j] + gc.mu_xx[i][j-1])/2.0;
                    double mu_xr_tmp = (gc.mu_xr[i][j] + gc.mu_xr[i][j-1])/2.0;
                    
                    double mu_inv_rr_tmp = mu_xx_tmp/(mu_para_tmp*mu_perp_tmp);
                    double mu_inv_xr_tmp = -mu_xr_tmp/(mu_para_tmp*mu_perp_tmp);
                    //------------------------------------
    
                    //values at cell interface
                    //------------------------------------
                    double rhoe_tmp = (gc.rhoe[i][j] + gc.rhoe[i][j-1])/2.0;
                    double dphidr = (gc.phi[i][j] - gc.phi[i][j-1])/pm.dr;
                    double dpedr = ph::Boltz*(gc.rhoe[i][j]*gc.Te[i][j] - gc.rhoe[i][j-1]*gc.Te[i][j-1])/pm.dr;
                    //------------------------------------

                    //BC - set adjacent flag (adjacent = 0, others = 1)
                    //adjacent "cell" to the wall for setting zero-flux on the wall
                    //------------------------------------
                    //left
                    double bLx_wall  = double(i!=gc.i_flc_bl[0][0] || iblock != 0); //z0
                    bLx_wall = bLx_wall*double(i!=gc.i_flc_bl[2][0] || j>gc.j_flc_bl[0][0]-1); //z1
                    bLx_wall = bLx_wall*double(i!=gc.i_flc_bl[3][0] || j>gc.j_flc_bl[2][0]-1); //z2
                    bLx_wall = bLx_wall*double(i!=gc.i_flc_bl[4][0] || j<gc.j_flc_bl[1][0]); //z4
                    double bLx  = bLx_wall;

                    //left-lower
                    double bLxD_wall  = double(i!=gc.i_flc_bl[0][0] || iblock != 0); //z0
                    bLxD_wall = bLxD_wall*double(i!=gc.i_flc_bl[2][0] || j-1>gc.j_flc_bl[0][0]-1); //z1
                    bLxD_wall = bLxD_wall*double(i!=gc.i_flc_bl[3][0] || j-1>gc.j_flc_bl[2][0]-1); //z2
                    bLxD_wall = bLxD_wall*double(i!=gc.i_flc_bl[4][0] || j-1<gc.j_flc_bl[1][0]); //z4
                    double bLxD  = bLxD_wall;

                    //right
                    double bRx_wall  = double(i!=gc.i_flc_bl[1][1] || j<gc.j_flc_bl[1][0]); //z3
                    bRx_wall = bRx_wall*double(i!=gc.i_flc_bl[4][1]); //z5
                    double bRx  = bRx_wall;

                    //right-lower
                    double bRxD_wall  = double(i!=gc.i_flc_bl[1][1] || j-1<gc.j_flc_bl[1][0]); //z3
                    bRxD_wall = bRxD_wall*double(i!=gc.i_flc_bl[4][1]); //z5
                    double bRxD  = bRxD_wall;
                    //------------------------------------

                    //interpolation of transverse flux component (should be old value)
                    //------------------------------------
                    double nUex_tmp = (gc.r[j]*nUex_old_tmp[i+1][j]*bRx + gc.r[j]*nUex_old_tmp[i][j]*bLx 
                        + gc.r[j-1]*nUex_old_tmp[i+1][j-1]*bRxD + gc.r[j-1]*nUex_old_tmp[i][j-1]*bLxD)/(4.0*r_tmp);
                    //------------------------------------

                    //update flux
                    //------------------------------------
                    gr.nUer[i][j] =  nUer_old_tmp[i][j] 
                        - (mu_inv_xr_tmp*nUex_tmp + mu_inv_rr_tmp*nUer_old_tmp[i][j])*mu_star_tmp 
                        + mu_star_tmp*(rhoe_tmp*dphidr - dpedr/ph::e0);
                    //------------------------------------
  
                }
            }
        }
        //------------------------------------

        //lower-wall-BC (x0)
        //------------------------------------
        for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[0][1];i++){
            int j=gc.j_flc_bl[0][0];
            double vth = sqrt(8.0*ph::Boltz*gc.Te[i][j]/(M_PI*pm.masse));
            double delta_phi = fmax(-Er_old[i][j]*pm.dr/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerN_wall_Lr = -0.25*alpha_wflux*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Lr = -pm.coefIISEE_ts*gc.rhoi[i][j]*gr.Uir[i][j] //IISEE
                            -(coefEISEE_eb + coefEISEE_rd + coefEISEE_ts)*gc.rhoe_old[i][j]*GammaPerN_wall_Lr //EISEE
                            -pm.coefMISEE_ts*gr.rhoUmr[i][j]; //MISEE
            gr.nUer[i][j] = gc.rhoe[i][j]*GammaPerN_wall_Lr + Gamma_wallSEE_Lr;
        }
        //------------------------------------

        //lower-wall-BC (x2)
        //------------------------------------
        for (int i=gc.i_flc_bl[2][0];i<=gc.i_flc_bl[2][1];i++){
            int j=gc.j_flc_bl[2][0];
            double vth = sqrt(8.0*ph::Boltz*gc.Te[i][j]/(M_PI*pm.masse));
            double delta_phi = fmax(-Er_old[i][j]*pm.dr/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerN_wall_Lr = -0.25*alpha_wflux*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Lr = -pm.coefIISEE_ts*gc.rhoi[i][j]*gr.Uir[i][j] //IISEE
                            -(coefEISEE_eb + coefEISEE_rd + coefEISEE_ts)*gc.rhoe_old[i][j]*GammaPerN_wall_Lr //EISEE
                            -pm.coefMISEE_ts*gr.rhoUmr[i][j]; //MISEE
            gr.nUer[i][j] = gc.rhoe[i][j]*GammaPerN_wall_Lr + Gamma_wallSEE_Lr;
        }
        //------------------------------------

        //upper-wall-BC (x1)
        //------------------------------------
        for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[1][1];i++){
            int j=gc.j_flc_bl[0][1];
            double vth = sqrt(8.0*ph::Boltz*gc.Te[i][j]/(M_PI*pm.masse));
            double delta_phi =  fmax(Er_old[i][j+1]*pm.dr/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerN_wall_Rr = 0.25*alpha_wflux*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Rr = -pm.coefIISEE_ts*gc.rhoi[i][j]*gr.Uir[i][j+1] //IISEE
                            -(coefEISEE_eb + coefEISEE_rd + coefEISEE_ts)*gc.rhoe_old[i][j]*GammaPerN_wall_Rr //EISEE
                            -pm.coefMISEE_ts*gr.rhoUmr[i][j+1]; //MISEE
            gr.nUer[i][j+1] = gc.rhoe[i][j]*GammaPerN_wall_Rr + Gamma_wallSEE_Rr;
        }
        //------------------------------------

        //upper-wall-BC (x4)
        //------------------------------------
        for (int i=gc.i_flc_bl[1][1]+1;i<=gc.i_flc_bl[3][1];i++){
            int j=gc.j_flc_bl[3][1];
            double vth = sqrt(8.0*ph::Boltz*gc.Te[i][j]/(M_PI*pm.masse));
            double delta_phi =  fmax(Er_old[i][j+1]*pm.dr/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerN_wall_Rr = 0.25*alpha_wflux*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Rr = -pm.coefIISEE_ts*gc.rhoi[i][j]*gr.Uir[i][j+1] //IISEE
                            -(coefEISEE_eb + coefEISEE_rd + coefEISEE_ts)*gc.rhoe_old[i][j]*GammaPerN_wall_Rr //EISEE
                            -pm.coefMISEE_ts*gr.rhoUmr[i][j+1]; //MISEE
            gr.nUer[i][j+1] = gc.rhoe[i][j]*GammaPerN_wall_Rr + Gamma_wallSEE_Rr;
        }
        //------------------------------------

        //upper-open-BC (x5)
        //------------------------------------
        for (int i=gc.i_flc_bl[4][0];i<=gc.i_flc_bl[4][1];i++){
            int j=gc.j_flc_bl[4][1];
            double rhoUir_Rr = gc.rhoi[i][j]*gr.Uir[i][j+1];
            double Gamma_open_Rr = rhoUir_Rr;
            gr.nUer[i][j+1] =  Gamma_open_Rr;
        }
        //------------------------------------

        //Make coefficient for Poisson equation for 5 point stencil
        //==============================================
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
            for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
                for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                    //radial position at upper and lower cell interface
                    //------------------------------------
                    double rR = (gc.r[j]+gc.r[j+1])/2.0;
                    double rL = (gc.r[j]+gc.r[j-1])/2.0;
                    //------------------------------------

                    //normalized radial position at upper and lower cell interface
                    //------------------------------------
                    double qR = (gc.r[j]+gc.r[j+1])/2.0/gc.r[j];
                    double qL = (gc.r[j]+gc.r[j-1])/2.0/gc.r[j];
                    //------------------------------------

                    //moblity at cell interface
                    //------------------------------------
                    double mu_star_Lx = (gc.mu_perp[i][j] + gc.mu_perp[i-1][j])/2.0*c_mu;
                    double mu_star_Rx = (gc.mu_perp[i][j] + gc.mu_perp[i+1][j])/2.0*c_mu;
                    double mu_star_Lr = (gc.mu_perp[i][j] + gc.mu_perp[i][j-1])/2.0*c_mu;
                    double mu_star_Rr = (gc.mu_perp[i][j] + gc.mu_perp[i][j+1])/2.0*c_mu;
                    //------------------------------------

                    //electron density at cell interface
                    //------------------------------------
                    double rhoe_Lx = (gc.rhoe[i][j] + gc.rhoe[i-1][j])/2.0;
                    double rhoe_Rx = (gc.rhoe[i][j] + gc.rhoe[i+1][j])/2.0;
                    double rhoe_Lr = (gc.rhoe[i][j] + gc.rhoe[i][j-1])/2.0;
                    double rhoe_Rr = (gc.rhoe[i][j] + gc.rhoe[i][j+1])/2.0;
                    //------------------------------------

                    //BC - set adjacent flag (adjacent = 0, others = 1)
                    //adjacent "cell" to the wall, open, or central axis BC
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
                    //==============================================

                    //coefficient for Poisson equation for 5 point stencil
                    //==============================================
                    //vertical and horizental 5 stencil
                    //------------------------------------
                    aP[i][j] =  (    (ph::eps0 + ph::e0*pm.dt*rhoe_Rx*mu_star_Rx*bRx) 
                                 +   (ph::eps0 + ph::e0*pm.dt*rhoe_Lx*mu_star_Lx*bLx)
                                )*pm.dr*pm.dr 
                              + ( qR*(ph::eps0 + ph::e0*pm.dt*rhoe_Rr*mu_star_Rr*bRr) 
                                + qL*(ph::eps0 + ph::e0*pm.dt*rhoe_Lr*mu_star_Lr*bLr)
                                )*pm.dx*pm.dx;

                    aE[i][j] =    (ph::eps0 + ph::e0*pm.dt*rhoe_Rx*mu_star_Rx*bRx)*pm.dr*pm.dr;
                    aW[i][j] =    (ph::eps0 + ph::e0*pm.dt*rhoe_Lx*mu_star_Lx*bLx)*pm.dr*pm.dr;
                    aN[i][j] = qR*(ph::eps0 + ph::e0*pm.dt*rhoe_Rr*mu_star_Rr*bRr)*pm.dx*pm.dx;
                    aS[i][j] = qL*(ph::eps0 + ph::e0*pm.dt*rhoe_Lr*mu_star_Lr*bLr)*pm.dx*pm.dx;
                    //------------------------------------

                    //RHS
                    //------------------------------------
                    double nabla_nUe_tmp = (gx.nUex[i+1][j]-gx.nUex[i][j])/pm.dx + (qR*gr.nUer[i][j+1]-qL*gr.nUer[i][j])/pm.dr;
                    double lapPhi =  -(gx.Ex[i+1][j] - gx.Ex[i][j])/pm.dx - (qR*gr.Er[i][j+1] - qL*gr.Er[i][j])/pm.dr;
                    
                    b[i][j] = (  ph::e0*gc.rhoi[i][j] - ph::e0*gc.rhoe[i][j] //rhoeはまだ更新されていないのでrhoe_oldでなくrhoeを使う
                               - ph::e0*pm.dt*gc.rate_ionize[i][j] 
                               + ph::eps0*lapPhi + ph::e0*pm.dt*nabla_nUe_tmp
                              )*pm.dx*pm.dx*pm.dr*pm.dr;
                    //------------------------------------
                }
            }
        }
        //==============================================

        //Set boundary condition
        //==============================================
        //left-wall-BC (z0) potential = 0
        //------------------------------------
        for (int j=gc.j_flc_bl[0][0];j<=gc.j_flc_bl[0][1];j++){
            int i=gc.i_flc_bl[0][0];
            //modify coefficient
            //------------------------------------
            aP[i][j] = aP[i][j] + aW[i][j];
            aW[i][j] = 0.0;
            //------------------------------------
        }
        //------------------------------------

        //left-wall-BC (z1) potential = 0
        //------------------------------------
        for (int j=gc.j_flc_bl[2][0];j<=gc.j_flc_bl[0][0]-1;j++){
            int i=gc.i_flc_bl[2][0];
            aP[i][j] = aP[i][j] + aW[i][j];
            aW[i][j] = 0.0;
        }
        //------------------------------------

        //left-wall-BC (z2) potential = 0
        //------------------------------------
        for (int j=gc.j_flc_bl[3][0];j<=gc.j_flc_bl[2][0]-1;j++){
            int i=gc.i_flc_bl[3][0];
            aP[i][j] = aP[i][j] + aW[i][j];
            aW[i][j] = 0.0;
        }
        //------------------------------------

        //left-wall-BC (z4) potential = 0
        //------------------------------------
        for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[4][1];j++){
            int i=gc.i_flc_bl[4][0];
            aP[i][j] = aP[i][j] + aW[i][j];
            aW[i][j] = 0.0;
        }
        //------------------------------------

        //right-wall-BC (z3) continuity potential = 0
        //------------------------------------
        for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[1][1];j++){
            int i=gc.i_flc_bl[1][1];
            aP[i][j] = aP[i][j] + aE[i][j];
            aE[i][j] = 0.0;
        }
        //------------------------------------

        //right-wall-BC (z5) continuity potential = pm.V_bias
        //------------------------------------
        for (int j=gc.j_flc_bl[4][0];j<=gc.j_flc_bl[4][1];j++){
            int i=gc.i_flc_bl[4][1];
            aP[i][j] = aP[i][j] + aE[i][j];
            aE[i][j] = 0.0;
        }
        //------------------------------------

        //lower-wall-BC (x0) potential = 0
        //------------------------------------
        for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[0][1];i++){
            int j=gc.j_flc_bl[0][0];
            aP[i][j] = aP[i][j] + aS[i][j];
            aS[i][j] = 0.0;
        }
        //------------------------------------

        //lower-wall-BC (x2) potential = 0
        //------------------------------------
        for (int i=gc.i_flc_bl[2][0];i<=gc.i_flc_bl[2][1];i++){
            int j=gc.j_flc_bl[2][0];
            aP[i][j] = aP[i][j] + aS[i][j];
            aS[i][j] = 0.0;
        }
        //------------------------------------

        //upper-wall-BC (x1) potential = 0
        //------------------------------------
        for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[1][1];i++){
            int j=gc.j_flc_bl[0][1];
            aP[i][j] = aP[i][j] + aN[i][j];
            aN[i][j] = 0.0;
        }
        //------------------------------------

        //upper-wall-BC (x4) potential = 0
        //------------------------------------
        for (int i=gc.i_flc_bl[1][1]+1;i<=gc.i_flc_bl[3][1];i++){
            int j=gc.j_flc_bl[3][1];
            aP[i][j] = aP[i][j] + aN[i][j];
            aN[i][j] = 0.0;
        }
        //------------------------------------

        //upper-open-BC (x5) grad(phi) = 0
        //------------------------------------
        for (int i=gc.i_flc_bl[4][0];i<=gc.i_flc_bl[4][1];i++){
            int j=gc.j_flc_bl[4][1];
            aP[i][j] = aP[i][j] - aN[i][j];
            aN[i][j] = 0.0;
        }
        //------------------------------------
        //==============================================

        //Solver
        //------------------------------------
        if(pm.flag_iter_phi == 0){
            //SOR
            //------------------------------------
            double alpha = 1.9;
            solver_SOR(aP,aE,aW,aN,aS,b,gc.i_flc_bl,gc.j_flc_bl,pm.n_bl-1,alpha,pm.maxITR_SOR_phi,pm.error_cnv_SOR_phi,0,dphi);
            //------------------------------------
        }else{
            //Semi-coarse Multi-Grid by Hypre
            //------------------------------------
            solver_SMG(aP,aE,aW,aN,aS,b,gc.i_flc_bl,gc.j_flc_bl,pm.n_bl-1,pm.maxITR_SOR_phi,pm.error_cnv_SOR_phi,0,dphi);
            //------------------------------------
        }
        //------------------------------------

        //Reflect BCs for ghost cell
        //==============================================
        
        //set BCs for PC method
        //=====================================================================
        //left-wall-BC (z0) potential = 0
        //------------------------------------
        for (int j=gc.j_flc_bl[0][0];j<=gc.j_flc_bl[0][1];j++){
            int i=gc.i_flc_bl[0][0];
            dphi[i-1][j] = - dphi[i][j];
            gx.Ex[i][j] = gx.Ex[i][j] - (dphi[i][j]-dphi[i-1][j])/pm.dx;
        }
        //------------------------------------

        //left-wall-BC (z1) potential = 0
        //------------------------------------
        for (int j=gc.j_flc_bl[2][0];j<=gc.j_flc_bl[0][0]-1;j++){
            int i=gc.i_flc_bl[2][0];
            dphi[i-1][j] = - dphi[i][j];
            gx.Ex[i][j] = gx.Ex[i][j] - (dphi[i][j]-dphi[i-1][j])/pm.dx;
        }
        //------------------------------------

        //left-wall-BC (z2) potential = 0
        //------------------------------------
        for (int j=gc.j_flc_bl[3][0];j<=gc.j_flc_bl[2][0]-1;j++){
            int i=gc.i_flc_bl[3][0];
            dphi[i-1][j] = - dphi[i][j];
            gx.Ex[i][j] = gx.Ex[i][j] - (dphi[i][j]-dphi[i-1][j])/pm.dx;
        }
        //------------------------------------

        //left-wall-BC (z4) potential = 0
        //------------------------------------
        for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[4][1];j++){
            int i=gc.i_flc_bl[4][0];
            dphi[i-1][j] = - dphi[i][j];
            gx.Ex[i][j] = gx.Ex[i][j] - (dphi[i][j]-dphi[i-1][j])/pm.dx;
        }
        //------------------------------------

        //right-wall-BC (z3) potential = 0
        //------------------------------------
        for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[1][1];j++){
            int i=gc.i_flc_bl[1][1];
            dphi[i+1][j] = - dphi[i][j];
            gx.Ex[i+1][j] = gx.Ex[i+1][j] - (dphi[i+1][j]-dphi[i][j])/pm.dx;
        }
        //------------------------------------

        //right-wall-BC (z5) potential = pm.V_bias
        //------------------------------------
        for (int j=gc.j_flc_bl[4][0];j<=gc.j_flc_bl[4][1];j++){
            int i=gc.i_flc_bl[4][1];
            dphi[i+1][j] = - dphi[i][j];
            gx.Ex[i+1][j] = gx.Ex[i+1][j] - (dphi[i+1][j]-dphi[i][j])/pm.dx;
        }
        //------------------------------------


        //lower-wall-BC (x0) potential = 0
        //------------------------------------
        for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[0][1];i++){
            int j=gc.j_flc_bl[0][0];
            dphi[i][j-1] = - dphi[i][j];
            gr.Er[i][j] = gr.Er[i][j] - (dphi[i][j]-dphi[i][j-1])/pm.dr;
        }
        //------------------------------------

        //lower-wall-BC (x2) potential = 0
        //------------------------------------
        for (int i=gc.i_flc_bl[2][0];i<=gc.i_flc_bl[2][1];i++){
            int j=gc.j_flc_bl[2][0];
            dphi[i][j-1] = - dphi[i][j];
            gr.Er[i][j] = gr.Er[i][j] - (dphi[i][j]-dphi[i][j-1])/pm.dr;
        }
        //------------------------------------

        //upper-wall-BC (x1) potential = 0
        //------------------------------------
        for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[1][1];i++){
            int j=gc.j_flc_bl[0][1];
            dphi[i][j+1] = - dphi[i][j];
            gr.Er[i][j+1] = gr.Er[i][j+1] - (dphi[i][j+1]-dphi[i][j])/pm.dr;
        }
        //------------------------------------

        //upper-wall-BC (x4) potential = 0
        //------------------------------------
        for (int i=gc.i_flc_bl[1][1]+1;i<=gc.i_flc_bl[3][1];i++){
            int j=gc.j_flc_bl[3][1];
            dphi[i][j+1] = - dphi[i][j];
            gr.Er[i][j+1] = gr.Er[i][j+1] - (dphi[i][j+1]-dphi[i][j])/pm.dr;
        }
        //------------------------------------

        //upper-open-BC (x5) grad(phi) = 0
        //------------------------------------
        for (int i=gc.i_flc_bl[4][0];i<=gc.i_flc_bl[4][1];i++){
            int j=gc.j_flc_bl[4][1];
            dphi[i][j+1] = dphi[i][j];
            gr.Er[i][j+1] = 0.0;
        }
        //------------------------------------
        //=====================================================================

        //Correction
        //==============================================
        //for phi
        //------------------------------------
        for(int i=0;i<=pm.ni+1;i++){
            for(int j=0;j<=pm.nj+1;j++){
                gc.phi[i][j] = gc.phi[i][j] + dphi[i][j];
            }
        }
        //------------------------------------

        //Update Ex
        //------------------------------------
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){
            for (int i=gx.i_flx_bl[iblock][0];i<=gx.i_flx_bl[iblock][1];i++){
                for (int j=gx.j_flx_bl[iblock][0];j<=gx.j_flx_bl[iblock][1];j++){
                    gx.Ex[i][j] = - (gc.phi[i][j]-gc.phi[i-1][j])/pm.dx;
                }
            }
        }
        //------------------------------------

        //Update Er
        //------------------------------------
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){
            for (int i=gr.i_flr_bl[iblock][0];i<=gr.i_flr_bl[iblock][1];i++){
                for (int j=gr.j_flr_bl[iblock][0];j<=gr.j_flr_bl[iblock][1];j++){
                    gr.Er[i][j] = - (gc.phi[i][j]-gc.phi[i][j-1])/pm.dr;
                }
            }
        }
        //------------------------------------

        //for rhoUex
        //------------------------------------
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
            for (int i=gx.i_flx_bl[iblock][0];i<=gx.i_flx_bl[iblock][1];i++){
                for (int j=gx.j_flx_bl[iblock][0];j<=gx.j_flx_bl[iblock][1];j++){
                    
                    //values at cell interface
                    //------------------------------------
                    double mu_star_tmp = (gc.mu_perp[i][j] + gc.mu_perp[i-1][j])/2.0*c_mu;
                    double rhoe_tmp = (gc.rhoe[i][j] + gc.rhoe[i-1][j])/2.0;
                    double ddphidx = (dphi[i][j] - dphi[i-1][j])/pm.dx;
                    //------------------------------------
                    
                    //calculate flux correction
                    //------------------------------------
                    dnUex[i][j] = mu_star_tmp*rhoe_tmp*ddphidx;
                    //------------------------------------

                    //Update
                    //------------------------------------
                    gx.nUex[i][j] = gx.nUex[i][j] + dnUex[i][j];
                    //------------------------------------
                }
            }
        }
        //------------------------------------

        //for rhoUer
        //------------------------------------
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){
            for (int i=gr.i_flr_bl[iblock][0];i<=gr.i_flr_bl[iblock][1];i++){ 
                for (int j=gr.j_flr_bl[iblock][0];j<=gr.j_flr_bl[iblock][1];j++){
                    
                    //values at cell interface
                    //------------------------------------
                    double mu_star_tmp = (gc.mu_perp[i][j] + gc.mu_perp[i][j-1])/2.0*c_mu;
                    double rhoe_tmp = (gc.rhoe[i][j] + gc.rhoe[i][j-1])/2.0;
                    double ddphidr = (dphi[i][j] - dphi[i][j-1])/pm.dr;
                    //------------------------------------

                    //calculate flux correction
                    //------------------------------------
                    dnUer[i][j] =  mu_star_tmp*rhoe_tmp*ddphidr;
                    //------------------------------------

                    //Update
                    //------------------------------------
                    gr.nUer[i][j] = gr.nUer[i][j] + dnUer[i][j];
                    //------------------------------------
                }
            }
        }
        //------------------------------------
        
        
        //==============================================
        

        //Caluculate residuals
        //------------------------------------
        double error_phi2 = 0.0;
        double error_nUex2 = 0.0;
        double error_nUer2 = 0.0;
        calcRes(error_phi2, gc.phi,  phi_old_tmp,  gc.i_flc_bl, gc.j_flc_bl, pm.n_bl-1,1e-100,flag_error2);
        calcRes(error_nUex2,gx.nUex, nUex_old_tmp, gx.i_flx_bl, gx.j_flx_bl, pm.n_bl-1,1e-100,flag_error2);
        calcRes(error_nUer2,gr.nUer, nUer_old_tmp, gr.i_flr_bl, gr.j_flr_bl, pm.n_bl-1,1e-100,flag_error2);
        
        if(pm.flag_chk == 1) output_residual_for_PC("phi",pm.gtime, pm.itime, error_phi2, error_nUex2, error_nUer2, pm.itime_PC_phi, ncount);
        
        double error_global = fmax(fmax(error_phi2,error_nUex2),error_nUer2);
        
        if(ncount % ndiv_out == 0){
            std::cout << "ncount = "<<ncount<<std::endl;
            std::cout << "error_phi = " << error_phi2
                <<  " error_nUex = " << error_nUex2
                <<  " error_nUer= " << error_nUer2
                <<  " error_max = " << fmax(fmax(error_phi2,error_nUex2),error_nUer2)
                << std::endl;
        }

        if(error_global < pm.error_cnv_HES_phi || ncount >= pm.maxITR_HES_phi){
            flag_end = 1;
        }

    } while (flag_end == 0);


    //Caluculate residuals
    //------------------------------------
    calcRes(pm.error_phi, gc.phi,gc.phi_old,gc.i_flc_bl, gc.j_flc_bl, pm.n_bl-1, 1e-100, pm.flag_error);
    calcRes(pm.error_nUex, gx.nUex, gx.nUex_old, gx.i_flx_bl, gx.j_flx_bl, pm.n_bl-1, 1e-100, pm.flag_error);
    calcRes(pm.error_nUer, gr.nUer, gr.nUer_old, gr.i_flr_bl, gr.j_flr_bl, pm.n_bl-1, 1e-100, pm.flag_error);
    //------------------------------------

}

//*****************************************************************
//**                                                             **
//**           void solve_rhoe_wdTe_wSEE_PC()                    **
//**                                                             **
//*****************************************************************
void FluidModule::solve_rhoe_wdTe_wSEE_PC(Params &pm,GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr, SeeVec &se){

    std::vector<std::vector<double> > aW(pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //係数
    std::vector<std::vector<double> > aE(pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //係数
    std::vector<std::vector<double> > aS(pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //係数
    std::vector<std::vector<double> > aN(pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //係数
    std::vector<std::vector<double> > aP(pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //係数
    std::vector<std::vector<double> > b (pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //ソース
    std::vector<std::vector<double> > rhoe_new(pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //ソース
    std::vector<std::vector<double> > rhoe_old_tmp(pm.ni+2,std::vector<double>(pm.nj+2,0.0));
    std::vector<std::vector<double> > rhoUex_old_tmp(pm.ni+2,std::vector<double>(pm.nj+2,0.0));
    std::vector<std::vector<double> > rhoUer_old_tmp(pm.ni+2,std::vector<double>(pm.nj+2,0.0));
    std::vector<std::vector<double> > drhoe(pm.ni+2,std::vector<double>(pm.nj+2,0.0));

    //Numerical treatment for exponential overflow and zero-division
    //------------------------------------
    double lim_exp = 1.0e2; //expの制限
    double EPS = 1.0e-3; //ゼロ割回避用の数字
    //------------------------------------

    //Acceleration factor for PC method
    //------------------------------------
    double c_mu = 1.0;
    //------------------------------------


    //convergence criteria
    //------------------------------------
    int flag_error2 = pm.flag_error; //擬似ステップのエラー収束判定 0:max, 1:RMS
    //------------------------------------

    //output duration
    //------------------------------------
    int ndiv_out = 100;
    //------------------------------------

    //scaling factor of wall flux
    //------------------------------------
    double alpha_wflux = 1.0;
    //------------------------------------
    
    //Reserve the past values
    //------------------------------------
    for (int i=0;i<=pm.ni+1;i++){ 
        for (int j=0;j<=pm.nj+1;j++){ 
            gc.rhoe_old[i][j] = gc.rhoe[i][j];
            gx.rhoUex_old[i][j] = gx.rhoUex[i][j];
            gr.rhoUer_old[i][j] = gr.rhoUer[i][j];
        }
    }
    //------------------------------------

    //Initial prediction of correction
    //------------------------------------
    for (int i=0;i<=pm.ni+1;i++){ 
        for (int j=0;j<=pm.nj+1;j++){ 
            drhoe[i][j] = gc.rhoe[i][j] - gc.rhoe_old[i][j];
        }
    }
    //------------------------------------

    //PC method iteration
    //------------------------------------
    int ncount = 0; //number of iteration of PC method
    int flag_end = 0;
    do{
        
        ncount ++;
        pm.itime_PC_rhoe++;

        //Reserve previous-step values
        //------------------------------------
        for (int i=0;i<=pm.ni+1;i++){ 
            for (int j=0;j<=pm.nj+1;j++){ 
                rhoe_old_tmp[i][j] = gc.rhoe[i][j];
                rhoUex_old_tmp[i][j] = gx.rhoUex[i][j];
                rhoUer_old_tmp[i][j] = gr.rhoUer[i][j];
            }
        }
        //------------------------------------

        //Predict flux
        //==============================================
        //for rhoUex
        //------------------------------------
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
            for (int i=gx.i_flx_bl[iblock][0];i<=gx.i_flx_bl[iblock][1];i++){
                for (int j=gx.j_flx_bl[iblock][0];j<=gx.j_flx_bl[iblock][1];j++){
                    
                    //radial position at upper and lower cell interface
                    //------------------------------------
                    double rL = (gc.r[j]+gc.r[j-1])/2.0;
                    double rR = (gc.r[j]+gc.r[j+1])/2.0;
                    //------------------------------------

                    //moblity at cell interface
                    //------------------------------------
                    double mu_para_tmp = (gc.mu_para[i][j] + gc.mu_para[i-1][j])/2.0;
                    double mu_perp_tmp = (gc.mu_perp[i][j] + gc.mu_perp[i-1][j])/2.0;
                    
                    double mu_star_tmp = mu_perp_tmp*c_mu;
                    //------------------------------------

                    //calculate inverse moblity tensor at cell interface
                    //------------------------------------
                    double mu_rr_tmp = (gc.mu_rr[i][j] + gc.mu_rr[i-1][j])/2.0;
                    double mu_xr_tmp = (gc.mu_xr[i][j] + gc.mu_xr[i-1][j])/2.0;

                    double mu_inv_xx_tmp = mu_rr_tmp/(mu_para_tmp*mu_perp_tmp);
                    double mu_inv_xr_tmp = -mu_xr_tmp/(mu_para_tmp*mu_perp_tmp);
                    //------------------------------------

                    //values at cell interface for Sharfetter-Gummmel scheme
                    //------------------------------------
                    double Te_tmp = (gc.Te[i][j] + gc.Te[i-1][j])/2.0; //electron temperature
                    double dTedx_tmp = (gc.Te[i][j] - gc.Te[i-1][j])/pm.dx; //electron temperature gradient
                    double vd_tmp = mu_star_tmp*(-gx.Ex[i][j] - ph::Boltz/ph::e0*dTedx_tmp) + EPS; //electron drift velocity
                    double D_tmp = ph::Boltz*Te_tmp/ph::e0*mu_star_tmp; //electron diffusion coefficient
                    double z_tmp = vd_tmp*pm.dx/D_tmp; //Peclet number
                    //------------------------------------
                    
                    //BC - set adjacent flag (adjacent = 0, others = 1)
                    //adjacent "cell" to the wall for setting zero-flux on the wall
                    //------------------------------------
                    //lower
                    double bLr_wall  = double(j!=gc.j_flc_bl[0][0] || i>gc.i_flc_bl[0][1]); //x0
                    bLr_wall = bLr_wall*double(j!=gc.j_flc_bl[2][0] || i>gc.i_flc_bl[2][1]); //x2
                    double bLr_cen = double(j!=gc.j_flc_bl[3][0]); //x6
                    double bLr  = bLr_wall*bLr_cen;

                    //lower-left
                    double bLrD_wall  = double(j!=gc.j_flc_bl[0][0] || i-1>gc.i_flc_bl[0][1]); //x0
                    bLrD_wall = bLrD_wall*double(j!=gc.j_flc_bl[2][0] || i-1>gc.i_flc_bl[2][1]); //x2
                    double bLrD_cen = double(j!=gc.j_flc_bl[3][0]); //x6
                    double bLrD  = bLrD_wall*bLrD_cen;

                    //upper
                    double bRr_wall  = double(j!=gc.j_flc_bl[0][1] || iblock > 2); //x1
                    bRr_wall = bRr_wall*double(j!=gc.j_flc_bl[3][1] || i<gc.i_flc_bl[1][1]+1 || i>gc.i_flc_bl[3][1]); //x4
                    //double bRr_open  = (j!=gc.j_flc_bl[4][1] || iblock != 4); //x5
                    double bRr = bRr_wall;

                    //upper-left
                    double bRrD_wall  = double(j!=gc.j_flc_bl[0][1] || iblock > 2); //x1
                    bRrD_wall = bRrD_wall*double(j!=gc.j_flc_bl[3][1] || i-1<gc.i_flc_bl[1][1]+1 || i-1>gc.i_flc_bl[3][1]); //x4
                    //double bRrD_open  = (j!=gc.j_flc_bl[4][1] || iblock != 4); //x5
                    double bRrD = bRrD_wall;
                    //------------------------------------

                    //interpolation of transverse flux component (should be old value)
                    //------------------------------------
                    double rhoUer_tmp = (rR*rhoUer_old_tmp[i][j+1]*bRr + rL*rhoUer_old_tmp[i][j]*bLr
                        + rR*rhoUer_old_tmp[i-1][j+1]*bRrD + rL*rhoUer_old_tmp[i-1][j]*bLrD)/(4.0*gc.r[j]);
                    //------------------------------------

                    //update flux
                    //------------------------------------
                    gx.rhoUex[i][j] =
                          vd_tmp/(1.0 - exp(fmin(-z_tmp,lim_exp)))*gc.rhoe[i-1][j]
                        - vd_tmp/(exp(fmin( z_tmp,lim_exp)) - 1.0)*gc.rhoe[i][j] 
                        + rhoUex_old_tmp[i][j] 
                        - mu_star_tmp*(mu_inv_xx_tmp*rhoUex_old_tmp[i][j] + mu_inv_xr_tmp*rhoUer_tmp);
                    //------------------------------------
                }
            }     
        }
        //------------------------------------

        //left-wall-BC (z0)
        //------------------------------------
        for (int j=gc.j_flc_bl[0][0];j<=gc.j_flc_bl[0][1];j++){
            int i=gc.i_flc_bl[0][0];
            double vth = sqrt(8.0*ph::Boltz*gc.Te[i][j]/(M_PI*pm.masse));
            double delta_phi = fmax(-gx.Ex[i][j]*pm.dx/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerN_wall_Lx = -0.25*alpha_wflux*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Lx = -pm.coefIISEE_ts*gc.rhoi[i][j]*gx.Uix[i][j] //IISEE
                            -(coefEISEE_eb + coefEISEE_rd + coefEISEE_ts)*gc.rhoe_old[i][j]*GammaPerN_wall_Lx //EISEE
                            -pm.coefMISEE_ts*gx.rhoUmx[i][j]; //MISEE
            gx.rhoUex[i][j] = gc.rhoe[i][j]*GammaPerN_wall_Lx + Gamma_wallSEE_Lx;
            gx.rhoUex_wall[i][j] = gx.rhoUex[i][j];
        }
        //------------------------------------
    
        //left-wall-BC (z1)
        //------------------------------------
        for (int j=gc.j_flc_bl[2][0];j<=gc.j_flc_bl[0][0]-1;j++){
            int i=gc.i_flc_bl[2][0];
            double vth = sqrt(8.0*ph::Boltz*gc.Te[i][j]/(M_PI*pm.masse));
            double delta_phi = fmax(-gx.Ex[i][j]*pm.dx/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerN_wall_Lx = -0.25*alpha_wflux*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Lx = -pm.coefIISEE_ts*gc.rhoi[i][j]*gx.Uix[i][j] //IISEE
                            -(coefEISEE_eb + coefEISEE_rd + coefEISEE_ts)*gc.rhoe_old[i][j]*GammaPerN_wall_Lx //EISEE
                            -pm.coefMISEE_ts*gx.rhoUmx[i][j]; //MISEE
            gx.rhoUex[i][j] = gc.rhoe[i][j]*GammaPerN_wall_Lx + Gamma_wallSEE_Lx;
            gx.rhoUex_wall[i][j] = gx.rhoUex[i][j];
        }
        //------------------------------------
    
        //left-wall-BC (z2)
        //------------------------------------
        for (int j=gc.j_flc_bl[3][0];j<=gc.j_flc_bl[2][0]-1;j++){
            int i=gc.i_flc_bl[3][0];
            double vth = sqrt(8.0*ph::Boltz*gc.Te[i][j]/(M_PI*pm.masse));
            double delta_phi = fmax(-gx.Ex[i][j]*pm.dx/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerN_wall_Lx = -0.25*alpha_wflux*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Lx = -pm.coefIISEE_ts*gc.rhoi[i][j]*gx.Uix[i][j] //IISEE
                            -(coefEISEE_eb + coefEISEE_rd + coefEISEE_ts)*gc.rhoe_old[i][j]*GammaPerN_wall_Lx //EISEE
                            -pm.coefMISEE_ts*gx.rhoUmx[i][j]; //MISEE
            gx.rhoUex[i][j] = gc.rhoe[i][j]*GammaPerN_wall_Lx + Gamma_wallSEE_Lx;
            gx.rhoUex_wall[i][j] = gx.rhoUex[i][j];
        }
        //------------------------------------
    
        //left-wall-BC (z4)
        //------------------------------------
        for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[4][1];j++){
            int i=gc.i_flc_bl[4][0];
            double vth = sqrt(8.0*ph::Boltz*gc.Te[i][j]/(M_PI*pm.masse));
            double delta_phi = fmax(-gx.Ex[i][j]*pm.dx/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerN_wall_Lx = -0.25*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Lx = -pm.coefIISEE_ts*gc.rhoi[i][j]*gx.Uix[i][j] //IISEE
                            -(coefEISEE_eb + coefEISEE_rd + coefEISEE_ts)*gc.rhoe_old[i][j]*GammaPerN_wall_Lx //EISEE
                            -pm.coefMISEE_ts*gx.rhoUmx[i][j]; //MISEE
            gx.rhoUex[i][j] = gc.rhoe[i][j]*GammaPerN_wall_Lx + Gamma_wallSEE_Lx;
            gx.rhoUex_wall[i][j] = gx.rhoUex[i][j];
        }
        //------------------------------------
    
        //right-wall-BC (z3)
        //------------------------------------
        for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[1][1];j++){
            int i=gc.i_flc_bl[1][1];
            double vth = sqrt(8.0*ph::Boltz*gc.Te[i][j]/(M_PI*pm.masse));
            double delta_phi =  fmax(gx.Ex[i+1][j]*pm.dx/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerN_wall_Rx = 0.25*alpha_wflux*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Rx = -pm.coefIISEE_ts*gc.rhoi[i][j]*gx.Uix[i+1][j] //IISEE
                            -(coefEISEE_eb + coefEISEE_rd + coefEISEE_ts)*gc.rhoe_old[i][j]*GammaPerN_wall_Rx //EISEE
                            -pm.coefMISEE_ts*gx.rhoUmx[i+1][j]; //MISEE
            gx.rhoUex[i+1][j] =  gc.rhoe[i][j]*GammaPerN_wall_Rx + Gamma_wallSEE_Rx;
            gx.rhoUex_wall[i+1][j] = gx.rhoUex[i+1][j];
        }
        //------------------------------------
    
        //right-wall-BC (z5)
        //------------------------------------
        for (int j=gc.j_flc_bl[4][0];j<=gc.j_flc_bl[4][1];j++){
            int i=gc.i_flc_bl[4][1];
            double vth = sqrt(8.0*ph::Boltz*gc.Te[i][j]/(M_PI*pm.masse));
            double delta_phi =  fmax(gx.Ex[i+1][j]*pm.dx/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerN_wall_Rx = 0.25*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Rx = -pm.coefIISEE_ts*gc.rhoi[i][j]*gx.Uix[i+1][j] //IISEE
                            -(coefEISEE_eb + coefEISEE_rd + coefEISEE_ts)*gc.rhoe_old[i][j]*GammaPerN_wall_Rx //EISEE
                            -pm.coefMISEE_ts*gx.rhoUmx[i+1][j]; //MISEE
            gx.rhoUex[i+1][j] =  gc.rhoe[i][j]*GammaPerN_wall_Rx + Gamma_wallSEE_Rx;
            gx.rhoUex_wall[i+1][j] = gx.rhoUex[i+1][j];
        }
        //=====================================================================================

        //for rhoUer
        //------------------------------------
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){
            for (int i=gr.i_flr_bl[iblock][0];i<=gr.i_flr_bl[iblock][1];i++){ 
                for (int j=gr.j_flr_bl[iblock][0];j<=gr.j_flr_bl[iblock][1];j++){

                    //radial position at cell interface
                    //------------------------------------
                    double r_tmp = (gc.r[j] + gc.r[j-1])/2.0;
                    //------------------------------------

                    //moblity at cell interface
                    //------------------------------------
                    double mu_para_tmp = (gc.mu_para[i][j] + gc.mu_para[i][j-1])/2.0;
                    double mu_perp_tmp = (gc.mu_perp[i][j] + gc.mu_perp[i][j-1])/2.0;

                    double mu_star_tmp = mu_perp_tmp*c_mu;
                    //------------------------------------

                    //calculate inverse moblity tensor at cell interface
                    //------------------------------------
                    double mu_xx_tmp = (gc.mu_xx[i][j] + gc.mu_xx[i][j-1])/2.0;
                    double mu_xr_tmp = (gc.mu_xr[i][j] + gc.mu_xr[i][j-1])/2.0;
                    
                    double mu_inv_rr_tmp = mu_xx_tmp/(mu_para_tmp*mu_perp_tmp);
                    double mu_inv_xr_tmp = -mu_xr_tmp/(mu_para_tmp*mu_perp_tmp);
                    //------------------------------------
    
                    //values at cell interface for Sharfetter-Gummmel scheme
                    //------------------------------------
                    double Te_tmp = (gc.Te[i][j] + gc.Te[i][j-1])/2.0; //electron temperature
                    double dTedr_tmp = (gc.Te[i][j] - gc.Te[i][j-1])/pm.dr; //electron temperature gradient
                    double vd_tmp = mu_star_tmp*(-gr.Er[i][j] - ph::Boltz/ph::e0*dTedr_tmp) + EPS; //electron drift velocity
                    double D_tmp = ph::Boltz*Te_tmp/ph::e0*mu_star_tmp; //electron diffusion coefficient
                    double z_tmp = vd_tmp*pm.dr/D_tmp; //Peclet number
                    //------------------------------------

                    //BC - set adjacent flag (adjacent = 0, others = 1)
                    //adjacent "cell" to the wall for setting zero-flux on the wall
                    //------------------------------------
                    //left
                    double bLx_wall  = double(i!=gc.i_flc_bl[0][0] || iblock != 0); //z0
                    bLx_wall = bLx_wall*double(i!=gc.i_flc_bl[2][0] || j>gc.j_flc_bl[0][0]-1); //z1
                    bLx_wall = bLx_wall*double(i!=gc.i_flc_bl[3][0] || j>gc.j_flc_bl[2][0]-1); //z2
                    bLx_wall = bLx_wall*double(i!=gc.i_flc_bl[4][0] || j<gc.j_flc_bl[1][0]); //z4
                    double bLx  = bLx_wall;

                    //left-lower
                    double bLxD_wall  = double(i!=gc.i_flc_bl[0][0] || iblock != 0); //z0
                    bLxD_wall = bLxD_wall*double(i!=gc.i_flc_bl[2][0] || j-1>gc.j_flc_bl[0][0]-1); //z1
                    bLxD_wall = bLxD_wall*double(i!=gc.i_flc_bl[3][0] || j-1>gc.j_flc_bl[2][0]-1); //z2
                    bLxD_wall = bLxD_wall*double(i!=gc.i_flc_bl[4][0] || j-1<gc.j_flc_bl[1][0]); //z4
                    double bLxD  = bLxD_wall;

                    //right
                    double bRx_wall  = double(i!=gc.i_flc_bl[1][1] || j<gc.j_flc_bl[1][0]); //z3
                    bRx_wall = bRx_wall*double(i!=gc.i_flc_bl[4][1]); //z5
                    double bRx  = bRx_wall;

                    //right-lower
                    double bRxD_wall  = double(i!=gc.i_flc_bl[1][1] || j-1<gc.j_flc_bl[1][0]); //z3
                    bRxD_wall = bRxD_wall*double(i!=gc.i_flc_bl[4][1]); //z5
                    double bRxD  = bRxD_wall;
                    //------------------------------------

                    //interpolation of transverse flux component (should be old value)
                    //------------------------------------
                    double rhoUex_tmp = (gc.r[j]*rhoUex_old_tmp[i+1][j]*bRx + gc.r[j]*rhoUex_old_tmp[i][j]*bLx 
                        + gc.r[j-1]*rhoUex_old_tmp[i+1][j-1]*bRxD + gc.r[j-1]*rhoUex_old_tmp[i][j-1]*bLxD)/(4.0*r_tmp);
                    //------------------------------------

                    //update flux
                    //------------------------------------
                    gr.rhoUer[i][j] =  
                          vd_tmp/(1.0 - exp(fmin(-z_tmp,lim_exp)))*gc.rhoe[i][j-1]
                        - vd_tmp/(exp(fmin( z_tmp,lim_exp)) - 1.0)*gc.rhoe[i][j]
                        + rhoUer_old_tmp[i][j]
                        - mu_star_tmp*(mu_inv_xr_tmp*rhoUex_tmp + mu_inv_rr_tmp*rhoUer_old_tmp[i][j]);
                    //------------------------------------
                }
            }
        }

        //lower-wall-BC (x0)
        //------------------------------------
        for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[0][1];i++){
            int j=gc.j_flc_bl[0][0];
            double vth = sqrt(8.0*ph::Boltz*gc.Te[i][j]/(M_PI*pm.masse));
            double delta_phi = fmax(-gr.Er[i][j]*pm.dr/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerN_wall_Lr = -0.25*alpha_wflux*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Lr = -pm.coefIISEE_ts*gc.rhoi[i][j]*gr.Uir[i][j] //IISEE
                            -(coefEISEE_eb + coefEISEE_rd + coefEISEE_ts)*gc.rhoe_old[i][j]*GammaPerN_wall_Lr //EISEE
                            -pm.coefMISEE_ts*gr.rhoUmr[i][j]; //MISEE
            gr.rhoUer[i][j] = gc.rhoe[i][j]*GammaPerN_wall_Lr + Gamma_wallSEE_Lr;
            gr.rhoUer_wall[i][j] = gr.rhoUer[i][j];
        }
        //------------------------------------

        //lower-wall-BC (x2)
        //------------------------------------
        for (int i=gc.i_flc_bl[2][0];i<=gc.i_flc_bl[2][1];i++){
            int j=gc.j_flc_bl[2][0];
            double vth = sqrt(8.0*ph::Boltz*gc.Te[i][j]/(M_PI*pm.masse));
            double delta_phi = fmax(-gr.Er[i][j]*pm.dr/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerN_wall_Lr = -0.25*alpha_wflux*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te[i][j] + 1e-100),lim_exp));
            
            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Lr = -pm.coefIISEE_ts*gc.rhoi[i][j]*gr.Uir[i][j] //IISEE
                            -(coefEISEE_eb + coefEISEE_rd + coefEISEE_ts)*gc.rhoe_old[i][j]*GammaPerN_wall_Lr //EISEE
                            -pm.coefMISEE_ts*gr.rhoUmr[i][j]; //MISEE
            gr.rhoUer[i][j] = gc.rhoe[i][j]*GammaPerN_wall_Lr + Gamma_wallSEE_Lr;
            gr.rhoUer_wall[i][j] = gr.rhoUer[i][j];
        }
        //------------------------------------

        //upper-wall-BC (x1)
        //------------------------------------
        for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[1][1];i++){
            int j=gc.j_flc_bl[0][1];
            double vth = sqrt(8.0*ph::Boltz*gc.Te[i][j]/(M_PI*pm.masse));
            double delta_phi =  fmax(gr.Er[i][j+1]*pm.dr/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerN_wall_Rr = 0.25*alpha_wflux*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Rr = -pm.coefIISEE_ts*gc.rhoi[i][j]*gr.Uir[i][j+1] //IISEE
                            -(coefEISEE_eb + coefEISEE_rd + coefEISEE_ts)*gc.rhoe_old[i][j]*GammaPerN_wall_Rr //EISEE
                            -pm.coefMISEE_ts*gr.rhoUmr[i][j+1]; //MISEE
            gr.rhoUer[i][j+1] = gc.rhoe[i][j]*GammaPerN_wall_Rr + Gamma_wallSEE_Rr;
            gr.rhoUer_wall[i][j+1] = gr.rhoUer[i][j+1];
        }
        //------------------------------------

        //upper-wall-BC (x4)
        //------------------------------------
        for (int i=gc.i_flc_bl[1][1]+1;i<=gc.i_flc_bl[3][1];i++){
            int j=gc.j_flc_bl[3][1];
            double vth = sqrt(8.0*ph::Boltz*gc.Te[i][j]/(M_PI*pm.masse));
            double delta_phi =  fmax(gr.Er[i][j+1]*pm.dr/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerN_wall_Rr = 0.25*alpha_wflux*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Rr = -pm.coefIISEE_ts*gc.rhoi[i][j]*gr.Uir[i][j+1] //IISEE
                            -(coefEISEE_eb + coefEISEE_rd + coefEISEE_ts)*gc.rhoe_old[i][j]*GammaPerN_wall_Rr //EISEE
                            -pm.coefMISEE_ts*gr.rhoUmr[i][j+1]; //MISEE
            gr.rhoUer[i][j+1] = gc.rhoe[i][j]*GammaPerN_wall_Rr + Gamma_wallSEE_Rr;
            gr.rhoUer_wall[i][j+1] = gr.rhoUer[i][j+1];
        }
        //------------------------------------

        //upper-open-BC (x5)
        //------------------------------------
        for (int i=gc.i_flc_bl[4][0];i<=gc.i_flc_bl[4][1];i++){
            int j=gc.j_flc_bl[4][1];
            double rhoUir_Rr = gc.rhoi[i][j]*gr.Uir[i][j+1];
            double Gamma_open_Rr = rhoUir_Rr;
            gr.rhoUer[i][j+1] =  Gamma_open_Rr;
            gr.rhoUer_wall[i][j+1] = 0.0;
        }
        //------------------------------------

        //Make coefficient for electron continuity equation for 5 point stencil
        //------------------------------------
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
            for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
                for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                    //radial position at upper and lower cell interface
                    //------------------------------------
                    double rR = (gc.r[j]+gc.r[j+1])/2.0;
                    double rL = (gc.r[j]+gc.r[j-1])/2.0;
                    //------------------------------------

                    //normalized radial position at upper and lower cell interface
                    //------------------------------------
                    double qR = (gc.r[j]+gc.r[j+1])/2.0/gc.r[j];
                    double qL = (gc.r[j]+gc.r[j-1])/2.0/gc.r[j];
                    //------------------------------------
                    
                    //moblity at cell interface
                    //------------------------------------
                    double mu_star_Lx = (gc.mu_perp[i][j]+gc.mu_perp[i-1][j])/2.0*c_mu;
                    double mu_star_Rx = (gc.mu_perp[i][j]+gc.mu_perp[i+1][j])/2.0*c_mu;
                    double mu_star_Lr = (gc.mu_perp[i][j]+gc.mu_perp[i][j-1])/2.0*c_mu;
                    double mu_star_Rr = (gc.mu_perp[i][j]+gc.mu_perp[i][j+1])/2.0*c_mu;
                    //------------------------------------

                    //electron temperature at cell interface
                    //------------------------------------
                    double Te_Lx = (gc.Te[i][j]+gc.Te[i-1][j])/2.0;
                    double Te_Rx = (gc.Te[i][j]+gc.Te[i+1][j])/2.0;
                    double Te_Lr = (gc.Te[i][j]+gc.Te[i][j-1])/2.0;
                    double Te_Rr = (gc.Te[i][j]+gc.Te[i][j+1])/2.0;
                    //------------------------------------

                    //electron diffusion coefficient at cell interface
                    //------------------------------------
                    double D_Lx = ph::Boltz*Te_Lx/ph::e0*mu_star_Lx;
                    double D_Rx = ph::Boltz*Te_Rx/ph::e0*mu_star_Rx;
                    double D_Lr = ph::Boltz*Te_Lr/ph::e0*mu_star_Lr;
                    double D_Rr = ph::Boltz*Te_Rr/ph::e0*mu_star_Rr;
                    //------------------------------------

                    //electron temperature gradient at cell interface
                    //------------------------------------
                    double dTedx_Lx = (gc.Te[i  ][j  ] - gc.Te[i-1][j  ])/pm.dx;
                    double dTedx_Rx = (gc.Te[i+1][j  ] - gc.Te[i  ][j  ])/pm.dx;
                    double dTedr_Lr = (gc.Te[i  ][j  ] - gc.Te[i  ][j-1])/pm.dr;
                    double dTedr_Rr = (gc.Te[i  ][j+1] - gc.Te[i  ][j  ])/pm.dr;
                    //------------------------------------

                    //electron drift velocity at cell interface
                    //------------------------------------
                    double vd_Lx = mu_star_Lx*(-gx.Ex[i  ][j  ] - ph::Boltz/ph::e0*dTedx_Lx) + EPS;
                    double vd_Rx = mu_star_Rx*(-gx.Ex[i+1][j  ] - ph::Boltz/ph::e0*dTedx_Rx) + EPS;
                    double vd_Lr = mu_star_Lr*(-gr.Er[i  ][j  ] - ph::Boltz/ph::e0*dTedr_Lr) + EPS;
                    double vd_Rr = mu_star_Rr*(-gr.Er[i  ][j+1] - ph::Boltz/ph::e0*dTedr_Rr) + EPS;
                    //------------------------------------

                    //Peclet number at cell interface
                    //------------------------------------
                    double z_Lx = vd_Lx*pm.dx/D_Lx;
                    double z_Rx = vd_Rx*pm.dx/D_Rx;
                    double z_Lr = vd_Lr*pm.dr/D_Lr;
                    double z_Rr = vd_Rr*pm.dr/D_Rr;
                    //------------------------------------

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

                    //coefficient for electron continuity equation for 5 point stencil
                    //==============================================
                    //vertical and horizental 5 stencil
                    //------------------------------------
                    aP[i][j] =    vd_Rx/(1.0 - exp(fmin(-z_Rx,lim_exp)))*bRx*pm.dr
                        +    vd_Lx/(exp(fmin( z_Lx,lim_exp)) - 1.0)*bLx*pm.dr
                        + qR*vd_Rr/(1.0 - exp(fmin(-z_Rr,lim_exp)))*bRr*pm.dx
                        + qL*vd_Lr/(exp(fmin( z_Lr,lim_exp)) - 1.0)*bLr*pm.dx
                        + pm.dx*pm.dr/pm.dt; //time-term
                    aE[i][j] =    vd_Rx/(exp(fmin( z_Rx,lim_exp)) - 1.0)*bRx*pm.dr;
                    aW[i][j] =    vd_Lx/(1.0 - exp(fmin(-z_Lx,lim_exp)))*bLx*pm.dr;
                    aN[i][j] = qR*vd_Rr/(exp(fmin( z_Rr,lim_exp)) - 1.0)*bRr*pm.dx;
                    aS[i][j] = qL*vd_Lr/(1.0 - exp(fmin(-z_Lr,lim_exp)))*bLr*pm.dx;
                    //------------------------------------

                    //RHS
                    //------------------------------------
                    double nabla_rhoUe_tmp = (gx.rhoUex[i+1][j]-gx.rhoUex[i][j])/pm.dx + (qR*gr.rhoUer[i][j+1]-qL*gr.rhoUer[i][j])/pm.dr;
                    b[i][j] = (-nabla_rhoUe_tmp 
                        + gc.rate_ionize[i][j]
                        - (gc.rhoe[i][j] - gc.rhoe_old[i][j])/pm.dt)*pm.dx*pm.dr; //time-term
                    //------------------------------------
                    //==============================================

                }
            }
        }
        //------------------------------------

        //Solver
        //------------------------------------
        if(pm.flag_iter_rhoe == 0){
            //SOR
            //------------------------------------
            double alpha = 1.9;
            solver_SOR(aP,aE,aW,aN,aS,b,gc.i_flc_bl,gc.j_flc_bl,pm.n_bl-1,alpha,pm.maxITR_SOR_rhoe,pm.error_cnv_SOR_rhoe,0,drhoe);
            //------------------------------------
        }else{
            //Semi-coarse Multi-Grid by Hypre
            //------------------------------------
            solver_SMG(aP,aE,aW,aN,aS,b,gc.i_flc_bl,gc.j_flc_bl,pm.n_bl-1,pm.maxITR_SOR_rhoe,pm.error_cnv_SOR_rhoe,0,drhoe);
            //------------------------------------
        }
        //------------------------------------

        //Correction
        //==============================================
        //for rhoUex
        //------------------------------------
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
            for (int i=gx.i_flx_bl[iblock][0];i<=gx.i_flx_bl[iblock][1];i++){
                for (int j=gx.j_flx_bl[iblock][0];j<=gx.j_flx_bl[iblock][1];j++){
                    
                    //values at cell interface for Sharfetter-Gummmel scheme
                    //------------------------------------
                    double mu_star_tmp = (gc.mu_perp[i][j] + gc.mu_perp[i-1][j])/2.0*c_mu;
                    double dTedx_tmp = (gc.Te[i][j] - gc.Te[i-1][j])/pm.dx;
                    double vd_tmp = mu_star_tmp*(-gx.Ex[i][j] - ph::Boltz/ph::e0*dTedx_tmp) + EPS;
                    double Te_tmp = (gc.Te[i][j] + gc.Te[i-1][j])/2.0;
                    double D_tmp = ph::Boltz*Te_tmp/ph::e0*mu_star_tmp;
                    double z_tmp = vd_tmp*pm.dx/D_tmp;
                    //------------------------------------

                    //calculate flux correction
                    //------------------------------------
                    double drhoUex = vd_tmp/(1.0 - exp(fmin(-z_tmp,lim_exp)))*drhoe[i-1][j]
                     - vd_tmp/(exp(fmin( z_tmp,lim_exp)) - 1.0)*drhoe[i][j];
                    //------------------------------------
                    
                    //Update
                    //------------------------------------
                    gx.rhoUex[i][j] = gx.rhoUex[i][j] + drhoUex;
                    //------------------------------------
                }
            }
        }
        //------------------------------------

        //for rhoUer
        //------------------------------------
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){
            for (int i=gr.i_flr_bl[iblock][0];i<=gr.i_flr_bl[iblock][1];i++){ 
                for (int j=gr.j_flr_bl[iblock][0];j<=gr.j_flr_bl[iblock][1];j++){
                    
                    //values at cell interface for Sharfetter-Gummmel scheme
                    //------------------------------------
                    double mu_star_tmp = (gc.mu_perp[i][j] + gc.mu_perp[i][j-1])/2.0*c_mu;
                    double dTedr_tmp = (gc.Te[i  ][j  ] - gc.Te[i  ][j-1])/pm.dr;
                    double vd_tmp = mu_star_tmp*(-gr.Er[i][j] - ph::Boltz/ph::e0*dTedr_tmp) + EPS;
                    double Te_tmp = (gc.Te[i][j] + gc.Te[i][j-1])/2.0;
                    double D_tmp = ph::Boltz*Te_tmp/ph::e0*mu_star_tmp;
                    double z_tmp = vd_tmp*pm.dr/D_tmp;
                    //------------------------------------

                    //calculate flux correction
                    //------------------------------------
                    double drhoUer =  vd_tmp/(1.0 - exp(fmin(-z_tmp,lim_exp)))*drhoe[i][j-1]
                     - vd_tmp/(exp(fmin( z_tmp,lim_exp)) - 1.0)*drhoe[i][j];
                    //------------------------------------
                    
                    //Update
                    //------------------------------------
                    gr.rhoUer[i][j] = gr.rhoUer[i][j] + drhoUer;
                    //------------------------------------
                }
            }
        }
        //------------------------------------

        //for rhoe
        //------------------------------------
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
            for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
                for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                    double rhoe_new = gc.rhoe[i][j] + drhoe[i][j];
                    if(rhoe_new > 0.0){
                        gc.rhoe[i][j] = rhoe_new;
                    }
                }
            }
        }
        //------------------------------------
        //==============================================

        //Caluculate residuals
        //------------------------------------
        double error_rhoe2 = 0.0;
        double error_rhoUex2 = 0.0;
        double error_rhoUer2 = 0.0;
        calcRes(error_rhoe2, gc.rhoe,  rhoe_old_tmp,  gc.i_flc_bl, gc.j_flc_bl, pm.n_bl-1,1e-100,flag_error2);
        calcRes(error_rhoUex2,gx.rhoUex, rhoUex_old_tmp, gx.i_flx_bl, gx.j_flx_bl, pm.n_bl-1,1e-100,flag_error2);
        calcRes(error_rhoUer2,gr.rhoUer, rhoUer_old_tmp, gr.i_flr_bl, gr.j_flr_bl, pm.n_bl-1,1e-100,flag_error2);
        
        if(pm.flag_chk == 1) output_residual_for_PC("rhoe",pm.gtime, pm.itime, error_rhoe2, error_rhoUex2, error_rhoUer2, pm.itime_PC_rhoe, ncount);
        
        double error_global = fmax(fmax(error_rhoe2,error_rhoUex2),error_rhoUer2);
        
        if(ncount % ndiv_out == 0){
            std::cout << "ncount = "<<ncount<<std::endl;
            std::cout << "error_rhoe = " << error_rhoe2
                <<  " error_rhoUex = " << error_rhoUex2
                <<  " error_rhoUer= " << error_rhoUer2
                <<  " error_max = " << fmax(fmax(error_rhoe2,error_rhoUex2),error_rhoUer2)
                << std::endl;
        }

        if(error_global < pm.error_cnv_HES_rhoe || ncount >= pm.maxITR_HES_rhoe){
            flag_end = 1;
        }

    } while (flag_end == 0);

    //calculation of particle balance
    //------------------------------------
    if(pm.itime%pm.ndiv_fout == 0 || pm.itime == pm.ntime){
        
        //volume integrated values
        //------------------------------------
        double dNedt_sum         = 0.0;
        double nabla_rhoUe_sum   = 0.0;
        double rate_d_ionz_sum   = 0.0;
        double rate_s_ionz_sum   = 0.0;
        //------------------------------------
    
        //calculate volume integrated of each term
        //------------------------------------
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
            for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
                for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                    
                    //radial position at upper and lower cell interface
                    //------------------------------------
                    double rR = (gc.r[j]+gc.r[j+1])/2.0;
                    double rL = (gc.r[j]+gc.r[j-1])/2.0;
                    //------------------------------------

                    //normalized radial position at upper and lower cell interface
                    //------------------------------------
                    double qR = (gc.r[j]+gc.r[j+1])/2.0/gc.r[j];
                    double qL = (gc.r[j]+gc.r[j-1])/2.0/gc.r[j];
                    //------------------------------------

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

                    double nabla_rhoUe = (gx.rhoUex[i+1][j]-gx.rhoUex[i][j])/pm.dx + (qR*gr.rhoUer[i][j+1]-qL*gr.rhoUer[i][j])/pm.dr;
                    
                    double volume = 2.0*M_PI*gc.r[j]*pm.dr*pm.dx;
                    
                    dNedt_sum         = dNedt_sum         + volume*(gc.rhoe[i][j] - gc.rhoe_old[i][j])/pm.dt;
                    nabla_rhoUe_sum   = nabla_rhoUe_sum   + volume*nabla_rhoUe;
                    rate_d_ionz_sum   = rate_d_ionz_sum   + volume*gc.rhoe_old[i][j]*gc.nu_ionz[i][j];
                    rate_s_ionz_sum   = rate_s_ionz_sum   + volume*gc.rhoe_old[i][j]*gc.nu_ionzStep[i][j];
                }
            }
        }

        //output particle balance for electron
        //------------------------------------
        std::string char1="results/particle_balance_electron";
        std::string char2=std::to_string(pm.nOut);
        std::string char_csv=".csv";
        std::ofstream outputfile8(char1+char2+char_csv);
        outputfile8 << std::setprecision(std::numeric_limits<double>::max_digits10) << std::scientific;
        outputfile8<< "dN/dt (pcl/s), loss_to_bnd (pcl/s), direct-iz (pcl/s), sw-iz (pcl/s), sum (pcl/s)"<< std::endl;
        outputfile8<< dNedt_sum << "," << -nabla_rhoUe_sum << "," << rate_d_ionz_sum << ","<< rate_s_ionz_sum 
            << ","<< -dNedt_sum + (-nabla_rhoUe_sum + rate_d_ionz_sum + rate_s_ionz_sum)<< std::endl;
        outputfile8.close();
        //------------------------------------
    }
    //------------------------------------

    //set values of ghost cells for ES-potential calculation
    //=====================================================================================
    
    //left-wall-BC (chamber)
    //------------------------------------
    for (int j=gc.j_flc_bl[0][0];j<=gc.j_flc_bl[0][1];j++){
        int i = gc.i_flc_bl[0][0];
        gc.rhoe[i-1][j] = 0.0;
    }
    //------------------------------------

    //left-wall-BC (magnet)
    //------------------------------------
    for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[0][0]-1;j++){
        int i = gc.i_flc_bl[1][0];
        gc.rhoe[i-1][j] = 0.0;
    }
    //------------------------------------

    //lower-wall-BC (antenna)
    //------------------------------------
    for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[0][1];i++){ 
        int j = gc.j_flc_bl[0][0];
        gc.rhoe[i][j-1] = 0.0;
    }

    //right-open-BC
    //------------------------------------
    for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[1][1];j++){
        int i = gc.i_flc_bl[1][1];
        gc.rhoe[i+1][j] = gc.rhoe[i][j];
    }
    //------------------------------------

    //upper-open-BC
    //------------------------------------
    for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[1][1];i++){ 
        int j = gc.j_flc_bl[0][1];
        gc.rhoe[i][j+1] = gc.rhoe[i][j];
    }
    //------------------------------------
    //=====================================================================================

    /*
    //Caluculate Current
    //------------------------------------
    Ii_Anode = 0.0;
    Ie_Anode = 0.0;
    I_Anode = 0.0;

    for (int j=gc.j_flc_bl[4][0];j<=gc.j_flc_bl[4][1];j++){
        int i=gc.i_flc_bl[4][0];
        double vth = sqrt(8.0*ph::Boltz*pm.Ti/(M_PI*pm.massi));
        double um = sqrt(2.0*ph::Boltz*pm.Ti/pm.massi);
        double utx_R = gx.Uix[i+1][j]/(um + 1e-100);

        double rhoUix_tmp = gc.rhoi[i][j]*0.25*vth*(exp(-utx_R*utx_R) + utx_R*sqrt(M_PI)*(erf(utx_R) + 1.0));
        
        double delta_Ii = rhoUix_tmp  *2.0*M_PI*gc.r[j]*pm.dr*( ph::e0);
        double delta_Ie = gx.rhoUex[i][j]*2.0*M_PI*gc.r[j]*pm.dr*(-ph::e0);

        Ii_Anode = Ii_Anode + delta_Ii;
        Ie_Anode = Ie_Anode + delta_Ie;
        I_Anode =  I_Anode + delta_Ii + delta_Ie;
    }
    //------------------------------------
    */

    //Caluculate Current
    //------------------------------------
    pm.Ii_Anode = 0.0;
    pm.Ie_Anode = 0.0;
    pm.I_Anode = 0.0;

    for (int j=gc.j_flc_bl[4][0];j<=gc.j_flc_bl[4][1];j++){
        int i=gc.i_flc_bl[4][1];
        double vth = sqrt(8.0*ph::Boltz*pm.Ti/(M_PI*pm.massi));
        double um = sqrt(2.0*ph::Boltz*pm.Ti/pm.massi);
        double utx_R = gx.Uix[i+1][j]/(um + 1e-100);

        double rhoUix_tmp = gc.rhoi[i][j]*0.25*vth*(exp(-utx_R*utx_R) + utx_R*sqrt(M_PI)*(erf(utx_R) + 1.0));
        
        double delta_Ii = rhoUix_tmp  *2.0*M_PI*gc.r[j]*pm.dr*( ph::e0);
        double delta_Ie = gx.rhoUex[i+1][j]*2.0*M_PI*gc.r[j]*pm.dr*(-ph::e0);

        pm.Ii_Anode = pm.Ii_Anode + delta_Ii;
        pm.Ie_Anode = pm.Ie_Anode + delta_Ie;
        pm.I_Anode =  pm.I_Anode + delta_Ii + delta_Ie;
    }
    //------------------------------------

    /*
    //Caluculate Current
    //------------------------------------
    Ii_Nozzle = 0.0;
    Ie_Nozzle = 0.0;
    I_Nozzle = 0.0;

    for (int j=gc.j_flc_bl[4][0];j<=gc.j_flc_bl[4][1];j++){
        int i=gc.i_flc_bl[4][0];
        double vth = sqrt(8.0*ph::Boltz*pm.Ti/(M_PI*pm.massi));
        double um = sqrt(2.0*ph::Boltz*pm.Ti/pm.massi);
        double utx_R = gx.Uix[i+1][j]/(um + 1e-100);
        

        double rhoUix_tmp = gc.rhoi[i][j]*0.25*vth*(exp(-utx_R*utx_R) + utx_R*sqrt(M_PI)*(erf(utx_R) + 1.0));
        
        double delta_Ii = rhoUix_tmp  *2.0*M_PI*gc.r[j]*pm.dr*( ph::e0);
        double delta_Ie = gx.rhoUex[i][j]*2.0*M_PI*gc.r[j]*pm.dr*(-ph::e0);

        Ii_Anode = Ii_Anode + delta_Ii;
        Ie_Anode = Ie_Anode + delta_Ie;
        I_Anode =  I_Anode + delta_Ii + delta_Ie;
    }
        */
    //------------------------------------

    //Caluculate residuals
    //------------------------------------
    calcRes(pm.error_rhoe, gc.rhoe,gc.rhoe_old,gc.i_flc_bl, gc.j_flc_bl, pm.n_bl-1, 1e-100, pm.flag_error);
    calcRes(pm.error_rhoUex, gx.rhoUex, gx.rhoUex_old, gx.i_flx_bl, gx.j_flx_bl, pm.n_bl-1, 1e-100, pm.flag_error);
    calcRes(pm.error_rhoUer, gr.rhoUer, gr.rhoUer_old, gr.i_flr_bl, gr.j_flr_bl, pm.n_bl-1, 1e-100, pm.flag_error);
    //------------------------------------

}
//*****************************************************************
//**                                                             **
//**           void solve_Te_wdTe_wSEE_PC()                      **
//**                                                             **
//*****************************************************************
void FluidModule::solve_Te_wdTe_wSEE_PC(Params &pm,GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr, SeeVec &se){

    std::vector<std::vector<double> > aW(pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //係数
    std::vector<std::vector<double> > aE(pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //係数
    std::vector<std::vector<double> > aS(pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //係数
    std::vector<std::vector<double> > aN(pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //係数
    std::vector<std::vector<double> > aP(pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //係数
    std::vector<std::vector<double> > b (pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //ソース
    std::vector<std::vector<double> > rhoeps_old_tmp(pm.ni+2,std::vector<double>(pm.nj+2,0.0));
    std::vector<std::vector<double> > Gx_old_tmp(pm.ni+2,std::vector<double>(pm.nj+2,0.0));
    std::vector<std::vector<double> > Gr_old_tmp(pm.ni+2,std::vector<double>(pm.nj+2,0.0));
    std::vector<std::vector<double> > drhoeps(pm.ni+2,std::vector<double>(pm.nj+2,0.0));

    //Numerical treatment for exponential overflow and zero-division
    //------------------------------------
    double lim_exp = 1.0e2; //expの制限
    double EPS = 1.0e-3; //ゼロ割回避用の数字
    //------------------------------------

    //Acceleration factor for PC method
    //------------------------------------
    double c_mu = 1.0;
    //------------------------------------

    //convergence criteria
    //------------------------------------
    int flag_error2 = pm.flag_error; //擬似ステップのエラー収束判定 0:max, 1:RMS
    //------------------------------------

    //output duration
    //------------------------------------
    int ndiv_out = 100;
    //------------------------------------

    //scaling factor of wall flux
    //------------------------------------
    double alpha_wflux = 1.0;
    //------------------------------------
    
    //Reserve the past values
    //------------------------------------
    for (int i=0;i<=pm.ni+1;i++){ 
        for (int j=0;j<=pm.nj+1;j++){ 
            gc.Te_old[i][j] = gc.Te[i][j];
            gc.rhoeps_old[i][j] = gc.rhoeps[i][j];
            gx.Gx_old[i][j] = gx.Gx[i][j];
            gr.Gr_old[i][j] = gr.Gr[i][j];
        }
    }
    //------------------------------------

    //Initial prediction of correction
    //------------------------------------
    for (int i=0;i<=pm.ni+1;i++){ 
        for (int j=0;j<=pm.nj+1;j++){ 
            drhoeps[i][j] = gc.rhoeps[i][j] - gc.rhoeps_old[i][j];
        }
    }
    //------------------------------------
    
    //PC method iteration
    //------------------------------------
    int ncount = 0; //number of iteration of PC method
    int flag_end = 0;
    do{
        
        ncount ++;
        pm.itime_PC_rhoeps++;

        //Reserve previous-step values
        //------------------------------------
        for (int i=0;i<=pm.ni+1;i++){ 
            for (int j=0;j<=pm.nj+1;j++){ 
                rhoeps_old_tmp[i][j] = gc.rhoeps[i][j];
                Gx_old_tmp[i][j] = gx.Gx[i][j];
                Gr_old_tmp[i][j] = gr.Gr[i][j];
            }
        }
        //------------------------------------

        //Predict flux
        //==============================================
        //for rhoUex
        //------------------------------------
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
            for (int i=gx.i_flx_bl[iblock][0];i<=gx.i_flx_bl[iblock][1];i++){
                for (int j=gx.j_flx_bl[iblock][0];j<=gx.j_flx_bl[iblock][1];j++){
                
                    //radial position at upper and lower cell interface
                    //------------------------------------
                    double rL = (gc.r[j]+gc.r[j-1])/2.0;
                    double rR = (gc.r[j]+gc.r[j+1])/2.0;
                    //------------------------------------

                    //moblity at cell interface
                    //------------------------------------
                    double mu_para_tmp = (gc.mu_para[i][j] + gc.mu_para[i-1][j])/2.0;
                    double mu_perp_tmp = (gc.mu_perp[i][j] + gc.mu_perp[i-1][j])/2.0;
                    double mu_star_tmp = mu_perp_tmp*c_mu;
                    //------------------------------------

                    //calculate inverse moblity tensor at cell interface
                    //------------------------------------
                    double mu_r_tmp = (gc.mu_rr[i][j] + gc.mu_rr[i-1][j])/2.0;
                    double mu_c_tmp = (gc.mu_xr[i][j] + gc.mu_xr[i-1][j])/2.0;

                    double mu_x_inv_tmp =  mu_r_tmp/(mu_para_tmp*mu_perp_tmp);
                    double mu_c_inv_tmp = -mu_c_tmp/(mu_para_tmp*mu_perp_tmp);
                    //------------------------------------

                    //values at cell interface for Sharfetter-Gummmel scheme
                    //------------------------------------
                    double dTedx_tmp = (gc.Te_old[i][j] - gc.Te_old[i-1][j])/pm.dx; //electron temperature
                    double vd_tmp = mu_star_tmp*(-gx.Ex[i][j] - ph::Boltz/ph::e0*dTedx_tmp) + EPS; //electron temperature gradient
                    double Te_tmp = (gc.Te_old[i][j] + gc.Te_old[i-1][j])/2.0; //electron drift velocity
                    double D_tmp = ph::Boltz*Te_tmp/ph::e0*mu_star_tmp; //electron diffusion coefficient
                    double z_tmp = vd_tmp*pm.dx/D_tmp; //Peclet number
                    //------------------------------------

                    //BC - set adjacent flag (adjacent = 0, others = 1)
                    //adjacent "cell" to the wall for setting zero-flux on the wall
                    //------------------------------------
                    //lower
                    double bLr_wall  = double(j!=gc.j_flc_bl[0][0] || i>gc.i_flc_bl[0][1]); //x0
                    bLr_wall = bLr_wall*double(j!=gc.j_flc_bl[2][0] || i>gc.i_flc_bl[2][1]); //x2
                    double bLr_cen = double(j!=gc.j_flc_bl[3][0]); //x6
                    double bLr  = bLr_wall*bLr_cen;

                    //lower-left
                    double bLrD_wall  = double(j!=gc.j_flc_bl[0][0] || i-1>gc.i_flc_bl[0][1]); //x0
                    bLrD_wall = bLrD_wall*double(j!=gc.j_flc_bl[2][0] || i-1>gc.i_flc_bl[2][1]); //x2
                    double bLrD_cen = double(j!=gc.j_flc_bl[3][0]); //x6
                    double bLrD  = bLrD_wall*bLrD_cen;

                    //upper
                    double bRr_wall  = double(j!=gc.j_flc_bl[0][1] || iblock > 2); //x1
                    bRr_wall = bRr_wall*double(j!=gc.j_flc_bl[3][1] || i<gc.i_flc_bl[1][1]+1 || i>gc.i_flc_bl[3][1]); //x4
                    //double bRr_open  = (j!=gc.j_flc_bl[4][1] || iblock != 4); //x5
                    double bRr = bRr_wall;

                    //upper-left
                    double bRrD_wall  = double(j!=gc.j_flc_bl[0][1] || iblock > 2); //x1
                    bRrD_wall = bRrD_wall*double(j!=gc.j_flc_bl[3][1] || i-1<gc.i_flc_bl[1][1]+1 || i-1>gc.i_flc_bl[3][1]); //x4
                    //double bRrD_open  = (j!=gc.j_flc_bl[4][1] || iblock != 4); //x5
                    double bRrD = bRrD_wall;
                    //------------------------------------

                    //interpolation of transverse flux component (should be old value)
                    //------------------------------------
                    double Gr_tmp = (rR*Gr_old_tmp[i][j+1]*bRr + rL*Gr_old_tmp[i][j]*bLr 
                        + rR*Gr_old_tmp[i-1][j+1]*bRrD + rL*Gr_old_tmp[i-1][j]*bLrD)/(4.0*gc.r[j]);
                    //------------------------------------

                    //update flux
                    //------------------------------------
                    gx.Gx[i][j] =
                        ( vd_tmp/(1.0 - exp(fmin(-z_tmp,lim_exp)))*gc.rhoeps[i-1][j]
                        - vd_tmp/(exp(fmin( z_tmp,lim_exp)) - 1.0)*gc.rhoeps[i][j]
                        )*5.0/3.0 
                        + Gx_old_tmp[i][j] 
                        - mu_star_tmp*(mu_x_inv_tmp*Gx_old_tmp[i][j] + mu_c_inv_tmp*Gr_tmp);
                    //------------------------------------
                }
            }     
        }
        //------------------------------------

        
        //left-wall-BC (z0)
        //------------------------------------
        for (int j=gc.j_flc_bl[0][0];j<=gc.j_flc_bl[0][1];j++){
            int i=gc.i_flc_bl[0][0];
            double vth = sqrt(8.0*ph::Boltz*gc.Te_old[i][j]/(M_PI*pm.masse));
            double delta_phi = fmax(-gx.Ex[i][j]*pm.dx/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerNeps_wall_Lx = -vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te_old[i][j] + 1e-100),lim_exp))/3.0*alpha_wflux;
            double eps_emitSEE = 2.0*ph::e0*pm.Te_emitSEE_eV;
            double GammaPerN_wall_Lx = -0.25*alpha_wflux*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te_old[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Lx = - eps_emitSEE*pm.coefIISEE_ts*gc.rhoi[i][j]*gx.Uix[i][j] //IISEE
                            - (coefEISEE_eb + pm.ratioEngy_EISEE_rd*coefEISEE_rd)*gc.rhoeps_old[i][j]*GammaPerNeps_wall_Lx //EISEE
                            - eps_emitSEE*coefEISEE_ts*gc.rhoe_old[i][j]*GammaPerN_wall_Lx //EISEE
                            - eps_emitSEE*pm.coefMISEE_ts*gx.rhoUmx[i][j]; //MISEE

            gx.Gx[i][j] = gc.rhoeps[i][j]*GammaPerNeps_wall_Lx + Gamma_wallSEE_Lx;
            gx.Gx_wall[i][j] = gx.Gx[i][j];
        }
        //------------------------------------

        //left-wall-BC (z1)
        //------------------------------------
        for (int j=gc.j_flc_bl[2][0];j<=gc.j_flc_bl[0][0]-1;j++){
            int i=gc.i_flc_bl[2][0];
            double vth = sqrt(8.0*ph::Boltz*gc.Te_old[i][j]/(M_PI*pm.masse));
            double delta_phi = fmax(-gx.Ex[i][j]*pm.dx/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerNeps_wall_Lx = -vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te_old[i][j] + 1e-100),lim_exp))/3.0*alpha_wflux;
            double eps_emitSEE = 2.0*ph::e0*pm.Te_emitSEE_eV;
            double GammaPerN_wall_Lx = -0.25*alpha_wflux*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te_old[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Lx = - eps_emitSEE*pm.coefIISEE_ts*gc.rhoi[i][j]*gx.Uix[i][j] //IISEE
                            - (coefEISEE_eb + pm.ratioEngy_EISEE_rd*coefEISEE_rd)*gc.rhoeps_old[i][j]*GammaPerNeps_wall_Lx //EISEE
                            - eps_emitSEE*coefEISEE_ts*gc.rhoe_old[i][j]*GammaPerN_wall_Lx //EISEE
                            - eps_emitSEE*pm.coefMISEE_ts*gx.rhoUmx[i][j]; //MISEE

            gx.Gx[i][j] = gc.rhoeps[i][j]*GammaPerNeps_wall_Lx + Gamma_wallSEE_Lx;
            gx.Gx_wall[i][j] = gx.Gx[i][j];
        }
        //------------------------------------

        //left-wall-BC (z2)
        //------------------------------------
        for (int j=gc.j_flc_bl[3][0];j<=gc.j_flc_bl[2][0]-1;j++){
            int i=gc.i_flc_bl[3][0];
            double vth = sqrt(8.0*ph::Boltz*gc.Te_old[i][j]/(M_PI*pm.masse));
            double delta_phi = fmax(-gx.Ex[i][j]*pm.dx/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerNeps_wall_Lx = -vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te_old[i][j] + 1e-100),lim_exp))/3.0*alpha_wflux;
            double eps_emitSEE = 2.0*ph::e0*pm.Te_emitSEE_eV;
            double GammaPerN_wall_Lx = -0.25*alpha_wflux*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te_old[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Lx = - eps_emitSEE*pm.coefIISEE_ts*gc.rhoi[i][j]*gx.Uix[i][j] //IISEE
                            - (coefEISEE_eb + pm.ratioEngy_EISEE_rd*coefEISEE_rd)*gc.rhoeps_old[i][j]*GammaPerNeps_wall_Lx //EISEE
                            - eps_emitSEE*coefEISEE_ts*gc.rhoe_old[i][j]*GammaPerN_wall_Lx //EISEE
                            - eps_emitSEE*pm.coefMISEE_ts*gx.rhoUmx[i][j]; //MISEE

            gx.Gx[i][j] = gc.rhoeps[i][j]*GammaPerNeps_wall_Lx + Gamma_wallSEE_Lx;
            gx.Gx_wall[i][j] = gx.Gx[i][j];
        }
        //------------------------------------

        //left-wall-BC (z4)
        //------------------------------------
        for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[4][1];j++){
            int i=gc.i_flc_bl[4][0];
            double vth = sqrt(8.0*ph::Boltz*gc.Te_old[i][j]/(M_PI*pm.masse));
            double delta_phi = fmax(-gx.Ex[i][j]*pm.dx/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerNeps_wall_Lx = -vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te_old[i][j] + 1e-100),lim_exp))/3.0;
            double eps_emitSEE = 2.0*ph::e0*pm.Te_emitSEE_eV;
            double GammaPerN_wall_Lx = -0.25*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te_old[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Lx = - eps_emitSEE*pm.coefIISEE_ts*gc.rhoi[i][j]*gx.Uix[i][j] //IISEE
                            - (coefEISEE_eb + pm.ratioEngy_EISEE_rd*coefEISEE_rd)*gc.rhoeps_old[i][j]*GammaPerNeps_wall_Lx //EISEE
                            - eps_emitSEE*coefEISEE_ts*gc.rhoe_old[i][j]*GammaPerN_wall_Lx //EISEE
                            - eps_emitSEE*pm.coefMISEE_ts*gx.rhoUmx[i][j]; //MISEE

            gx.Gx[i][j] = gc.rhoeps[i][j]*GammaPerNeps_wall_Lx + Gamma_wallSEE_Lx;
            gx.Gx_wall[i][j] = gx.Gx[i][j];
        }
        //------------------------------------

        //right-wall-BC (z3)
        //------------------------------------
        for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[1][1];j++){
            int i=gc.i_flc_bl[1][1];
            double vth = sqrt(8.0*ph::Boltz*gc.Te_old[i][j]/(M_PI*pm.masse));
            double delta_phi =  fmax(gx.Ex[i+1][j]*pm.dx/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerNeps_wall_Rx = vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te_old[i][j] + 1e-100),lim_exp))/3.0*alpha_wflux;
            double eps_emitSEE = 2.0*ph::e0*pm.Te_emitSEE_eV;
            double GammaPerN_wall_Rx = 0.25*alpha_wflux*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te_old[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Rx = - eps_emitSEE*pm.coefIISEE_ts*gc.rhoi[i][j]*gx.Uix[i+1][j] //IISEE
                - (coefEISEE_eb + pm.ratioEngy_EISEE_rd*coefEISEE_rd)*gc.rhoeps_old[i][j]*GammaPerNeps_wall_Rx //EISEE
                - eps_emitSEE*coefEISEE_ts*gc.rhoe_old[i][j]*GammaPerN_wall_Rx //EISEE
                - eps_emitSEE*pm.coefMISEE_ts*gx.rhoUmx[i+1][j]; //MISEE

            gx.Gx[i+1][j] = gc.rhoeps[i][j]*GammaPerNeps_wall_Rx + Gamma_wallSEE_Rx;
            gx.Gx_wall[i+1][j] = gx.Gx[i+1][j];
        }
        //------------------------------------

        //right-wall-BC (z5)
        //------------------------------------
        for (int j=gc.j_flc_bl[4][0];j<=gc.j_flc_bl[4][1];j++){
            int i=gc.i_flc_bl[4][1];
            double vth = sqrt(8.0*ph::Boltz*gc.Te_old[i][j]/(M_PI*pm.masse));
            double delta_phi =  fmax(gx.Ex[i+1][j]*pm.dx/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerNeps_wall_Rx = vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te_old[i][j] + 1e-100),lim_exp))/3.0;
            double eps_emitSEE = 2.0*ph::e0*pm.Te_emitSEE_eV;
            double GammaPerN_wall_Rx = 0.25*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te_old[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Rx = - eps_emitSEE*pm.coefIISEE_ts*gc.rhoi[i][j]*gx.Uix[i+1][j] //IISEE
                - (coefEISEE_eb + pm.ratioEngy_EISEE_rd*coefEISEE_rd)*gc.rhoeps_old[i][j]*GammaPerNeps_wall_Rx //EISEE
                - eps_emitSEE*coefEISEE_ts*gc.rhoe_old[i][j]*GammaPerN_wall_Rx //EISEE
                - eps_emitSEE*pm.coefMISEE_ts*gx.rhoUmx[i+1][j]; //MISEE

            gx.Gx[i+1][j] = gc.rhoeps[i][j]*GammaPerNeps_wall_Rx + Gamma_wallSEE_Rx;
            gx.Gx_wall[i+1][j] = gx.Gx[i+1][j];
        }
        //------------------------------------

        //for rhoUer
        //------------------------------------
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){
            for (int i=gr.i_flr_bl[iblock][0];i<=gr.i_flr_bl[iblock][1];i++){ 
                for (int j=gr.j_flr_bl[iblock][0];j<=gr.j_flr_bl[iblock][1];j++){
                    
                    //radial position at cell interface
                    //------------------------------------
                    double r_tmp = (gc.r[j] + gc.r[j-1])/2.0;
                    //------------------------------------

                    //moblity at cell interface
                    //------------------------------------
                    double mu_para_tmp = (gc.mu_para[i][j] + gc.mu_para[i][j-1])/2.0;
                    double mu_perp_tmp = (gc.mu_perp[i][j] + gc.mu_perp[i][j-1])/2.0;

                    double mu_star_tmp = mu_perp_tmp*c_mu;
                    //------------------------------------

                    //calculate inverse moblity tensor at cell interface
                    //------------------------------------
                    double mu_x_tmp = (gc.mu_xx[i][j] + gc.mu_xx[i][j-1])/2.0;
                    double mu_c_tmp = (gc.mu_xr[i][j] + gc.mu_xr[i][j-1])/2.0;

                    double mu_r_inv_tmp = mu_x_tmp/(mu_para_tmp*mu_perp_tmp);
                    double mu_c_inv_tmp = -mu_c_tmp/(mu_para_tmp*mu_perp_tmp);
                    //------------------------------------
    
                    //values at cell interface for Sharfetter-Gummmel scheme
                    //------------------------------------
                    double dTedr_tmp = (gc.Te_old[i][j] - gc.Te_old[i][j-1])/pm.dr; //electron temperature
                    double vd_tmp = mu_star_tmp*(-gr.Er[i][j] - ph::Boltz/ph::e0*dTedr_tmp) + EPS; //electron temperature gradient
                    double Te_tmp = (gc.Te_old[i][j] + gc.Te_old[i][j-1])/2.0; //electron drift velocity
                    double D_tmp = ph::Boltz*Te_tmp/ph::e0*mu_star_tmp; //electron diffusion coefficient
                    double z_tmp = vd_tmp*pm.dr/D_tmp; //Peclet number
                    //------------------------------------

                    //BC - set adjacent flag (adjacent = 0, others = 1)
                    //adjacent "cell" to the wall for setting zero-flux on the wall
                    //------------------------------------
                    //left
                    double bLx_wall  = double(i!=gc.i_flc_bl[0][0] || iblock != 0); //z0
                    bLx_wall = bLx_wall*double(i!=gc.i_flc_bl[2][0] || j>gc.j_flc_bl[0][0]-1); //z1
                    bLx_wall = bLx_wall*double(i!=gc.i_flc_bl[3][0] || j>gc.j_flc_bl[2][0]-1); //z2
                    bLx_wall = bLx_wall*double(i!=gc.i_flc_bl[4][0] || j<gc.j_flc_bl[1][0]); //z4
                    double bLx  = bLx_wall;

                    //left-lower
                    double bLxD_wall  = double(i!=gc.i_flc_bl[0][0] || iblock != 0); //z0
                    bLxD_wall = bLxD_wall*double(i!=gc.i_flc_bl[2][0] || j-1>gc.j_flc_bl[0][0]-1); //z1
                    bLxD_wall = bLxD_wall*double(i!=gc.i_flc_bl[3][0] || j-1>gc.j_flc_bl[2][0]-1); //z2
                    bLxD_wall = bLxD_wall*double(i!=gc.i_flc_bl[4][0] || j-1<gc.j_flc_bl[1][0]); //z4
                    double bLxD  = bLxD_wall;

                    //right
                    double bRx_wall  = double(i!=gc.i_flc_bl[1][1] || j<gc.j_flc_bl[1][0]); //z3
                    bRx_wall = bRx_wall*double(i!=gc.i_flc_bl[4][1]); //z5
                    double bRx  = bRx_wall;

                    //right-lower
                    double bRxD_wall  = double(i!=gc.i_flc_bl[1][1] || j-1<gc.j_flc_bl[1][0]); //z3
                    bRxD_wall = bRxD_wall*double(i!=gc.i_flc_bl[4][1]); //z5
                    double bRxD  = bRxD_wall;
                    //------------------------------------

                    //interpolation of transverse flux component (should be old value)
                    //------------------------------------
                    double Gx_tmp = (gc.r[j]*Gx_old_tmp[i+1][j]*bRx + gc.r[j]*Gx_old_tmp[i][j]*bLx
                        + gc.r[j-1]*Gx_old_tmp[i+1][j-1]*bRxD + gc.r[j-1]*Gx_old_tmp[i][j-1]*bLxD)/(4.0*r_tmp);
                    //------------------------------------

                    //update flux
                    //------------------------------------
                    gr.Gr[i][j] =  
                        ( vd_tmp/(1.0 - exp(fmin(-z_tmp,lim_exp)))*gc.rhoeps[i][j-1]
                        - vd_tmp/(exp(fmin( z_tmp,lim_exp)) - 1.0)*gc.rhoeps[i][j]
                        )*5.0/3.0
                        + Gr_old_tmp[i][j] 
                        - mu_star_tmp*(mu_c_inv_tmp*Gx_tmp + mu_r_inv_tmp*Gr_old_tmp[i][j]);
                    //------------------------------------
                }
            }
        }
        //------------------------------------

        //lower-wall-BC (x0)
        //------------------------------------
        for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[0][1];i++){
            int j=gc.j_flc_bl[0][0];
            double vth = sqrt(8.0*ph::Boltz*gc.Te_old[i][j]/(M_PI*pm.masse));
            double delta_phi = fmax(-gr.Er[i][j]*pm.dr/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerNeps_wall_Lr = -vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te_old[i][j] + 1e-100),lim_exp))/3.0*alpha_wflux;
            double eps_emitSEE = 2.0*ph::e0*pm.Te_emitSEE_eV;
            double GammaPerN_wall_Lr = -0.25*alpha_wflux*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te_old[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Lr = - eps_emitSEE*pm.coefIISEE_ts*gc.rhoi[i][j]*gr.Uir[i][j] //IISEE
                - (coefEISEE_eb + pm.ratioEngy_EISEE_rd*coefEISEE_rd)*gc.rhoeps_old[i][j]*GammaPerNeps_wall_Lr //EISEE
                - eps_emitSEE*coefEISEE_ts*gc.rhoe_old[i][j]*GammaPerN_wall_Lr //EISEE
                - eps_emitSEE*pm.coefMISEE_ts*gr.rhoUmr[i][j]; //MISEE
             
            gr.Gr[i][j] = gc.rhoeps[i][j]*GammaPerNeps_wall_Lr + Gamma_wallSEE_Lr;
            gr.Gr_wall[i][j] = gr.Gr[i][j];
        }
        //------------------------------------
    
        //lower-wall-BC (x2)
        //------------------------------------
        for (int i=gc.i_flc_bl[2][0];i<=gc.i_flc_bl[2][1];i++){
            int j=gc.j_flc_bl[2][0];
            double vth = sqrt(8.0*ph::Boltz*gc.Te_old[i][j]/(M_PI*pm.masse));
            double delta_phi = fmax(-gr.Er[i][j]*pm.dr/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerNeps_wall_Lr = -vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te_old[i][j] + 1e-100),lim_exp))/3.0*alpha_wflux;
            double eps_emitSEE = 2.0*ph::e0*pm.Te_emitSEE_eV;
            double GammaPerN_wall_Lr = -0.25*alpha_wflux*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te_old[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Lr = - eps_emitSEE*pm.coefIISEE_ts*gc.rhoi[i][j]*gr.Uir[i][j] //IISEE
                - (coefEISEE_eb + pm.ratioEngy_EISEE_rd*coefEISEE_rd)*gc.rhoeps_old[i][j]*GammaPerNeps_wall_Lr //EISEE
                - eps_emitSEE*coefEISEE_ts*gc.rhoe_old[i][j]*GammaPerN_wall_Lr //EISEE
                - eps_emitSEE*pm.coefMISEE_ts*gr.rhoUmr[i][j]; //MISEE
        
            gr.Gr[i][j] = gc.rhoeps[i][j]*GammaPerNeps_wall_Lr + Gamma_wallSEE_Lr;
        }
        //------------------------------------
    
        //upper-wall-BC (x1)
        //------------------------------------
        for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[1][1];i++){
            int j=gc.j_flc_bl[0][1];
            double vth = sqrt(8.0*ph::Boltz*gc.Te_old[i][j]/(M_PI*pm.masse));
            double delta_phi =  fmax(gr.Er[i][j+1]*pm.dr/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerNeps_wall_Rr = vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te_old[i][j] + 1e-100),lim_exp))/3.0*alpha_wflux;
            double eps_emitSEE = 2.0*ph::e0*pm.Te_emitSEE_eV;
            double GammaPerN_wall_Rr = 0.25*alpha_wflux*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te_old[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Rr = - eps_emitSEE*pm.coefIISEE_ts*gc.rhoi[i][j]*gr.Uir[i][j+1] //IISEE
                - (coefEISEE_eb + pm.ratioEngy_EISEE_rd*coefEISEE_rd)*gc.rhoeps_old[i][j]*GammaPerNeps_wall_Rr //EISEE
                - eps_emitSEE*coefEISEE_ts*gc.rhoe_old[i][j]*GammaPerN_wall_Rr //EISEE
                - eps_emitSEE*pm.coefMISEE_ts*gr.rhoUmr[i][j+1]; //MISEE
        
            gr.Gr[i][j+1] = gc.rhoeps[i][j]*GammaPerNeps_wall_Rr + Gamma_wallSEE_Rr;
            gr.Gr_wall[i][j+1] = gr.Gr[i][j+1];
        }
        //------------------------------------
    
        //upper-wall-BC (x4)
        //------------------------------------
        for (int i=gc.i_flc_bl[1][1]+1;i<=gc.i_flc_bl[3][1];i++){
            int j=gc.j_flc_bl[3][1];
            double vth = sqrt(8.0*ph::Boltz*gc.Te_old[i][j]/(M_PI*pm.masse));
            double delta_phi =  fmax(gr.Er[i][j+1]*pm.dr/2.0,0.0); //adjacent cell-center potential from wall
            double GammaPerNeps_wall_Rr = vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te_old[i][j] + 1e-100),lim_exp))/3.0*alpha_wflux;
            double eps_emitSEE = 2.0*ph::e0*pm.Te_emitSEE_eV;
            double GammaPerN_wall_Rr = 0.25*alpha_wflux*vth*exp(fmin(-ph::e0*delta_phi/(ph::Boltz*gc.Te_old[i][j] + 1e-100),lim_exp));

            double coefEISEE_eb = func_Te(se.coefEISEE_eb_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_rd = func_Te(se.coefEISEE_rd_table,gc.Te[i][j],pm.dTe_SEE);
            double coefEISEE_ts = func_Te(se.coefEISEE_ts_table,gc.Te[i][j],pm.dTe_SEE);

            double Gamma_wallSEE_Rr = - eps_emitSEE*pm.coefIISEE_ts*gc.rhoi[i][j]*gr.Uir[i][j+1] //IISEE
                - (coefEISEE_eb + pm.ratioEngy_EISEE_rd*coefEISEE_rd)*gc.rhoeps_old[i][j]*GammaPerNeps_wall_Rr //EISEE
                - eps_emitSEE*coefEISEE_ts*gc.rhoe_old[i][j]*GammaPerN_wall_Rr //EISEE
                - eps_emitSEE*pm.coefMISEE_ts*gr.rhoUmr[i][j+1]; //MISEE
        
            gr.Gr[i][j+1] = gc.rhoeps[i][j]*GammaPerNeps_wall_Rr + Gamma_wallSEE_Rr;
            gr.Gr_wall[i][j+1] = gr.Gr[i][j+1];
        }
        //------------------------------------
    
        //upper-open-BC (x5)
        //------------------------------------
        for (int i=gc.i_flc_bl[4][0];i<=gc.i_flc_bl[4][1];i++){
            int j=gc.j_flc_bl[4][1];
            double Gamma_open_Rr = 5.0/2.0*ph::Boltz*gc.Te_old[i][j]*fmax(gr.rhoUer[i][j+1],0.0);
            gr.Gr[i][j+1] = Gamma_open_Rr;
            gr.Gr_wall[i][j+1] = 0.0;
        }
        //------------------------------------

        //Make coefficient for electron energy conservation for 5 point stencil
        //------------------------------------
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
            for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
                for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                    //radial position at upper and lower cell interface
                    //------------------------------------
                    double rR = (gc.r[j]+gc.r[j+1])/2.0;
                    double rL = (gc.r[j]+gc.r[j-1])/2.0;
                    //------------------------------------

                    //normalized radial position at upper and lower cell interface
                    //------------------------------------
                    double qR = (gc.r[j]+gc.r[j+1])/2.0/gc.r[j];
                    double qL = (gc.r[j]+gc.r[j-1])/2.0/gc.r[j];
                    //------------------------------------

                    //moblity at cell interface
                    //------------------------------------
                    double mu_star_Lx = (gc.mu_perp[i][j]+gc.mu_perp[i-1][j])/2.0*c_mu;
                    double mu_star_Rx = (gc.mu_perp[i][j]+gc.mu_perp[i+1][j])/2.0*c_mu;
                    double mu_star_Lr = (gc.mu_perp[i][j]+gc.mu_perp[i][j-1])/2.0*c_mu;
                    double mu_star_Rr = (gc.mu_perp[i][j]+gc.mu_perp[i][j+1])/2.0*c_mu;
                    //------------------------------------

                    //electron temperature at cell interface
                    //------------------------------------
                    double Te_Lx = (gc.Te_old[i][j]+gc.Te_old[i-1][j])/2.0;
                    double Te_Rx = (gc.Te_old[i][j]+gc.Te_old[i+1][j])/2.0;
                    double Te_Lr = (gc.Te_old[i][j]+gc.Te_old[i][j-1])/2.0;
                    double Te_Rr = (gc.Te_old[i][j]+gc.Te_old[i][j+1])/2.0;
                    //------------------------------------

                    //electron diffusion coefficient at cell interface
                    //------------------------------------
                    double D_Lx = ph::Boltz*Te_Lx/ph::e0*mu_star_Lx;
                    double D_Rx = ph::Boltz*Te_Rx/ph::e0*mu_star_Rx;
                    double D_Lr = ph::Boltz*Te_Lr/ph::e0*mu_star_Lr;
                    double D_Rr = ph::Boltz*Te_Rr/ph::e0*mu_star_Rr;
                    //------------------------------------

                    //electron temperature gradient at cell interface
                    //------------------------------------
                    double dTedx_Lx = (gc.Te_old[i  ][j  ] - gc.Te_old[i-1][j  ])/pm.dx;
                    double dTedx_Rx = (gc.Te_old[i+1][j  ] - gc.Te_old[i  ][j  ])/pm.dx;
                    double dTedr_Lr = (gc.Te_old[i  ][j  ] - gc.Te_old[i  ][j-1])/pm.dr;
                    double dTedr_Rr = (gc.Te_old[i  ][j+1] - gc.Te_old[i  ][j  ])/pm.dr;
                    //------------------------------------

                    //electron drift velocity at cell interface
                    //------------------------------------
                    double vd_Lx = mu_star_Lx*(-gx.Ex[i  ][j  ] - ph::Boltz/ph::e0*dTedx_Lx) + EPS;
                    double vd_Rx = mu_star_Rx*(-gx.Ex[i+1][j  ] - ph::Boltz/ph::e0*dTedx_Rx) + EPS;
                    double vd_Lr = mu_star_Lr*(-gr.Er[i  ][j  ] - ph::Boltz/ph::e0*dTedr_Lr) + EPS;
                    double vd_Rr = mu_star_Rr*(-gr.Er[i  ][j+1] - ph::Boltz/ph::e0*dTedr_Rr) + EPS;
                    //------------------------------------

                    //Peclet number at cell interface
                    //------------------------------------
                    double z_Lx = vd_Lx*pm.dx/D_Lx;
                    double z_Rx = vd_Rx*pm.dx/D_Rx;
                    double z_Lr = vd_Lr*pm.dr/D_Lr;
                    double z_Rr = vd_Rr*pm.dr/D_Rr;
                    //------------------------------------

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

                    //coefficient for electron energy conservation for 5 point stencil
                    //==============================================
                    //vertical and horizental 5 stencil
                    //------------------------------------
                    aP[i][j] =    vd_Rx/(1.0 - exp(fmin(-z_Rx,lim_exp)))*bRx*pm.dr 
                             +    vd_Lx/(exp(fmin( z_Lx,lim_exp)) - 1.0)*bLx*pm.dr
                             + qR*vd_Rr/(1.0 - exp(fmin(-z_Rr,lim_exp)))*bRr*pm.dx
                             + qL*vd_Lr/(exp(fmin( z_Lr,lim_exp)) - 1.0)*bLr*pm.dx
                             //+ 3.0/5.0*delossdeps[i][j]*pm.dx*pm.dr //part of implicit source term prediction
                             + 3.0/5.0*pm.dx*pm.dr/pm.dt; //time-term
                    aE[i][j] =    vd_Rx/(exp(fmin( z_Rx,lim_exp)) - 1.0)*bRx*pm.dr;
                    aW[i][j] =    vd_Lx/(1.0 - exp(fmin(-z_Lx,lim_exp)))*bLx*pm.dr;
                    aN[i][j] = qR*vd_Rr/(exp(fmin( z_Rr,lim_exp)) - 1.0)*bRr*pm.dx;
                    aS[i][j] = qL*vd_Lr/(1.0 - exp(fmin(-z_Lr,lim_exp)))*bLr*pm.dx;
                    //------------------------------------

                    //RHS
                    //------------------------------------
                    double nabla_G = (gx.Gx[i+1][j]-gx.Gx[i][j])/pm.dx + (qR*gr.Gr[i][j+1]-qL*gr.Gr[i][j])/pm.dr;
                    //double rhoUex_tmp = (gx.rhoUex[i+1][j]*bRx_wall + gx.rhoUex[i][j]*bLx_wall)/2.0;
                    //double rhoUer_tmp = (rR*gr.rhoUer[i][j+1]*bRr_wall + rL*gr.rhoUer[i][j]*bLr_wall)/(2.0*gc.r[j]);
                    //double Ex_tmp = (gx.Ex[i+1][j] + gx.Ex[i][j])/2.0;
                    //double Er_tmp = (rR*gr.Er[i][j+1] + rL*gr.Er[i][j])/(2.0*gc.r[j]);
                    double JdotE_x = -ph::e0*(gx.rhoUex[i+1][j]*gx.Ex[i+1][j]*bRx_wall + gx.rhoUex[i][j]*gx.Ex[i][j]*bLx_wall)/2.0;
                    double JdotE_r = -ph::e0*(rR*gr.rhoUer[i][j+1]*gr.Er[i][j+1]*bRr_wall + rL*gr.rhoUer[i][j]*gr.Er[i][j]*bLr_wall)/(2.0*gc.r[j]);

                    double RHS1 = (-gc.rhoeps[i][j] + gc.rhoeps_old[i][j])/pm.dt;
                    //double RHS2 = -ph::e0*(Ex_tmp*rhoUex_tmp + Er_tmp*rhoUer_tmp);
                    double RHS2 = JdotE_x + JdotE_r;
                    double RHS3 = gc.Pabs[i][j] 
                        - gc.rhoe[i][j]*gc.rate_eloss[i][j];
                        //- delossdeps[i][j]*gc.rhoeps[i][j]
                        //+ gc.rhoe[i][j]*(1.5*ph::Boltz*gc.Te[i][j])*delossdeps[i][j]; //part of implicit source term prediction
                        
                    double RHS4 = -nabla_G;

                    b[i][j] = 3.0/5.0*(RHS1 + RHS2 + RHS3 + RHS4)*pm.dx*pm.dr;
                    //------------------------------------
                }
            }
        }
        //------------------------------------

        //solve
        if(pm.flag_iter_rhoeps==0){
            //SOR
            //------------------------------------
            double alpha = 1.9;
            solver_SOR(aP,aE,aW,aN,aS,b,gc.i_flc_bl,gc.j_flc_bl,pm.n_bl-1,alpha,pm.maxITR_SOR_rhoeps,pm.error_cnv_SOR_rhoeps,0,drhoeps);
            //------------------------------------
        }else{
            //------------------------------------
            solver_SMG(aP,aE,aW,aN,aS,b,gc.i_flc_bl,gc.j_flc_bl,pm.n_bl-1,pm.maxITR_SOR_rhoeps,pm.error_cnv_SOR_rhoeps,0,drhoeps);
            //------------------------------------
        }

        //Correction
        //==============================================
        //for Gx
        //------------------------------------
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
            for (int i=gx.i_flx_bl[iblock][0];i<=gx.i_flx_bl[iblock][1];i++){
                for (int j=gx.j_flx_bl[iblock][0];j<=gx.j_flx_bl[iblock][1];j++){
                    
                    //values at cell interface for Sharfetter-Gummmel scheme
                    //------------------------------------
                    double mu_star_tmp = (gc.mu_perp[i][j] + gc.mu_perp[i-1][j])/2.0*c_mu;
                    double dTedx_tmp = (gc.Te_old[i][j] - gc.Te_old[i-1][j])/pm.dx;
                    double vd_tmp = mu_star_tmp*(-gx.Ex[i][j] - ph::Boltz/ph::e0*dTedx_tmp) + EPS;
                    double Te_tmp = (gc.Te_old[i][j] + gc.Te_old[i-1][j])/2.0;
                    double D_tmp = ph::Boltz*Te_tmp/ph::e0*mu_star_tmp;
                    double z_tmp = vd_tmp*pm.dx/D_tmp;
                    //------------------------------------
                    
                    //calculate flux correction
                    //------------------------------------
                    double dGx = (vd_tmp/(1.0 - exp(fmin(-z_tmp,lim_exp)))*drhoeps[i-1][j]
                        - vd_tmp/(exp(fmin( z_tmp,lim_exp)) - 1.0)*drhoeps[i][j]
                        )*5.0/3.0;
                    //------------------------------------
                    
                    //Update
                    //------------------------------------
                    gx.Gx[i][j] = gx.Gx[i][j] + dGx;
                    //------------------------------------
                }
            }   
        }
        //------------------------------------

        //for Gr
        //------------------------------------
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){
            for (int i=gr.i_flr_bl[iblock][0];i<=gr.i_flr_bl[iblock][1];i++){ 
                for (int j=gr.j_flr_bl[iblock][0];j<=gr.j_flr_bl[iblock][1];j++){
                    
                    //values at cell interface for Sharfetter-Gummmel scheme
                    //------------------------------------
                    double mu_star_tmp = (gc.mu_perp[i][j] + gc.mu_perp[i][j-1])/2.0*c_mu;
                    double dTedr_tmp = (gc.Te_old[i  ][j  ] - gc.Te_old[i  ][j-1])/pm.dr;
                    double vd_tmp = mu_star_tmp*(-gr.Er[i][j] - ph::Boltz/ph::e0*dTedr_tmp) + EPS;
                    double Te_tmp = (gc.Te_old[i][j] + gc.Te_old[i][j-1])/2.0;
                    double D_tmp = ph::Boltz*Te_tmp/ph::e0*mu_star_tmp;
                    double z_tmp = vd_tmp*pm.dr/D_tmp;
                    //------------------------------------

                    //calculate flux correction
                    //------------------------------------
                    double dGr =  (vd_tmp/(1.0 - exp(fmin(-z_tmp,lim_exp)))*drhoeps[i][j-1]
                        - vd_tmp/(exp(fmin( z_tmp,lim_exp)) - 1.0)*drhoeps[i][j]
                        )*5.0/3.0;
                    //------------------------------------
                    
                    //Update
                    //------------------------------------
                    gr.Gr[i][j] = gr.Gr[i][j] + dGr;
                    //------------------------------------
                }
            }
        }
        //------------------------------------
        
        //for rhoeps
        //------------------------------------
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
            for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
                for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                    //gc.rhoeps[i][j] = fmax(gc.rhoeps[i][j] + drhoeps[i][j],1e-6);
                    double rhoeps_new = gc.rhoeps[i][j] + drhoeps[i][j];
                    if(rhoeps_new > 0.0){
                        gc.rhoeps[i][j] = rhoeps_new;
                    }
                }
            }
        }
        //------------------------------------
        //==============================================
       
        //Caluculate residuals
        //------------------------------------
        double error_rhoeps2 = 0.0;
        double error_Gx2 = 0.0;
        double error_Gr2 = 0.0;
        calcRes(error_rhoeps2, gc.rhoeps,  rhoeps_old_tmp,  gc.i_flc_bl, gc.j_flc_bl, pm.n_bl-1,1e-100,flag_error2);
        calcRes(error_Gx2,gx.Gx, Gx_old_tmp, gx.i_flx_bl, gx.j_flx_bl, pm.n_bl-1,1e-100,flag_error2);
        calcRes(error_Gr2,gr.Gr, Gr_old_tmp, gr.i_flr_bl, gr.j_flr_bl, pm.n_bl-1,1e-100,flag_error2);
        
        if(pm.flag_chk == 1) output_residual_for_PC("rhoeps",pm.gtime, pm.itime, error_rhoeps2, error_Gx2, error_Gr2, pm.itime_PC_rhoeps, ncount);
        
        double error_global = fmax(fmax(error_rhoeps2,error_Gx2),error_Gr2);
        
        if(ncount % ndiv_out == 0){
            std::cout << "ncount = "<<ncount<<std::endl;
            std::cout << "error_rhoeps = " << error_rhoeps2
                <<  " error_Gx = " << error_Gx2
                <<  " error_Gr= " << error_Gr2
                <<  " error_max = " << fmax(fmax(error_rhoeps2,error_Gx2),error_Gr2)
                << std::endl;
        }

        if(error_global < pm.error_cnv_HES_rhoeps || ncount >= pm.maxITR_HES_rhoeps){
            flag_end = 1;
        }
    
    } while (flag_end == 0);

    
    //calculation of energy balance
    //------------------------------------
    if(pm.itime%2000 == 0 || pm.itime == pm.ntime){
        
        //volume integrated values
        //------------------------------------
        double dEdt_sum          = 0.0;
        double nabla_G_sum       = 0.0;
        double rate_eloss_sum    = 0.0;
        double Joule_sum         = 0.0;
        double input_sum         = 0.0;
        //------------------------------------

        //calculate volume integrated of each term
        //------------------------------------
        for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
            for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
                for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                    
                    //radial position at upper and lower cell interface
                    //------------------------------------
                    double rR = (gc.r[j]+gc.r[j+1])/2.0;
                    double rL = (gc.r[j]+gc.r[j-1])/2.0;
                    //------------------------------------

                    //normalized radial position at upper and lower cell interface
                    //------------------------------------
                    double qR = (gc.r[j]+gc.r[j+1])/2.0/gc.r[j];
                    double qL = (gc.r[j]+gc.r[j-1])/2.0/gc.r[j];
                    //------------------------------------

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

                    //RHS
                    //------------------------------------
                    double nabla_G = (gx.Gx[i+1][j]-gx.Gx[i][j])/pm.dx + (qR*gr.Gr[i][j+1]-qL*gr.Gr[i][j])/pm.dr;
                    //double rhoUex_tmp = (gx.rhoUex[i+1][j]*bRx_wall + gx.rhoUex[i][j]*bLx_wall)/2.0;
                    //double rhoUer_tmp = (rR*gr.rhoUer[i][j+1]*bRr_wall + rL*gr.rhoUer[i][j]*bLr_wall)/(2.0*gc.r[j]);
                    //double Ex_tmp = (gx.Ex[i+1][j] + gx.Ex[i][j])/2.0;
                    //double Er_tmp = (rR*gr.Er[i][j+1] + rL*gr.Er[i][j])/(2.0*gc.r[j]);
                    double JdotE_x = -ph::e0*(gx.rhoUex[i+1][j]*gx.Ex[i+1][j]*bRx_wall + gx.rhoUex[i][j]*gx.Ex[i][j]*bLx_wall)/2.0;
                    double JdotE_r = -ph::e0*(rR*gr.rhoUer[i][j+1]*gr.Er[i][j+1]*bRr_wall + rL*gr.rhoUer[i][j]*gr.Er[i][j]*bLr_wall)/(2.0*gc.r[j]);
                    double RHS2 = JdotE_x + JdotE_r;

                    //double P_ES = (-ph::e0)*(Ex_tmp*rhoUex_tmp + Er_tmp*rhoUer_tmp);
                    double volume = 2.0*M_PI*gc.r[j]*pm.dr*pm.dx;
                    
                    dEdt_sum          = dEdt_sum          + volume*(gc.rhoeps[i][j] - gc.rhoeps_old[i][j])/pm.dt;
                    nabla_G_sum       = nabla_G_sum       + volume*nabla_G;
                    rate_eloss_sum    = rate_eloss_sum    + volume*(-gc.rhoe[i][j]*gc.rate_eloss[i][j]);
                    Joule_sum         = Joule_sum         + volume*RHS2;
                    input_sum         = input_sum         + volume*gc.Pabs[i][j];
                }
            }
        }

        //output energy balance for electron
        //------------------------------------
        std::string char1="results/energy_balance_electron";
        std::string char2=std::to_string(pm.nOut);
        std::string char_csv=".csv";
        std::ofstream outputfile11(char1+char2+char_csv);
        outputfile11 << std::setprecision(std::numeric_limits<double>::max_digits10) << std::scientific;
        outputfile11<< "dE/dt (J/s), loss_to_bnd (J/s), collision_loss (J/s), Joule (J/s), input (J/s), sum (J/s)"<< std::endl;
        outputfile11<< dEdt_sum << "," << -nabla_G_sum << "," << rate_eloss_sum 
            << "," << Joule_sum << "," << input_sum 
            << ","<< -dEdt_sum + (-nabla_G_sum + rate_eloss_sum + Joule_sum + input_sum)<<std::endl;
        outputfile11.close();
        //------------------------------------
    }

    
    //std::cout << "ncount_Te = "<<ncount<<std::endl;

    //Don’t update unless the value is positive.
    //------------------------------------
    for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
        for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
            for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                if(gc.rhoeps[i][j]<0.0){
                    gc.rhoeps[i][j] = fmax(gc.rhoeps_old[i][j],0.0);
                }
            }
        }
    }
    //------------------------------------

    //calculate maximum electron density to avoid zero-division in Te calc.
    //------------------------------------
    double rhoe_max = 0.0;
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

    //calculate electron temperature
    //------------------------------------
    double Te_max = 1000.0*ph::e0/ph::Boltz; //limit maximum
    double ratio = 0.001; //ratio of additional value to maximum
    for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
        for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
            for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                gc.Te[i][j] = fmin(gc.rhoeps[i][j]/(3.0/2.0*(gc.rhoe[i][j] + rhoe_max*ratio)*ph::Boltz + 1e-100),Te_max);
            }
        }
    }
    //------------------------------------

    //set values of ghost cells for transport coefficient
    //=====================================================================================
    //left-wall-BC (chamber)
    //------------------------------------
    for (int j=gc.j_flc_bl[0][0];j<=gc.j_flc_bl[0][1];j++){
        int i = gc.i_flc_bl[0][0];
        gc.Te[i-1][j] = gc.Te[i][j];
    }
    //------------------------------------

    //left-wall-BC (magnet)
    //------------------------------------
    for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[0][0]-1;j++){
        int i = gc.i_flc_bl[1][0];
        gc.Te[i-1][j] = gc.Te[i][j];
    }
    //------------------------------------

    //lower-wall-BC (antenna)
    //------------------------------------
    for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[0][1];i++){ 
        int j = gc.j_flc_bl[0][0];
        gc.Te[i][j-1] = gc.Te[i][j];
    }
    //------------------------------------

    //right-open-BC
    //------------------------------------
    for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[1][1];j++){
        int i = gc.i_flc_bl[1][1];
        gc.Te[i+1][j] = gc.Te[i][j];
    }
    //------------------------------------

    //upper-open-BC
    //------------------------------------
    for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[1][1];i++){ 
        int j = gc.j_flc_bl[0][1];
        gc.Te[i][j+1] = gc.Te[i][j];
    }
    //------------------------------------

    //Caluculate residuals
    //------------------------------------
    calcRes(pm.error_rhoeps, gc.rhoeps,gc.rhoeps_old,gc.i_flc_bl, gc.j_flc_bl, pm.n_bl-1, 1e-100, pm.flag_error);
    calcRes(pm.error_Gx, gx.Gx, gx.Gx_old, gx.i_flx_bl, gx.j_flx_bl, pm.n_bl-1, 1e-100, pm.flag_error);
    calcRes(pm.error_Gr, gr.Gr, gr.Gr_old, gr.i_flr_bl, gr.j_flr_bl, pm.n_bl-1, 1e-100, pm.flag_error);
    //------------------------------------

}
//*****************************************************************
//**                                                             **
//**           void update_rhon()                                **
//**                                                             **
//*****************************************************************
void FluidModule::update_rhon(Params &pm,GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr){


    std::vector<std::vector<double> > aW(pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //係数
    std::vector<std::vector<double> > aE(pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //係数
    std::vector<std::vector<double> > aS(pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //係数
    std::vector<std::vector<double> > aN(pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //係数
    std::vector<std::vector<double> > aP(pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //係数
    std::vector<std::vector<double> > b (pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //ソース

    //time step for metastable neutral
    //------------------------------------
    double dt_n = pm.dt*double(pm.ndt_n);
    //------------------------------------
    
    double vth = sqrt(8.0*ph::Boltz*pm.Tn/(M_PI*pm.massi));
    double sigma = M_PI*(2.0*pm.ri)*(2.0*pm.ri); //Attention：有効半径なので\sigma = \pi*d^2
    //double nu_tmp = sqrt(2.0)*rhon*sigma*vth;
    //double Dn = ph::Boltz*pm.Tn/(pm.massi*nu_tmp);
    double DnN = 2.0/(3.0*sigma)*sqrt(ph::Boltz*pm.Tn/(M_PI*pm.massi));

    //double di_tmp = 2.0*ri;
    //double DnN2 = 3.0/(8.0*di_tmp*di_tmp)*sqrt(ph::Boltz*pm.Tn/(M_PI*pm.massi));
    //std::cout << "DnN = " <<DnN << " DnN2 = " << DnN2 << std::endl;

    //std::cout << "DnN = " << DnN <<std::endl;

    //double Dn = pm.DmN/rhon_ini;

    //過去の値を保存 (ゴーストセルも含む)
    for (int i=0;i<=pm.ni+1;i++){ 
        for (int j=0;j<=pm.nj+1;j++){ 
            gc.rhon_old[i][j] = gc.rhon[i][j];
            gx.rhoUnx_old[i][j] = gx.rhoUnx[i][j];
            gr.rhoUnr_old[i][j] = gr.rhoUnr[i][j];
        }
    }

    //係数作成
    for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
        for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
            for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                double rL = (gc.r[j]+gc.r[j-1])/2.0;
                double rR = (gc.r[j]+gc.r[j+1])/2.0;
                double qL = (gc.r[j]+gc.r[j-1])/2.0/gc.r[j];
                double qR = (gc.r[j]+gc.r[j+1])/2.0/gc.r[j];

                double rhon_Lx = (gc.rhon_old[i][j]+gc.rhon_old[i-1][j])/2.0;
                double rhon_Rx = (gc.rhon_old[i][j]+gc.rhon_old[i+1][j])/2.0;
                double rhon_Lr = (gc.rhon_old[i][j]+gc.rhon_old[i][j-1])/2.0;
                double rhon_Rr = (gc.rhon_old[i][j]+gc.rhon_old[i][j+1])/2.0;

                double Dn_Lx = DnN/(rhon_Lx + 1e-100);
                double Dn_Rx = DnN/(rhon_Rx + 1e-100);
                double Dn_Lr = DnN/(rhon_Lr + 1e-100);
                double Dn_Rr = DnN/(rhon_Rr + 1e-100);

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
                double bRr_wall  = double(j!=gc.j_flc_bl[0][1] || iblock > 2 || gc.x[i] <= pm.width_neutIn + pm.x6); //x1 not inlet
                double bRr_in  = double(j!=gc.j_flc_bl[0][1] || iblock > 2 || gc.x[i] > pm.width_neutIn + pm.x6); //x1 inlet
                bRr_wall = bRr_wall*double(j!=gc.j_flc_bl[3][1] || i<gc.i_flc_bl[1][1]+1 || i>gc.i_flc_bl[3][1]); //x4
                double bRr_open  = (j!=gc.j_flc_bl[4][1] || iblock != 4); //x5
                double bRr = bRr_wall*bRr_in*bRr_open;
                //------------------------------------

                //wall-flux divide by density (for implicit method)
                //------------------------------------
                //double GammaPerN_open_Rr = 0.0;
                //if(bRr_open == 0){
                //    double vth = sqrt(8.0*ph::Boltz*gc.Te[i][j]/(M_PI*pm.masse));
                //    GammaPerN_open_Rr = 0.25*vth;
                //}
                //------------------------------------

                //boundary-flux
                //------------------------------------
                double Gamma_wall_Lx = 0.0;
                if(bLx_wall == 0){
                    Gamma_wall_Lx = -gx.rhoUix_wall[i][j] - gx.rhoUmx_wall[i][j];
                }
                double Gamma_wall_Rx = 0.0;
                if(bRx_wall == 0){
                    Gamma_wall_Rx = -gx.rhoUix_wall[i+1][j] - gx.rhoUmx_wall[i+1][j];
                }
                double Gamma_wall_Lr = 0.0;
                if(bLr_wall == 0){
                    Gamma_wall_Lr = -gr.rhoUir_wall[i][j]   - gr.rhoUmr_wall[i][j];
                }
                double Gamma_wall_Rr = 0.0;
                if(bRr_wall == 0){
                    Gamma_wall_Rr = -gr.rhoUir_wall[i][j+1] - gr.rhoUmr_wall[i][j+1];
                }
                double Gamma_in_Rr = 0.0;
                if(bRr_in== 0){
                    Gamma_in_Rr = -gr.rhoUir_wall[i][j+1] - gr.rhoUmr_wall[i][j+1] - pm.fn_In;
                    //std::cout << Gamma_in_Rr << ","<<pm.fn_In << std::endl;
                }
                double Gamma_open_Rr = 0.0;
                if(bRr_open== 0){
                    Gamma_open_Rr = fmax(rL/rR*gr.rhoUnr_old[i][j],0.0);
                }
                //------------------------------------
   
                aP[i][j] = (Dn_Rx*bRx + Dn_Lx*bLx)*pm.dr*pm.dr 
                    + (qR*Dn_Rr*bRr + qL*Dn_Lr*bLr)*pm.dx*pm.dx 
                    //+    0.25*vth*pm.dr*pm.dx*pm.dr/Dn*nbLx //陰的指定
                    //+    0.25*vth*pm.dr*pm.dx*pm.dr/Dn*nbRx //陰的指定
                    //+ qL*0.25*vth*pm.dx*pm.dx*pm.dr/Dn*nbLr //陰的指定
                    //+ qR*GammaPerN_open_Rr*pm.dx*pm.dx*pm.dr //陰的指定
                    + pm.dx*pm.dx*pm.dr*pm.dr/dt_n; //時間項;
                aE[i][j] = Dn_Rx*pm.dr*pm.dr*bRx;
                aW[i][j] = Dn_Lx*pm.dr*pm.dr*bLx;
                aN[i][j] = qR*Dn_Rr*pm.dx*pm.dx*bRr;
                aS[i][j] = qL*Dn_Lr*pm.dx*pm.dx*bLr;

                //double rhoUnx_Rx = gx.rhoUnx_old[i][j];
                //double rhoUnx_Rr = gr.rhoUnr_old[i][j]*rL/rR;

                //metaからのstepwise-excitationとsuper-elasticのde-excitationによる生成
                double G = (gc.nu_excStep[i][j] + gc.nu_super[i][j])*gc.rhoe[i][j]; 
                
                //metaに行ってしまう分と「direct」-ionizationによる分
                double L = (gc.nu_excMeta[i][j] + 0.5*gc.nu_exc[i][j] + gc.nu_ionz[i][j])*gc.rhoe[i][j];

                b[i][j] = gc.rhon_old[i][j]/dt_n*pm.dx*pm.dx*pm.dr*pm.dr + (G - L)*pm.dx*pm.dx*pm.dr*pm.dr
                    + (
                        Gamma_wall_Lx //左壁境界1,2 //陽的指定
                      )*pm.dr*pm.dx*pm.dr
                    - (
                        Gamma_wall_Rx //右壁境界
                      )*pm.dr*pm.dx*pm.dr
                    + qL*(
                        Gamma_wall_Lr //下壁境界
                      )*pm.dx*pm.dx*pm.dr
                    - qR*(
                        Gamma_wall_Rr //上壁境界
                        + Gamma_open_Rr //上壁境界
                        + Gamma_in_Rr
                      )*pm.dx*pm.dx*pm.dr;
            }
        }
    }

    if(pm.flag_iter_phi==0){
        double alpha = 1.9;
        solver_SOR(aP,aE,aW,aN,aS,b,gc.i_flc_bl,gc.j_flc_bl,pm.n_bl-1,alpha,pm.maxITR_SOR_phi,pm.error_cnv_SOR_phi,0,gc.rhon);
    }else{
        solver_SMG(aP,aE,aW,aN,aS,b,gc.i_flc_bl,gc.j_flc_bl,pm.n_bl-1,pm.maxITR_SOR_phi,pm.error_cnv_SOR_phi,0,gc.rhon);
    }

    //rhoUnx更新
    for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
        for (int i=gx.i_flx_bl[iblock][0];i<=gx.i_flx_bl[iblock][1];i++){ 
            for (int j=gx.j_flx_bl[iblock][0];j<=gx.j_flx_bl[iblock][1];j++){
                double rhon_Lx = (gc.rhon_old[i][j]+gc.rhon_old[i-1][j])/2.0;
                double Dn_Lx = DnN/(rhon_Lx + 1e-100);

                gx.rhoUnx[i][j] = -Dn_Lx*(gc.rhon[i][j] - gc.rhon[i-1][j])/pm.dx;
            }
        }
    }

    //rhoUnr更新
    for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
        for (int i=gr.i_flr_bl[iblock][0];i<=gr.i_flr_bl[iblock][1];i++){ 
            for (int j=gr.j_flr_bl[iblock][0];j<=gr.j_flr_bl[iblock][1];j++){
                double rhon_Lr = (gc.rhon_old[i][j]+gc.rhon_old[i][j-1])/2.0;
                double Dn_Lr = DnN/(rhon_Lr + 1e-100);

                gr.rhoUnr[i][j] = -Dn_Lr*(gc.rhon[i][j] - gc.rhon[i][j-1])/pm.dr;
            }
        }   
    }


    //左 壁 z0
    for (int j=gc.j_flc_bl[0][0];j<=gc.j_flc_bl[0][1];j++){
        int i=gc.i_flc_bl[0][0];
        //gx.rhoUnx[i][j] = -0.25*gc.rhon[i][j]*vth;
        //gx.rhoUnx[i][j] = -gc.rhoi[i][j]*gx.Uix[i][j]; //Recombination
        double Gamma_wall_Lx = -gx.rhoUix_wall[i][j]   - gx.rhoUmx_wall[i][j];
        gx.rhoUnx[i][j] = Gamma_wall_Lx;
        gx.rhoUnx_wall[i][j] = gx.rhoUnx[i][j];
    }
    //左 壁 z1
    for (int j=gc.j_flc_bl[2][0];j<=gc.j_flc_bl[0][0]-1;j++){
        int i=gc.i_flc_bl[2][0];
        //gx.rhoUnx[i][j] = -0.25*gc.rhon[i][j]*vth;
        //gx.rhoUnx[i][j] = -gc.rhoi[i][j]*gx.Uix[i][j]; //Recombination
        double Gamma_wall_Lx = -gx.rhoUix_wall[i][j]   - gx.rhoUmx_wall[i][j];
        gx.rhoUnx[i][j] = Gamma_wall_Lx;
        gx.rhoUnx_wall[i][j] = gx.rhoUnx[i][j];
    }
    //左 壁 z2
    for (int j=gc.j_flc_bl[3][0];j<=gc.j_flc_bl[2][0]-1;j++){
        int i=gc.i_flc_bl[3][0];
        //gx.rhoUnx[i][j] = -0.25*gc.rhon[i][j]*vth;
        //gx.rhoUnx[i][j] = -gc.rhoi[i][j]*gx.Uix[i][j]; //Recombination
        double Gamma_wall_Lx = -gx.rhoUix_wall[i][j]   - gx.rhoUmx_wall[i][j];
        gx.rhoUnx[i][j] = Gamma_wall_Lx;
        gx.rhoUnx_wall[i][j] = gx.rhoUnx[i][j];
    }
    //左 壁 z4
    for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[4][1];j++){
        int i=gc.i_flc_bl[4][0];
        //gx.rhoUnx[i][j] = -0.25*gc.rhon[i][j]*vth;
        //gx.rhoUnx[i][j] = -gc.rhoi[i][j]*gx.Uix[i][j]; //Recombination
        double Gamma_wall_Lx = -gx.rhoUix_wall[i][j]   - gx.rhoUmx_wall[i][j];
        gx.rhoUnx[i][j] = Gamma_wall_Lx;
        gx.rhoUnx_wall[i][j] = gx.rhoUnx[i][j];
    }
    //右 壁 z3
    for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[1][1];j++){
        int i=gc.i_flc_bl[1][1];
        //gx.rhoUnx[i+1][j] = 0.25*gc.rhon[i][j]*vth;
        double Gamma_wall_Rx = -gx.rhoUix_wall[i+1][j] - gx.rhoUmx_wall[i+1][j];
        gx.rhoUnx[i+1][j] = Gamma_wall_Rx;
        gx.rhoUnx_wall[i+1][j] = gx.rhoUnx[i+1][j];
    }
    //右 壁 z5
    for (int j=gc.j_flc_bl[4][0];j<=gc.j_flc_bl[4][1];j++){
        int i=gc.i_flc_bl[4][1];
        //gx.rhoUnx[i+1][j] = 0.25*gc.rhon[i][j]*vth;
        double Gamma_wall_Rx = -gx.rhoUix_wall[i+1][j] - gx.rhoUmx_wall[i+1][j];
        gx.rhoUnx[i+1][j] = Gamma_wall_Rx;
        gx.rhoUnx_wall[i+1][j] = gx.rhoUnx[i+1][j];
    }
    //下 壁 x0
    for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[0][1];i++){
        int j=gc.j_flc_bl[0][0];
        //gr.rhoUnr[i][j] = -0.25*gc.rhon[i][j]*vth;
        //gr.rhoUnr[i][j] = -gc.rhoi[i][j]*gr.Uir[i][j]; //Recombination
        double Gamma_wall_Lr = -gr.rhoUir_wall[i][j]   - gr.rhoUmr_wall[i][j];
        gr.rhoUnr[i][j] = Gamma_wall_Lr;
        gr.rhoUnr_wall[i][j] = gr.rhoUnr[i][j];
    }
    //下 壁 x2
    for (int i=gc.i_flc_bl[2][0];i<=gc.i_flc_bl[2][1];i++){
        int j=gc.j_flc_bl[2][0];
        //gr.rhoUnr[i][j] = -0.25*gc.rhon[i][j]*vth;
        //gr.rhoUnr[i][j] = -gc.rhoi[i][j]*gr.Uir[i][j]; //Recombination
        double Gamma_wall_Lr = -gr.rhoUir_wall[i][j]   - gr.rhoUmr_wall[i][j];
        gr.rhoUnr[i][j] = Gamma_wall_Lr;
        gr.rhoUnr_wall[i][j] = gr.rhoUnr[i][j];
    }
    //上 壁 x1
    for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[1][1];i++){
        int j=gc.j_flc_bl[0][1];
        if(gc.x[i] < pm.width_neutIn + pm.x6){ //流入部
            //gr.rhoUnr[i][j+1] = -pm.fn_In;
            double Gamma_in_Rr = -gr.rhoUir_wall[i][j+1] - gr.rhoUmr_wall[i][j+1] - pm.fn_In;
            gr.rhoUnr[i][j+1] = Gamma_in_Rr;
        }else{ //ただの壁
            //gr.rhoUnr[i][j+1] = 0.25*gc.rhon[i][j]*vth;
            double Gamma_wall_Rr = -gr.rhoUir_wall[i][j+1] - gr.rhoUmr_wall[i][j+1];
            gr.rhoUnr[i][j+1] = Gamma_wall_Rr;
            gr.rhoUnr_wall[i][j+1] = gr.rhoUnr[i][j+1];
        }
    }
    //上 壁 x4
    for (int i=gc.i_flc_bl[1][1]+1;i<=gc.i_flc_bl[3][1];i++){
        int j=gc.j_flc_bl[3][1];
        //gr.rhoUnr[i][j+1] = 0.25*gc.rhon[i][j]*vth;
        //gr.rhoUnr[i][j+1] = -gc.rhoi[i][j]*gr.Uir[i][j+1]; //Recombination
        double Gamma_wall_Rr = -gr.rhoUir_wall[i][j+1] - gr.rhoUmr_wall[i][j+1];
        gr.rhoUnr[i][j+1] = Gamma_wall_Rr;
        gr.rhoUnr_wall[i][j+1] = gr.rhoUnr[i][j+1];
    }
    //上 開放 x5
    for (int i=gc.i_flc_bl[4][0];i<=gc.i_flc_bl[4][1];i++){
        int j=gc.j_flc_bl[4][1];
        //continuous
        double rL = (gc.r[j]+gc.r[j-1])/2.0;
        double rR = (gc.r[j]+gc.r[j+1])/2.0;
        double Gamma_open_Rr = fmax(rL/rR*gr.rhoUnr_old[i][j],0.0);
        gr.rhoUnr[i][j+1] = Gamma_open_Rr;
        gr.rhoUnr_wall[i][j+1] = 0.0;

        //double vth = sqrt(8.0*ph::Boltz*gc.Te[i][j]/(M_PI*pm.masse));
        //double GammaPerN_open_Rr = 0.25*vth;
        //gr.rhoUnr[i][j+1] = gc.rhon[i][j]*GammaPerN_open_Rr ;
        
        //std::cout << Gamma_open_Rr << std::endl;
    }

    //calculation of particle balance
    //------------------------------------
    if(pm.itime%pm.ndiv_fout == 0 || pm.itime == pm.ntime){
        
        //volume integrated values
        //------------------------------------
        double dNndt_sum           = 0.0;
        double nabla_rhoUn_sum     = 0.0;
        double rate_exc_m_to_r_sum = 0.0;
        double rate_super_sum      = 0.0;
        double rate_exc_g_to_m_sum = 0.0;
        double rate_exc_g_to_h_sum = 0.0;
        double rate_d_ionz_sum     = 0.0;
        //------------------------------------
    
        //calculate volume integrated of each term
        //------------------------------------
        for (int iblock=0;iblock<2;iblock++){ 
            for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
                for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                    
                    //radial position at upper and lower cell interface
                    //------------------------------------
                    double rR = (gc.r[j]+gc.r[j+1])/2.0;
                    double rL = (gc.r[j]+gc.r[j-1])/2.0;
                    //------------------------------------

                    //normalized radial position at upper and lower cell interface
                    //------------------------------------
                    double qR = (gc.r[j]+gc.r[j+1])/2.0/gc.r[j];
                    double qL = (gc.r[j]+gc.r[j-1])/2.0/gc.r[j];
                    //------------------------------------

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
                    double bRr_wall  = double(j!=gc.j_flc_bl[0][1] || iblock > 2 || gc.x[i] <= pm.width_neutIn + pm.x6); //x1 not inlet
                    double bRr_in  = double(j!=gc.j_flc_bl[0][1] || iblock > 2 || gc.x[i] > pm.width_neutIn + pm.x6); //x1 inlet
                    bRr_wall = bRr_wall*double(j!=gc.j_flc_bl[3][1] || i<gc.i_flc_bl[1][1]+1 || i>gc.i_flc_bl[3][1]); //x4
                    double bRr_open  = (j!=gc.j_flc_bl[4][1] || iblock != 4); //x5
                    double bRr = bRr_wall*bRr_in*bRr_open;
                    //------------------------------------

                    double Gamma_in_Rr = 0.0;
                    if(bRr_in== 0){
                        Gamma_in_Rr = -gr.rhoUir_wall[i][j+1] - gr.rhoUmr_wall[i][j+1] - pm.fn_In;
                        //std::cout << Gamma_in_Rr << ","<<pm.fn_In << std::endl;
                    }

                    double nabla_rhoUn = (gx.rhoUnx[i+1][j]*bRx_wall - gx.rhoUnx[i][j]*bLx_wall)/pm.dx 
                                       + (qR*gr.rhoUnr[i][j+1]*bRr_wall - qL*gr.rhoUnr[i][j]*bLr_wall)/pm.dr
                                       + (gx.rhoUnx_wall[i+1][j] - gx.rhoUnx_wall[i][j])/pm.dx 
                                       + (qR*gr.rhoUnr_wall[i][j+1] - qL*gr.rhoUnr_wall[i][j])/pm.dr;

                    double volume = 2.0*M_PI*gc.r[j]*pm.dr*pm.dx;
                    
                    dNndt_sum         = dNndt_sum         + volume*(gc.rhon[i][j] - gc.rhon_old[i][j])/pm.dt;
                    nabla_rhoUn_sum   = nabla_rhoUn_sum   + volume*nabla_rhoUn;

                    //generation
                    rate_exc_m_to_r_sum = rate_exc_m_to_r_sum  + volume*gc.rhoe[i][j]*gc.nu_excStep[i][j];
                    rate_super_sum      = rate_super_sum       + volume*gc.rhoe[i][j]*gc.nu_super[i][j];
                    rate_exc_g_to_m_sum = rate_exc_g_to_m_sum  - volume*gc.rhoe[i][j]*gc.nu_excMeta[i][j];
                    rate_exc_g_to_h_sum = rate_exc_g_to_h_sum  - volume*gc.rhoe[i][j]*0.5*gc.nu_exc[i][j];
                    rate_d_ionz_sum     = rate_d_ionz_sum      - volume*gc.rhoe[i][j]*gc.nu_ionz[i][j];
                }
            }
        }

        //output particle balance for metastable
        //------------------------------------
        std::string char1="results/particle_balance_ground";
        std::string char2=std::to_string(pm.nOut);
        std::string char_csv=".csv";
        std::ofstream outputfile10(char1+char2+char_csv);
        outputfile10 << std::setprecision(std::numeric_limits<double>::max_digits10) << std::scientific;
        outputfile10<< "dN/dt (pcl/s), loss_to_bnd (pcl/s), sw-ex (pcl/s), super (pcl/s), ex-to-m (pcl/s), ex-to-h (pcl/s), direct-iz (pcl/s), sum (pcl/s)"<< std::endl;
        outputfile10 << dNndt_sum 
            << "," << -nabla_rhoUn_sum 
            << "," << rate_exc_m_to_r_sum  
            << "," << rate_super_sum 
            << "," << rate_exc_g_to_m_sum 
            << "," << rate_exc_g_to_h_sum
            << "," << rate_d_ionz_sum
            << "," << -dNndt_sum 
                + (
                    - nabla_rhoUn_sum
                    + rate_exc_m_to_r_sum 
                    + rate_super_sum 
                    + rate_exc_g_to_m_sum 
                    + rate_exc_g_to_h_sum 
                    + rate_d_ionz_sum
                )
            << std::endl;
        outputfile10.close();
        //------------------------------------
    }
    //------------------------------------

    //ゴーストセル以外で負になっているところがあれば戻す
    int ncount = 0;
    for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
        for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
            for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                if(gc.rhon[i][j]<0.0){
                    //std::cout << "Warning negative density at (" << i <<","<<j<<") ne = " << gc.rhon[i][j]<< std::endl; 
                    ncount ++;
                    gc.rhon[i][j] = fmax(gc.rhon_old[i][j],0.0);
                }
            }
        }
    }
    if(ncount!=0){
        std::cout << "negative rhon count = " << ncount << std::endl;
    }

    //Caluculate residuals
    //------------------------------------
    calcRes(pm.error_rhon, gc.rhon,gc.rhon_old,gc.i_flc_bl, gc.j_flc_bl, pm.n_bl-1, 1e-100, pm.flag_error);
    //------------------------------------

}

//*****************************************************************
//**                                                             **
//**           void update_rhom()                                **
//**                                                             **
//*****************************************************************
void FluidModule::update_rhom(Params &pm,GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr){


    std::vector<std::vector<double> > aW(pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //係数
    std::vector<std::vector<double> > aE(pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //係数
    std::vector<std::vector<double> > aS(pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //係数
    std::vector<std::vector<double> > aN(pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //係数
    std::vector<std::vector<double> > aP(pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //係数
    std::vector<std::vector<double> > b (pm.ni+2,std::vector<double>(pm.nj+2,0.0)); //ソース

    //time step for metastable neutral
    //------------------------------------
    double dt_m = pm.dt*double(pm.ndt_m);
    //------------------------------------

    //common parameter
    //------------------------------------
    double vth = sqrt(8.0*ph::Boltz*pm.Tn/(M_PI*pm.massi)); //thermal speed
    //------------------------------------

    //Reserve the past values
    //------------------------------------
    for (int i=0;i<=pm.ni+1;i++){ 
        for (int j=0;j<=pm.nj+1;j++){ 
            gc.rhom_old[i][j] = gc.rhom[i][j];
            gx.rhoUmx_old[i][j] = gx.rhoUmx[i][j];
            gr.rhoUmr_old[i][j] = gr.rhoUmr[i][j];
        }
    }
    //------------------------------------

    //Make coefficient for continuity equation for 5 point stencil
    //=================================================================================
    for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
        for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
            for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                //radial position at upper and lower cell interface
                //------------------------------------
                double rL = (gc.r[j]+gc.r[j-1])/2.0;
                double rR = (gc.r[j]+gc.r[j+1])/2.0;
                //------------------------------------
                
                //normalized radial position at upper and lower cell interface
                //------------------------------------
                double qL = (gc.r[j]+gc.r[j-1])/2.0/gc.r[j];
                double qR = (gc.r[j]+gc.r[j+1])/2.0/gc.r[j];
                //------------------------------------

                //diffusion coefficient at cell interface
                //------------------------------------
                double rhon_Lx = (gc.rhon[i][j]+gc.rhon[i-1][j])/2.0;
                double rhon_Rx = (gc.rhon[i][j]+gc.rhon[i+1][j])/2.0;
                double rhon_Lr = (gc.rhon[i][j]+gc.rhon[i][j-1])/2.0;
                double rhon_Rr = (gc.rhon[i][j]+gc.rhon[i][j+1])/2.0;

                double Dm_Lx = pm.DmN/(rhon_Lx + 1e-100);
                double Dm_Rx = pm.DmN/(rhon_Rx + 1e-100);
                double Dm_Lr = pm.DmN/(rhon_Lr + 1e-100);
                double Dm_Rr = pm.DmN/(rhon_Rr + 1e-100);
                //------------------------------------

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
                double bRr_wall  = double(j!=gc.j_flc_bl[0][1] || iblock > 2); //x1 not inlet
                bRr_wall = bRr_wall*double(j!=gc.j_flc_bl[3][1] || i<gc.i_flc_bl[1][1]+1 || i>gc.i_flc_bl[3][1]); //x4
                double bRr_open  = (j!=gc.j_flc_bl[4][1] || iblock != 4); //x5
                double bRr = bRr_wall*bRr_open;
                //------------------------------------

                //BC - normal flux component to the wall or open boundaries
                //*****************************************************
                //wall-flux divide by density (for implicit method)
                //------------------------------------
                double GammaPerNeps_wall_Lx = 0.0;
                double GammaPerNeps_wall_Rx = 0.0;
                double GammaPerNeps_wall_Lr = 0.0;
                double GammaPerNeps_wall_Rr = 0.0;
                if(bLx_wall == 0){
                    GammaPerNeps_wall_Lx = -0.25*vth;
                }
                if(bRx_wall == 0){
                    GammaPerNeps_wall_Rx =  0.25*vth;
                }
                if(bLr_wall == 0){
                    GammaPerNeps_wall_Lr = -0.25*vth;
                }
                if(bRr_wall == 0){
                    GammaPerNeps_wall_Rr =  0.25*vth;
                }
                //------------------------------------

                //open-flux (equal to ion)
                //------------------------------------
                double Gamma_open_Rr = 0.0;
                if(bRr_open == 0){
                    Gamma_open_Rr = rL/rR*gr.rhoUmr_old[i][j];
                }
                //*****************************************************


                double um = sqrt(2.0*ph::Boltz*pm.Tn/pm.massi);

                aP[i][j] = (Dm_Rx*bRx + Dm_Lx*bLx)*pm.dr*pm.dr 
                    + (qR*Dm_Rr*bRr + qL*Dm_Lr*bLr)*pm.dx*pm.dx 
                    -    GammaPerNeps_wall_Lx*pm.dr*pm.dx*pm.dr //陰的指定
                    +    GammaPerNeps_wall_Rx*pm.dr*pm.dx*pm.dr //陰的指定
                    - qL*GammaPerNeps_wall_Lr*pm.dx*pm.dx*pm.dr //陰的指定
                    + qR*GammaPerNeps_wall_Rr*pm.dx*pm.dx*pm.dr //陰的指定
                    + pm.dx*pm.dx*pm.dr*pm.dr/dt_m; //時間項;
                aE[i][j] =    Dm_Rx*pm.dr*pm.dr*bRx;
                aW[i][j] =    Dm_Lx*pm.dr*pm.dr*bLx;
                aN[i][j] = qR*Dm_Rr*pm.dx*pm.dx*bRr;
                aS[i][j] = qL*Dm_Lr*pm.dx*pm.dx*bLr;

                //ソース項
                double G = gc.rhoe[i][j]*(gc.nu_excMeta[i][j] + 0.5*gc.nu_exc[i][j]); //生成　直接励起と upperの半分がmetaに最終的に落ちてくる
                double L = gc.rhoe[i][j]*(gc.nu_ionzStep[i][j] + gc.nu_excStep[i][j] + gc.nu_super[i][j]); //損失 (stepwise Ionizationとstepwise Excitation)

                b[i][j] = gc.rhom_old[i][j]/dt_m*pm.dx*pm.dx*pm.dr*pm.dr + (G - L)*pm.dx*pm.dx*pm.dr*pm.dr
                    + (0.0
                          //rhoUmx_Lx*nbLx //左壁境界1,2 //陽的指定
                      )*pm.dr*pm.dx*pm.dr
                    - (0.0
                         //rhoUmx_Rx*nbRx //右壁境界
                      )*pm.dr*pm.dx*pm.dr
                    + qL*(0.0
                          //rhoUmr_Lr*nbLr //下壁境界　 //陽的指定
                      )*pm.dx*pm.dx*pm.dr
                    - qR*(
                        Gamma_open_Rr //上壁境界
                      )*pm.dx*pm.dx*pm.dr;
            }
        }
    }
    //=================================================================================

    if(pm.flag_iter_phi==0){
        double alpha = 1.9;
        solver_SOR(aP,aE,aW,aN,aS,b,gc.i_flc_bl,gc.j_flc_bl,pm.n_bl-1,alpha,pm.maxITR_SOR_phi,pm.error_cnv_SOR_phi,0,gc.rhom);
    }else{
        solver_SMG(aP,aE,aW,aN,aS,b,gc.i_flc_bl,gc.j_flc_bl,pm.n_bl-1,pm.maxITR_SOR_phi,pm.error_cnv_SOR_phi,0,gc.rhom);
    }

    //ゴーストセル以外で負になっているところがあれば戻す
    int ncount = 0;
    for (int iblock=0;iblock<2;iblock++){ 
        for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
            for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                if(gc.rhom[i][j]<0.0){
                    //std::cout << "Warning negative density at (" << i <<","<<j<<") ne = " << gc.rhom[i][j]<< std::endl; 
                    ncount ++;
                    gc.rhom[i][j] = fmax(gc.rhom_old[i][j],0.0);
                }
            }
        }
    }
    if(ncount!=0){
        std::cout << "negative rhom count = " << ncount << std::endl;
    }

    //rhoUmx更新
    for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
        for (int i=gx.i_flx_bl[iblock][0];i<=gx.i_flx_bl[iblock][1];i++){ 
            for (int j=gx.j_flx_bl[iblock][0];j<=gx.j_flx_bl[iblock][1];j++){
                double rhon_Lx = (gc.rhon[i][j]+gc.rhon[i-1][j])/2.0;
                double Dm_Lx = pm.DmN/(rhon_Lx + 1e-100);

                gx.rhoUmx[i][j] = -Dm_Lx*(gc.rhom[i][j] - gc.rhom[i-1][j])/pm.dx;
            }
        }
    }

    //rhoUmr更新
    for (int iblock=0;iblock<pm.n_bl-1;iblock++){ 
        for (int i=gr.i_flr_bl[iblock][0];i<=gr.i_flr_bl[iblock][1];i++){ 
            for (int j=gr.j_flr_bl[iblock][0];j<=gr.j_flr_bl[iblock][1];j++){
                double rhon_Lr = (gc.rhon[i][j]+gc.rhon[i][j-1])/2.0;
                double Dm_Lr = pm.DmN/(rhon_Lr + 1e-100);
                
                gr.rhoUmr[i][j] = -Dm_Lr*(gc.rhom[i][j] - gc.rhom[i][j-1])/pm.dr;
            }
        }   
    }

    //左 壁 z0
    for (int j=gc.j_flc_bl[0][0];j<=gc.j_flc_bl[0][1];j++){
        int i=gc.i_flc_bl[0][0];
        gx.rhoUmx[i][j] = -0.25*gc.rhom[i][j]*vth;
        gx.rhoUmx_wall[i][j] = gx.rhoUmx[i][j];
    }
    //左 壁 z1
    for (int j=gc.j_flc_bl[2][0];j<=gc.j_flc_bl[0][0]-1;j++){
        int i=gc.i_flc_bl[2][0];
        gx.rhoUmx[i][j] = -0.25*gc.rhom[i][j]*vth;
        gx.rhoUmx_wall[i][j] = gx.rhoUmx[i][j];
    }
    //左 壁 z2
    for (int j=gc.j_flc_bl[3][0];j<=gc.j_flc_bl[2][0]-1;j++){
        int i=gc.i_flc_bl[3][0];
        gx.rhoUmx[i][j] = -0.25*gc.rhom[i][j]*vth;
        gx.rhoUmx_wall[i][j] = gx.rhoUmx[i][j];
    }
    //左 壁 z4
    for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[4][1];j++){
        int i=gc.i_flc_bl[4][0];
        gx.rhoUmx[i][j] = -0.25*gc.rhom[i][j]*vth;
        gx.rhoUmx_wall[i][j] = gx.rhoUmx[i][j];
    }
    //右 壁 z3
    for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[1][1];j++){
        int i=gc.i_flc_bl[1][1];
        gx.rhoUmx[i+1][j] = 0.25*gc.rhom[i][j]*vth;
        gx.rhoUmx_wall[i+1][j] = gx.rhoUmx[i+1][j];
    }
    //右 壁 z5
    for (int j=gc.j_flc_bl[4][0];j<=gc.j_flc_bl[4][1];j++){
        int i=gc.i_flc_bl[4][1];
        gx.rhoUmx[i+1][j] = 0.25*gc.rhom[i][j]*vth;
        gx.rhoUmx_wall[i+1][j] = gx.rhoUmx[i+1][j];
    }
    //下 壁 x0
    for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[0][1];i++){
        int j=gc.j_flc_bl[0][0];
        gr.rhoUmr[i][j] = -0.25*gc.rhom[i][j]*vth;
        gr.rhoUmr_wall[i][j] = gr.rhoUmr[i][j];
    }
    //下 壁 x2
    for (int i=gc.i_flc_bl[2][0];i<=gc.i_flc_bl[2][1];i++){
        int j=gc.j_flc_bl[2][0];
        gr.rhoUmr[i][j] = -0.25*gc.rhom[i][j]*vth;
        gr.rhoUmr_wall[i][j] = gr.rhoUmr[i][j];
    }

    //上 壁 x1
    for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[1][1];i++){
        int j=gc.j_flc_bl[0][1];
        gr.rhoUmr[i][j+1] = 0.25*gc.rhom[i][j]*vth;
        gr.rhoUmr_wall[i][j+1] = gr.rhoUmr[i][j+1];
    }
    //上 壁 x4
    for (int i=gc.i_flc_bl[1][1]+1;i<=gc.i_flc_bl[3][1];i++){
        int j=gc.j_flc_bl[3][1];
        gr.rhoUmr[i][j+1] = 0.25*gc.rhom[i][j]*vth;
        gr.rhoUmr_wall[i][j+1] = gr.rhoUmr[i][j+1];
    }
    //上 開放 x5
    for (int i=gc.i_flc_bl[4][0];i<=gc.i_flc_bl[4][1];i++){
        int j=gc.j_flc_bl[4][1];
        //continuous
        double rL = (gc.r[j]+gc.r[j-1])/2.0;
        double rR = (gc.r[j]+gc.r[j+1])/2.0;
        double rhoUmr_Rr = gr.rhoUmr_old[i][j]*rL/rR;
        gr.rhoUmr[i][j+1] = rhoUmr_Rr;
        gr.rhoUmr_wall[i][j+1] = 0.0;
    }

    
    //calculation of particle balance
    //------------------------------------
    if(pm.itime%pm.ndiv_fout == 0 || pm.itime == pm.ntime){
        
        //volume integrated values
        //------------------------------------
        double dNmdt_sum         = 0.0;
        double nabla_rhoUm_sum   = 0.0;
        double rate_exc_to_h_sum = 0.0;
        double rate_super_sum    = 0.0;
        double rate_exc_to_m_sum = 0.0;
        double rate_drop_sum     = 0.0;
        double rate_s_ionz_sum   = 0.0;
        //------------------------------------
    
        //calculate volume integrated of each term
        //------------------------------------
        for (int iblock=0;iblock<2;iblock++){ 
            for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
                for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                    
                    //radial position at upper and lower cell interface
                    //------------------------------------
                    double rR = (gc.r[j]+gc.r[j+1])/2.0;
                    double rL = (gc.r[j]+gc.r[j-1])/2.0;
                    //------------------------------------

                    //normalized radial position at upper and lower cell interface
                    //------------------------------------
                    double qR = (gc.r[j]+gc.r[j+1])/2.0/gc.r[j];
                    double qL = (gc.r[j]+gc.r[j-1])/2.0/gc.r[j];
                    //------------------------------------

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
                    double bRr_wall  = double(j!=gc.j_flc_bl[0][1] || iblock > 2); //x1 not inlet
                    bRr_wall = bRr_wall*double(j!=gc.j_flc_bl[3][1] || i<gc.i_flc_bl[1][1]+1 || i>gc.i_flc_bl[3][1]); //x4
                    double bRr_open  = (j!=gc.j_flc_bl[4][1] || iblock != 4); //x5
                    double bRr = bRr_wall*bRr_open;
                    //------------------------------------

                    double nabla_rhoUm = (gx.rhoUmx[i+1][j]*bRx_wall - gx.rhoUmx[i][j]*bLx_wall)/pm.dx 
                                       + (qR*gr.rhoUmr[i][j+1]*bRr_wall - qL*gr.rhoUmr[i][j]*bLr_wall)/pm.dr
                                       + (gx.rhoUmx_wall[i+1][j] -gx.rhoUmx_wall[i][j])/pm.dx 
                                       + (qR*gr.rhoUmr_wall[i][j+1] - qL*gr.rhoUmr_wall[i][j])/pm.dr;

                    double volume = 2.0*M_PI*gc.r[j]*pm.dr*pm.dx;
                    
                    dNmdt_sum         = dNmdt_sum         + volume*(gc.rhom[i][j] - gc.rhom_old[i][j])/pm.dt;
                    nabla_rhoUm_sum   = nabla_rhoUm_sum   + volume*nabla_rhoUm;
                    rate_s_ionz_sum   = rate_s_ionz_sum   + volume*gc.rhoe[i][j]*gc.nu_ionzStep[i][j];
                    rate_exc_to_h_sum = rate_exc_to_h_sum + volume*gc.rhoe[i][j]*gc.nu_excStep[i][j];
                    rate_super_sum    = rate_super_sum    + volume*gc.rhoe[i][j]*gc.nu_super[i][j];
                    rate_exc_to_m_sum = rate_exc_to_m_sum + volume*gc.rhoe[i][j]*gc.nu_excMeta[i][j];
                    rate_drop_sum     = rate_drop_sum     + volume*gc.rhoe[i][j]*0.5*gc.nu_exc[i][j];
                }
            }
        }

        //output particle balance for metastable
        //------------------------------------
        std::string char1="results/particle_balance_meta";
        std::string char2=std::to_string(pm.nOut);
        std::string char_csv=".csv";
        std::ofstream outputfile10(char1+char2+char_csv);
        outputfile10 << std::setprecision(std::numeric_limits<double>::max_digits10) << std::scientific;
        outputfile10<< "dN/dt (pcl/s), loss_to_bnd (pcl/s), sw-iz (pcl/s), sw-ex (pcl/s), super (pcl/s), direct-ex (pcl/s), drop (pcl/s), sum (pcl/s)"<< std::endl;
        outputfile10<< dNmdt_sum << ","<< -nabla_rhoUm_sum <<  ","<< -rate_s_ionz_sum  
            <<"," << -rate_exc_to_h_sum << "," << -rate_super_sum << "," <<  rate_exc_to_m_sum << "," <<  rate_drop_sum
            << "," << -dNmdt_sum + (-nabla_rhoUm_sum - rate_s_ionz_sum - rate_exc_to_h_sum - rate_super_sum + rate_exc_to_m_sum + rate_drop_sum)<< std::endl;
        outputfile10.close();
        //------------------------------------
    }
    //------------------------------------
    
    //Caluculate residuals
    //------------------------------------
    calcRes(pm.error_rhom, gc.rhom,gc.rhom_old,gc.i_flc_bl, gc.j_flc_bl, pm.n_bl-1, 1e-100, pm.flag_error);
    //------------------------------------
}

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

                //cout << Mache[i][j] << "," <<nu_ano_IAT[i][j] << ","<<0.01*omega_pi*ud/cs*Te[i][j]/Ti << endl;
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
                + double(pm.flag_Bohm)*gc.nu_ano[i][j] 
                + double(pm.flag_Sagdeev)*gc.nu_ano_IAT[i][j];
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
