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

    if(pm.icon_adp_dt == 1){ //アダプティブdt (dtを計算)
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
    if(pm.icon_iter_Ui == 0){
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
    if(pm.icon_iter_Ui == 0){
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
    if(pm.icon_iter_Ui == 0){
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
    if(pm.icon_iter_rhoi==0){
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
        gx.rhoUix_wall[i][j] = gc.rhoi[i][j]*GammaPerN_wall_Lx;
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
    //=====================================================================================
    //ion density
    //------------------------------------
    pm.error_rhoi = 0.0;
    if(pm.icon_error == 0){ //normalized maximum
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
    if(pm.icon_error == 0){ //normalized maximum
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
    if(pm.icon_error == 0){ //normalized maximum
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
    if(pm.icon_error == 0){ //normalized maximum
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
