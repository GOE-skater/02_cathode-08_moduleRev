
#pragma once
#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <time.h>
#include <vector>
#include <sstream>
#include <iomanip>
#include <limits>

#include "constants.hpp"
#include "params.hpp"
#include "arrays.hpp"
#include "miniFuncs.hpp"

using namespace std;

//*****************************************************************
//**                                                             **
//**           class FluidModule                                 **
//**                                                             **
//*****************************************************************
class OutputModule{
    private:

    public:
        void output_phase(Params &pm, GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr);
        void output(Params &pm, GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr, BolsigVec &bo);
};

//*****************************************************************
//**                                                             **
//**           void output                                       **
//**                                                             **
//*****************************************************************
void OutputModule::output_phase(Params &pm, GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr){
    
    cout<<"output!"<<endl;
    string char_cnv,char1,char2,char_csv;
    int ndeg = 18;

    int ncount = 0;
    double Emag_rep = 0;
    for(int i=1;i<=pm.ni;i++){
        for(int j=1;j<=pm.nj;j++){
            
            double rL = (gc.r[j-1] + gc.r[j])/2.0;
            double rR = (gc.r[j+1] + gc.r[j])/2.0;
            
            complex<double> Ex_tmp = (gx.E1x[i][j]+gx.E1x[i+1][j])/2.0;
            complex<double> Er_tmp = (rR*gr.E1r[i][j+1] + rL*gr.E1r[i][j])/(2.0*gc.r[j]);
            complex<double> Ep_tmp =  gc.E1p[i][j];
            double Emag = norm(Ex_tmp) + norm(Er_tmp) + norm(Ep_tmp);
            Emag_rep = Emag_rep + Emag;
            ncount ++;
        }
    }
    Emag_rep = sqrt(Emag_rep/double(ncount));

    for(int n=0;n<ndeg;n++){
        double deg = 2.0*M_PI/ndeg*n;
        char1="results/result_deg";
        char2=to_string(n);
        char_csv=".csv";

        ofstream outputfile1(char1+char2+char_csv);
        //outputfile1<<"x,rho,Ui,rhoUix,Ue,rhoUex,rhoUex_E,rhoUex_D,E,phi,psi,rho_th,U_th,E_th,phi_th,psi_th,rate_ionize,divi,dive1,dive2,dive3, div_poisson, LHS, RHS,zero" << endl;
        //outputfile1<<"x,rhon,rho,rhoe(test),Uix,rhoUix,Uex,Uey,rhoUex,rhoUey,Te,heat_flux,Bz,E,phi,nu_m,nu_en,nu_ei,nu_ionz,nu_exc,nu_wall,nu_ano,mue_parae,mue_perp,Halle,jd,Id,rate_eloss,rate_ionize,divi,dive,divn,engy_LHS,engy_LHS1,engy_LHS2,engy_LHS3,engy_LHS4,engy_RHS,engy_RHS1,engy_RHS2,engy_RHS3,zero" << endl;
        outputfile1<<"i,j,x,r,Ex(R),Ex(I),|Ex|,arg(Ex),Er(R),Er(I),|Er|,arg(Er),Ep(R),Ep(I),|Ep|,arg(Ep),Er_anly(R),Er_anly(I),|Er_anly|,arg(Er_anly),Er_err(R),Er_err(I),|Er_err|,arg(Er_err),sigma_xx(R),sigma_xx(I),|sigma_xx|,arg(sigma_xx),sigma_rr(R),sigma_rr(I),|sigma_rr|,arg(sigma_rr),sigma_pp(R),sigma_pp(I),|sigma_pp|,arg(sigma_pp),sigma_xr(R),sigma_xr(I),|sigma_xr|,arg(sigma_xr),sigma_xp(R),sigma_xp(I),|sigma_xp|,arg(sigma_xp),sigma_rp(R),sigma_rp(I),|sigma_rp|,arg(sigma_rp),J1x(R),J1x(I),|J1x|,arg(J1x),J1r(R),J1r(I),|J1r|,arg(J1r),J1p(R),J1p(I),|J1p|,arg(J1p),J1x_exc(R),J1x_exc(I),|J1x_exc|,arg(J1x_exc),J1r_exc(R),J1r_exc(I),|J1r_exc|,arg(J1r_exc),J1p_exc(R),J1p_exc(I),|J1p_exc|,arg(J1p_exc),divE,divJ,div(E+J),jdgBnd_Ex,jdgBnd_Er,jdgBnd_Ep,Pabs,Pabs(x),Pabs(r),Pabs(p),rho,Te,nu_m1,Bx,Br,Bmag,Ap,epsr,zero" << endl;
        complex<double> EPS(1e-100,1e-100);

        
        for(int i=1;i<=pm.ni;i++){
            for(int j=1;j<=pm.nj;j++){
                
                double rL = (gc.r[j-1] + gc.r[j])/2.0;
                double rR = (gc.r[j+1] + gc.r[j])/2.0;
                
                complex<double> Ex_tmp = (gx.E1x[i][j]+gx.E1x[i+1][j])/2.0*exp(complex<double>(0,deg))*double(gc.jdgBnd_Ep[i][j]);
                complex<double> Er_tmp = (rR*gr.E1r[i][j+1] + rL*gr.E1r[i][j])/(2.0*gc.r[j])*exp(complex<double>(0,deg))*double(gc.jdgBnd_Ep[i][j]);
                complex<double> Ep_tmp =  gc.E1p[i][j]*exp(complex<double>(0,deg))*double(gc.jdgBnd_Ep[i][j]);
                complex<double> J1x_tmp = (gx.J1x[i][j]+gx.J1x[i+1][j])/2.0*exp(complex<double>(0,deg))*double(gc.jdgBnd_Ep[i][j]);
                complex<double> J1r_tmp = (rR*gr.J1r[i][j+1] + rL*gr.J1r[i][j])/(2.0*gc.r[j])*exp(complex<double>(0,deg))*double(gc.jdgBnd_Ep[i][j]);
                complex<double> J1p_tmp =  gc.J1p[i][j]*exp(complex<double>(0,deg))*double(gc.jdgBnd_Ep[i][j]);

                complex<double> J1x_exc_tmp = (gx.J1x_exc[i][j]+gx.J1x_exc[i+1][j])/2.0*exp(complex<double>(0,deg))*double(gc.jdgBnd_Ep[i][j]);
                complex<double> J1r_exc_tmp = (rR*gr.J1r_exc[i][j+1] + rL*gr.J1r_exc[i][j])/(2.0*gc.r[j])*exp(complex<double>(0,deg))*double(gc.jdgBnd_Ep[i][j]);
                complex<double> J1p_exc_tmp =  gc.J1p_exc[i][j]*exp(complex<double>(0,deg))*double(gc.jdgBnd_Ep[i][j]);
                
                complex<double> divE = (
                    (rR*gr.E1r[i][j+1]-rL*gr.E1r[i][j])/(pm.dr*gc.r[j]) + (gx.E1x[i+1][j]-gx.E1x[i][j])/pm.dx
                )*(pm.dx/Emag_rep)*double(gc.jdgBnd_Ep[i][j]);
                complex<double> divJ = (
                    (rR*gr.J1r[i][j+1]-rL*gr.J1r[i][j])/(pm.dr*gc.r[j]) + (gx.J1x[i+1][j]-gx.J1x[i][j])/pm.dx
                )*(pm.dx/(Emag_rep*pm.omegam*ph::eps0))*double(gc.jdgBnd_Ep[i][j]);
                complex<double> divEJ = (
                    divJ + complex<double>(0.0,pm.omegam*ph::eps0)*divE
                )*(pm.dx/(Emag_rep*pm.omegam*ph::eps0))*double(gc.jdgBnd_Ep[i][j]);

                //cout << "Emag_rep/dx = " << Emag_rep/pm.dx << endl;

                //double aveE = (abs(gr.E1r[i][j+1]) + abs(gr.E1r[i][j]) + abs(gx.E1x[i+1][j]) + abs(gx.E1x[i][j]))/4.0;
                //double aveE = sqrt(real(Ex_tmp*conj(Ex_tmp) + Er_tmp*conj(Er_tmp)));
                //double aveE = 1.0/pm.dx; //5000
                //divE = divE/((aveE+1e-100)/pm.dx);


                double nu_m_tmp = gc.nu_m1[i][j];
                double rho_tmp  = gc.rhoe[i][j];
                double Bx_tmp   = gc.Bx[i][j];
                double Br_tmp   = gc.Br[i][j];

                complex<double> nu_cmp(nu_m_tmp,pm.omegam);
                complex<double> Hallcmpx = ph::e0*Bx_tmp/(pm.masse*nu_cmp + EPS);
                complex<double> Hallcmpr = ph::e0*Br_tmp/(pm.masse*nu_cmp + EPS);
                complex<double> HallcmpMag2 =  Hallcmpx*Hallcmpx + Hallcmpr*Hallcmpr;
                complex<double> coef = rho_tmp*ph::e0*ph::e0/(pm.masse*nu_cmp + EPS)/(1.0 + HallcmpMag2 + EPS);
                complex<double> sigmaexr = coef*(Hallcmpx*Hallcmpr);
                complex<double> sigmaexp = coef*(Hallcmpr);
                complex<double> sigmaexx = coef*(1.0 + Hallcmpx*Hallcmpx);
                complex<double> sigmaerr = coef*(1.0 + Hallcmpr*Hallcmpr);
                complex<double> sigmaerp = coef*(-Hallcmpx);
                complex<double> sigmaerx = coef*(Hallcmpx*Hallcmpr);
                complex<double> sigmaepr = coef*(Hallcmpx);
                complex<double> sigmaepp = coef*(1.0);
                complex<double> sigmaepx = coef*(-Hallcmpr);
                
                double Pabsx_tmp = (gx.Pabsx[i][j] + gx.Pabsx[i+1][j])/2.0;
                double Pabsr_tmp = (gr.Pabsr[i][j] + gr.Pabsr[i][j+1])/2.0;
                
                complex<double> Er_anly(0.0,0.0);
                complex<double> Er_err(0.0,0.0);
                
                if(pm.flag_impTest == 1){
                    double c0 = sqrt(1.0/(ph::eps0*gc.epsr[i][j]*ph::mu0));
                    Er_anly = pm.J1r_exc/gc.r[j]*1e3/log(pm.r2/pm.r1)*exp(complex<double>(0,-pm.omegam/c0*gc.x[i]+deg))*double(gc.jdgBnd_Ep[i][j]);
                }
                Er_err = Er_tmp - Er_anly;
                
                //cout << nu_cmp<< " , "<< 1.0/nu_cmp << endl;

                //complex<double> Jx_tmp = sigmaexr*Er_tmp + sigmaexp*Ep_tmp + sigmaexx*Ex_tmp;
                //complex<double> Jr_tmp = sigmaerr*Er_tmp + sigmaerp*Ep_tmp + sigmaerx*Ex_tmp;
                //complex<double> Jp_tmp = sigmaepr*Er_tmp + sigmaepp*Ep_tmp + sigmaepx*Ex_tmp;

                outputfile1<< i << ","<< j << "," << gc.x[i]<< ","<< gc.r[j]
                    << ","<< real(Ex_tmp) << "," << imag(Ex_tmp) << ","<< abs(Ex_tmp) << ","<< arg(Ex_tmp)
                    << ","<< real(Er_tmp) << "," << imag(Er_tmp) << ","<< abs(Er_tmp) << ","<< arg(Er_tmp)
                    << ","<< real(Ep_tmp) << "," << imag(Ep_tmp) << ","<< abs(Ep_tmp) << ","<< arg(Ep_tmp)
                    << ","<< real(Er_anly) << "," << imag(Er_anly) << ","<< abs(Er_anly) << ","<< arg(Er_anly)
                    << ","<< real(Er_err) << "," << imag(Er_err) << ","<< abs(Er_err) << ","<< arg(Er_err)
                    << ","<< real(sigmaexx) << "," << imag(sigmaexx) << ","<< abs(sigmaexx) << ","<< arg(sigmaexx)
                    << ","<< real(sigmaerr) << "," << imag(sigmaerr) << ","<< abs(sigmaerr) << ","<< arg(sigmaerr)
                    << ","<< real(sigmaepp) << "," << imag(sigmaepp) << ","<< abs(sigmaepp) << ","<< arg(sigmaepp)
                    << ","<< real(sigmaexr) << "," << imag(sigmaexr) << ","<< abs(sigmaexr) << ","<< arg(sigmaexr)
                    << ","<< real(sigmaexp) << "," << imag(sigmaexp) << ","<< abs(sigmaexp) << ","<< arg(sigmaexp)
                    << ","<< real(sigmaerp) << "," << imag(sigmaerp) << ","<< abs(sigmaerp) << ","<< arg(sigmaerp)
                    << ","<< real(J1x_tmp) << "," << imag(J1x_tmp) << ","<< abs(J1x_tmp) << "," << arg(J1x_tmp)
                    << ","<< real(J1r_tmp) << "," << imag(J1r_tmp) << ","<< abs(J1r_tmp) << "," << arg(J1r_tmp)
                    << ","<< real(J1p_tmp) << "," << imag(J1p_tmp) << ","<< abs(J1p_tmp) << "," << arg(J1p_tmp)
                    << ","<< real(J1x_exc_tmp) << "," << imag(J1x_exc_tmp) << ","<< abs(J1x_exc_tmp) << "," << arg(J1x_exc_tmp)
                    << ","<< real(J1r_exc_tmp) << "," << imag(J1r_exc_tmp) << ","<< abs(J1r_exc_tmp) << "," << arg(J1r_exc_tmp)
                    << ","<< real(J1p_exc_tmp) << "," << imag(J1p_exc_tmp) << ","<< abs(J1p_exc_tmp) << "," << arg(J1p_exc_tmp)
                    << ","<< abs(divE)<< ","<< abs(divJ)<< ","<< abs(divEJ)
                    << "," << gx.jdgBnd_Ex[i][j]  << "," << gr.jdgBnd_Er[i][j]  << "," << gc.jdgBnd_Ep[i][j]
                    << ","<< gc.Pabs[i][j]*gc.jdgBnd_flc[i][j]<< ","<< Pabsx_tmp*gc.jdgBnd_flc[i][j]<< ","<< Pabsr_tmp*gc.jdgBnd_flc[i][j]<< ","<< gc.Pabsp[i][j]*gc.jdgBnd_flc[i][j]
                    << ","<< gc.rhoe[i][j]*gc.jdgBnd_flc[i][j] << ","<< gc.Te[i][j]*ph::Boltz/ph::e0*gc.jdgBnd_flc[i][j]<< ","<< gc.nu_m1[i][j]*gc.jdgBnd_flc[i][j]
                    << ","<< gc.Bx[i][j]*gc.jdgBnd_flc[i][j] << ","<< gc.Br[i][j]*gc.jdgBnd_flc[i][j]<< ","<< sqrt(gc.Bx[i][j]*gc.Bx[i][j]+gc.Br[i][j]*gc.Br[i][j])*gc.jdgBnd_flc[i][j]<< ","<< gc.Ap[i][j]
                    << ","<< gc.epsr[i][j]*gc.jdgBnd_Ep[i][j]
                    << "," << 0.0<< endl;
            }
        }

        outputfile1.close();
    }
}

//*****************************************************************
//**                                                             **
//**           void output                                       **
//**                                                             **
//*****************************************************************
void OutputModule::output(Params &pm, GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr, BolsigVec &bo){
    
    //calculate maximum electron density to avoid zero-division
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

    //time step
    //------------------------------------
    double dt_i = pm.dt*double(pm.ndt_i);
    double dt_m = pm.dt*double(pm.ndt_m);
    double dt_n = pm.dt*double(pm.ndt_n);
    //------------------------------------

    string char_cnv,char1,char2,char_csv;

    char1="results/result";
    char2=to_string(pm.nOut);
    char_csv=".csv";

    cout<<"output!"<<endl;

    ofstream outputfile1(char1+char2+char_csv);
    //ofstream outputfile1_1(char1+"_tmp"+char2+char_csv);
    //outputfile1<<"x,rho,Ui,rhoUix,Ue,rhoUex,rhoUex_E,rhoUex_D,E,phi,psi,rho_th,U_th,E_th,phi_th,psi_th,rate_ionize,divi,dive1,dive2,dive3, div_poisson, LHS, RHS,zero" << endl;
    //outputfile1<<"x,rhon,rho,rhoe(test),Uix,rhoUix,Uex,Uey,rhoUex,rhoUey,Te,heat_flux,Bz,E,phi,nu_m,nu_en,nu_ei,nu_ionz,nu_exc,nu_wall,nu_ano,mue_parae,mue_perp,Halle,jd,Id,rate_eloss,rate_ionize,divi,dive,divn,engy_LHS,engy_LHS1,engy_LHS2,engy_LHS3,engy_LHS4,engy_RHS,engy_RHS1,engy_RHS2,engy_RHS3,zero" << endl;
    outputfile1<<"i,j,x,r,rhoi,rhoe,Uix,Uir,Uip,Ui_mag,rhoUix,rhoUir,rhoUip,Te,Te_eff,rhoeps,Gx,Gr,phi,Ex,Er,rhoUex,rhoUer,rhoUep,Uex,Uer,Uep,Ue_mag,rate_ionize,nu_m,nu_m1,nu_ionz,nu_ionzStep,nu_excStep,nu_exc,nu_excMeta,nui_m,mu_para,mu_perp,mu_x,mu_r,mu_c,Ap,Bx,Br,Bmag,Omega,Omegax,Omegar,Pabs,rate_eloss,Loss,P_ES,Pnet,nabla_rhoUe,nabla_rhoUa,nabla_rhoUi,nabla_rhoUm,rate_rhom,rateG_rhom,rateL_rhom,nabla_G,nabla_Ga,rhom,rhoUmx,rhoUmr,Umx,Umr,Debye,phi_analytical,|Ex_mw|,|Er_mw|,|Ep_mw|,rate_ionize_direct,rate_ionize_stepwise,rhon,d|B|/dx,d|B|/dr,rate_eloss_n,rate_eloss_m,distECR,nu_eff,Jx_exc,Jr_exc,Jp_exc,scx,scr,rhoq_sc,nUex,nUer,nabla_nUe,rhoUnx,rhoUnr,nabla_rhoUn,rate_rhon,nu_super,deltaB,nu_ano,nu_ano_IAT,jdgBnd_flc,Mache,zero" << endl;
    //outputfile1_1 << "i,j,x,r,rhom,rhoUmx,rhoUmr,nabla_rhoUm,rhon,rhoUnx,rhoUnr,nabla_rhoUn,zero" << endl;
    
    for(int i=1;i<=pm.ni;i++){
        for(int j=1;j<=pm.nj;j++){
            double rL = (gc.r[j]+gc.r[j-1])/2.0;
            double rR = (gc.r[j]+gc.r[j+1])/2.0;
            double qR = (gc.r[j]+gc.r[j+1])/2.0/gc.r[j];
            double qL = (gc.r[j]+gc.r[j-1])/2.0/gc.r[j];

            //BC - set flag for adjacent to the wall, open, or central axis (adjacent = 0, others = 1)
            //------------------------------------
            //left
            //double bLx_wall  = double(i!=gc.i_flc_bl[0][0] || iblock != 0); //z0
            double bLx_wall  = double(i!=gc.i_flc_bl[0][0]  || j< gc.j_flc_bl[0][0]); //z0
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
            //double bRr_wall  = double(j!=gc.j_flc_bl[0][1] || iblock > 2); //x1
            double bRr_wall  = double(j!=gc.j_flc_bl[0][1] || i>gc.i_flc_bl[1][1]); //x1
            bRr_wall = bRr_wall*double(j!=gc.j_flc_bl[3][1] || i<gc.i_flc_bl[1][1]+1 || i>gc.i_flc_bl[3][1]); //x4
            //w/o x5
            double bRr = bRr_wall;
            //------------------------------------


            double Bmag = sqrt(gc.Bx[i][j]*gc.Bx[i][j] + gc.Br[i][j]*gc.Br[i][j]);
            double absB_Lx = sqrt(gc.Bx[i-1][j]*gc.Bx[i-1][j] + gc.Br[i-1][j]*gc.Br[i-1][j]);
            double absB_Rx = sqrt(gc.Bx[i+1][j]*gc.Bx[i+1][j] + gc.Br[i+1][j]*gc.Br[i+1][j]);
            double absB_Lr = sqrt(gc.Bx[i][j-1]*gc.Bx[i][j-1] + gc.Br[i][j-1]*gc.Br[i][j-1]);
            double absB_Rr = sqrt(gc.Bx[i][j+1]*gc.Bx[i][j+1] + gc.Br[i][j+1]*gc.Br[i][j+1]);
            double dBdx = (absB_Rx-absB_Lx)/(2.0*pm.dx);
            double dBdr = (absB_Rr-absB_Lr)/(2.0*pm.dr);

            //double nu_ano = ph::e0*Bmag/pm.masse/16.0;
            double Omega = ph::e0*Bmag/(pm.masse*(gc.nu_m[i][j])+1e-100);
            double Omega_max = 100;
            double ratio = 1.0;//fmin(Omega_max/Omega,1.0);
            Omega = Omega*ratio;
            double Omega_x = ph::e0*gc.Bx[i][j]/(pm.masse*gc.nu_m[i][j]+1e-100)*ratio;
            double Omega_r = ph::e0*gc.Br[i][j]/(pm.masse*gc.nu_m[i][j]+1e-100)*ratio;

            double nablaB = sqrt(pow((gc.Bx[i+1][j] - gc.Bx[i-1][j])/(2.0*pm.dx),2) 
                + pow((gc.Br[i+1][j] - gc.Br[i-1][j])/(2.0*pm.dr),2));

            double deltaB = Bmag/(nablaB+1e-100);


            double Uix_tmp = (gx.Uix[i][j] + gx.Uix[i+1][j])/2.0;
            double Uir_tmp = (rL*gr.Uir[i][j] + rR*gr.Uir[i][j+1])/(2.0*gc.r[j]);
            double Uip_tmp = gc.Uip[i][j];

            double rhoUix_tmp = (gx.rhoUix[i][j] + gx.rhoUix[i+1][j])/2.0;
            double rhoUir_tmp = (rL*gr.rhoUir[i][j] + rR*gr.rhoUir[i][j+1])/(2.0*gc.r[j]);
            double rhoUip_tmp = gc.rhoi[i][j]*gc.Uip[i][j];

            double Ui_mag = sqrt(Uix_tmp*Uix_tmp + Uir_tmp*Uir_tmp);
            
            double rhoUex_tmp = (gx.rhoUex[i][j] + gx.rhoUex[i+1][j])/2.0;
            double rhoUer_tmp = (rL*gr.rhoUer[i][j] + rR*gr.rhoUer[i][j+1])/(2.0*gc.r[j]);
            double rhoUep_tmp = rhoUer_tmp*Omega_x - rhoUex_tmp*Omega_r;
            double Uex_tmp = rhoUex_tmp/(gc.rhoe[i][j]+ratio_rhoe*rhoe_max);
            double Uer_tmp = rhoUer_tmp/(gc.rhoe[i][j]+ratio_rhoe*rhoe_max);
            double Uep_tmp = rhoUep_tmp/(gc.rhoe[i][j]+ratio_rhoe*rhoe_max);

            double Ue_mag = sqrt(Uex_tmp*Uex_tmp + Uer_tmp*Uer_tmp);

            double Te_eff = gc.Te[i][j] + 1.0/(3.0*ph::Boltz)*pm.masse*Ue_mag*Ue_mag;

            double rhoUmx_tmp = (gx.rhoUmx[i][j] + gx.rhoUmx[i+1][j])/2.0;
            double rhoUmr_tmp = (rL*gr.rhoUmr[i][j] + rR*gr.rhoUmr[i][j+1])/(2.0*gc.r[j]);

            double rhoUnx_tmp = (gx.rhoUnx[i][j] + gx.rhoUnx[i+1][j])/2.0;
            double rhoUnr_tmp = (rL*gr.rhoUnr[i][j] + rR*gr.rhoUnr[i][j+1])/(2.0*gc.r[j]);

            double Ex_tmp = (gx.Ex[i][j] + gx.Ex[i+1][j])/2.0;
            double Er_tmp = (rL*gr.Er[i][j] + rR*gr.Er[i][j+1])/(2.0*gc.r[j]);

            double rhoi_Lx = (gc.rhoi[i][j] + gc.rhoi[i-1][j])/2.0;
            double rhoi_Rx = (gc.rhoi[i][j] + gc.rhoi[i+1][j])/2.0;
            double rhoi_Lr = (gc.rhoi[i][j] + gc.rhoi[i][j-1])/(2.0);
            double rhoi_Rr = (gc.rhoi[i][j] + gc.rhoi[i][j+1])/(2.0);

            //if(j==1) cout << gr.rhoUer[i][j-1] << endl;

            //double nabla_rhoUe = (gx.rhoUex[i+1][j]-gx.rhoUex[i][j])/pm.dx + (qR*gr.rhoUer[i][j+1]-qL*gr.rhoUer[i][j])/pm.dr;
            double nabla_rhoUe = (gx.rhoUex[i+1][j] - gx.rhoUex[i][j])/pm.dx 
                + (qR*gr.rhoUer[i][j+1] - qL*gr.rhoUer[i][j])/pm.dr;

            //double nabla_rhoUi = (rhoi_Rx*gx.Uix[i+1][j]-rhoi_Lx*gx.Uix[i][j])/pm.dx + (qR*rhoi_Rr*gr.Uir[i][j+1]-qL*rhoi_Lr*gr.Uir[i][j])/pm.dr;
            double nabla_rhoUi = (gx.rhoUix[i+1][j]*bRx_wall - gx.rhoUix[i][j]*bLx_wall)/pm.dx 
                                       + (qR*gr.rhoUir[i][j+1]*bRr_wall -qL*gr.rhoUir[i][j]*bLr_wall)/pm.dr
                                       + (gx.rhoUix_wall[i+1][j] - gx.rhoUix_wall[i][j])/pm.dx 
                                       + (qR*gr.rhoUir_wall[i][j+1] - qL*gr.rhoUir_wall[i][j])/pm.dr;

            double nabla_rhoUm = (gx.rhoUmx[i+1][j]-gx.rhoUmx[i][j])/pm.dx + (qR*gr.rhoUmr[i][j+1]-qL*gr.rhoUmr[i][j])/pm.dr;
            double nabla_rhoUn = (gx.rhoUnx[i+1][j]-gx.rhoUnx[i][j])/pm.dx + (qR*gr.rhoUnr[i][j+1]-qL*gr.rhoUnr[i][j])/pm.dr;

            double rate_rhon = (gc.nu_excStep[i][j] + gc.nu_super[i][j])*gc.rhoe[i][j] - (gc.nu_excMeta[i][j] + 0.5*gc.nu_exc[i][j] + gc.nu_ionz[i][j])*gc.rhoe[i][j];
            double lapPhi =  -(gx.Ex[i+1][j] - gx.Ex[i][j])/pm.dx - (qR*gr.Er[i][j+1] - qL*gr.Er[i][j])/pm.dr;
            double nabla_nUe = (gx.nUex[i+1][j] - gx.nUex[i][j])/pm.dx + (qR*gr.nUer[i][j+1] - qL*gr.nUer[i][j])/pm.dr;
            double nUex_tmp = (gx.nUex[i][j] + gx.nUex[i+1][j])/2.0;
            double nUer_tmp = (rL*gr.nUer[i][j] + rR*gr.nUer[i][j+1])/(2.0*gc.r[j]);

            double Gm = gc.rhoe[i][j]*(gc.nu_excMeta[i][j] + 0.5*gc.nu_exc[i][j]); //基底からの励起 + 上の準位からのdrop
            double Lm = gc.rhoe[i][j]*(gc.nu_ionzStep[i][j] + gc.nu_excStep[i][j] + gc.nu_super[i][j]);

            double Gn = (gc.nu_excStep[i][j] + gc.nu_super[i][j])*gc.rhoe[i][j]; 
            double Ln = (gc.nu_excMeta[i][j] + 0.5*gc.nu_exc[i][j] + gc.nu_ionz[i][j])*gc.rhoe[i][j];
            
            double nabla_G = (gx.Gx[i+1][j]-gx.Gx[i][j])/pm.dx + (qR*gr.Gr[i][j+1]-qL*gr.Gr[i][j])/pm.dr;
            double RHS2 =  ph::e0*(Ex_tmp*rhoUex_tmp + Er_tmp*rhoUer_tmp);

            double P_ES = -RHS2;

            double Ex_mw_tmp = abs(gx.E1x[i][j] + gx.E1x[i+1][j])/2.0;
            double Er_mw_tmp = abs(rL*gr.E1r[i][j+1] + rR*gr.E1r[i][j])/(2.0*gc.r[j]);
            double Ep_mw_tmp = abs(gc.E1p[i][j]);

            double Jx_exc_tmp = abs(gx.J1x[i][j] + gx.J1x[i+1][j])/2.0;
            double Jr_exc_tmp = abs(rL*gr.J1r[i][j+1] + rR*gr.J1r[i][j])/(2.0*gc.r[j]);
            double Jp_exc_tmp = abs(gc.J1p[i][j]);

            double ut = sqrt(2.0*ph::Boltz*(gc.Te[i][j]+1.0*ph::e0/ph::Boltz)/pm.masse);
            double Mache_tmp = sqrt(pow(Uex_tmp,2) + pow(Uer_tmp,2))/(ut + 100.0);

            //表面電荷の電荷密度換算量
            double rhoq_sc = (gx.scx[i][j]*2.0*M_PI*gc.r[j]*pm.dr + gr.scr[i][j]*2.0*M_PI*rR*pm.dx)/(2.0*M_PI*gc.r[j]*pm.dx*pm.dr); 

            outputfile1<< i << ","<< j << "," << gc.x[i]<< ","<< gc.r[j]
                << "," << gc.rhoi[i][j]*gc.jdgBnd_flc[i][j] << "," << gc.rhoe[i][j] *gc.jdgBnd_flc[i][j]
                << ","<< Uix_tmp*gc.jdgBnd_flc[i][j] << ","<< Uir_tmp*gc.jdgBnd_flc[i][j]<< ","<< Uip_tmp*gc.jdgBnd_flc[i][j]<< "," << Ui_mag*gc.jdgBnd_flc[i][j]
                << ","<< rhoUix_tmp*gc.jdgBnd_flc[i][j] << ","<< rhoUir_tmp*gc.jdgBnd_flc[i][j]<< ","<< rhoUip_tmp*gc.jdgBnd_flc[i][j]
                << ","<< gc.Te[i][j]*ph::Boltz/ph::e0*gc.jdgBnd_flc[i][j]
                << ","<< Te_eff*ph::Boltz/ph::e0*gc.jdgBnd_flc[i][j]
                << ","<< gc.rhoeps[i][j]*gc.jdgBnd_flc[i][j]
                << ","<< gx.Gx[i][j]*gc.jdgBnd_flc[i][j]<< ","<< gr.Gr[i][j]*gc.jdgBnd_flc[i][j]
                << ","<< gc.phi[i][j]*gc.jdgBnd_flc[i][j]<< ","<< gx.Ex[i][j]*gc.jdgBnd_flc[i][j]<< ","<< gr.Er[i][j]*gc.jdgBnd_flc[i][j]
                << ","<< rhoUex_tmp*gc.jdgBnd_flc[i][j]<< ","<< rhoUer_tmp*gc.jdgBnd_flc[i][j]<< ","<< rhoUep_tmp*gc.jdgBnd_flc[i][j]
                << ","<< Uex_tmp*gc.jdgBnd_flc[i][j]<< ","<< Uer_tmp*gc.jdgBnd_flc[i][j]<< ","<< Uep_tmp*gc.jdgBnd_flc[i][j]<< ","<< Ue_mag*gc.jdgBnd_flc[i][j]
                << ","<< gc.rate_ionize[i][j]*gc.jdgBnd_flc[i][j]<< ","<< gc.nu_m[i][j]*gc.jdgBnd_flc[i][j]<< ","<< gc.nu_m1[i][j]*gc.jdgBnd_flc[i][j]
                << ","<< gc.nu_ionz[i][j]*gc.jdgBnd_flc[i][j] << ","<< gc.nu_ionzStep[i][j]*gc.jdgBnd_flc[i][j]<< ","<< gc.nu_excStep[i][j]*gc.jdgBnd_flc[i][j] 
                << ","<< gc.nu_exc[i][j] *gc.jdgBnd_flc[i][j]<< ","<< gc.nu_excMeta[i][j] *gc.jdgBnd_flc[i][j]<< ","<< gc.nui_m[i][j] *gc.jdgBnd_flc[i][j]
                << ","<< gc.mu_para[i][j]*gc.jdgBnd_flc[i][j]<< ","<< gc.mu_perp[i][j]*gc.jdgBnd_flc[i][j]
                << ","<< gc.mu_xx[i][j]*gc.jdgBnd_flc[i][j]<< ","<< gc.mu_rr[i][j]*gc.jdgBnd_flc[i][j]<< ","<< gc.mu_xr[i][j]*gc.jdgBnd_flc[i][j]
                << ","<< gc.Ap[i][j]*gc.jdgBnd_flc[i][j]<< ","<< gc.Bx[i][j]*gc.jdgBnd_flc[i][j]<< ","<< gc.Br[i][j]*gc.jdgBnd_flc[i][j]<< ","<< Bmag*gc.jdgBnd_flc[i][j]
                << ","<< Omega*gc.jdgBnd_flc[i][j]<< ","<< Omega_x*gc.jdgBnd_flc[i][j]<< ","<< Omega_r*gc.jdgBnd_flc[i][j]
                << ","<< gc.Pabs[i][j]*gc.jdgBnd_flc[i][j]
                << ","<< gc.rate_eloss[i][j]*gc.jdgBnd_flc[i][j]
                << ","<< gc.rhoe[i][j]*gc.rate_eloss[i][j]*gc.jdgBnd_flc[i][j]
                << ","<< P_ES*gc.jdgBnd_flc[i][j]
                << ","<< (gc.Pabs[i][j] - gc.rhoe[i][j]*gc.rate_eloss[i][j] + P_ES)*gc.jdgBnd_flc[i][j]
                << ","<< (nabla_rhoUe + (gc.rhoe[i][j] - gc.rhoe_old[i][j])/pm.dt)*gc.jdgBnd_flc[i][j]
                << ","<< (nabla_rhoUe)*gc.jdgBnd_flc[i][j]
                << ","<< (nabla_rhoUi + (gc.rhoi[i][j] - gc.rhoi_old[i][j])/dt_i)*gc.jdgBnd_flc[i][j]
                << ","<< (nabla_rhoUm + (gc.rhom[i][j] - gc.rhom_old[i][j])/dt_m)*gc.jdgBnd_flc[i][j]
                << ","<< (Gm - Lm)*gc.jdgBnd_flc[i][j]
                << ","<< (Gm)*gc.jdgBnd_flc[i][j]
                << ","<< (Lm)*gc.jdgBnd_flc[i][j]
                << ","<< (nabla_G + (gc.rhoeps[i][j] - gc.rhoeps_old[i][j])/pm.dt - P_ES + gc.rhoe[i][j]*gc.rate_eloss[i][j])*gc.jdgBnd_flc[i][j]
                << ","<< (nabla_G  + gc.rhoe[i][j]*gc.rate_eloss[i][j])*gc.jdgBnd_flc[i][j]
                << "," << gc.rhom[i][j]*gc.jdgBnd_flc[i][j]<< "," << rhoUmx_tmp*gc.jdgBnd_flc[i][j]<< "," << rhoUmr_tmp*gc.jdgBnd_flc[i][j]
                << "," << rhoUmx_tmp/(gc.rhom[i][j]+1e-100)*gc.jdgBnd_flc[i][j]<< "," << rhoUmr_tmp/(gc.rhom[i][j]+1e-100)*gc.jdgBnd_flc[i][j]
                << "," << sqrt(ph::eps0*ph::Boltz*gc.Te[i][j]/((gc.rhoe[i][j] + 1)*ph::e0*ph::e0))*gc.jdgBnd_flc[i][j]
                << "," << ph::Boltz*gc.Te[i][j]/ph::e0*(0.5*log(pm.massi/(2.0*M_PI*pm.masse)) + log(fabs(rhoi_Lr*gr.Uir[i][j]/(gr.rhoUer[i][j] + 1e-100))))*gc.jdgBnd_flc[i][j]
                << "," << Ex_mw_tmp*double(gc.jdgBnd_Ep[i][j])<< "," << Er_mw_tmp*double(gc.jdgBnd_Ep[i][j])<< "," << Ep_mw_tmp*double(gc.jdgBnd_Ep[i][j])
                << "," << gc.rhoe[i][j]*gc.nu_ionz[i][j]*gc.jdgBnd_flc[i][j] << ","<<gc.rhoe[i][j]*gc.nu_ionzStep[i][j]*gc.jdgBnd_flc[i][j]
                << "," << gc.rhon[i][j]*gc.jdgBnd_flc[i][j]
                << "," << dBdx*gc.jdgBnd_flc[i][j] << ","<<dBdr*gc.jdgBnd_flc[i][j]
                << "," << gc.rhon[i][j]*func_Te(bo.rate_eloss_n_N,gc.Te[i][j],pm.dTe)*gc.jdgBnd_flc[i][j]<< "," << gc.rhom[i][j]*func_Te(bo.rate_eloss_m_N,gc.Te[i][j],pm.dTe)*gc.jdgBnd_flc[i][j]
                << "," << gc.distECR[i][j]*gc.jdgBnd_flc[i][j]
                << "," << pm.nu_eff/(1.0 + 4.0*pow(gc.distECR[i][j]/pm.deltaECR,2))*gc.jdgBnd_flc[i][j]
                << "," << Jx_exc_tmp<< "," << Jr_exc_tmp<< "," << Jp_exc_tmp
                << "," << gx.scx[i][j]*gc.jdgBnd_flc[i][j]
                << "," << gr.scr[i][j]*gc.jdgBnd_flc[i][j]
                << "," << rhoq_sc*gc.jdgBnd_flc[i][j]
                << ","<< nUex_tmp*gc.jdgBnd_flc[i][j]<< ","<< nUer_tmp*gc.jdgBnd_flc[i][j]
                << ","<< (nabla_nUe + (gc.rhoi[i][j] - gc.rhoe_old[i][j])/pm.dt + ph::eps0/(ph::e0*pm.dt)*lapPhi)*gc.jdgBnd_flc[i][j]
                << "," << rhoUnx_tmp*gc.jdgBnd_flc[i][j]<< "," << rhoUnr_tmp*gc.jdgBnd_flc[i][j]
                << ","<< (nabla_rhoUn + (gc.rhon[i][j] - gc.rhon_old[i][j])/dt_n)*gc.jdgBnd_flc[i][j]
                << ","<< (Gn - Ln)*gc.jdgBnd_flc[i][j]
                << ","<< gc.nu_super[i][j]*gc.jdgBnd_flc[i][j]
                << ","<< deltaB*gc.jdgBnd_flc[i][j]
                << ","<< gc.nu_ano[i][j]*gc.jdgBnd_flc[i][j]
                << ","<< gc.nu_ano_IAT[i][j]*gc.jdgBnd_flc[i][j]
                << ","<< gc.jdgBnd_flc[i][j]
                << ","<< Mache_tmp*gc.jdgBnd_flc[i][j]
                << "," << 0.0<< endl;
        }
    }

    outputfile1.close();
    //outputfile1_1.close();

    char1="results/restart";
    string char_num = + "_ni" + to_string(pm.ni) + "_nj" + to_string(pm.nj)+ "_";
    char2=to_string(pm.nOut);
    char_csv=".csv";

    ofstream outputfile2(char1+char2+char_csv);

    outputfile2<<"i,j,x,r,rhoi,rhoi_old,Uix,Uix_old,Uir,Uir_old,Uip,Uip_old,rhoe,rhoe_old,rhoUex,rhoUex_old,rhoUer,rhoUer_old,Te,Te_old,Gx,Gx_old,Gr,Gr_old,phi,phi_old,rhom,rhom_old,rhoUmx,rhoUmr,rate_ionize,scx,scr,nUex,nUex_old,nUer,nUer_old,rhon,rhon_old,rhoUnx,rhoUnr,rhoeps,rhoeps_old,rhoUix_wall,rhoUir_wall,rhoUex_wall,rhoUer_wall,rhoUmx_wall,rhoUmr_wall,rhoUnx_wall,rhoUnr_wall,Gx_wall,Gr_wall,zero"<< endl;
    outputfile2 << setprecision(numeric_limits<double>::max_digits10) << scientific;
    
    for(int i=0;i<pm.ni+2;i++){
        for(int j=0;j<pm.nj+2;j++){
            outputfile2<< i << ","<< j << "," << gc.x[i]<< ","<< gc.r[j]
                << "," << truncate(gc.rhoi[i][j]) << "," << truncate(gc.rhoi_old[i][j]) 
                << "," << truncate(gx.Uix[i][j])  << "," << truncate(gx.Uix_old[i][j])  
                << "," << truncate(gr.Uir[i][j])<< "," << truncate(gr.Uir_old[i][j])
                << "," << truncate(gc.Uip[i][j])<< "," << truncate(gc.Uip_old[i][j])
                << "," << truncate(gc.rhoe[i][j]) << "," << truncate(gc.rhoe_old[i][j]) 
                << "," << truncate(gx.rhoUex[i][j]) << "," << truncate(gx.rhoUex_old[i][j]) 
                << "," << truncate(gr.rhoUer[i][j])<< "," << truncate(gr.rhoUer_old[i][j])
                << "," << truncate(gc.Te[i][j])*ph::Boltz/ph::e0 << "," << truncate(gc.Te_old[i][j])*ph::Boltz/ph::e0
                << "," << truncate(gx.Gx[i][j]) << "," << truncate(gx.Gx_old[i][j]) 
                << "," << truncate(gr.Gr[i][j])<< "," << truncate(gr.Gr_old[i][j])
                << "," << truncate(gc.phi[i][j])<< "," << truncate(gc.phi_old[i][j])
                << "," << truncate(gc.rhom[i][j])<< "," << truncate(gc.rhom_old[i][j])
                << "," << truncate(gx.rhoUmx[i][j])<< "," << truncate(gr.rhoUmr[i][j])
                << "," << truncate(gc.rate_ionize[i][j])
                << "," << truncate(gx.scx[i][j])
                << "," << truncate(gr.scr[i][j])
                << "," << truncate(gx.nUex[i][j]) << "," << truncate(gx.nUex_old[i][j]) 
                << "," << truncate(gr.nUer[i][j]) << "," << truncate(gr.nUer_old[i][j]) 
                << "," << truncate(gc.rhon[i][j]) << "," << truncate(gc.rhon_old[i][j]) 
                << "," << truncate(gx.rhoUnx[i][j])      << "," << truncate(gr.rhoUnr[i][j]) 
                << "," << truncate(gc.rhoeps[i][j])      << "," << truncate(gc.rhoeps_old[i][j]) 
                << "," << truncate(gx.rhoUix_wall[i][j]) << "," << truncate(gr.rhoUir_wall[i][j]) 
                << "," << truncate(gx.rhoUex_wall[i][j]) << "," << truncate(gr.rhoUer_wall[i][j]) 
                << "," << truncate(gx.rhoUmx_wall[i][j]) << "," << truncate(gr.rhoUmr_wall[i][j]) 
                << "," << truncate(gx.rhoUnx_wall[i][j]) << "," << truncate(gr.rhoUnr_wall[i][j]) 
                << "," << truncate(gx.Gx_wall[i][j])     << "," << truncate(gr.Gr_wall[i][j]) 
                << "," << 0.0<< endl;
                //if(i == gc.i_flc_bl[2][0]) cout << "output, j = "<<j << " , " << truncate(gx.rhoUix_wall[i][j]) << ","<<gx.rhoUix_wall[i][j]<<endl;
        }
    }
    outputfile2.close();

    //File operation
    //===============================================
    //setup for sum of current per wall
    //------------------------------------
    char1="results/current_sum";
    ofstream outputfile3(char1+char2+char_csv);
    outputfile3 << setprecision(numeric_limits<double>::max_digits10) << scientific;

    pm.Ii_Anode = 0.0;
    pm.Ie_Anode = 0.0;
    pm.I_Anode = 0.0;

    //wall current
    //------------------------------------
    //ion
    double Ii_z0_leftInWall = 0.0;
    double Ii_z1_leftAntBaseWall = 0.0;
    double Ii_z2_leftAntTopWall = 0.0;
    double Ii_z3_rightInWall = 0.0;
    double Ii_z4_leftOutWall = 0.0;
    double Ii_z5_rightAnodeWall = 0.0;
    double Ii_x0_botInWall = 0.0;
    double Ii_x1_topInWall = 0.0;
    double Ii_x2_botAntSideWall = 0.0;
    double Ii_x4_topOrfWall = 0.0;
    double Ii_x5_topOpen = 0.0;

    //electron
    double Ie_z0_leftInWall = 0.0;
    double Ie_z1_leftAntBaseWall = 0.0;
    double Ie_z2_leftAntTopWall = 0.0;
    double Ie_z3_rightInWall = 0.0;
    double Ie_z4_leftOutWall = 0.0;
    double Ie_z5_rightAnodeWall = 0.0;
    double Ie_x0_botInWall = 0.0;
    double Ie_x1_topInWall = 0.0;
    double Ie_x2_botAntSideWall = 0.0;
    double Ie_x4_topOrfWall = 0.0;
    double Ie_x5_topOpen = 0.0;

    //metastable
    double Im_z0_leftInWall = 0.0;
    double Im_z1_leftAntBaseWall = 0.0;
    double Im_z2_leftAntTopWall = 0.0;
    double Im_z3_rightInWall = 0.0;
    double Im_z4_leftOutWall = 0.0;
    double Im_z5_rightAnodeWall = 0.0;
    double Im_x0_botInWall = 0.0;
    double Im_x1_topInWall = 0.0;
    double Im_x2_botAntSideWall = 0.0;
    double Im_x4_topOrfWall = 0.0;
    double Im_x5_topOpen = 0.0;

    //grand neutral
    double In_z0_leftInWall = 0.0;
    double In_z1_leftAntBaseWall = 0.0;
    double In_z2_leftAntTopWall = 0.0;
    double In_z3_rightInWall = 0.0;
    double In_z4_leftOutWall = 0.0;
    double In_z5_rightAnodeWall = 0.0;
    double In_x0_botInWall = 0.0;
    double In_x1_topInWall = 0.0;
    double In_x2_botAntSideWall = 0.0;
    double In_x4_topOrfWall = 0.0;
    double In_x5_topOpen = 0.0;
    //------------------------------------

    //wall energy flux
    //------------------------------------
    //ion
    double Wi_z0_leftInWall = 0.0;
    double Wi_z1_leftAntBaseWall = 0.0;
    double Wi_z2_leftAntTopWall = 0.0;
    double Wi_z3_rightInWall = 0.0;
    double Wi_z4_leftOutWall = 0.0;
    double Wi_z5_rightAnodeWall = 0.0;
    double Wi_x0_botInWall = 0.0;
    double Wi_x1_topInWall = 0.0;
    double Wi_x2_botAntSideWall = 0.0;
    double Wi_x4_topOrfWall = 0.0;
    double Wi_x5_topOpen = 0.0;

    //electron
    double We_z0_leftInWall = 0.0;
    double We_z1_leftAntBaseWall = 0.0;
    double We_z2_leftAntTopWall = 0.0;
    double We_z3_rightInWall = 0.0;
    double We_z4_leftOutWall = 0.0;
    double We_z5_rightAnodeWall = 0.0;
    double We_x0_botInWall = 0.0;
    double We_x1_topInWall = 0.0;
    double We_x2_botAntSideWall = 0.0;
    double We_x4_topOrfWall = 0.0;
    double We_x5_topOpen = 0.0;
    //------------------------------------


    //setup for boundary current and electric field
    //------------------------------------
    char1="results/bnd_I_and_E";
    ofstream outputfile4(char1+char2+char_csv);
    outputfile4 << setprecision(numeric_limits<double>::max_digits10) << scientific;
    outputfile4 << "i,j,x,r, Ii (mA), Ie (mA), I (mA), ji (mA/mm^2), je (mA/mm^2), j (mA/mm^2), En, zero"<< endl;
    //------------------------------------

    //setup for boundary velocity for kinetic flux
    //------------------------------------
    //char1="results/bnd_KF_shift";
    //ofstream outputfile6(char1+char2+char_csv);
    //outputfile6 << setprecision(numeric_limits<double>::max_digits10) << scientific;
    //outputfile6 << "i,j,x,r, uten, Uten, uten_cc,uten_cc, utin, Utin, utin_cc, Utin_cc, utmn, vth2, omega_ce, Larmor, LarmorFactor, angle, sin(angle), LarmorFactor*cos(angle),coef,coef_mirror,zero"<< endl;
    //------------------------------------
    
    //setup for sum of current per wall
    //------------------------------------
    char1="results/energy_flux_sum";
    ofstream outputfile7(char1+char2+char_csv);
    outputfile7 << setprecision(numeric_limits<double>::max_digits10) << scientific;
    //------------------------------------
    //===============================================

    //Caluculate wall information
    //===============================================
    //left-wall-BC (z0)
    //------------------------------------
    for (int j=gc.j_flc_bl[0][0];j<=gc.j_flc_bl[0][1];j++){
        int i=gc.i_flc_bl[0][0];
        double x_tmp = (gc.x[i] + gc.x[i-1])/2.0;

        double En = (-1.0)*gx.Ex[i][j];

        double ji = (-1.0)*( ph::e0)*gx.rhoUix_wall[i][j];
        double je = (-1.0)*(-ph::e0)*gx.rhoUex_wall[i][j];

        double Ii = ji*2.0*M_PI*gc.r[j]*pm.dr;
        double Ie = je*2.0*M_PI*gc.r[j]*pm.dr;
        
        outputfile4<< i << ","<< j << "," << x_tmp << "," << gc.r[j] 
            << ","<< Ii*1000 << ","<< Ie*1000<< ","<< (Ii+Ie)*1000 
            << ","<< ji/1000 << ","<< je/1000<< ","<< (ji+je)/1000 
            << ","<< En << ","<< 0 << endl;

        Ii_z0_leftInWall = Ii_z0_leftInWall - ( ph::e0)*gx.rhoUix_wall[i][j]*2.0*M_PI*gc.r[j]*pm.dr;
        Ie_z0_leftInWall = Ie_z0_leftInWall - (-ph::e0)*gx.rhoUex_wall[i][j]*2.0*M_PI*gc.r[j]*pm.dr;
        Im_z0_leftInWall = Im_z0_leftInWall -       gx.rhoUmx_wall[i][j]*2.0*M_PI*gc.r[j]*pm.dr;
        In_z0_leftInWall = In_z0_leftInWall -       gx.rhoUnx_wall[i][j]*2.0*M_PI*gc.r[j]*pm.dr;

        double engyi = 0.5*pm.massi*(gx.Uix[i][j]*gx.Uix[i][j] + gr.Uir[i][j]*gr.Uir[i][j] + gc.Uip[i][j]*gc.Uip[i][j]) + 5.0/2.0*ph::Boltz*pm.Ti;
        Wi_z0_leftInWall = Wi_z0_leftInWall - engyi*gx.rhoUix_wall[i][j]*2.0*M_PI*gc.r[j]*pm.dr;
        We_z0_leftInWall = We_z0_leftInWall - gx.Gx_wall[i][j]*2.0*M_PI*gc.r[j]*pm.dr;
    }

    //left-wall-BC (z1)
    //------------------------------------
    for (int j=gc.j_flc_bl[2][0];j<=gc.j_flc_bl[0][0]-1;j++){
        int i=gc.i_flc_bl[2][0];
        double x_tmp = (gc.x[i] + gc.x[i-1])/2.0;

        double En = (-1.0)*gx.Ex[i][j];

        double ji = (-1.0)*( ph::e0)*gx.rhoUix_wall[i][j];
        double je = (-1.0)*(-ph::e0)*gx.rhoUex_wall[i][j];

        double Ii = ji*2.0*M_PI*gc.r[j]*pm.dr;
        double Ie = je*2.0*M_PI*gc.r[j]*pm.dr;
        
        outputfile4<< i << ","<< j << "," << x_tmp << "," << gc.r[j] 
            << ","<< Ii*1000 << ","<< Ie*1000<< ","<< (Ii+Ie)*1000 
            << ","<< ji/1000 << ","<< je/1000<< ","<< (ji+je)/1000 
            << ","<< En << ","<< 0 << endl;

        Ii_z1_leftAntBaseWall = Ii_z1_leftAntBaseWall - ( ph::e0)*gx.rhoUix_wall[i][j]*2.0*M_PI*gc.r[j]*pm.dr;
        Ie_z1_leftAntBaseWall = Ie_z1_leftAntBaseWall - (-ph::e0)*gx.rhoUex_wall[i][j]*2.0*M_PI*gc.r[j]*pm.dr;
        Im_z1_leftAntBaseWall = Im_z1_leftAntBaseWall -       gx.rhoUmx_wall[i][j]*2.0*M_PI*gc.r[j]*pm.dr;
        In_z1_leftAntBaseWall = In_z1_leftAntBaseWall -       gx.rhoUnx_wall[i][j]*2.0*M_PI*gc.r[j]*pm.dr;

        double engyi = 0.5*pm.massi*(gx.Uix[i][j]*gx.Uix[i][j] + gr.Uir[i][j]*gr.Uir[i][j] + gc.Uip[i][j]*gc.Uip[i][j]) + 5.0/2.0*ph::Boltz*pm.Ti;
        Wi_z1_leftAntBaseWall = Wi_z1_leftAntBaseWall - engyi*gx.rhoUix_wall[i][j]*2.0*M_PI*gc.r[j]*pm.dr;
        We_z1_leftAntBaseWall = We_z1_leftAntBaseWall - gx.Gx_wall[i][j]*2.0*M_PI*gc.r[j]*pm.dr;
    }
    //------------------------------------

    //left-wall-BC (z2)
    //------------------------------------
    for (int j=gc.j_flc_bl[3][0];j<=gc.j_flc_bl[2][0]-1;j++){
        int i=gc.i_flc_bl[3][0];
        double x_tmp = (gc.x[i] + gc.x[i-1])/2.0;

        double En = (-1.0)*gx.Ex[i][j];

        double ji = (-1.0)*( ph::e0)*gx.rhoUix_wall[i][j];
        double je = (-1.0)*(-ph::e0)*gx.rhoUex_wall[i][j];

        double Ii = ji*2.0*M_PI*gc.r[j]*pm.dr;
        double Ie = je*2.0*M_PI*gc.r[j]*pm.dr;
        
        outputfile4<< i << ","<< j << "," << x_tmp << "," << gc.r[j] 
            << ","<< Ii*1000 << ","<< Ie*1000<< ","<< (Ii+Ie)*1000 
            << ","<< ji/1000 << ","<< je/1000<< ","<< (ji+je)/1000 
            << ","<< En << ","<< 0 << endl;

        Ii_z2_leftAntTopWall = Ii_z2_leftAntTopWall - ( ph::e0)*gx.rhoUix_wall[i][j]*2.0*M_PI*gc.r[j]*pm.dr;
        Ie_z2_leftAntTopWall = Ie_z2_leftAntTopWall - (-ph::e0)*gx.rhoUex_wall[i][j]*2.0*M_PI*gc.r[j]*pm.dr;
        Im_z2_leftAntTopWall = Im_z2_leftAntTopWall -       gx.rhoUmx_wall[i][j]*2.0*M_PI*gc.r[j]*pm.dr;
        In_z2_leftAntTopWall = In_z2_leftAntTopWall -       gx.rhoUnx_wall[i][j]*2.0*M_PI*gc.r[j]*pm.dr;

        double engyi = 0.5*pm.massi*(gx.Uix[i][j]*gx.Uix[i][j] + gr.Uir[i][j]*gr.Uir[i][j] + gc.Uip[i][j]*gc.Uip[i][j]) + 5.0/2.0*ph::Boltz*pm.Ti;
        Wi_z2_leftAntTopWall = Wi_z2_leftAntTopWall - engyi*gx.rhoUix_wall[i][j]*2.0*M_PI*gc.r[j]*pm.dr;
        We_z2_leftAntTopWall = We_z2_leftAntTopWall - gx.Gx_wall[i][j]*2.0*M_PI*gc.r[j]*pm.dr;
    }
    //------------------------------------

    //left-wall-BC (z4)
    //------------------------------------
    for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[4][1];j++){
        int i=gc.i_flc_bl[4][0];
        double x_tmp = (gc.x[i] + gc.x[i-1])/2.0;

        double En = (-1.0)*gx.Ex[i][j];

        double ji = (-1.0)*( ph::e0)*gx.rhoUix_wall[i][j];
        double je = (-1.0)*(-ph::e0)*gx.rhoUex_wall[i][j];

        double Ii = ji*2.0*M_PI*gc.r[j]*pm.dr;
        double Ie = je*2.0*M_PI*gc.r[j]*pm.dr;
        
        outputfile4<< i << ","<< j << "," << x_tmp << "," << gc.r[j] 
            << ","<< Ii*1000 << ","<< Ie*1000<< ","<< (Ii+Ie)*1000 
            << ","<< ji/1000 << ","<< je/1000<< ","<< (ji+je)/1000 
            << ","<< En << ","<< 0 << endl;

        Ii_z4_leftOutWall = Ii_z4_leftOutWall - ( ph::e0)*gx.rhoUix_wall[i][j]*2.0*M_PI*gc.r[j]*pm.dr;
        Ie_z4_leftOutWall = Ie_z4_leftOutWall - (-ph::e0)*gx.rhoUex_wall[i][j]*2.0*M_PI*gc.r[j]*pm.dr;
        Im_z4_leftOutWall = Im_z4_leftOutWall -       gx.rhoUmx_wall[i][j]*2.0*M_PI*gc.r[j]*pm.dr;
        In_z4_leftOutWall = In_z4_leftOutWall -       gx.rhoUnx_wall[i][j]*2.0*M_PI*gc.r[j]*pm.dr;

        double engyi = 0.5*pm.massi*(gx.Uix[i][j]*gx.Uix[i][j] + gr.Uir[i][j]*gr.Uir[i][j] + gc.Uip[i][j]*gc.Uip[i][j]) + 5.0/2.0*ph::Boltz*pm.Ti;
        Wi_z4_leftOutWall = Wi_z4_leftOutWall - engyi*gx.rhoUix_wall[i][j]*2.0*M_PI*gc.r[j]*pm.dr;
        We_z4_leftOutWall = We_z4_leftOutWall - gx.Gx_wall[i][j]*2.0*M_PI*gc.r[j]*pm.dr;
    }
    //------------------------------------

    //right-wall-BC (z3)
    //------------------------------------
    for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[1][1];j++){
        int i=gc.i_flc_bl[1][1];

        double x_tmp = (gc.x[i] + gc.x[i+1])/2.0;

        double En = (+1.0)*gx.Ex[i+1][j];
        //cout << gx.Ex[i+1][j] << ","<<Ex_old[i+1][j] << endl;

        double ji = (+1.0)*( ph::e0)*gx.rhoUix_wall[i+1][j];
        double je = (+1.0)*(-ph::e0)*gx.rhoUex_wall[i+1][j];
        double Ii = ji*2.0*M_PI*gc.r[j]*pm.dr;
        double Ie = je*2.0*M_PI*gc.r[j]*pm.dr;

        outputfile4<< i+1 << ","<< j << "," << x_tmp << "," << gc.r[j] 
            << ","<< Ii*1000 << ","<< Ie*1000<< ","<< (Ii+Ie)*1000 
            << ","<< ji/1000 << ","<< je/1000<< ","<< (ji+je)/1000 
            << ","<< En << ","<< 0 << endl;

        Ii_z3_rightInWall = Ii_z3_rightInWall + ( ph::e0)*gx.rhoUix_wall[i+1][j]*2.0*M_PI*gc.r[j]*pm.dr;
        Ie_z3_rightInWall = Ie_z3_rightInWall + (-ph::e0)*gx.rhoUex_wall[i+1][j]*2.0*M_PI*gc.r[j]*pm.dr;
        Im_z3_rightInWall = Im_z3_rightInWall +       gx.rhoUmx_wall[i+1][j]*2.0*M_PI*gc.r[j]*pm.dr;
        In_z3_rightInWall = In_z3_rightInWall +       gx.rhoUnx_wall[i+1][j]*2.0*M_PI*gc.r[j]*pm.dr;

        double engyi = 0.5*pm.massi*(gx.Uix[i+1][j]*gx.Uix[i+1][j] + gr.Uir[i][j]*gr.Uir[i][j] + gc.Uip[i][j]*gc.Uip[i][j]) + 5.0/2.0*ph::Boltz*pm.Ti;
        Wi_z3_rightInWall = Wi_z3_rightInWall + engyi*gx.rhoUix_wall[i+1][j]*2.0*M_PI*gc.r[j]*pm.dr;
        We_z3_rightInWall = We_z3_rightInWall + gx.Gx_wall[i+1][j]*2.0*M_PI*gc.r[j]*pm.dr;
    }
    //------------------------------------

    //right-wall-BC (z5)
    //------------------------------------
    for (int j=gc.j_flc_bl[4][0];j<=gc.j_flc_bl[4][1];j++){
        int i=gc.i_flc_bl[4][1];

        double x_tmp = (gc.x[i] + gc.x[i+1])/2.0;

        double En = (+1.0)*gx.Ex[i+1][j];
        //cout << gx.Ex[i+1][j] << ","<<Ex_old[i+1][j] << endl;

        double ji = (+1.0)*( ph::e0)*gx.rhoUix_wall[i+1][j];
        double je = (+1.0)*(-ph::e0)*gx.rhoUex_wall[i+1][j];
        double Ii = ji*2.0*M_PI*gc.r[j]*pm.dr;
        double Ie = je*2.0*M_PI*gc.r[j]*pm.dr;

        outputfile4<< i+1 << ","<< j << "," << x_tmp << "," << gc.r[j] 
            << ","<< Ii*1000 << ","<< Ie*1000<< ","<< (Ii+Ie)*1000 
            << ","<< ji/1000 << ","<< je/1000<< ","<< (ji+je)/1000 
            << ","<< En << ","<< 0 << endl;

        Ii_z5_rightAnodeWall = Ii_z5_rightAnodeWall + ( ph::e0)*gx.rhoUix_wall[i+1][j]*2.0*M_PI*gc.r[j]*pm.dr;
        Ie_z5_rightAnodeWall = Ie_z5_rightAnodeWall + (-ph::e0)*gx.rhoUex_wall[i+1][j]*2.0*M_PI*gc.r[j]*pm.dr;
        Im_z5_rightAnodeWall = Im_z5_rightAnodeWall +       gx.rhoUmx_wall[i+1][j]*2.0*M_PI*gc.r[j]*pm.dr;
        In_z5_rightAnodeWall = In_z5_rightAnodeWall +       gx.rhoUnx_wall[i+1][j]*2.0*M_PI*gc.r[j]*pm.dr;

        double engyi = 0.5*pm.massi*(gx.Uix[i+1][j]*gx.Uix[i+1][j] + gr.Uir[i][j]*gr.Uir[i][j] + gc.Uip[i][j]*gc.Uip[i][j]) + 5.0/2.0*ph::Boltz*pm.Ti;
        Wi_z5_rightAnodeWall = Wi_z5_rightAnodeWall + engyi*gx.rhoUix_wall[i+1][j]*2.0*M_PI*gc.r[j]*pm.dr;
        We_z5_rightAnodeWall = We_z5_rightAnodeWall + gx.Gx_wall[i+1][j]*2.0*M_PI*gc.r[j]*pm.dr;
    }
    //------------------------------------
    
    //lower-wall-BC (x0)
    //------------------------------------
    for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[0][1];i++){
        int j=gc.j_flc_bl[0][0];

        double r_tmp = (gc.r[j]+gc.r[j-1])/2.0;

        double En = (-1.0)*gr.Er[i][j];

        double ji = (-1.0)*( ph::e0)*gr.rhoUir_wall[i][j];
        double je = (-1.0)*(-ph::e0)*gr.rhoUer_wall[i][j];
        double Ii = ji*2.0*M_PI*r_tmp*pm.dx;
        double Ie = je*2.0*M_PI*r_tmp*pm.dx;
        
        outputfile4<< i << ","<< j << "," << gc.x[i] << "," << r_tmp 
            << ","<< Ii*1000 << ","<< Ie*1000<< ","<< (Ii+Ie)*1000 
            << ","<< ji/1000 << ","<< je/1000<< ","<< (ji+je)/1000 
            << ","<< En << ","<< 0 << endl;

        Ii_x0_botInWall = Ii_x0_botInWall - ( ph::e0)*gr.rhoUir_wall[i][j]*2.0*M_PI*r_tmp*pm.dx;
        Ie_x0_botInWall = Ie_x0_botInWall - (-ph::e0)*gr.rhoUer_wall[i][j]*2.0*M_PI*r_tmp*pm.dx;
        Im_x0_botInWall = Im_x0_botInWall -       gr.rhoUmr_wall[i][j]*2.0*M_PI*r_tmp*pm.dx;
        In_x0_botInWall = In_x0_botInWall -       gr.rhoUnr_wall[i][j]*2.0*M_PI*r_tmp*pm.dx;
        
        double engyi = 0.5*pm.massi*(gx.Uix[i][j]*gx.Uix[i][j] + gr.Uir[i][j]*gr.Uir[i][j] + gc.Uip[i][j]*gc.Uip[i][j]) + 5.0/2.0*ph::Boltz*pm.Ti;
        Wi_x0_botInWall = Wi_x0_botInWall - engyi*gr.rhoUir_wall[i][j]*2.0*M_PI*r_tmp*pm.dx;
        We_x0_botInWall = We_x0_botInWall - gr.Gr_wall[i][j]*2.0*M_PI*r_tmp*pm.dx;
    }
    //------------------------------------

    //lower-wall-BC (x2)
    //------------------------------------
    for (int i=gc.i_flc_bl[2][0];i<=gc.i_flc_bl[2][1];i++){
        int j=gc.j_flc_bl[2][0];

         double r_tmp = (gc.r[j]+gc.r[j-1])/2.0;

        double En = (-1.0)*gr.Er[i][j];

        double ji = (-1.0)*( ph::e0)*gr.rhoUir_wall[i][j];
        double je = (-1.0)*(-ph::e0)*gr.rhoUer_wall[i][j];
        double Ii = ji*2.0*M_PI*r_tmp*pm.dx;
        double Ie = je*2.0*M_PI*r_tmp*pm.dx;
        
        outputfile4<< i << ","<< j << "," << gc.x[i] << "," << r_tmp 
            << ","<< Ii*1000 << ","<< Ie*1000<< ","<< (Ii+Ie)*1000 
            << ","<< ji/1000 << ","<< je/1000<< ","<< (ji+je)/1000 
            << ","<< En << ","<< 0 << endl;

        Ii_x2_botAntSideWall = Ii_x2_botAntSideWall - ( ph::e0)*gr.rhoUir_wall[i][j]*2.0*M_PI*r_tmp*pm.dx;
        Ie_x2_botAntSideWall = Ie_x2_botAntSideWall - (-ph::e0)*gr.rhoUer_wall[i][j]*2.0*M_PI*r_tmp*pm.dx;
        Im_x2_botAntSideWall = Im_x2_botAntSideWall -       gr.rhoUmr_wall[i][j]*2.0*M_PI*r_tmp*pm.dx;
        In_x2_botAntSideWall = In_x2_botAntSideWall -       gr.rhoUnr_wall[i][j]*2.0*M_PI*r_tmp*pm.dx;
        
        double engyi = 0.5*pm.massi*(gx.Uix[i][j]*gx.Uix[i][j] + gr.Uir[i][j]*gr.Uir[i][j] + gc.Uip[i][j]*gc.Uip[i][j]) + 5.0/2.0*ph::Boltz*pm.Ti;
        Wi_x2_botAntSideWall = Wi_x2_botAntSideWall - engyi*gr.rhoUir_wall[i][j]*2.0*M_PI*r_tmp*pm.dx;
        We_x2_botAntSideWall = We_x2_botAntSideWall - gr.Gr_wall[i][j]*2.0*M_PI*r_tmp*pm.dx;
    }
    //------------------------------------

    //upper-wall-BC (x1)
    //------------------------------------
    for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[1][1];i++){
        int j=gc.j_flc_bl[0][1];

        double r_tmp = (gc.r[j]+gc.r[j+1])/2.0;

        double En = (+1.0)*gr.Er[i][j+1];

        double ji = (+1.0)*( ph::e0)*gr.rhoUir_wall[i][j+1];
        double je = (+1.0)*(-ph::e0)*gr.rhoUer_wall[i][j+1];
        double Ii = ji*2.0*M_PI*r_tmp*pm.dx;
        double Ie = je*2.0*M_PI*r_tmp*pm.dx;
        
        outputfile4<< i << ","<< j+1 << "," << gc.x[i] << "," << r_tmp 
            << ","<< Ii*1000 << ","<< Ie*1000<< ","<< (Ii+Ie)*1000 
            << ","<< ji/1000 << ","<< je/1000<< ","<< (ji+je)/1000 
            << ","<< En << ","<< 0 << endl;

        Ii_x1_topInWall = Ii_x1_topInWall + ( ph::e0)*gr.rhoUir_wall[i][j+1]*2.0*M_PI*r_tmp*pm.dx;
        Ie_x1_topInWall = Ie_x1_topInWall + (-ph::e0)*gr.rhoUer_wall[i][j+1]*2.0*M_PI*r_tmp*pm.dx;
        Im_x1_topInWall = Im_x1_topInWall +       gr.rhoUmr_wall[i][j+1]*2.0*M_PI*r_tmp*pm.dx;
        In_x1_topInWall = In_x1_topInWall +       gr.rhoUnr_wall[i][j+1]*2.0*M_PI*r_tmp*pm.dx;

        double engyi = 0.5*pm.massi*(gx.Uix[i][j]*gx.Uix[i][j] + gr.Uir[i][j+1]*gr.Uir[i][j+1] + gc.Uip[i][j]*gc.Uip[i][j]) + 5.0/2.0*ph::Boltz*pm.Ti;
        Wi_x1_topInWall = Wi_x1_topInWall + engyi*gr.rhoUir_wall[i][j+1]*2.0*M_PI*r_tmp*pm.dx;
        We_x1_topInWall = We_x1_topInWall + gr.Gr_wall[i][j+1]*2.0*M_PI*r_tmp*pm.dx;
    }
    //------------------------------------


    //upper-wall-BC (x4)
    //------------------------------------
    for (int i=gc.i_flc_bl[1][1]+1;i<=gc.i_flc_bl[3][1];i++){
        int j=gc.j_flc_bl[3][1];

        double r_tmp = (gc.r[j]+gc.r[j+1])/2.0;

        double En = (+1.0)*gr.Er[i][j+1];

        double ji = (+1.0)*( ph::e0)*gr.rhoUir_wall[i][j+1];
        double je = (+1.0)*(-ph::e0)*gr.rhoUer_wall[i][j+1];
        double Ii = ji*2.0*M_PI*r_tmp*pm.dx;
        double Ie = je*2.0*M_PI*r_tmp*pm.dx;
        
        outputfile4<< i << ","<< j+1 << "," << gc.x[i] << "," << r_tmp 
            << ","<< Ii*1000 << ","<< Ie*1000<< ","<< (Ii+Ie)*1000 
            << ","<< ji/1000 << ","<< je/1000<< ","<< (ji+je)/1000 
            << ","<< En << ","<< 0 << endl;

        Ii_x4_topOrfWall = Ii_x4_topOrfWall + ( ph::e0)*gr.rhoUir_wall[i][j+1]*2.0*M_PI*r_tmp*pm.dx;
        Ie_x4_topOrfWall = Ie_x4_topOrfWall + (-ph::e0)*gr.rhoUer_wall[i][j+1]*2.0*M_PI*r_tmp*pm.dx;
        Im_x4_topOrfWall = Im_x4_topOrfWall +       gr.rhoUmr_wall[i][j+1]*2.0*M_PI*r_tmp*pm.dx;
        In_x4_topOrfWall = In_x4_topOrfWall +       gr.rhoUnr_wall[i][j+1]*2.0*M_PI*r_tmp*pm.dx;

        double engyi = 0.5*pm.massi*(gx.Uix[i][j]*gx.Uix[i][j] + gr.Uir[i][j+1]*gr.Uir[i][j+1] + gc.Uip[i][j]*gc.Uip[i][j]) + 5.0/2.0*ph::Boltz*pm.Ti;
        Wi_x4_topOrfWall = Wi_x4_topOrfWall + engyi*gr.rhoUir_wall[i][j+1]*2.0*M_PI*r_tmp*pm.dx;
        We_x4_topOrfWall = We_x4_topOrfWall + gr.Gr_wall[i][j+1]*2.0*M_PI*r_tmp*pm.dx;
    }


    //upper-open-BC (x5)
    //------------------------------------
    for (int i=gc.i_flc_bl[4][0];i<=gc.i_flc_bl[4][1];i++){
        int j=gc.j_flc_bl[4][1];

        double r_tmp = (gc.r[j]+gc.r[j+1])/2.0;

        double En = (+1.0)*gr.Er[i][j+1];

        double rhoUir_Rr = gc.rhoi[i][j]*gr.Uir[i][j+1];
        double ji = (+1.0)*( ph::e0)*rhoUir_Rr;
        double je = (+1.0)*(-ph::e0)*gr.rhoUer[i][j+1];
        double Ii = ji*2.0*M_PI*r_tmp*pm.dx;
        double Ie = je*2.0*M_PI*r_tmp*pm.dx;

        outputfile4<< i << ","<< j+1 << "," << gc.x[i] << "," << r_tmp 
            << ","<< Ii*1000 << ","<< Ie*1000<< ","<< (Ii+Ie)*1000 
            << ","<< ji/1000 << ","<< je/1000<< ","<< (ji+je)/1000 
            << ","<< En << ","<< 0 << endl;

        Ii_x5_topOpen = Ii_x5_topOpen + ( ph::e0)*rhoUir_Rr*2.0*M_PI*r_tmp*pm.dx;
        Ie_x5_topOpen = Ie_x5_topOpen + (-ph::e0)*gr.rhoUer[i][j+1]*2.0*M_PI*r_tmp*pm.dx;
        Im_x5_topOpen = Im_x5_topOpen +       gr.rhoUmr[i][j+1]*2.0*M_PI*r_tmp*pm.dx;
        In_x5_topOpen = In_x5_topOpen +       gr.rhoUnr[i][j+1]*2.0*M_PI*r_tmp*pm.dx;

        double engyi = 0.5*pm.massi*(gx.Uix[i][j]*gx.Uix[i][j] + gr.Uir[i][j+1]*gr.Uir[i][j+1] + gc.Uip[i][j]*gc.Uip[i][j]) + 5.0/2.0*ph::Boltz*pm.Ti;
        Wi_x5_topOpen = Wi_x5_topOpen + engyi*rhoUir_Rr*2.0*M_PI*r_tmp*pm.dx;
        We_x5_topOpen = We_x5_topOpen + gr.Gr[i][j+1]*2.0*M_PI*r_tmp*pm.dx;
    }
    //------------------------------------

    outputfile4.close();
    //outputfile6.close();

    /*
    //------------------------------------
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
    */

    pm.Ii_Anode = Ii_z5_rightAnodeWall;
    pm.Ie_Anode = Ie_z5_rightAnodeWall;
    pm.I_Anode =  pm.Ii_Anode + pm.Ie_Anode;

    //Caluculate nozzle Current
    //------------------------------------
    pm.Ii_Nozzle = 0.0;
    pm.Ie_Nozzle = 0.0;
    pm.I_Nozzle = 0.0;
    double Im_Nozzle = 0.0;
    double In_Nozzle = 0.0;

    for (int j=gc.j_flc_bl[3][0];j<=gc.j_flc_bl[3][1];j++){
        int i=gc.i_flc_bl[1][1];
        double rhoi_tmp = (gc.rhoi[i][j] + gc.rhoi[i+1][j])/2.0;
        double rhoUix_tmp = rhoi_tmp*gx.Uix[i+1][j];
        
        double delta_Ii = rhoUix_tmp  *2.0*M_PI*gc.r[j]*pm.dr*( ph::e0);
        double delta_Ie = gx.rhoUex[i][j]*2.0*M_PI*gc.r[j]*pm.dr*(-ph::e0);
        double delta_Im = gx.rhoUmx[i][j]*2.0*M_PI*gc.r[j]*pm.dr;
        double delta_In = gx.rhoUnx[i][j]*2.0*M_PI*gc.r[j]*pm.dr;

        pm.Ii_Nozzle = pm.Ii_Nozzle + delta_Ii;
        pm.Ie_Nozzle = pm.Ie_Nozzle + delta_Ie;
        pm.I_Nozzle =  pm.I_Nozzle + delta_Ii + delta_Ie;

        Im_Nozzle = Im_Nozzle + delta_Im;
        In_Nozzle = In_Nozzle + delta_In;
    }
    //------------------------------------

    //output current sum
    //------------------------------------
    double Ii_sum = Ii_z0_leftInWall + Ii_z1_leftAntBaseWall + Ii_z2_leftAntTopWall + Ii_z3_rightInWall + Ii_z4_leftOutWall + Ii_z5_rightAnodeWall 
            + Ii_x0_botInWall+ Ii_x1_topInWall + Ii_x2_botAntSideWall + Ii_x4_topOrfWall + Ii_x5_topOpen;
    double Ie_sum = Ie_z0_leftInWall + Ie_z1_leftAntBaseWall + Ie_z2_leftAntTopWall + Ie_z3_rightInWall + Ie_z4_leftOutWall + Ie_z5_rightAnodeWall 
            + Ie_x0_botInWall+ Ie_x1_topInWall + Ie_x2_botAntSideWall + Ie_x4_topOrfWall + Ie_x5_topOpen;
    double Im_sum = Im_z0_leftInWall + Im_z1_leftAntBaseWall + Im_z2_leftAntTopWall + Im_z3_rightInWall + Im_z4_leftOutWall + Im_z5_rightAnodeWall 
            + Im_x0_botInWall+ Im_x1_topInWall + Im_x2_botAntSideWall + Im_x4_topOrfWall + Im_x5_topOpen;
    double In_sum = In_z0_leftInWall + In_z1_leftAntBaseWall + In_z2_leftAntTopWall + In_z3_rightInWall + In_z4_leftOutWall + In_z5_rightAnodeWall 
            + In_x0_botInWall+ In_x1_topInWall + In_x2_botAntSideWall + In_x4_topOrfWall + In_x5_topOpen;
    
    outputfile3<< "No, wall, Ii (mA), Ie (mA), I (mA), Ndoti (pcl/s), Ndote (pcl/s), Ndotm (pcl/s), Ndotn (pcl/s)"<< endl;
    outputfile3<< "z0 , leftInWall, " << Ii_z0_leftInWall*1000 << "," << Ie_z0_leftInWall*1000 << "," << (Ii_z0_leftInWall + Ie_z0_leftInWall)*1000 << "," << 1.0/ph::e0*Ii_z0_leftInWall << "," << 1.0/(-ph::e0)*Ie_z0_leftInWall << "," << Im_z0_leftInWall << "," << In_z0_leftInWall  << endl;
    outputfile3<< "z1 , leftAntBaseWall, " << Ii_z1_leftAntBaseWall*1000 << "," << Ie_z1_leftAntBaseWall*1000 << "," << (Ii_z1_leftAntBaseWall + Ie_z1_leftAntBaseWall)*1000 << "," << 1.0/ph::e0*Ii_z1_leftAntBaseWall << "," << 1.0/(-ph::e0)*Ie_z1_leftAntBaseWall << "," << Im_z1_leftAntBaseWall << "," << In_z1_leftAntBaseWall  << endl;
    outputfile3<< "z2 , leftAntTopWall, " << Ii_z2_leftAntTopWall*1000 << "," << Ie_z2_leftAntTopWall*1000 << "," << (Ii_z2_leftAntTopWall + Ie_z2_leftAntTopWall)*1000 << "," << 1.0/ph::e0*Ii_z2_leftAntTopWall << "," << 1.0/(-ph::e0)*Ie_z2_leftAntTopWall << "," << Im_z2_leftAntTopWall << "," << In_z2_leftAntTopWall  << endl;
    outputfile3<< "z3 , rightInWall, " << Ii_z3_rightInWall*1000 << "," << Ie_z3_rightInWall*1000 << "," << (Ii_z3_rightInWall + Ie_z3_rightInWall)*1000 << "," << 1.0/ph::e0*Ii_z3_rightInWall << "," << 1.0/(-ph::e0)*Ie_z3_rightInWall << "," << Im_z3_rightInWall << "," << In_z3_rightInWall  << endl;
    outputfile3<< "z4 , leftOutWall, " << Ii_z4_leftOutWall*1000 << "," << Ie_z4_leftOutWall*1000 << "," << (Ii_z4_leftOutWall + Ie_z4_leftOutWall)*1000 << "," << 1.0/ph::e0*Ii_z4_leftOutWall << "," << 1.0/(-ph::e0)*Ie_z4_leftOutWall << "," << Im_z4_leftOutWall << "," << In_z4_leftOutWall  << endl;
    outputfile3<< "z5 , rightAnodeWall, " << Ii_z5_rightAnodeWall*1000 << "," << Ie_z5_rightAnodeWall*1000 << "," << (Ii_z5_rightAnodeWall + Ie_z5_rightAnodeWall)*1000 << "," << 1.0/ph::e0*Ii_z5_rightAnodeWall << "," << 1.0/(-ph::e0)*Ie_z5_rightAnodeWall << "," << Im_z5_rightAnodeWall << "," << In_z5_rightAnodeWall  << endl;
    outputfile3<< "x0 , botInWall, " << Ii_x0_botInWall*1000 << "," << Ie_x0_botInWall*1000 << "," << (Ii_x0_botInWall + Ie_x0_botInWall)*1000 << "," << 1.0/ph::e0*Ii_x0_botInWall << "," << 1.0/(-ph::e0)*Ie_x0_botInWall << "," << Im_x0_botInWall << "," << In_x0_botInWall  << endl;
    outputfile3<< "x1 , topInWall, " << Ii_x1_topInWall*1000 << "," << Ie_x1_topInWall*1000 << "," << (Ii_x1_topInWall + Ie_x1_topInWall)*1000 << "," << 1.0/ph::e0*Ii_x1_topInWall << "," << 1.0/(-ph::e0)*Ie_x1_topInWall << "," << Im_x1_topInWall << "," << In_x1_topInWall  << endl;
    outputfile3<< "x2 , botAntSideWall, " << Ii_x2_botAntSideWall*1000 << "," << Ie_x2_botAntSideWall*1000 << "," << (Ii_x2_botAntSideWall + Ie_x2_botAntSideWall)*1000 << "," << 1.0/ph::e0*Ii_x2_botAntSideWall << "," << 1.0/(-ph::e0)*Ie_x2_botAntSideWall << "," << Im_x2_botAntSideWall << "," << In_x2_botAntSideWall  << endl;
    outputfile3<< "x4 , topOrfWall, " << Ii_x4_topOrfWall*1000 << "," << Ie_x4_topOrfWall*1000 << "," << (Ii_x4_topOrfWall + Ie_x4_topOrfWall)*1000 << "," << 1.0/ph::e0*Ii_x4_topOrfWall << "," << 1.0/(-ph::e0)*Ie_x4_topOrfWall << "," << Im_x4_topOrfWall << "," << In_x4_topOrfWall  << endl;
    outputfile3<< "x5 , topOpen, " << Ii_x5_topOpen*1000 << "," << Ie_x5_topOpen*1000 << "," << (Ii_x5_topOpen + Ie_x5_topOpen)*1000 << "," << 1.0/ph::e0*Ii_x5_topOpen << "," << 1.0/(-ph::e0)*Ie_x5_topOpen << "," << Im_x5_topOpen << "," << In_x5_topOpen  << endl;
    outputfile3<< "*  , Nozzle , " << pm.Ii_Nozzle*1000 << "," << pm.Ie_Nozzle*1000 << "," << (pm.Ii_Nozzle + pm.Ie_Nozzle)*1000 << "," << 1.0/ph::e0*pm.Ii_Nozzle << "," << 1.0/(-ph::e0)*pm.Ie_Nozzle << "," << Im_Nozzle << "," << In_Nozzle  << endl;
    outputfile3<< "*  , sum, " << Ii_sum*1000 << "," << Ie_sum*1000 << "," << (Ii_sum + Ie_sum)*1000 << "," << 1.0/ph::e0*Ii_sum << "," << 1.0/(-ph::e0)*Ie_sum << "," << Im_sum << "," << In_sum  << endl;
    outputfile3.close();
    //------------------------------------

    //output energy sum
    //------------------------------------
    double Wi_sum = Wi_z0_leftInWall + Wi_z1_leftAntBaseWall + Wi_z3_rightInWall + Wi_z4_leftOutWall + Wi_z5_rightAnodeWall 
            + Wi_x0_botInWall+ Wi_x1_topInWall + Wi_x2_botAntSideWall + Wi_x4_topOrfWall + Wi_x5_topOpen;
    double We_sum = We_z0_leftInWall + We_z1_leftAntBaseWall + We_z3_rightInWall + We_z4_leftOutWall + We_z5_rightAnodeWall 
            + We_x0_botInWall+ We_x1_topInWall + We_x2_botAntSideWall + We_x4_topOrfWall + We_x5_topOpen;

    outputfile7<< "No, wall, Wi (J/s), We (J/s), W (J/s)"<< endl;
    outputfile7<< "z0 , leftInWall, "      << Wi_z0_leftInWall      << "," << We_z0_leftInWall      << "," << (Wi_z0_leftInWall      + We_z0_leftInWall)       << endl;
    outputfile7<< "z1 , leftAntBaseWall, " << Wi_z1_leftAntBaseWall << "," << We_z1_leftAntBaseWall << "," << (Wi_z1_leftAntBaseWall + We_z1_leftAntBaseWall)  << endl;
    outputfile7<< "z2 , leftAntTopWall, "  << Wi_z2_leftAntTopWall  << "," << We_z2_leftAntTopWall  << "," << (Wi_z2_leftAntTopWall  + We_z2_leftAntTopWall)   << endl;
    outputfile7<< "z3 , rightInWall, "     << Wi_z3_rightInWall     << "," << We_z3_rightInWall     << "," << (Wi_z3_rightInWall     + We_z3_rightInWall)      << endl;
    outputfile7<< "z4 , leftOutWall, "     << Wi_z4_leftOutWall     << "," << We_z4_leftOutWall     << "," << (Wi_z4_leftOutWall     + We_z4_leftOutWall)      << endl;
    outputfile7<< "z5 , rightAnodeWall, "  << Wi_z5_rightAnodeWall  << "," << We_z5_rightAnodeWall  << "," << (Wi_z5_rightAnodeWall  + We_z5_rightAnodeWall)   << endl;
    outputfile7<< "x0 , botInWall, "       << Wi_x0_botInWall       << "," << We_x0_botInWall       << "," << (Wi_x0_botInWall       + We_x0_botInWall)        << endl;
    outputfile7<< "x1 , topInWall, "       << Wi_x1_topInWall       << "," << We_x1_topInWall       << "," << (Wi_x1_topInWall       + We_x1_topInWall)        << endl;
    outputfile7<< "x2 , botAntSideWall, "  << Wi_x2_botAntSideWall  << "," << We_x2_botAntSideWall  << "," << (Wi_x2_botAntSideWall  + We_x2_botAntSideWall)   << endl;
    outputfile7<< "x4 , topOrfWall, "      << Wi_x4_topOrfWall      << "," << We_x4_topOrfWall      << "," << (Wi_x4_topOrfWall      + We_x4_topOrfWall)       << endl;
    outputfile7<< "x5 , topOpen, "         << Wi_x5_topOpen         << "," << We_x5_topOpen         << "," << (Wi_x5_topOpen         + We_x5_topOpen)          << endl;
    outputfile7<< "*  , sum,"              << Wi_sum                << "," << We_sum                << "," << (Wi_sum                + We_sum)                 << endl;
    outputfile7.close();
    //------------------------------------

    double direct_ionz_sum = 0.0;
    double stepwise_ionz_sum = 0.0;
    for (int iblock=0;iblock<2;iblock++){ 
        for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
            for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                direct_ionz_sum = direct_ionz_sum + gc.rhoe[i][j]*gc.nu_ionz[i][j]*(pm.dx*pm.dr*2.0*M_PI*gc.r[j]);
                stepwise_ionz_sum = stepwise_ionz_sum + gc.rhoe[i][j]*gc.nu_ionzStep[i][j]*(pm.dx*pm.dr*2.0*M_PI*gc.r[j]);
            }
        }
    }

    cout << "direct_ionz_sum = " << direct_ionz_sum  
        << " stepwise_ionz_sum = " << stepwise_ionz_sum << " ratio = " << stepwise_ionz_sum/(stepwise_ionz_sum+direct_ionz_sum)
        << endl;

    pm.nOut = pm.nOut + 1;
}

