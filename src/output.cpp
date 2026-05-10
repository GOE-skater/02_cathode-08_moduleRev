//*****************************************************************
//*****************************************************************
//**                                                             **
//**          plasma_check                                       **
//**          coded by Ryo Shirakawa.                            **
//*****************************************************************
//*****************************************************************

#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <time.h>
#include <vector>
#include "common.hpp"
#include "functions.hpp"
#include <sstream>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <complex>
#include <iomanip>
#include <limits>

//*****************************************************************
//**                                                             **
//**           void output                                       **
//**                                                             **
//*****************************************************************
void output(){
    
    //calculate maximum electron density to avoid zero-division
    //------------------------------------
    double rhoe_max = 0.0;
    double ratio_rhoe = 0.001;
    for (int iblock=0;iblock<n_bl-1;iblock++){ 
        for (int i=i_flc_bl[iblock][0];i<=i_flc_bl[iblock][1];i++){ 
            for (int j=j_flc_bl[iblock][0];j<=j_flc_bl[iblock][1];j++){
                if(rhoe_max < rhoe[i][j]){
                    rhoe_max = rhoe[i][j];
                }
            }
        }
    }
    //------------------------------------

    //time step
    //------------------------------------
    double dt_i = dt*double(ndt_i);
    double dt_m = dt*double(ndt_m);
    double dt_n = dt*double(ndt_n);
    //------------------------------------

    std::string char_cnv,char1,char2,char_csv;

    char1="results/result";
    char2=std::to_string(nOut);
    char_csv=".csv";

    std::cout<<"output!"<<std::endl;

    std::ofstream outputfile1(char1+char2+char_csv);
    //std::ofstream outputfile1_1(char1+"_tmp"+char2+char_csv);
    //outputfile1<<"x,rho,Ui,rhoUix,Ue,rhoUex,rhoUex_E,rhoUex_D,E,phi,psi,rho_th,U_th,E_th,phi_th,psi_th,rate_ionize,divi,dive1,dive2,dive3, div_poisson, LHS, RHS,zero" << std::endl;
    //outputfile1<<"x,rhon,rho,rhoe(test),Uix,rhoUix,Uex,Uey,rhoUex,rhoUey,Te,heat_flux,Bz,E,phi,nu_m,nu_en,nu_ei,nu_ionz,nu_exc,nu_wall,nu_ano,mue_parae,mue_perp,Halle,jd,Id,rate_eloss,rate_ionize,divi,dive,divn,engy_LHS,engy_LHS1,engy_LHS2,engy_LHS3,engy_LHS4,engy_RHS,engy_RHS1,engy_RHS2,engy_RHS3,zero" << std::endl;
    outputfile1<<"i,j,x,r,rhoi,rhoe,Uix,Uir,Uip,Ui_mag,Te,Te_eff,rhoeps,qex,qer,Gx,Gr,phi,Ex,Er,rhoUex,rhoUer,rhoUep,Uex,Uer,Uep,Ue_mag,rate_ionize,nu_m,nu_m1,nu_ionz,nu_ionzStep,nu_excStep,nu_exc,nu_excMeta,nui_m,mu_para,mu_perp,mu_x,mu_r,mu_c,kappa_para,kappa_perp,kappa_x,kappa_r,kappa_c,cx,cr,sx,sr,Ap,Bx,Br,Bmag,Omega,Omegax,Omegar,Pabs,rate_eloss,Loss,P_ES,Pnet,nabla_rhoUe,nabla_rhoUa,nabla_rhoUi,nabla_rhoUm,rate_rhom,rateG_rhom,rateL_rhom,nabla_G,nabla_Ga,rhom,rhoUmx,rhoUmr,Umx,Umr,limit,Debye,phi_analytical,|Ex_mw|,|Er_mw|,|Ep_mw|,rate_ionize_direct,rate_ionize_stepwise,rhon,d|B|/dx,d|B|/dr,rate_eloss_n,rate_eloss_m,distECR,nu_eff,Jx_exc,Jr_exc,Jp_exc,scx,scr,rhoq_sc,nUex,nUer,nabla_nUe,rhoUnx,rhoUnr,nabla_rhoUn,rate_rhon,nu_super,deltaB,nu_ano,nu_ano_IAT,jdgBnd_flc,Mache,zero" << std::endl;
    //outputfile1_1 << "i,j,x,r,rhom,rhoUmx,rhoUmr,nabla_rhoUm,rhon,rhoUnx,rhoUnr,nabla_rhoUn,zero" << std::endl;
    
    for(int i=1;i<=ni;i++){
        for(int j=1;j<=nj;j++){
            double rL = (r[j]+r[j-1])/2.0;
            double rR = (r[j]+r[j+1])/2.0;
            double qR = (r[j]+r[j+1])/2.0/r[j];
            double qL = (r[j]+r[j-1])/2.0/r[j];

            double Bmag = sqrt(Bx[i][j]*Bx[i][j] + Br[i][j]*Br[i][j]);
            double absB_Lx = sqrt(Bx[i-1][j]*Bx[i-1][j] + Br[i-1][j]*Br[i-1][j]);
            double absB_Rx = sqrt(Bx[i+1][j]*Bx[i+1][j] + Br[i+1][j]*Br[i+1][j]);
            double absB_Lr = sqrt(Bx[i][j-1]*Bx[i][j-1] + Br[i][j-1]*Br[i][j-1]);
            double absB_Rr = sqrt(Bx[i][j+1]*Bx[i][j+1] + Br[i][j+1]*Br[i][j+1]);
            double dBdx = (absB_Rx-absB_Lx)/(2.0*dx);
            double dBdr = (absB_Rr-absB_Lr)/(2.0*dr);

            //double nu_ano = e0*Bmag/masse/16.0;
            double Omega = e0*Bmag/(masse*(nu_m[i][j])+1e-100);
            double Omega_max = 100;
            double ratio = 1.0;//fmin(Omega_max/Omega,1.0);
            Omega = Omega*ratio;
            double Omega_x = e0*Bx[i][j]/(masse*nu_m[i][j]+1e-100)*ratio;
            double Omega_r = e0*Br[i][j]/(masse*nu_m[i][j]+1e-100)*ratio;

            double nablaB = sqrt(pow((Bx[i+1][j] - Bx[i-1][j])/(2.0*dx),2) 
                + pow((Br[i+1][j] - Br[i-1][j])/(2.0*dr),2));

            double deltaB = Bmag/(nablaB+1e-100);


            double Uix_tmp = (Uix[i][j] + Uix[i+1][j])/2.0;
            double Uir_tmp = (rL*Uir[i][j] + rR*Uir[i][j+1])/(2.0*r[j]);
            double Uip_tmp = Uip[i][j];

            double Ui_mag = sqrt(Uix_tmp*Uix_tmp + Uir_tmp*Uir_tmp);
            
            double rhoUex_tmp = (rhoUex[i][j] + rhoUex[i+1][j])/2.0;
            double rhoUer_tmp = (rL*rhoUer[i][j] + rR*rhoUer[i][j+1])/(2.0*r[j]);
            double rhoUep_tmp = rhoUer_tmp*Omega_x - rhoUex_tmp*Omega_r;
            double Uex_tmp = rhoUex_tmp/(rhoe[i][j]+ratio_rhoe*rhoe_max);
            double Uer_tmp = rhoUer_tmp/(rhoe[i][j]+ratio_rhoe*rhoe_max);
            double Uep_tmp = rhoUep_tmp/(rhoe[i][j]+ratio_rhoe*rhoe_max);

            

            double Ue_mag = sqrt(Uex_tmp*Uex_tmp + Uer_tmp*Uer_tmp);

            double Te_eff = Te[i][j] + 1.0/(3.0*Boltz)*masse*Ue_mag*Ue_mag;

            double rhoUmx_tmp = (rhoUmx[i][j] + rhoUmx[i+1][j])/2.0;
            double rhoUmr_tmp = (rL*rhoUmr[i][j] + rR*rhoUmr[i][j+1])/(2.0*r[j]);

            double rhoUnx_tmp = (rhoUnx[i][j] + rhoUnx[i+1][j])/2.0;
            double rhoUnr_tmp = (rL*rhoUnr[i][j] + rR*rhoUnr[i][j+1])/(2.0*r[j]);

            double Ex_tmp = (Ex[i][j] + Ex[i+1][j])/2.0;
            double Er_tmp = (rL*Er[i][j] + rR*Er[i][j+1])/(2.0*r[j]);

            double rhoi_Lx = (rhoi[i][j] + rhoi[i-1][j])/2.0;
            double rhoi_Rx = (rhoi[i][j] + rhoi[i+1][j])/2.0;
            double rhoi_Lr = (rhoi[i][j] + rhoi[i][j-1])/(2.0);
            double rhoi_Rr = (rhoi[i][j] + rhoi[i][j+1])/(2.0);

            //if(j==1) std::cout << rhoUer[i][j-1] << std::endl;

            double nabla_rhoUe = (rhoUex[i+1][j]-rhoUex[i][j])/dx + (qR*rhoUer[i][j+1]-qL*rhoUer[i][j])/dr;

            double nabla_rhoUi = (rhoi_Rx*Uix[i+1][j]-rhoi_Lx*Uix[i][j])/dx + (qR*rhoi_Rr*Uir[i][j+1]-qL*rhoi_Lr*Uir[i][j])/dr;
            double nabla_rhoUm = (rhoUmx[i+1][j]-rhoUmx[i][j])/dx + (qR*rhoUmr[i][j+1]-qL*rhoUmr[i][j])/dr;
            double nabla_rhoUn = (rhoUnx[i+1][j]-rhoUnx[i][j])/dx + (qR*rhoUnr[i][j+1]-qL*rhoUnr[i][j])/dr;

            double lapPhi =  -(Ex[i+1][j] - Ex[i][j])/dx - (qR*Er[i][j+1] - qL*Er[i][j])/dr;
            double nabla_nUe = (nUex[i+1][j]-nUex[i][j])/dx + (qR*nUer[i][j+1]-qL*nUer[i][j])/dr;
            double nUex_tmp = (nUex[i][j] + nUex[i+1][j])/2.0;
            double nUer_tmp = (rL*nUer[i][j] + rR*nUer[i][j+1])/(2.0*r[j]);

            double Gm = rhoe[i][j]*(nu_excMeta[i][j] + 0.5*nu_exc[i][j]); //基底からの励起 + 上の準位からのdrop
            double Lm = rhoe[i][j]*(nu_ionzStep[i][j] + nu_excStep[i][j] + nu_super[i][j]);

            double Gn = (nu_excStep[i][j] + nu_super[i][j])*rhoe[i][j]; 
            double Ln = (nu_excMeta[i][j] + 0.5*nu_exc[i][j] + nu_ionz[i][j])*rhoe[i][j];
            
            double nabla_G = (Gx[i+1][j]-Gx[i][j])/dx + (qR*Gr[i][j+1]-qL*Gr[i][j])/dr;
            double RHS2 =  e0*(Ex_tmp*rhoUex_tmp + Er_tmp*rhoUer_tmp);

            double P_ES = -RHS2;

            double Ex_mw_tmp = std::abs(E1x[i][j] + E1x[i+1][j])/2.0;
            double Er_mw_tmp = std::abs(rL*E1r[i][j+1] + rR*E1r[i][j])/(2.0*r[j]);
            double Ep_mw_tmp = std::abs(E1p[i][j]);

            double Jx_exc_tmp = std::abs(J1x[i][j] + J1x[i+1][j])/2.0;
            double Jr_exc_tmp = std::abs(rL*J1r[i][j+1] + rR*J1r[i][j])/(2.0*r[j]);
            double Jp_exc_tmp = std::abs(J1p[i][j]);

            double ut = sqrt(2.0*Boltz*(Te[i][j]+1.0*e0/Boltz)/masse);
            double Mache_tmp = sqrt(pow(Uex_tmp,2) + pow(Uer_tmp,2))/(ut + 100.0);

            //表面電荷の電荷密度換算量
            double rhoq_sc = (scx[i][j]*2.0*M_PI*r[j]*dr + scr[i][j]*2.0*M_PI*rR*dx)/(2.0*M_PI*r[j]*dx*dr); 

            outputfile1<< i << ","<< j << "," << x[i]<< ","<< r[j]
                << "," << rhoi[i][j]*jdgBnd_flc[i][j] << "," << rhoe[i][j] *jdgBnd_flc[i][j]
                << ","<< Uix_tmp*jdgBnd_flc[i][j] << ","<< Uir_tmp*jdgBnd_flc[i][j]<< ","<< Uip_tmp*jdgBnd_flc[i][j]<< "," << Ui_mag*jdgBnd_flc[i][j]
                << ","<< Te[i][j]*Boltz/e0*jdgBnd_flc[i][j]
                << ","<< Te_eff*Boltz/e0*jdgBnd_flc[i][j]
                << ","<< rhoeps[i][j]*jdgBnd_flc[i][j]
                << ","<< qex[i][j]*jdgBnd_flc[i][j]<< ","<< qer[i][j]*jdgBnd_flc[i][j]
                << ","<< Gx[i][j]*jdgBnd_flc[i][j]<< ","<< Gr[i][j]*jdgBnd_flc[i][j]
                << ","<< phi[i][j]*jdgBnd_flc[i][j]<< ","<< Ex[i][j]*jdgBnd_flc[i][j]<< ","<< Er[i][j]*jdgBnd_flc[i][j]
                << ","<< rhoUex_tmp*jdgBnd_flc[i][j]<< ","<< rhoUer_tmp*jdgBnd_flc[i][j]<< ","<< rhoUep_tmp*jdgBnd_flc[i][j]
                << ","<< Uex_tmp*jdgBnd_flc[i][j]<< ","<< Uer_tmp*jdgBnd_flc[i][j]<< ","<< Uep_tmp*jdgBnd_flc[i][j]<< ","<< Ue_mag*jdgBnd_flc[i][j]
                << ","<< rate_ionize[i][j]*jdgBnd_flc[i][j]<< ","<< nu_m[i][j]*jdgBnd_flc[i][j]<< ","<< nu_m1[i][j]*jdgBnd_flc[i][j]
                << ","<< nu_ionz[i][j]*jdgBnd_flc[i][j] << ","<< nu_ionzStep[i][j]*jdgBnd_flc[i][j]<< ","<< nu_excStep[i][j]*jdgBnd_flc[i][j] 
                << ","<< nu_exc[i][j] *jdgBnd_flc[i][j]<< ","<< nu_excMeta[i][j] *jdgBnd_flc[i][j]<< ","<< nui_m[i][j] *jdgBnd_flc[i][j]
                << ","<< mu_para[i][j]*jdgBnd_flc[i][j]<< ","<< mu_perp[i][j]*jdgBnd_flc[i][j]
                << ","<< mu_xx[i][j]*jdgBnd_flc[i][j]<< ","<< mu_rr[i][j]*jdgBnd_flc[i][j]<< ","<< mu_xr[i][j]*jdgBnd_flc[i][j]
                << ","<< kappa_para[i][j]*jdgBnd_flc[i][j]<< ","<< kappa_perp[i][j]*jdgBnd_flc[i][j]
                << ","<< kappa_x[i][j]*jdgBnd_flc[i][j]<< ","<< kappa_r[i][j]*jdgBnd_flc[i][j]<< ","<< kappa_c[i][j]*jdgBnd_flc[i][j]
                << ","<< cx[i][j]*jdgBnd_flc[i][j]<< ","<< cr[i][j]*jdgBnd_flc[i][j]<< ","<< sx[i][j]*jdgBnd_flc[i][j]<< ","<< sr[i][j]*jdgBnd_flc[i][j]
                << ","<< Ap[i][j]*jdgBnd_flc[i][j]<< ","<< Bx[i][j]*jdgBnd_flc[i][j]<< ","<< Br[i][j]*jdgBnd_flc[i][j]<< ","<< Bmag*jdgBnd_flc[i][j]
                << ","<< Omega*jdgBnd_flc[i][j]<< ","<< Omega_x*jdgBnd_flc[i][j]<< ","<< Omega_r*jdgBnd_flc[i][j]
                << ","<< Pabs[i][j]*jdgBnd_flc[i][j]
                << ","<< rate_eloss[i][j]*jdgBnd_flc[i][j]
                << ","<< rhoe[i][j]*rate_eloss[i][j]*jdgBnd_flc[i][j]
                << ","<< P_ES*jdgBnd_flc[i][j]
                << ","<< (Pabs[i][j] - rhoe[i][j]*rate_eloss[i][j] + P_ES)*jdgBnd_flc[i][j]
                << ","<< (nabla_rhoUe + (rhoe[i][j] - rhoe_old[i][j])/dt)*jdgBnd_flc[i][j]
                << ","<< (nabla_rhoUe)*jdgBnd_flc[i][j]
                << ","<< (nabla_rhoUi + (rhoi[i][j] - rhoi_old[i][j])/dt_i)*jdgBnd_flc[i][j]
                << ","<< (nabla_rhoUm + (rhom[i][j] - rhom_old[i][j])/dt_m)*jdgBnd_flc[i][j]
                << ","<< (Gm - Lm)*jdgBnd_flc[i][j]
                << ","<< (Gm)*jdgBnd_flc[i][j]
                << ","<< (Lm)*jdgBnd_flc[i][j]
                << ","<< (nabla_G + (rhoeps[i][j] - rhoeps_old[i][j])/dt - P_ES + rhoe[i][j]*rate_eloss[i][j])*jdgBnd_flc[i][j]
                << ","<< (nabla_G  + rhoe[i][j]*rate_eloss[i][j])*jdgBnd_flc[i][j]
                << "," << rhom[i][j]*jdgBnd_flc[i][j]<< "," << rhoUmx_tmp*jdgBnd_flc[i][j]<< "," << rhoUmr_tmp*jdgBnd_flc[i][j]
                << "," << rhoUmx_tmp/(rhom[i][j]+1e-100)*jdgBnd_flc[i][j]<< "," << rhoUmr_tmp/(rhom[i][j]+1e-100)*jdgBnd_flc[i][j]
                << "," << limit[i][j]*jdgBnd_flc[i][j]
                << "," << sqrt(eps0*Boltz*Te[i][j]/((rhoe[i][j] + 1)*e0*e0))*jdgBnd_flc[i][j]
                << "," << Boltz*Te[i][j]/e0*(0.5*log(massi/(2.0*M_PI*masse)) + log(fabs(rhoi_Lr*Uir[i][j]/(rhoUer[i][j] + 1e-100))))*jdgBnd_flc[i][j]
                << "," << Ex_mw_tmp*double(jdgBnd_Ep[i][j])<< "," << Er_mw_tmp*double(jdgBnd_Ep[i][j])<< "," << Ep_mw_tmp*double(jdgBnd_Ep[i][j])
                << "," << rhoe[i][j]*nu_ionz[i][j]*jdgBnd_flc[i][j] << ","<<rhoe[i][j]*nu_ionzStep[i][j]*jdgBnd_flc[i][j]
                << "," << rhon[i][j]*jdgBnd_flc[i][j]
                << "," << dBdx*jdgBnd_flc[i][j] << ","<<dBdr*jdgBnd_flc[i][j]
                << "," << rhon[i][j]*func_Te(rate_eloss_n_N,Te[i][j],dTe)*jdgBnd_flc[i][j]<< "," << rhom[i][j]*func_Te(rate_eloss_m_N,Te[i][j],dTe)*jdgBnd_flc[i][j]
                << "," << distECR[i][j]*jdgBnd_flc[i][j]
                << "," << nu_eff/(1.0 + 4.0*pow(distECR[i][j]/deltaECR,2))*jdgBnd_flc[i][j]
                << "," << Jx_exc_tmp<< "," << Jr_exc_tmp<< "," << Jp_exc_tmp
                << "," << scx[i][j]*jdgBnd_flc[i][j]
                << "," << scr[i][j]*jdgBnd_flc[i][j]
                << "," << rhoq_sc*jdgBnd_flc[i][j]
                << ","<< nUex_tmp*jdgBnd_flc[i][j]<< ","<< nUer_tmp*jdgBnd_flc[i][j]
                << ","<< (nabla_nUe + (rhoi[i][j] - rhoe_old[i][j])/dt + eps0/(e0*dt)*lapPhi)*jdgBnd_flc[i][j]
                << "," << rhoUnx_tmp*jdgBnd_flc[i][j]<< "," << rhoUnr_tmp*jdgBnd_flc[i][j]
                << ","<< (nabla_rhoUn + (rhon[i][j] - rhon_old[i][j])/dt_n)*jdgBnd_flc[i][j]
                << ","<< (Gn - Ln)*jdgBnd_flc[i][j]
                << ","<< nu_super[i][j]*jdgBnd_flc[i][j]
                << ","<< deltaB*jdgBnd_flc[i][j]
                << ","<< nu_ano[i][j]*jdgBnd_flc[i][j]
                << ","<< nu_ano_IAT[i][j]*jdgBnd_flc[i][j]
                << ","<< jdgBnd_flc[i][j]
                << ","<< Mache_tmp*jdgBnd_flc[i][j]
                << "," << 0.0<< std::endl;
        }
    }

    outputfile1.close();
    //outputfile1_1.close();

    char1="results/restart";
    std::string char_num = + "_ni" + std::to_string(ni) + "_nj" + std::to_string(nj)+ "_";
    char2=std::to_string(nOut);
    char_csv=".csv";

    std::ofstream outputfile2(char1+char2+char_csv);

    outputfile2<<"i,j,x,r,rhoi,rhoi_old,Uix,Uix_old,Uir,Uir_old,Uip,Uip_old,rhoe,rhoe_old,rhoUex,rhoUex_old,rhoUer,rhoUer_old,Te,Te_old,Gx,Gx_old,Gr,Gr_old,phi,phi_old,rhom,rhom_old,rhoUmx,rhoUmr,rate_ionize,scx,scr,nUex,nUex_old,nUer,nUer_old,rhon,rhon_old,rhoUnx,rhoUnr,rhoeps,rhoeps_old,zero"<< std::endl;
    outputfile2 << std::setprecision(std::numeric_limits<double>::max_digits10) << std::scientific;
    
    for(int i=0;i<ni+2;i++){
        for(int j=0;j<nj+2;j++){
            outputfile2<< i << ","<< j << "," << x[i]<< ","<< r[j]
                << "," << truncate(rhoi[i][j]) << "," << truncate(rhoi_old[i][j]) 
                << "," << truncate(Uix[i][j])  << "," << truncate(Uix_old[i][j])  
                << "," << truncate(Uir[i][j])<< "," << truncate(Uir_old[i][j])
                << "," << truncate(Uip[i][j])<< "," << truncate(Uip_old[i][j])
                << "," << truncate(rhoe[i][j]) << "," << truncate(rhoe_old[i][j]) 
                << "," << truncate(rhoUex[i][j]) << "," << truncate(rhoUex_old[i][j]) 
                << "," << truncate(rhoUer[i][j])<< "," << truncate(rhoUer_old[i][j])
                << "," << truncate(Te[i][j])*Boltz/e0 << "," << truncate(Te_old[i][j])*Boltz/e0
                << "," << truncate(Gx[i][j]) << "," << truncate(Gx_old[i][j]) 
                << "," << truncate(Gr[i][j])<< "," << truncate(Gr_old[i][j])
                << "," << truncate(phi[i][j])<< "," << truncate(phi_old[i][j])
                << "," << truncate(rhom[i][j])<< "," << truncate(rhom_old[i][j])
                << "," << truncate(rhoUmx[i][j])<< "," << truncate(rhoUmr[i][j])
                << "," << truncate(rate_ionize[i][j])
                << "," << truncate(scx[i][j])
                << "," << truncate(scr[i][j])
                << "," << truncate(nUex[i][j]) << "," << truncate(nUex_old[i][j]) 
                << "," << truncate(nUer[i][j]) << "," << truncate(nUer_old[i][j]) 
                << "," << truncate(rhon[i][j]) << "," << truncate(rhon_old[i][j]) 
                << "," << truncate(rhoUnx[i][j]) << "," << truncate(rhoUnr[i][j]) 
                << "," << truncate(rhoeps[i][j]) << "," << truncate(rhoeps_old[i][j]) 
                << "," << 0.0<< std::endl;
        }
    }
    outputfile2.close();

    //File operation
    //===============================================
    //setup for sum of current per wall
    //------------------------------------
    char1="results/current_sum";
    std::ofstream outputfile3(char1+char2+char_csv);
    outputfile3 << std::setprecision(std::numeric_limits<double>::max_digits10) << std::scientific;

    Ii_Anode = 0.0;
    Ie_Anode = 0.0;
    I_Anode = 0.0;

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
    std::ofstream outputfile4(char1+char2+char_csv);
    outputfile4 << std::setprecision(std::numeric_limits<double>::max_digits10) << std::scientific;
    outputfile4 << "i,j,x,r, Ii (mA), Ie (mA), I (mA), ji (mA/mm^2), je (mA/mm^2), j (mA/mm^2), En, zero"<< std::endl;
    //------------------------------------

    //setup for boundary velocity for kinetic flux
    //------------------------------------
    //char1="results/bnd_KF_shift";
    //std::ofstream outputfile6(char1+char2+char_csv);
    //outputfile6 << std::setprecision(std::numeric_limits<double>::max_digits10) << std::scientific;
    //outputfile6 << "i,j,x,r, uten, Uten, uten_cc,uten_cc, utin, Utin, utin_cc, Utin_cc, utmn, vth2, omega_ce, Larmor, LarmorFactor, angle, sin(angle), LarmorFactor*cos(angle),coef,coef_mirror,zero"<< std::endl;
    //------------------------------------
    
    //setup for sum of current per wall
    //------------------------------------
    char1="results/energy_flux_sum";
    std::ofstream outputfile7(char1+char2+char_csv);
    outputfile7 << std::setprecision(std::numeric_limits<double>::max_digits10) << std::scientific;
    //------------------------------------
    //===============================================

    //Caluculate wall information
    //===============================================
    //left-wall-BC (z0)
    //------------------------------------
    for (int j=j_flc_bl[0][0];j<=j_flc_bl[0][1];j++){
        int i=i_flc_bl[0][0];
        double x_tmp = (x[i] + x[i-1])/2.0;

        double En = (-1.0)*Ex[i][j];

        double ji = (-1.0)*( e0)*rhoUix_wall[i][j];
        double je = (-1.0)*(-e0)*rhoUex_wall[i][j];

        double Ii = ji*2.0*M_PI*r[j]*dr;
        double Ie = je*2.0*M_PI*r[j]*dr;
        
        outputfile4<< i << ","<< j << "," << x_tmp << "," << r[j] 
            << ","<< Ii*1000 << ","<< Ie*1000<< ","<< (Ii+Ie)*1000 
            << ","<< ji/1000 << ","<< je/1000<< ","<< (ji+je)/1000 
            << ","<< En << ","<< 0 << std::endl;

        Ii_z0_leftInWall = Ii_z0_leftInWall - ( e0)*rhoUix_wall[i][j]*2.0*M_PI*r[j]*dr;
        Ie_z0_leftInWall = Ie_z0_leftInWall - (-e0)*rhoUex_wall[i][j]*2.0*M_PI*r[j]*dr;
        Im_z0_leftInWall = Im_z0_leftInWall -       rhoUmx_wall[i][j]*2.0*M_PI*r[j]*dr;
        In_z0_leftInWall = In_z0_leftInWall -       rhoUnx_wall[i][j]*2.0*M_PI*r[j]*dr;

        double engyi = 0.5*massi*(Uix[i][j]*Uix[i][j] + Uir[i][j]*Uir[i][j] + Uip[i][j]*Uip[i][j]) + 5.0/2.0*Boltz*Ti;
        Wi_z0_leftInWall = Wi_z0_leftInWall - engyi*rhoUix_wall[i][j]*2.0*M_PI*r[j]*dr;
        We_z0_leftInWall = We_z0_leftInWall - Gx_wall[i][j]*2.0*M_PI*r[j]*dr;
    }

    //left-wall-BC (z1)
    //------------------------------------
    for (int j=j_flc_bl[2][0];j<=j_flc_bl[0][0]-1;j++){
        int i=i_flc_bl[2][0];
        double x_tmp = (x[i] + x[i-1])/2.0;

        double En = (-1.0)*Ex[i][j];

        double ji = (-1.0)*( e0)*rhoUix_wall[i][j];
        double je = (-1.0)*(-e0)*rhoUex_wall[i][j];

        double Ii = ji*2.0*M_PI*r[j]*dr;
        double Ie = je*2.0*M_PI*r[j]*dr;
        
        outputfile4<< i << ","<< j << "," << x_tmp << "," << r[j] 
            << ","<< Ii*1000 << ","<< Ie*1000<< ","<< (Ii+Ie)*1000 
            << ","<< ji/1000 << ","<< je/1000<< ","<< (ji+je)/1000 
            << ","<< En << ","<< 0 << std::endl;

        Ii_z1_leftAntBaseWall = Ii_z1_leftAntBaseWall - ( e0)*rhoUix_wall[i][j]*2.0*M_PI*r[j]*dr;
        Ie_z1_leftAntBaseWall = Ie_z1_leftAntBaseWall - (-e0)*rhoUex_wall[i][j]*2.0*M_PI*r[j]*dr;
        Im_z1_leftAntBaseWall = Im_z1_leftAntBaseWall -       rhoUmx_wall[i][j]*2.0*M_PI*r[j]*dr;
        In_z1_leftAntBaseWall = In_z1_leftAntBaseWall -       rhoUnx_wall[i][j]*2.0*M_PI*r[j]*dr;

        double engyi = 0.5*massi*(Uix[i][j]*Uix[i][j] + Uir[i][j]*Uir[i][j] + Uip[i][j]*Uip[i][j]) + 5.0/2.0*Boltz*Ti;
        Wi_z1_leftAntBaseWall = Wi_z1_leftAntBaseWall - engyi*rhoUix_wall[i][j]*2.0*M_PI*r[j]*dr;
        We_z1_leftAntBaseWall = We_z1_leftAntBaseWall - Gx_wall[i][j]*2.0*M_PI*r[j]*dr;
    }
    //------------------------------------

    //left-wall-BC (z2)
    //------------------------------------
    for (int j=j_flc_bl[3][0];j<=j_flc_bl[2][0]-1;j++){
        int i=i_flc_bl[3][0];
        double x_tmp = (x[i] + x[i-1])/2.0;

        double En = (-1.0)*Ex[i][j];

        double ji = (-1.0)*( e0)*rhoUix_wall[i][j];
        double je = (-1.0)*(-e0)*rhoUex_wall[i][j];

        double Ii = ji*2.0*M_PI*r[j]*dr;
        double Ie = je*2.0*M_PI*r[j]*dr;
        
        outputfile4<< i << ","<< j << "," << x_tmp << "," << r[j] 
            << ","<< Ii*1000 << ","<< Ie*1000<< ","<< (Ii+Ie)*1000 
            << ","<< ji/1000 << ","<< je/1000<< ","<< (ji+je)/1000 
            << ","<< En << ","<< 0 << std::endl;

        Ii_z2_leftAntTopWall = Ii_z2_leftAntTopWall - ( e0)*rhoUix_wall[i][j]*2.0*M_PI*r[j]*dr;
        Ie_z2_leftAntTopWall = Ie_z2_leftAntTopWall - (-e0)*rhoUex_wall[i][j]*2.0*M_PI*r[j]*dr;
        Im_z2_leftAntTopWall = Im_z2_leftAntTopWall -       rhoUmx_wall[i][j]*2.0*M_PI*r[j]*dr;
        In_z2_leftAntTopWall = In_z2_leftAntTopWall -       rhoUnx_wall[i][j]*2.0*M_PI*r[j]*dr;

        double engyi = 0.5*massi*(Uix[i][j]*Uix[i][j] + Uir[i][j]*Uir[i][j] + Uip[i][j]*Uip[i][j]) + 5.0/2.0*Boltz*Ti;
        Wi_z2_leftAntTopWall = Wi_z2_leftAntTopWall - engyi*rhoUix_wall[i][j]*2.0*M_PI*r[j]*dr;
        We_z2_leftAntTopWall = We_z2_leftAntTopWall - Gx_wall[i][j]*2.0*M_PI*r[j]*dr;
    }
    //------------------------------------

    //left-wall-BC (z4)
    //------------------------------------
    for (int j=j_flc_bl[1][0];j<=j_flc_bl[4][1];j++){
        int i=i_flc_bl[4][0];
        double x_tmp = (x[i] + x[i-1])/2.0;

        double En = (-1.0)*Ex[i][j];

        double ji = (-1.0)*( e0)*rhoUix_wall[i][j];
        double je = (-1.0)*(-e0)*rhoUex_wall[i][j];

        double Ii = ji*2.0*M_PI*r[j]*dr;
        double Ie = je*2.0*M_PI*r[j]*dr;
        
        outputfile4<< i << ","<< j << "," << x_tmp << "," << r[j] 
            << ","<< Ii*1000 << ","<< Ie*1000<< ","<< (Ii+Ie)*1000 
            << ","<< ji/1000 << ","<< je/1000<< ","<< (ji+je)/1000 
            << ","<< En << ","<< 0 << std::endl;

        Ii_z4_leftOutWall = Ii_z4_leftOutWall - ( e0)*rhoUix_wall[i][j]*2.0*M_PI*r[j]*dr;
        Ie_z4_leftOutWall = Ie_z4_leftOutWall - (-e0)*rhoUex_wall[i][j]*2.0*M_PI*r[j]*dr;
        Im_z4_leftOutWall = Im_z4_leftOutWall -       rhoUmx_wall[i][j]*2.0*M_PI*r[j]*dr;
        In_z4_leftOutWall = In_z4_leftOutWall -       rhoUnx_wall[i][j]*2.0*M_PI*r[j]*dr;

        double engyi = 0.5*massi*(Uix[i][j]*Uix[i][j] + Uir[i][j]*Uir[i][j] + Uip[i][j]*Uip[i][j]) + 5.0/2.0*Boltz*Ti;
        Wi_z4_leftOutWall = Wi_z4_leftOutWall - engyi*rhoUix_wall[i][j]*2.0*M_PI*r[j]*dr;
        We_z4_leftOutWall = We_z4_leftOutWall - Gx_wall[i][j]*2.0*M_PI*r[j]*dr;
    }
    //------------------------------------

    //right-wall-BC (z3)
    //------------------------------------
    for (int j=j_flc_bl[1][0];j<=j_flc_bl[1][1];j++){
        int i=i_flc_bl[1][1];

        double x_tmp = (x[i] + x[i+1])/2.0;

        double En = (+1.0)*Ex[i+1][j];
        //std::cout << Ex[i+1][j] << ","<<Ex_old[i+1][j] << std::endl;

        double ji = (+1.0)*( e0)*rhoUix_wall[i+1][j];
        double je = (+1.0)*(-e0)*rhoUex_wall[i+1][j];
        double Ii = ji*2.0*M_PI*r[j]*dr;
        double Ie = je*2.0*M_PI*r[j]*dr;

        outputfile4<< i+1 << ","<< j << "," << x_tmp << "," << r[j] 
            << ","<< Ii*1000 << ","<< Ie*1000<< ","<< (Ii+Ie)*1000 
            << ","<< ji/1000 << ","<< je/1000<< ","<< (ji+je)/1000 
            << ","<< En << ","<< 0 << std::endl;

        Ii_z3_rightInWall = Ii_z3_rightInWall + ( e0)*rhoUix_wall[i+1][j]*2.0*M_PI*r[j]*dr;
        Ie_z3_rightInWall = Ie_z3_rightInWall + (-e0)*rhoUex_wall[i+1][j]*2.0*M_PI*r[j]*dr;
        Im_z3_rightInWall = Im_z3_rightInWall +       rhoUmx_wall[i+1][j]*2.0*M_PI*r[j]*dr;
        In_z3_rightInWall = In_z3_rightInWall +       rhoUnx_wall[i+1][j]*2.0*M_PI*r[j]*dr;

        double engyi = 0.5*massi*(Uix[i+1][j]*Uix[i+1][j] + Uir[i][j]*Uir[i][j] + Uip[i][j]*Uip[i][j]) + 5.0/2.0*Boltz*Ti;
        Wi_z3_rightInWall = Wi_z3_rightInWall + engyi*rhoUix_wall[i+1][j]*2.0*M_PI*r[j]*dr;
        We_z3_rightInWall = We_z3_rightInWall + Gx_wall[i+1][j]*2.0*M_PI*r[j]*dr;
    }
    //------------------------------------

    //right-wall-BC (z5)
    //------------------------------------
    for (int j=j_flc_bl[4][0];j<=j_flc_bl[4][1];j++){
        int i=i_flc_bl[4][1];

        double x_tmp = (x[i] + x[i+1])/2.0;

        double En = (+1.0)*Ex[i+1][j];
        //std::cout << Ex[i+1][j] << ","<<Ex_old[i+1][j] << std::endl;

        double ji = (+1.0)*( e0)*rhoUix_wall[i+1][j];
        double je = (+1.0)*(-e0)*rhoUex_wall[i+1][j];
        double Ii = ji*2.0*M_PI*r[j]*dr;
        double Ie = je*2.0*M_PI*r[j]*dr;

        outputfile4<< i+1 << ","<< j << "," << x_tmp << "," << r[j] 
            << ","<< Ii*1000 << ","<< Ie*1000<< ","<< (Ii+Ie)*1000 
            << ","<< ji/1000 << ","<< je/1000<< ","<< (ji+je)/1000 
            << ","<< En << ","<< 0 << std::endl;

        Ii_z5_rightAnodeWall = Ii_z5_rightAnodeWall + ( e0)*rhoUix_wall[i+1][j]*2.0*M_PI*r[j]*dr;
        Ie_z5_rightAnodeWall = Ie_z5_rightAnodeWall + (-e0)*rhoUex_wall[i+1][j]*2.0*M_PI*r[j]*dr;
        Im_z5_rightAnodeWall = Im_z5_rightAnodeWall +       rhoUmx_wall[i+1][j]*2.0*M_PI*r[j]*dr;
        In_z5_rightAnodeWall = In_z5_rightAnodeWall +       rhoUnx_wall[i+1][j]*2.0*M_PI*r[j]*dr;

        double engyi = 0.5*massi*(Uix[i+1][j]*Uix[i+1][j] + Uir[i][j]*Uir[i][j] + Uip[i][j]*Uip[i][j]) + 5.0/2.0*Boltz*Ti;
        Wi_z5_rightAnodeWall = Wi_z5_rightAnodeWall + engyi*rhoUix_wall[i+1][j]*2.0*M_PI*r[j]*dr;
        We_z5_rightAnodeWall = We_z5_rightAnodeWall + Gx_wall[i+1][j]*2.0*M_PI*r[j]*dr;
    }
    //------------------------------------
    
    //lower-wall-BC (x0)
    //------------------------------------
    for (int i=i_flc_bl[0][0];i<=i_flc_bl[0][1];i++){
        int j=j_flc_bl[0][0];

        double r_tmp = (r[j]+r[j-1])/2.0;

        double En = (-1.0)*Er[i][j];

        double ji = (-1.0)*( e0)*rhoUir_wall[i][j];
        double je = (-1.0)*(-e0)*rhoUer_wall[i][j];
        double Ii = ji*2.0*M_PI*r_tmp*dx;
        double Ie = je*2.0*M_PI*r_tmp*dx;
        
        outputfile4<< i << ","<< j << "," << x[i] << "," << r_tmp 
            << ","<< Ii*1000 << ","<< Ie*1000<< ","<< (Ii+Ie)*1000 
            << ","<< ji/1000 << ","<< je/1000<< ","<< (ji+je)/1000 
            << ","<< En << ","<< 0 << std::endl;

        Ii_x0_botInWall = Ii_x0_botInWall - ( e0)*rhoUir_wall[i][j]*2.0*M_PI*r_tmp*dx;
        Ie_x0_botInWall = Ie_x0_botInWall - (-e0)*rhoUer_wall[i][j]*2.0*M_PI*r_tmp*dx;
        Im_x0_botInWall = Im_x0_botInWall -       rhoUmr_wall[i][j]*2.0*M_PI*r_tmp*dx;
        In_x0_botInWall = In_x0_botInWall -       rhoUnr_wall[i][j]*2.0*M_PI*r_tmp*dx;
        
        double engyi = 0.5*massi*(Uix[i][j]*Uix[i][j] + Uir[i][j]*Uir[i][j] + Uip[i][j]*Uip[i][j]) + 5.0/2.0*Boltz*Ti;
        Wi_x0_botInWall = Wi_x0_botInWall - engyi*rhoUir_wall[i][j]*2.0*M_PI*r_tmp*dx;
        We_x0_botInWall = We_x0_botInWall - Gr_wall[i][j]*2.0*M_PI*r_tmp*dx;
    }
    //------------------------------------

    //lower-wall-BC (x2)
    //------------------------------------
    for (int i=i_flc_bl[2][0];i<=i_flc_bl[2][1];i++){
        int j=j_flc_bl[2][0];

         double r_tmp = (r[j]+r[j-1])/2.0;

        double En = (-1.0)*Er[i][j];

        double ji = (-1.0)*( e0)*rhoUir_wall[i][j];
        double je = (-1.0)*(-e0)*rhoUer_wall[i][j];
        double Ii = ji*2.0*M_PI*r_tmp*dx;
        double Ie = je*2.0*M_PI*r_tmp*dx;
        
        outputfile4<< i << ","<< j << "," << x[i] << "," << r_tmp 
            << ","<< Ii*1000 << ","<< Ie*1000<< ","<< (Ii+Ie)*1000 
            << ","<< ji/1000 << ","<< je/1000<< ","<< (ji+je)/1000 
            << ","<< En << ","<< 0 << std::endl;

        Ii_x2_botAntSideWall = Ii_x2_botAntSideWall - ( e0)*rhoUir_wall[i][j]*2.0*M_PI*r_tmp*dx;
        Ie_x2_botAntSideWall = Ie_x2_botAntSideWall - (-e0)*rhoUer_wall[i][j]*2.0*M_PI*r_tmp*dx;
        Im_x2_botAntSideWall = Im_x2_botAntSideWall -       rhoUmr_wall[i][j]*2.0*M_PI*r_tmp*dx;
        In_x2_botAntSideWall = In_x2_botAntSideWall -       rhoUnr_wall[i][j]*2.0*M_PI*r_tmp*dx;
        
        double engyi = 0.5*massi*(Uix[i][j]*Uix[i][j] + Uir[i][j]*Uir[i][j] + Uip[i][j]*Uip[i][j]) + 5.0/2.0*Boltz*Ti;
        Wi_x2_botAntSideWall = Wi_x2_botAntSideWall - engyi*rhoUir_wall[i][j]*2.0*M_PI*r_tmp*dx;
        We_x2_botAntSideWall = We_x2_botAntSideWall - Gr_wall[i][j]*2.0*M_PI*r_tmp*dx;
    }
    //------------------------------------

    //upper-wall-BC (x1)
    //------------------------------------
    for (int i=i_flc_bl[0][0];i<=i_flc_bl[1][1];i++){
        int j=j_flc_bl[0][1];

        double r_tmp = (r[j]+r[j+1])/2.0;

        double En = (+1.0)*Er[i][j+1];

        double ji = (+1.0)*( e0)*rhoUir_wall[i][j+1];
        double je = (+1.0)*(-e0)*rhoUer_wall[i][j+1];
        double Ii = ji*2.0*M_PI*r_tmp*dx;
        double Ie = je*2.0*M_PI*r_tmp*dx;
        
        outputfile4<< i << ","<< j+1 << "," << x[i] << "," << r_tmp 
            << ","<< Ii*1000 << ","<< Ie*1000<< ","<< (Ii+Ie)*1000 
            << ","<< ji/1000 << ","<< je/1000<< ","<< (ji+je)/1000 
            << ","<< En << ","<< 0 << std::endl;

        Ii_x1_topInWall = Ii_x1_topInWall + ( e0)*rhoUir_wall[i][j+1]*2.0*M_PI*r_tmp*dx;
        Ie_x1_topInWall = Ie_x1_topInWall + (-e0)*rhoUer_wall[i][j+1]*2.0*M_PI*r_tmp*dx;
        Im_x1_topInWall = Im_x1_topInWall +       rhoUmr_wall[i][j+1]*2.0*M_PI*r_tmp*dx;
        In_x1_topInWall = In_x1_topInWall +       rhoUnr_wall[i][j+1]*2.0*M_PI*r_tmp*dx;

        double engyi = 0.5*massi*(Uix[i][j]*Uix[i][j] + Uir[i][j+1]*Uir[i][j+1] + Uip[i][j]*Uip[i][j]) + 5.0/2.0*Boltz*Ti;
        Wi_x1_topInWall = Wi_x1_topInWall + engyi*rhoUir_wall[i][j+1]*2.0*M_PI*r_tmp*dx;
        We_x1_topInWall = We_x1_topInWall + Gr_wall[i][j+1]*2.0*M_PI*r_tmp*dx;
    }
    //------------------------------------


    //upper-wall-BC (x4)
    //------------------------------------
    for (int i=i_flc_bl[1][1]+1;i<=i_flc_bl[3][1];i++){
        int j=j_flc_bl[3][1];

        double r_tmp = (r[j]+r[j+1])/2.0;

        double En = (+1.0)*Er[i][j+1];

        double ji = (+1.0)*( e0)*rhoUir_wall[i][j+1];
        double je = (+1.0)*(-e0)*rhoUer_wall[i][j+1];
        double Ii = ji*2.0*M_PI*r_tmp*dx;
        double Ie = je*2.0*M_PI*r_tmp*dx;
        
        outputfile4<< i << ","<< j+1 << "," << x[i] << "," << r_tmp 
            << ","<< Ii*1000 << ","<< Ie*1000<< ","<< (Ii+Ie)*1000 
            << ","<< ji/1000 << ","<< je/1000<< ","<< (ji+je)/1000 
            << ","<< En << ","<< 0 << std::endl;

        Ii_x4_topOrfWall = Ii_x4_topOrfWall + ( e0)*rhoUir_wall[i][j+1]*2.0*M_PI*r_tmp*dx;
        Ie_x4_topOrfWall = Ie_x4_topOrfWall + (-e0)*rhoUer_wall[i][j+1]*2.0*M_PI*r_tmp*dx;
        Im_x4_topOrfWall = Im_x4_topOrfWall +       rhoUmr_wall[i][j+1]*2.0*M_PI*r_tmp*dx;
        In_x4_topOrfWall = In_x4_topOrfWall +       rhoUnr_wall[i][j+1]*2.0*M_PI*r_tmp*dx;

        double engyi = 0.5*massi*(Uix[i][j]*Uix[i][j] + Uir[i][j+1]*Uir[i][j+1] + Uip[i][j]*Uip[i][j]) + 5.0/2.0*Boltz*Ti;
        Wi_x4_topOrfWall = Wi_x4_topOrfWall + engyi*rhoUir_wall[i][j+1]*2.0*M_PI*r_tmp*dx;
        We_x4_topOrfWall = We_x4_topOrfWall + Gr_wall[i][j+1]*2.0*M_PI*r_tmp*dx;
    }


    //upper-open-BC (x5)
    //------------------------------------
    for (int i=i_flc_bl[4][0];i<=i_flc_bl[4][1];i++){
        int j=j_flc_bl[4][1];

        double r_tmp = (r[j]+r[j+1])/2.0;

        double En = (+1.0)*Er[i][j+1];

        double rhoUir_Rr = rhoi[i][j]*Uir[i][j+1];
        double ji = (+1.0)*( e0)*rhoUir_Rr;
        double je = (+1.0)*(-e0)*rhoUer[i][j+1];
        double Ii = ji*2.0*M_PI*r_tmp*dx;
        double Ie = je*2.0*M_PI*r_tmp*dx;

        outputfile4<< i << ","<< j+1 << "," << x[i] << "," << r_tmp 
            << ","<< Ii*1000 << ","<< Ie*1000<< ","<< (Ii+Ie)*1000 
            << ","<< ji/1000 << ","<< je/1000<< ","<< (ji+je)/1000 
            << ","<< En << ","<< 0 << std::endl;

        Ii_x5_topOpen = Ii_x5_topOpen + ( e0)*rhoUir_Rr*2.0*M_PI*r_tmp*dx;
        Ie_x5_topOpen = Ie_x5_topOpen + (-e0)*rhoUer[i][j+1]*2.0*M_PI*r_tmp*dx;
        Im_x5_topOpen = Im_x5_topOpen +       rhoUmr[i][j+1]*2.0*M_PI*r_tmp*dx;
        In_x5_topOpen = In_x5_topOpen +       rhoUnr[i][j+1]*2.0*M_PI*r_tmp*dx;

        double engyi = 0.5*massi*(Uix[i][j]*Uix[i][j] + Uir[i][j+1]*Uir[i][j+1] + Uip[i][j]*Uip[i][j]) + 5.0/2.0*Boltz*Ti;
        Wi_x5_topOpen = Wi_x5_topOpen + engyi*rhoUir_Rr*2.0*M_PI*r_tmp*dx;
        We_x5_topOpen = We_x5_topOpen + Gr[i][j+1]*2.0*M_PI*r_tmp*dx;
    }
    //------------------------------------

    outputfile4.close();
    //outputfile6.close();

    /*
    //------------------------------------
    for (int j=j_flc_bl[4][0];j<=j_flc_bl[4][1];j++){
        int i=i_flc_bl[4][1];
        double vth = sqrt(8.0*Boltz*Ti/(M_PI*massi));
        double um = sqrt(2.0*Boltz*Ti/massi);
        double utx_R = Uix[i+1][j]/(um + 1e-100);

        double rhoUix_tmp = rhoi[i][j]*0.25*vth*(exp(-utx_R*utx_R) + utx_R*sqrt(M_PI)*(erf(utx_R) + 1.0));
        
        double delta_Ii = rhoUix_tmp  *2.0*M_PI*r[j]*dr*( e0);
        double delta_Ie = rhoUex[i+1][j]*2.0*M_PI*r[j]*dr*(-e0);

        Ii_Anode = Ii_Anode + delta_Ii;
        Ie_Anode = Ie_Anode + delta_Ie;
        I_Anode =  I_Anode + delta_Ii + delta_Ie;
    }
    //------------------------------------
    */

    Ii_Anode = Ii_z5_rightAnodeWall;
    Ie_Anode = Ie_z5_rightAnodeWall;
    I_Anode =  Ii_Anode + Ie_Anode;

    //Caluculate nozzle Current
    //------------------------------------
    Ii_Nozzle = 0.0;
    Ie_Nozzle = 0.0;
    I_Nozzle = 0.0;
    double Im_Nozzle = 0.0;
    double In_Nozzle = 0.0;

    for (int j=j_flc_bl[3][0];j<=j_flc_bl[3][1];j++){
        int i=i_flc_bl[1][1];
        double rhoi_tmp = (rhoi[i][j] + rhoi[i+1][j])/2.0;
        double rhoUix_tmp = rhoi_tmp*Uix[i+1][j];
        
        double delta_Ii = rhoUix_tmp  *2.0*M_PI*r[j]*dr*( e0);
        double delta_Ie = rhoUex[i][j]*2.0*M_PI*r[j]*dr*(-e0);
        double delta_Im = rhoUmx[i][j]*2.0*M_PI*r[j]*dr;
        double delta_In = rhoUnx[i][j]*2.0*M_PI*r[j]*dr;

        Ii_Nozzle = Ii_Nozzle + delta_Ii;
        Ie_Nozzle = Ie_Nozzle + delta_Ie;
        I_Nozzle =  I_Nozzle + delta_Ii + delta_Ie;

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
    
    outputfile3<< "No, wall, Ii (mA), Ie (mA), I (mA), Ndoti (pcl/s), Ndote (pcl/s), Ndotm (pcl/s), Ndotn (pcl/s)"<< std::endl;
    outputfile3<< "z0 , leftInWall, " << Ii_z0_leftInWall*1000 << "," << Ie_z0_leftInWall*1000 << "," << (Ii_z0_leftInWall + Ie_z0_leftInWall)*1000 << "," << 1.0/e0*Ii_z0_leftInWall << "," << 1.0/(-e0)*Ie_z0_leftInWall << "," << Im_z0_leftInWall << "," << In_z0_leftInWall  << std::endl;
    outputfile3<< "z1 , leftAntBaseWall, " << Ii_z1_leftAntBaseWall*1000 << "," << Ie_z1_leftAntBaseWall*1000 << "," << (Ii_z1_leftAntBaseWall + Ie_z1_leftAntBaseWall)*1000 << "," << 1.0/e0*Ii_z1_leftAntBaseWall << "," << 1.0/(-e0)*Ie_z1_leftAntBaseWall << "," << Im_z1_leftAntBaseWall << "," << In_z1_leftAntBaseWall  << std::endl;
    outputfile3<< "z2 , leftAntTopWall, " << Ii_z2_leftAntTopWall*1000 << "," << Ie_z2_leftAntTopWall*1000 << "," << (Ii_z2_leftAntTopWall + Ie_z2_leftAntTopWall)*1000 << "," << 1.0/e0*Ii_z2_leftAntTopWall << "," << 1.0/(-e0)*Ie_z2_leftAntTopWall << "," << Im_z2_leftAntTopWall << "," << In_z2_leftAntTopWall  << std::endl;
    outputfile3<< "z3 , rightInWall, " << Ii_z3_rightInWall*1000 << "," << Ie_z3_rightInWall*1000 << "," << (Ii_z3_rightInWall + Ie_z3_rightInWall)*1000 << "," << 1.0/e0*Ii_z3_rightInWall << "," << 1.0/(-e0)*Ie_z3_rightInWall << "," << Im_z3_rightInWall << "," << In_z3_rightInWall  << std::endl;
    outputfile3<< "z4 , leftOutWall, " << Ii_z4_leftOutWall*1000 << "," << Ie_z4_leftOutWall*1000 << "," << (Ii_z4_leftOutWall + Ie_z4_leftOutWall)*1000 << "," << 1.0/e0*Ii_z4_leftOutWall << "," << 1.0/(-e0)*Ie_z4_leftOutWall << "," << Im_z4_leftOutWall << "," << In_z4_leftOutWall  << std::endl;
    outputfile3<< "z5 , rightAnodeWall, " << Ii_z5_rightAnodeWall*1000 << "," << Ie_z5_rightAnodeWall*1000 << "," << (Ii_z5_rightAnodeWall + Ie_z5_rightAnodeWall)*1000 << "," << 1.0/e0*Ii_z5_rightAnodeWall << "," << 1.0/(-e0)*Ie_z5_rightAnodeWall << "," << Im_z5_rightAnodeWall << "," << In_z5_rightAnodeWall  << std::endl;
    outputfile3<< "x0 , botInWall, " << Ii_x0_botInWall*1000 << "," << Ie_x0_botInWall*1000 << "," << (Ii_x0_botInWall + Ie_x0_botInWall)*1000 << "," << 1.0/e0*Ii_x0_botInWall << "," << 1.0/(-e0)*Ie_x0_botInWall << "," << Im_x0_botInWall << "," << In_x0_botInWall  << std::endl;
    outputfile3<< "x1 , topInWall, " << Ii_x1_topInWall*1000 << "," << Ie_x1_topInWall*1000 << "," << (Ii_x1_topInWall + Ie_x1_topInWall)*1000 << "," << 1.0/e0*Ii_x1_topInWall << "," << 1.0/(-e0)*Ie_x1_topInWall << "," << Im_x1_topInWall << "," << In_x1_topInWall  << std::endl;
    outputfile3<< "x2 , botAntSideWall, " << Ii_x2_botAntSideWall*1000 << "," << Ie_x2_botAntSideWall*1000 << "," << (Ii_x2_botAntSideWall + Ie_x2_botAntSideWall)*1000 << "," << 1.0/e0*Ii_x2_botAntSideWall << "," << 1.0/(-e0)*Ie_x2_botAntSideWall << "," << Im_x2_botAntSideWall << "," << In_x2_botAntSideWall  << std::endl;
    outputfile3<< "x4 , topOrfWall, " << Ii_x4_topOrfWall*1000 << "," << Ie_x4_topOrfWall*1000 << "," << (Ii_x4_topOrfWall + Ie_x4_topOrfWall)*1000 << "," << 1.0/e0*Ii_x4_topOrfWall << "," << 1.0/(-e0)*Ie_x4_topOrfWall << "," << Im_x4_topOrfWall << "," << In_x4_topOrfWall  << std::endl;
    outputfile3<< "x5 , topOpen, " << Ii_x5_topOpen*1000 << "," << Ie_x5_topOpen*1000 << "," << (Ii_x5_topOpen + Ie_x5_topOpen)*1000 << "," << 1.0/e0*Ii_x5_topOpen << "," << 1.0/(-e0)*Ie_x5_topOpen << "," << Im_x5_topOpen << "," << In_x5_topOpen  << std::endl;
    outputfile3<< "*  , Nozzle , " << Ii_Nozzle*1000 << "," << Ie_Nozzle*1000 << "," << (Ii_Nozzle + Ie_Nozzle)*1000 << "," << 1.0/e0*Ii_Nozzle << "," << 1.0/(-e0)*Ie_Nozzle << "," << Im_Nozzle << "," << In_Nozzle  << std::endl;
    outputfile3<< "*  , sum, " << Ii_sum*1000 << "," << Ie_sum*1000 << "," << (Ii_sum + Ie_sum)*1000 << "," << 1.0/e0*Ii_sum << "," << 1.0/(-e0)*Ie_sum << "," << Im_sum << "," << In_sum  << std::endl;
    outputfile3.close();
    //------------------------------------

    //output energy sum
    //------------------------------------
    double Wi_sum = Wi_z0_leftInWall + Wi_z1_leftAntBaseWall + Wi_z3_rightInWall + Wi_z4_leftOutWall + Wi_z5_rightAnodeWall 
            + Wi_x0_botInWall+ Wi_x1_topInWall + Wi_x2_botAntSideWall + Wi_x4_topOrfWall + Wi_x5_topOpen;
    double We_sum = We_z0_leftInWall + We_z1_leftAntBaseWall + We_z3_rightInWall + We_z4_leftOutWall + We_z5_rightAnodeWall 
            + We_x0_botInWall+ We_x1_topInWall + We_x2_botAntSideWall + We_x4_topOrfWall + We_x5_topOpen;

    outputfile7<< "No, wall, Wi (J/s), We (J/s), W (J/s)"<< std::endl;
    outputfile7<< "z0 , leftInWall, "      << Wi_z0_leftInWall      << "," << We_z0_leftInWall      << "," << (Wi_z0_leftInWall      + We_z0_leftInWall)       << std::endl;
    outputfile7<< "z1 , leftAntBaseWall, " << Wi_z1_leftAntBaseWall << "," << We_z1_leftAntBaseWall << "," << (Wi_z1_leftAntBaseWall + We_z1_leftAntBaseWall)  << std::endl;
    outputfile7<< "z2 , leftAntTopWall, "  << Wi_z2_leftAntTopWall  << "," << We_z2_leftAntTopWall  << "," << (Wi_z2_leftAntTopWall  + We_z2_leftAntTopWall)   << std::endl;
    outputfile7<< "z3 , rightInWall, "     << Wi_z3_rightInWall     << "," << We_z3_rightInWall     << "," << (Wi_z3_rightInWall     + We_z3_rightInWall)      << std::endl;
    outputfile7<< "z4 , leftOutWall, "     << Wi_z4_leftOutWall     << "," << We_z4_leftOutWall     << "," << (Wi_z4_leftOutWall     + We_z4_leftOutWall)      << std::endl;
    outputfile7<< "z5 , rightAnodeWall, "  << Wi_z5_rightAnodeWall  << "," << We_z5_rightAnodeWall  << "," << (Wi_z5_rightAnodeWall  + We_z5_rightAnodeWall)   << std::endl;
    outputfile7<< "x0 , botInWall, "       << Wi_x0_botInWall       << "," << We_x0_botInWall       << "," << (Wi_x0_botInWall       + We_x0_botInWall)        << std::endl;
    outputfile7<< "x1 , topInWall, "       << Wi_x1_topInWall       << "," << We_x1_topInWall       << "," << (Wi_x1_topInWall       + We_x1_topInWall)        << std::endl;
    outputfile7<< "x2 , botAntSideWall, "  << Wi_x2_botAntSideWall  << "," << We_x2_botAntSideWall  << "," << (Wi_x2_botAntSideWall  + We_x2_botAntSideWall)   << std::endl;
    outputfile7<< "x4 , topOrfWall, "      << Wi_x4_topOrfWall      << "," << We_x4_topOrfWall      << "," << (Wi_x4_topOrfWall      + We_x4_topOrfWall)       << std::endl;
    outputfile7<< "x5 , topOpen, "         << Wi_x5_topOpen         << "," << We_x5_topOpen         << "," << (Wi_x5_topOpen         + We_x5_topOpen)          << std::endl;
    outputfile7<< "*  , sum,"              << Wi_sum                << "," << We_sum                << "," << (Wi_sum                + We_sum)                 << std::endl;
    outputfile7.close();
    //------------------------------------

    double direct_ionz_sum = 0.0;
    double stepwise_ionz_sum = 0.0;
    for (int iblock=0;iblock<2;iblock++){ 
        for (int i=i_flc_bl[iblock][0];i<=i_flc_bl[iblock][1];i++){ 
            for (int j=j_flc_bl[iblock][0];j<=j_flc_bl[iblock][1];j++){
                direct_ionz_sum = direct_ionz_sum + rhoe[i][j]*nu_ionz[i][j]*(dx*dr*2.0*M_PI*r[j]);
                stepwise_ionz_sum = stepwise_ionz_sum + rhoe[i][j]*nu_ionzStep[i][j]*(dx*dr*2.0*M_PI*r[j]);
            }
        }
    }

    std::cout << "direct_ionz_sum = " << direct_ionz_sum  
        << " stepwise_ionz_sum = " << stepwise_ionz_sum << " ratio = " << stepwise_ionz_sum/(stepwise_ionz_sum+direct_ionz_sum)
        << std::endl;

    nOut = nOut + 1;
}

//*****************************************************************
//**                                                             **
//**           void output                                       **
//**                                                             **
//*****************************************************************
void output_phase(){
    
    std::cout<<"output!"<<std::endl;
    std::string char_cnv,char1,char2,char_csv;
    int ndeg = 18;

    for(int n=0;n<ndeg;n++){
        double deg = 2.0*M_PI/ndeg*n;
        char1="results/result_deg";
        char2=std::to_string(n);
        char_csv=".csv";

        std::ofstream outputfile1(char1+char2+char_csv);
        //outputfile1<<"x,rho,Ui,rhoUix,Ue,rhoUex,rhoUex_E,rhoUex_D,E,phi,psi,rho_th,U_th,E_th,phi_th,psi_th,rate_ionize,divi,dive1,dive2,dive3, div_poisson, LHS, RHS,zero" << std::endl;
        //outputfile1<<"x,rhon,rho,rhoe(test),Uix,rhoUix,Uex,Uey,rhoUex,rhoUey,Te,heat_flux,Bz,E,phi,nu_m,nu_en,nu_ei,nu_ionz,nu_exc,nu_wall,nu_ano,mue_parae,mue_perp,Halle,jd,Id,rate_eloss,rate_ionize,divi,dive,divn,engy_LHS,engy_LHS1,engy_LHS2,engy_LHS3,engy_LHS4,engy_RHS,engy_RHS1,engy_RHS2,engy_RHS3,zero" << std::endl;
        outputfile1<<"i,j,x,r,Ex(R),Ex(I),|Ex|,arg(Ex),Er(R),Er(I),|Er|,arg(Er),Ep(R),Ep(I),|Ep|,arg(Ep),sigma_xx(R),sigma_xx(I),|sigma_xx|,arg(sigma_xx),sigma_rr(R),sigma_rr(I),|sigma_rr|,arg(sigma_rr),sigma_pp(R),sigma_pp(I),|sigma_pp|,arg(sigma_pp),sigma_xr(R),sigma_xr(I),|sigma_xr|,arg(sigma_xr),sigma_xp(R),sigma_xp(I),|sigma_xp|,arg(sigma_xp),sigma_rp(R),sigma_rp(I),|sigma_rp|,arg(sigma_rp),Jx(R),Jx(I),|Jx|,arg(Jx),Jr(R),Jr(I),|Jr|,arg(Jr),Jp(R),Jp(I),|Jp|,arg(Jp),Jx_exc(R),Jx_exc(I),|Jx_exc|,arg(Jx_exc),Jr_exc(R),Jr_exc(I),|Jr_exc|,arg(Jr_exc),Jp_exc(R),Jp_exc(I),|Jp_exc|,arg(Jp_exc),divE,jdgBnd_Ex,jdgBnd_Er,jdgBnd_Ep,Pabs,rho,Te,nu_m,Bx,Br,Bmag,Ap,epsr,zero" << std::endl;
        std::complex<double> EPS(1e-100,1e-100);
        for(int i=1;i<=ni;i++){
            for(int j=1;j<=nj;j++){
                std::complex<double> Ex_tmp = (E1x[i][j]+E1x[i+1][j])/2.0*exp(std::complex<double>(0,deg))*double(jdgBnd_Ep[i][j]);
                std::complex<double> Er_tmp = (E1r[i][j+1] + E1r[i][j])/2.0*exp(std::complex<double>(0,deg))*double(jdgBnd_Ep[i][j]);
                std::complex<double> Ep_tmp =  E1p[i][j]*exp(std::complex<double>(0,deg))*double(jdgBnd_Ep[i][j]);
                std::complex<double> Jx_exc_tmp = (J1x[i][j]+J1x[i+1][j])/2.0*exp(std::complex<double>(0,deg))*double(jdgBnd_Ep[i][j]);
                std::complex<double> Jr_exc_tmp = (J1r[i][j+1] + J1r[i][j])/2.0*exp(std::complex<double>(0,deg))*double(jdgBnd_Ep[i][j]);
                std::complex<double> Jp_exc_tmp =  J1p[i][j]*exp(std::complex<double>(0,deg))*double(jdgBnd_Ep[i][j]);
                double rL = (r[j-1] + r[j])/2.0;
                double rR = (r[j+1] + r[j])/2.0;
                std::complex<double> divE = ((rR*E1r[i][j+1]-rL*E1r[i][j])/(dr*r[j]) + (E1x[i+1][j]-E1x[i][j])/dx)*double(jdgBnd_Ep[i][j]);
                //double aveE = (std::abs(E1r[i][j+1]) + std::abs(E1r[i][j]) + std::abs(E1x[i+1][j]) + std::abs(E1x[i][j]))/4.0;
                //double aveE = sqrt(real(Ex_tmp*conj(Ex_tmp) + Er_tmp*conj(Er_tmp)));
                double aveE = 1.0/dx; //5000
                divE = divE/((aveE+1e-100)/dx);

                double nu_m_tmp = nu_m1[i][j];
                double rho_tmp  = rhoe[i][j];
                double Bx_tmp   = Bx[i][j];
                double Br_tmp   = Br[i][j];

                std::complex<double> nu_cmp(nu_m_tmp,omegam);
                std::complex<double> Hallcmpx = e0*Bx_tmp/(masse*nu_cmp + EPS);
                std::complex<double> Hallcmpr = e0*Br_tmp/(masse*nu_cmp + EPS);
                std::complex<double> HallcmpMag2 =  Hallcmpx*Hallcmpx + Hallcmpr*Hallcmpr;
                std::complex<double> coef = rho_tmp*e0*e0/(masse*nu_cmp + EPS)/(1.0 + HallcmpMag2 + EPS);
                std::complex<double> sigmaexr = coef*(Hallcmpx*Hallcmpr);
                std::complex<double> sigmaexp = coef*(Hallcmpr);
                std::complex<double> sigmaexx = coef*(1.0 + Hallcmpx*Hallcmpx);
                std::complex<double> sigmaerr = coef*(1.0 + Hallcmpr*Hallcmpr);
                std::complex<double> sigmaerp = coef*(-Hallcmpx);
                std::complex<double> sigmaerx = coef*(Hallcmpx*Hallcmpr);
                std::complex<double> sigmaepr = coef*(Hallcmpx);
                std::complex<double> sigmaepp = coef*(1.0);
                std::complex<double> sigmaepx = coef*(-Hallcmpr);
                
                //std::cout << nu_cmp<< " , "<< 1.0/nu_cmp << std::endl;


                std::complex<double> Jx_tmp = sigmaexr*Er_tmp + sigmaexp*Ep_tmp + sigmaexx*Ex_tmp;
                std::complex<double> Jr_tmp = sigmaerr*Er_tmp + sigmaerp*Ep_tmp + sigmaerx*Ex_tmp;
                std::complex<double> Jp_tmp = sigmaepr*Er_tmp + sigmaepp*Ep_tmp + sigmaepx*Ex_tmp;


                outputfile1<< i << ","<< j << "," << x[i]<< ","<< r[j]
                    << ","<< std::real(Ex_tmp) << "," << std::imag(Ex_tmp) << ","<< std::abs(Ex_tmp) << ","<< std::arg(Ex_tmp)
                    << ","<< std::real(Er_tmp) << "," << std::imag(Er_tmp) << ","<< std::abs(Er_tmp) << ","<< std::arg(Er_tmp)
                    << ","<< std::real(Ep_tmp) << "," << std::imag(Ep_tmp) << ","<< std::abs(Ep_tmp) << ","<< std::arg(Ep_tmp)
                    << ","<< std::real(sigmaexx) << "," << std::imag(sigmaexx) << ","<< std::abs(sigmaexx) << ","<< std::arg(sigmaexx)
                    << ","<< std::real(sigmaerr) << "," << std::imag(sigmaerr) << ","<< std::abs(sigmaerr) << ","<< std::arg(sigmaerr)
                    << ","<< std::real(sigmaepp) << "," << std::imag(sigmaepp) << ","<< std::abs(sigmaepp) << ","<< std::arg(sigmaepp)
                    << ","<< std::real(sigmaexr) << "," << std::imag(sigmaexr) << ","<< std::abs(sigmaexr) << ","<< std::arg(sigmaexr)
                    << ","<< std::real(sigmaexp) << "," << std::imag(sigmaexp) << ","<< std::abs(sigmaexp) << ","<< std::arg(sigmaexp)
                    << ","<< std::real(sigmaerp) << "," << std::imag(sigmaerp) << ","<< std::abs(sigmaerp) << ","<< std::arg(sigmaerp)
                    << ","<< std::real(Jx_tmp) << "," << std::imag(Jx_tmp) << ","<< std::abs(Jx_tmp) << ","<< std::arg(Jx_tmp)
                    << ","<< std::real(Jr_tmp) << "," << std::imag(Jr_tmp) << ","<< std::abs(Jr_tmp) << "," << std::arg(Jr_tmp)
                    << ","<< std::real(Jp_tmp) << "," << std::imag(Jp_tmp) << ","<< std::abs(Jp_tmp) << "," << std::arg(Jp_tmp)
                    << ","<< std::real(Jx_exc_tmp) << "," << std::imag(Jx_exc_tmp) << ","<< std::abs(Jx_exc_tmp) << "," << std::arg(Jx_exc_tmp)
                    << ","<< std::real(Jr_exc_tmp) << "," << std::imag(Jr_exc_tmp) << ","<< std::abs(Jr_exc_tmp) << "," << std::arg(Jr_exc_tmp)
                    << ","<< std::real(Jp_exc_tmp) << "," << std::imag(Jp_exc_tmp) << ","<< std::abs(Jp_exc_tmp) << "," << std::arg(Jp_exc_tmp)
                    << ","<< std::abs(divE)
                    << "," << jdgBnd_Ex[i][j]  << "," << jdgBnd_Er[i][j]  << "," << jdgBnd_Ep[i][j]
                    << ","<< Pabs[i][j]*jdgBnd_flc[i][j]
                    << ","<< rhoe[i][j]*jdgBnd_flc[i][j] << ","<< Te[i][j]*Boltz/e0*jdgBnd_flc[i][j]<< ","<< nu_m[i][j]*jdgBnd_flc[i][j]
                    << ","<< Bx[i][j]*jdgBnd_flc[i][j] << ","<< Br[i][j]*jdgBnd_flc[i][j]<< ","<< sqrt(Bx[i][j]*Bx[i][j]+Br[i][j]*Br[i][j])*jdgBnd_flc[i][j]<< ","<< Ap[i][j]
                    << ","<< epsr[i][j]*jdgBnd_Ep[i][j]
                    << "," << 0.0<< std::endl;
            }
        }

        outputfile1.close();
    }
}
