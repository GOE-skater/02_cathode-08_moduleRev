#pragma once
#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <time.h>
#include <vector>
#include <sstream>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <complex>

#include "constants.hpp"
#include "params.hpp"
#include "arrays.hpp"
#include "solverModule.hpp"

using namespace std;

//*****************************************************************
//**                                                             **
//**           class EmfieldModule                               **
//**                                                             **
//*****************************************************************
class EmfieldModule
{
    private:
    
    public:
        void solve_Microwave(Params &pm, GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr, GridK &gk, MicrowaveBC &mb);
        void update_energy_profile(Params &pm, GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr);
        
};
//*****************************************************************
//**                                                             **
//**           void solve_Microwave()                            **
//**                                                             **
//*****************************************************************
void EmfieldModule::solve_Microwave(Params &pm, GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr, GridK &gk, MicrowaveBC &mb){

    //テスト用 ※ 通常はコメントアウトすること
    //for (int i=0;i<pm.ni+2;i++){
    //    for (int j=0;j<pm.nj+2;j++){
    //        gc.rhoe[i][j] = 0.0;
    //        gc.nu_m1[i][j] = pm.nu_eff;
    //    }
    //}

    //Ex用 係数
    vector<vector<complex<double> > > aPx(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aWx(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aEx(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aNx(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aSx(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aNEx(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aNWx(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aSEx(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aSWx(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aEEx(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aWWx(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > bx(pm.ni+2,vector<complex<double> >  (pm.nj+2,0.0)); //係数

    //Er用 係数
    vector<vector<complex<double> > > aPr(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aWr(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aEr(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aNr(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aSr(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aNEr(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aNWr(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aSEr(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aSWr(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aNNr(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aSSr(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > br(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数

    //Ephi用 係数
    vector<vector<complex<double> > > aPp(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aWp(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aEp(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aNp(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aSp(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aNNp(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aSSp(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aEEp(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > aWWp(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数
    vector<vector<complex<double> > > bp(pm.ni+2,vector<complex<double> >(pm.nj+2,0.0)); //係数

    /******************** 係数作成 いらないところも含めて普通に全部やる方が速い ********************/
    complex<double> EPS(1e-100,1e-100); //微小複素数
    //for Ex
    for (int iblock=0;iblock<6;iblock++){ 
        for (int i=gx.i_flx_bl[iblock][0];i<=gx.i_flx_bl[iblock][1];i++){ 
            for (int j=gx.j_flx_bl[iblock][0];j<=gx.j_flx_bl[iblock][1];j++){
                double qL = (gc.r[j-1] + gc.r[j])/2.0/(gc.r[j]+1e-100);
                double qR = (gc.r[j+1] + gc.r[j])/2.0/(gc.r[j]+1e-100);
                
                double nu_m_tmp = (gc.nu_m1[i][j] + gc.nu_m1[i-1][j])/2.0;
                double rho_tmp  = (gc.rhoe[i][j]  + gc.rhoe[i-1][j] )/2.0;
                double Bx_tmp   = (gc.Bx[i][j]   + gc.Bx[i-1][j]  )/2.0;
                double Br_tmp   = (gc.Br[i][j]   + gc.Br[i-1][j]  )/2.0;

                complex<double> nu_cmp(nu_m_tmp,pm.omegam);
                complex<double> Hallcmpx = ph::e0*Bx_tmp/(pm.masse*nu_cmp + EPS);
                complex<double> Hallcmpr = ph::e0*Br_tmp/(pm.masse*nu_cmp + EPS);
                complex<double> HallcmpMag2 =  Hallcmpx*Hallcmpx + Hallcmpr*Hallcmpr;
                complex<double> coef = rho_tmp*ph::e0*ph::e0/(pm.masse*nu_cmp + EPS)/(1.0 + HallcmpMag2 + EPS);
                complex<double> sigmaexr = coef*(Hallcmpx*Hallcmpr);
                complex<double> sigmaexp = coef*(Hallcmpr);
                complex<double> sigmaexx = coef*(1.0 + Hallcmpx*Hallcmpx);
                complex<double> iomegaMu0(0.0,pm.omegam*ph::mu0);

                aPx[i][j] = 2.0*pm.dx*pm.dx + 2.0*pm.dr*pm.dr
                    - ph::eps0*gc.epsr[i][j]*ph::mu0*pm.omegam*pm.omegam*pm.dx*pm.dx*pm.dr*pm.dr + iomegaMu0*sigmaexx*pm.dx*pm.dx*pm.dr*pm.dr;
                aEx[i][j] = pm.dr*pm.dr;
                aWx[i][j] = pm.dr*pm.dr;
                aNx[i][j] = qR*pm.dx*pm.dx;
                aSx[i][j] = qL*pm.dx*pm.dx;
                bx[i][j] = complex<double>(0,-pm.omegam*ph::mu0*pm.dx*pm.dx*pm.dr*pm.dr)*gx.J1x[i][j];

                //Epの要素
                aEEx[i][j] = -iomegaMu0*sigmaexp*0.5*pm.dx*pm.dx*pm.dr*pm.dr;
                aWWx[i][j] = -iomegaMu0*sigmaexp*0.5*pm.dx*pm.dx*pm.dr*pm.dr;

                //Erの要素
                aNEx[i][j] = -iomegaMu0*sigmaexr*0.25*qR*pm.dx*pm.dx*pm.dr*pm.dr;
                aNWx[i][j] = -iomegaMu0*sigmaexr*0.25*qR*pm.dx*pm.dx*pm.dr*pm.dr;
                aSEx[i][j] = -iomegaMu0*sigmaexr*0.25*qL*pm.dx*pm.dx*pm.dr*pm.dr;
                aSWx[i][j] = -iomegaMu0*sigmaexr*0.25*qL*pm.dx*pm.dx*pm.dr*pm.dr;

                //cout << i <<","<< j 
                //    <<","<< aPx[i][j] << "," << aEx[i][j]<< "," << aWx[i][j]<< "," << aNx[i][j]<< "," << aNx[i][j]
                //    //<<","<< aEEx[i][j] << "," << aWWx[i][j]
                //    //<< "," << aNEx[i][j]<< "," << aNWx[i][j]
                //    //<< "," << aSEx[i][j]<< "," << aSWx[i][j]
                //    << endl;
            }
        }
    }

    //for Er
    for (int iblock=0;iblock<6;iblock++){ 
        for (int i=gr.i_flr_bl[iblock][0];i<=gr.i_flr_bl[iblock][1];i++){ 
            for (int j=gr.j_flr_bl[iblock][0];j<=gr.j_flr_bl[iblock][1];j++){
                double r_tmp = (gc.r[j-1] + gc.r[j])/2.0;
                double qL = gc.r[j-1]/(r_tmp+1e-100);
                double qR = gc.r[j]  /(r_tmp+1e-100);

                double nu_m_tmp = (gc.nu_m1[i][j] + gc.nu_m1[i][j-1])/2.0;
                double rho_tmp  = (gc.rhoe[i][j]  + gc.rhoe[i][j-1] )/2.0;
                double Bx_tmp   = (gc.r[j]*gc.Bx[i][j]   + gc.r[j-1]*gc.Bx[i][j-1]  )/(2.0*r_tmp);
                double Br_tmp   = (gc.r[j]*gc.Br[i][j]   + gc.r[j-1]*gc.Br[i][j-1]  )/(2.0*r_tmp);

                complex<double> nu_cmp(nu_m_tmp,pm.omegam);
                complex<double> Hallcmpx = ph::e0*Bx_tmp/(pm.masse*nu_cmp + EPS);
                complex<double> Hallcmpr = ph::e0*Br_tmp/(pm.masse*nu_cmp + EPS);
                complex<double> HallcmpMag2 =  Hallcmpx*Hallcmpx + Hallcmpr*Hallcmpr;
                complex<double> coef = rho_tmp*ph::e0*ph::e0/(pm.masse*nu_cmp + EPS)/(1.0 + HallcmpMag2 + EPS);
                complex<double> sigmaerr = coef*(1.0+Hallcmpr*Hallcmpr);
                complex<double> sigmaerp = coef*(-Hallcmpx);
                complex<double> sigmaerx = coef*(Hallcmpx*Hallcmpr);
                complex<double> iomegaMu0(0.0,pm.omegam*ph::mu0);

                aPr[i][j] = 2.0*pm.dx*pm.dx + 2.0*pm.dr*pm.dr + pm.dx*pm.dx*pm.dr*pm.dr/(r_tmp*r_tmp)
                    - ph::eps0*gc.epsr[i][j]*ph::mu0*pm.omegam*pm.omegam*pm.dx*pm.dx*pm.dr*pm.dr + iomegaMu0*sigmaerr*pm.dx*pm.dx*pm.dr*pm.dr;
                aEr[i][j] = pm.dr*pm.dr;
                aWr[i][j] = pm.dr*pm.dr;
                aNr[i][j] = qR*pm.dx*pm.dx;
                aSr[i][j] = qL*pm.dx*pm.dx;
                br[i][j] = complex<double>(0,-pm.omegam*ph::mu0*pm.dx*pm.dx*pm.dr*pm.dr)*gr.J1r[i][j];

                //Epの要素
                aNNr[i][j] = -iomegaMu0*sigmaerp*0.5*qR*pm.dx*pm.dx*pm.dr*pm.dr;
                aSSr[i][j] = -iomegaMu0*sigmaerp*0.5*qL*pm.dx*pm.dx*pm.dr*pm.dr;

                //Exの要素
                aNEr[i][j] = -iomegaMu0*sigmaerx*0.25*qR*pm.dx*pm.dx*pm.dr*pm.dr;
                aNWr[i][j] = -iomegaMu0*sigmaerx*0.25*qR*pm.dx*pm.dx*pm.dr*pm.dr;
                aSEr[i][j] = -iomegaMu0*sigmaerx*0.25*qL*pm.dx*pm.dx*pm.dr*pm.dr;
                aSWr[i][j] = -iomegaMu0*sigmaerx*0.25*qL*pm.dx*pm.dx*pm.dr*pm.dr;

                //std::cout << i << ","<< j 
                //    << ","<<aPr[i][j] << ","<<br[i][j] 
                //    << "," << 2.0*pm.dx*pm.dx 
                //    << "," << 2.0*pm.dr*pm.dr 
                //    << "," << pm.dx*pm.dx*pm.dr*pm.dr/(r_tmp*r_tmp)
                //    << std::endl;

                //std::cout << i << ","<< j 
                //    << ","<< br[i][j] << ","<< gr.J1r[i][j] 
                //    << "," << -pm.omegam*ph::mu0*pm.dx*pm.dx*pm.dr*pm.dr
                //    << "," << pm.omegam
                //    << "," << ph::mu0
                //    << "," << pm.dx
                //    << "," << pm.dr
                //    << std::endl;

            }
        }
    }

    //for Ephi
    for (int iblock=0;iblock<6;iblock++){ 
        for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
            for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                double qL = (gc.r[j-1] + gc.r[j])/2.0/(gc.r[j]+1e-100);
                double qR = (gc.r[j+1] + gc.r[j])/2.0/(gc.r[j]+1e-100);

                double nu_m_tmp = gc.nu_m1[i][j];
                double rho_tmp  = gc.rhoe[i][j];
                double Bx_tmp   = gc.Bx[i][j];
                double Br_tmp   = gc.Br[i][j];

                complex<double> nu_cmp(nu_m_tmp,pm.omegam);
                complex<double> Hallcmpx = ph::e0*Bx_tmp/(pm.masse*nu_cmp + EPS);
                complex<double> Hallcmpr = ph::e0*Br_tmp/(pm.masse*nu_cmp + EPS);
                complex<double> HallcmpMag2 =  Hallcmpx*Hallcmpx + Hallcmpr*Hallcmpr;
                complex<double> coef = rho_tmp*ph::e0*ph::e0/(pm.masse*nu_cmp + EPS)/(1.0 + HallcmpMag2 + EPS);
                complex<double> sigmaepr = coef*(Hallcmpx);
                complex<double> sigmaepp = coef*(1.0);
                complex<double> sigmaepx = coef*(-Hallcmpr);
                complex<double> iomegaMu0(0.0,pm.omegam*ph::mu0);

                aPp[i][j] = 2.0*pm.dx*pm.dx + 2.0*pm.dr*pm.dr + pm.dx*pm.dx*pm.dr*pm.dr/(gc.r[j]*gc.r[j])
                    - ph::eps0*gc.epsr[i][j]*ph::mu0*pm.omegam*pm.omegam*pm.dx*pm.dx*pm.dr*pm.dr + iomegaMu0*sigmaepp*pm.dx*pm.dx*pm.dr*pm.dr;
                aEp[i][j] = pm.dr*pm.dr;
                aWp[i][j] = pm.dr*pm.dr;
                aNp[i][j] = qR*pm.dx*pm.dx;
                aSp[i][j] = qL*pm.dx*pm.dx;
                bp[i][j] = complex<double>(0,-pm.omegam*ph::mu0*pm.dx*pm.dx*pm.dr*pm.dr)*gc.J1p[i][j];

                //Exの要素
                aEEp[i][j] = -iomegaMu0*sigmaepx*0.5*pm.dx*pm.dx*pm.dr*pm.dr;
                aWWp[i][j] = -iomegaMu0*sigmaepx*0.5*pm.dx*pm.dx*pm.dr*pm.dr;

                //Erの要素
                aNNp[i][j] = -iomegaMu0*sigmaepr*0.5*qR*pm.dx*pm.dx*pm.dr*pm.dr;
                aSSp[i][j] = -iomegaMu0*sigmaepr*0.5*qL*pm.dx*pm.dx*pm.dr*pm.dr;

            }
        }
    }

    /******************** 境界条件設定 (係数修正) ********************/
    //********* Ex境界条件 (左) *********
    for (int k=0;k<mb.iBndWx.size();k++){
        int i = mb.iBndWx[k];
        int j = mb.jBndWx[k];

        if(mb.sBndWx[k] == 0){ //ディリクレ
            aWx[i][j] = 0.0;
        }else if(mb.sBndWx[k] == 2){ //ガウス
            double qL = (gc.r[j-1] + gc.r[j])/2.0/(gc.r[j]+1e-100);
            double qR = (gc.r[j+1] + gc.r[j])/2.0/(gc.r[j]+1e-100);

            aPx[i][j] = aPx[i][j] - aWx[i][j];
            aNWx[i][j] = aNWx[i][j] + aWx[i][j]*pm.dx/pm.dr*qR;
            aSWx[i][j] = aSWx[i][j] - aWx[i][j]*pm.dx/pm.dr*qL;
            aWx[i][j] = 0.0;
        }else if(mb.sBndWx[k] == 3){ //凹角 下側壁
            double qL = (gc.r[j-1] + gc.r[j])/2.0/(gc.r[j]+1e-100);
            double qR = (gc.r[j+1] + gc.r[j])/2.0/(gc.r[j]+1e-100);
            double deno = 1.0 + pm.dx/pm.dr*qL;

            aPx[i][j] = aPx[i][j] - aWx[i][j]/deno;
            aNWx[i][j] = aNWx[i][j] + aWx[i][j]*(pm.dx/pm.dr*qR)/deno;
            aWx[i][j] = 0.0;
        }else if(mb.sBndWx[k] == 4){ //凹角 上側壁
            double qL = (gc.r[j-1] + gc.r[j])/2.0/(gc.r[j]+1e-100);
            double qR = (gc.r[j+1] + gc.r[j])/2.0/(gc.r[j]+1e-100);
            double deno = 1.0 + pm.dx/pm.dr*qR;

            aPx[i][j] = aPx[i][j] - aWx[i][j]/deno;
            aSWx[i][j] = aSWx[i][j] - aWx[i][j]*(pm.dx/pm.dr*qL)/deno;
            aWx[i][j] = 0.0;
        }
    }
    //********* Ex境界条件 (右) *********
    for (int k=0;k<mb.iBndEx.size();k++){
        int i = mb.iBndEx[k];
        int j = mb.jBndEx[k];

        if(mb.sBndEx[k] == 2){ //ガウス
            double qL = (gc.r[j-1] + gc.r[j])/2.0/(gc.r[j]+1e-100);
            double qR = (gc.r[j+1] + gc.r[j])/2.0/(gc.r[j]+1e-100);

            aPx[i][j] = aPx[i][j] - aEx[i][j];
            aNEx[i][j] = aNEx[i][j] - aEx[i][j]*pm.dx/pm.dr*qR;
            aSEx[i][j] = aSEx[i][j] + aEx[i][j]*pm.dx/pm.dr*qL;
            aEx[i][j] = 0.0;
        }else if(mb.sBndEx[k] == 4){ //凹角 上側壁
            double qL = (gc.r[j-1] + gc.r[j])/2.0/(gc.r[j]+1e-100);
            double qR = (gc.r[j+1] + gc.r[j])/2.0/(gc.r[j]+1e-100);
            double deno = 1.0 + pm.dx/pm.dr*qR;
            
            aPx[i][j] = aPx[i][j] - aEx[i][j]/deno;
            aSEx[i][j] = aSEx[i][j] + aEx[i][j]*(pm.dx/pm.dr*qL)/deno;
            aEx[i][j] = 0.0;
        }
    }
    //********* Ex境界条件 (下) *********
    for (int k=0;k<mb.iBndSx.size();k++){
        int i = mb.iBndSx[k];
        int j = mb.jBndSx[k];
        if(mb.sBndSx[k] == 0){ //ディリクレ
            aPx[i][j] = aPx[i][j] + gc.r[j]/gc.r[j-1]*aSx[i][j];
            aSx[i][j] = 0.0;
        }else if(mb.sBndSx[k] == 6){ //凸角 右側Open
            double rR = (gc.r[j+1] + gc.r[j])/2.0;
            double rL = (gc.r[j] + gc.r[j-1])/2.0;

            aPx[i][j] = aPx[i][j] - aSx[i][j]*( pm.dr/pm.dx*gc.r[j]/rL);
            aWx[i][j] = aWx[i][j] + aSx[i][j]*(-pm.dr/pm.dx*gc.r[j]/rL);
            aNWx[i][j] = aNWx[i][j] + aSx[i][j]*(rR/rL);
            aSx[i][j] = 0.0;
        }
    }
    //********* Ex境界条件 (上) *********
    for (int k=0;k<mb.iBndNx.size();k++){
        int i = mb.iBndNx[k];
        int j = mb.jBndNx[k];
        if(mb.sBndNx[k] == 0){ //ディリクレ
            aPx[i][j] = aPx[i][j] + gc.r[j]/gc.r[j+1]*aNx[i][j];
            aNx[i][j] = 0.0;
        }if(mb.sBndNx[k] == 1){ //開放
            double c0 = sqrt(1.0/(ph::eps0*gc.epsr[i][j]*ph::mu0));
            double r_tmp = (gc.r[j+1] + gc.r[j])/2.0;
            double qL = gc.r[j]/(r_tmp+1e-100);
            double qR = gc.r[j+1]  /(r_tmp+1e-100);
            double deno = pow(2.0*c0*sqrt(qR),2) + pow(pm.omegam*pm.dr*qR,2);
            complex<double> R = complex<double>
                (pow(2.0*c0,2)*sqrt(qL*qR) - pow(pm.omegam*pm.dr,2)*qL*qR
                ,-2.0*c0*pm.omegam*pm.dr*(qL*sqrt(qR)+qR*sqrt(qL)));
            R = R/(deno+1e-100);

            aPx[i][j] = aPx[i][j] - aNx[i][j]*R;
            aNx[i][j] = 0.0;
        }else if(mb.sBndNx[k] == 5){ //凸角 左側Open (*)
            double rR = (gc.r[j+1] + gc.r[j])/2.0;
            double rL = (gc.r[j] + gc.r[j-1])/2.0;

            aPx[i][j] = aPx[i][j] - aNx[i][j]*( pm.dr/pm.dx*gc.r[j]/rR);
            aEx[i][j] = aEx[i][j] + aNx[i][j]*(-pm.dr/pm.dx*gc.r[j]/rR);
            aSEx[i][j] = aSEx[i][j] + aNx[i][j]*( rL/rR);
            aNx[i][j] = 0.0;
        }else if(mb.sBndNx[k] == 6){ //凸角 右側Open
            double rR = (gc.r[j+1] + gc.r[j])/2.0;
            double rL = (gc.r[j] + gc.r[j-1])/2.0;

            aPx[i][j] = aPx[i][j] - aNx[i][j]*( pm.dr/pm.dx*gc.r[j]/rR);
            aWx[i][j] = aWx[i][j] + aNx[i][j]*(-pm.dr/pm.dx*gc.r[j]/rR);
            aSWx[i][j] = aSWx[i][j] + aNx[i][j]*(-rL/rR);
            aNx[i][j] = 0.0;
        }
    }

    //********* Er境界条件 (左) *********
    for (int k=0;k<mb.iBndWr.size();k++){
        int i = mb.iBndWr[k];
        int j = mb.jBndWr[k];

        if(mb.sBndWr[k] == 0){ //ディリクレ
            aPr[i][j] = aPr[i][j] + aWr[i][j];
            aWr[i][j] = 0.0;
        }else if(mb.sBndWr[k] == 1){ //開放
            double c0 = sqrt(1.0/(ph::eps0*gc.epsr[i][j]*ph::mu0));
            double deno = pow(2.0*c0,2) + pow(pm.omegam*pm.dx,2);
            complex<double> R = complex<double>(pow(2.0*c0,2) - pow(pm.omegam*pm.dx,2),-4.0*c0*pm.omegam*pm.dx);
            R = R/(deno+1e-100);
            
            aPr[i][j] = aPr[i][j] - aWr[i][j]*R;
            aWr[i][j] = 0.0;
        }else if(mb.sBndWr[k] == 5){ //凸角 下側Open
            double rR = (gc.r[j+1] + gc.r[j])/2.0;
            double rL = (gc.r[j] + gc.r[j-1])/2.0;
            
            aPr[i][j] = aPr[i][j] - aWr[i][j]*( pm.dx/pm.dr*rL/gc.r[j]);
            aNr[i][j] = aNr[i][j] + aWr[i][j]*(-pm.dx/pm.dr*rR/gc.r[j]);
            aNEr[i][j] = aNEr[i][j] + aWr[i][j]*(-1.0);
            aWr[i][j] = 0.0;
        }else if(mb.sBndWr[k] == 6){ //凸角 上側Open
            double rR = (gc.r[j] + gc.r[j-1])/2.0;
            double rL = (gc.r[j-1] + gc.r[j-2])/2.0;

            aPr[i][j] = aPr[i][j] - aWr[i][j]*( pm.dx/pm.dr*rR/gc.r[j-1]);
            aSr[i][j] = aSr[i][j] + aWr[i][j]*(-pm.dx/pm.dr*rL/gc.r[j-1]);
            aSEr[i][j] = aSEr[i][j] + aWr[i][j]*( 1.0);
            aWr[i][j] = 0.0;
        }
    }
    //********* Er境界条件 (右) *********
    for (int k=0;k<mb.iBndEr.size();k++){
        int i = mb.iBndEr[k];
        int j = mb.jBndEr[k];

        if(mb.sBndEr[k] == 0){ //ディリクレ
            aPr[i][j] = aPr[i][j] + aEr[i][j];
            aEr[i][j] = 0.0;
        }else if(mb.sBndEr[k] == 5){ //凸角 下側Open (*)
            double rR = (gc.r[j+1] + gc.r[j])/2.0;
            double rL = (gc.r[j] + gc.r[j-1])/2.0;

            aPr[i][j] = aPr[i][j] - aEr[i][j]*( pm.dx/pm.dr*rL/gc.r[j]);
            aNr[i][j] = aNr[i][j] + aEr[i][j]*(-pm.dx/pm.dr*rR/gc.r[j]);
            aNWr[i][j] = aNWr[i][j] + aEr[i][j]*( 1.0);
            aEr[i][j] = 0.0;
        }
    }
    //********* Er境界条件 (下) *********
    for (int k=0;k<mb.iBndSr.size();k++){
        int i = mb.iBndSr[k];
        int j = mb.jBndSr[k];

        if(mb.sBndSr[k] == 2){ //ガウス
            double rR = (gc.r[j] + gc.r[j-1])/2.0;
            double rL = (gc.r[j-1] + gc.r[j-2])/2.0;
            aPr[i][j] = aPr[i][j] - aSr[i][j]*rR/rL;
            aSEr[i][j] = aSEr[i][j] + aSr[i][j]*pm.dr/pm.dx*gc.r[j-1]/rL;
            aSWr[i][j] = aSWr[i][j] - aSr[i][j]*pm.dr/pm.dx*gc.r[j-1]/rL;
            aSr[i][j] = 0.0;
        }else if(mb.sBndSr[k] == 3){ //凹角 左側壁
            double rR = (gc.r[j] + gc.r[j-1])/2.0;
            double rL = (gc.r[j-1] + gc.r[j-2])/2.0;
            double deno = pm.dr/pm.dx*gc.r[j-1]/rL;
            if(gc.r[j] < gc.j_flc_bl[0][0]){ //マイクロ波の入り口部の特別処理
                deno = 1.0;
            }

            aPr[i][j] = aPr[i][j] - aSr[i][j]*rR/rL/deno;
            aNEr[i][j] = aNEr[i][j] + aSr[i][j]*pm.dr/pm.dx*gc.r[j-1]/rL/deno;
            aSr[i][j] = 0.0;
        }
    }
    //********* Er境界条件 (上) *********
    for (int k=0;k<mb.iBndNr.size();k++){
        int i = mb.iBndNr[k];
        int j = mb.jBndNr[k];

        if(mb.sBndNr[k] == 2){ //ガウス
            double rR = (gc.r[j+1] + gc.r[j])/2.0;
            double rL = (gc.r[j] + gc.r[j-1])/2.0;
            aPr[i][j] = aPr[i][j] - aNr[i][j]*rL/rR;
            aNEr[i][j] = aNEr[i][j] - aNr[i][j]*pm.dr/pm.dx*gc.r[j]/rR;
            aNWr[i][j] = aNWr[i][j] + aNr[i][j]*pm.dr/pm.dx*gc.r[j]/rR;
            aNr[i][j] = 0.0;
        }if(mb.sBndNr[k] == 3){ //凹角 左側壁
            double rR = (gc.r[j+1] + gc.r[j])/2.0;
            double rL = (gc.r[j] + gc.r[j-1])/2.0;
            double deno = pm.dr/pm.dx*gc.r[j]/rR;

            aPr[i][j] = aPr[i][j] - aNr[i][j]*rL/rR/deno;
            aSEr[i][j] = aSEr[i][j] - aNr[i][j]*pm.dr/pm.dx*gc.r[j]/rR/deno;
            aNr[i][j] = 0.0;
        }if(mb.sBndNr[k] == 4){ //凹角 右側壁
            double rR = (gc.r[j+1] + gc.r[j])/2.0;
            double rL = (gc.r[j] + gc.r[j-1])/2.0;
            double deno = pm.dr/pm.dx*gc.r[j]/rR;

            aPr[i][j] = aPr[i][j] - aNr[i][j]*rL/rR/deno;
            aSWr[i][j] = aSWr[i][j] + aNr[i][j]*pm.dr/pm.dx*gc.r[j]/rR/deno;
            aNr[i][j] = 0.0;
        }
    }

    //********* Ep境界条件 (左) *********
    for (int k=0;k<mb.iBndWp.size();k++){
        int i = mb.iBndWp[k];
        int j = mb.jBndWp[k];

        if(mb.sBndWp[k] == 0){ //ディリクレ
            aPp[i][j] = aPp[i][j] + aWp[i][j];
            aWp[i][j] = 0.0;
        }else if(mb.sBndWp[k] == 1){ //開放
            double c0 = sqrt(1.0/(ph::eps0*gc.epsr[i][j]*ph::mu0));
            double deno = pow(2.0*c0,2) + pow(pm.omegam*pm.dx,2);
            complex<double> R = complex<double>
                (pow(2.0*c0,2) - pow(pm.omegam*pm.dx,2),-4.0*c0*pm.omegam*pm.dx);
            R = R/(deno+1e-100);
            aPp[i][j] = aPp[i][j] - aWp[i][j]*R;
            aWp[i][j] = 0.0;
        }
    }
    //********* Ep境界条件 (右) *********
    for (int k=0;k<mb.iBndEp.size();k++){
        int i = mb.iBndEp[k];
        int j = mb.jBndEp[k];

        if(mb.sBndEp[k] == 0){ //ディリクレ
            aPp[i][j] = aPp[i][j] + aEp[i][j];
            aEp[i][j] = 0.0;
        }
    }
    //********* Ep境界条件 (下) *********
    for (int k=0;k<mb.iBndSp.size();k++){
        int i = mb.iBndSp[k];
        int j = mb.jBndSp[k];

        if(mb.sBndSp[k] == 0){ //ディリクレ
            aPp[i][j] = aPp[i][j] + gc.r[j]/gc.r[j-1]*aSp[i][j];
            aSp[i][j] = 0.0;
        }
    }
    //********* Ep境界条件 (上) *********
    for (int k=0;k<mb.iBndNp.size();k++){
        int i = mb.iBndNp[k];
        int j = mb.jBndNp[k];

        if(mb.sBndNp[k] == 0){ //ディリクレ
            aPp[i][j] = aPp[i][j] + gc.r[j]/gc.r[j+1]*aNp[i][j];
            aNp[i][j] = 0.0;
        }else if(mb.sBndNp[k] == 1){ //開放
            double c0 = sqrt(1.0/(ph::eps0*gc.epsr[i][j]*ph::mu0));
            double r_tmp = (gc.r[j+1] + gc.r[j])/2.0;
            double qL = gc.r[j]/(r_tmp+1e-100);
            double qR = gc.r[j+1]  /(r_tmp+1e-100);
            double deno = pow(2.0*c0*sqrt(qR),2) + pow(pm.omegam*pm.dr*qR,2);
            complex<double> R = complex<double>
                (pow(2.0*c0,2)*sqrt(qL*qR) - pow(pm.omegam*pm.dr,2)*qL*qR
                ,-2.0*c0*pm.omegam*pm.dr*(qL*sqrt(qR)+qR*sqrt(qL)));
            R = R/(deno+1e-100);

            aPp[i][j] = aPp[i][j] - aNp[i][j]*R;
            aNp[i][j] = 0.0;
        }
    }

    /******************** LUソルバーの準備 ********************/
    Eigen::SparseMatrix<complex<double> > A(pm.nk, pm.nk);
    A.reserve(Eigen::VectorXi::Constant(pm.nk,11)); //ここの数字を変えて帯域幅を確保する
    Eigen::VectorXcd b(pm.nk);
    Eigen::VectorXcd xv(pm.nk);

    //Ex-足し込み
    for (int k=0;k<pm.nkx;k++){

        int i = gk.ikx[k];
        int j = gk.jkx[k];

        double kE = gx.kx[i+1][j];
        double kW = gx.kx[i-1][j];
        double kN = gx.kx[i][j+1];
        double kS = gx.kx[i][j-1];

        double kNE = gr.kr[i][j+1];
        double kSE = gr.kr[i][j];
        double kNW = gr.kr[i-1][j+1];
        double kSW = gr.kr[i-1][j];

        double kEE = gc.kp[i][j];
        double kWW = gc.kp[i-1][j];

        //Ex用
        A.insert(k, k)      =   aPx[i][j];
        b[k]                =   bx[i][j];
        if(kE!=-1){
            A.insert(k, kE) = -aEx[i][j];
        }
        if(kW!=-1){
            A.insert(k, kW) = -aWx[i][j];
        }
        if(kN!=-1){
            A.insert(k, kN) = -aNx[i][j];
        }
        if(kS!=-1){
            A.insert(k, kS) = -aSx[i][j];
        }
        //Er用
        if(kNE!=-1){
            A.insert(k, pm.nkx+kNE) = -aNEx[i][j];
        }
        if(kSE!=-1){
            A.insert(k, pm.nkx+kSE) = -aSEx[i][j];
        }
        if(kNW!=-1){
            A.insert(k, pm.nkx+kNW) = -aNWx[i][j];
        }
        if(kSW!=-1){
            A.insert(k, pm.nkx+kSW) = -aSWx[i][j];
        }
        //Ep用
        if(kEE!=-1){
            A.insert(k, pm.nkx+pm.nkr+kEE) = -aEEx[i][j];
        }
        if(kWW!=-1){
            A.insert(k, pm.nkx+pm.nkr+kWW) = -aWWx[i][j];
        }

    }
    
    //Er-足し込み
    for (int k=0;k<pm.nkr;k++){
        int i = gk.ikr[k];
        int j = gk.jkr[k];

        double kE  = gr.kr[i+1][j];
        double kW  = gr.kr[i-1][j];
        double kN  = gr.kr[i][j+1];
        double kS  = gr.kr[i][j-1];

        double kNE = gx.kx[i+1][j];
        double kSE = gx.kx[i+1][j-1];
        double kNW = gx.kx[i][j];
        double kSW = gx.kx[i][j-1];

        double kNN = gc.kp[i][j];
        double kSS = gc.kp[i][j-1];

        //Er用
        A.insert(pm.nkx+k, pm.nkx+k)      =   aPr[i][j];
        b[pm.nkx+k]                    =   br[i][j];
        if(kN!=-1){
            A.insert(pm.nkx+k, pm.nkx+kN) = -aNr[i][j];
        }
        if(kS!=-1){
            A.insert(pm.nkx+k, pm.nkx+kS) = -aSr[i][j];
        }
        if(kE!=-1){
            A.insert(pm.nkx+k, pm.nkx+kE) = -aEr[i][j];
        }
        if(kW!=-1){
            A.insert(pm.nkx+k, pm.nkx+kW) = -aWr[i][j];
        }
        //Ex用
        if(kNE!=-1){
            A.insert(pm.nkx+k, kNE) = -aNEr[i][j];
        }
        if(kSE!=-1){
            A.insert(pm.nkx+k, kSE) = -aSEr[i][j];
        }
        if(kNW!=-1){
            A.insert(pm.nkx+k, kNW) = -aNWr[i][j];
        }
        if(kSW!=-1){
            A.insert(pm.nkx+k, kSW) = -aSWr[i][j];
        }
        //Ep用
        if(kNN!=-1){
            A.insert(pm.nkx+k, pm.nkx + pm.nkr+kNN) = -aNNr[i][j];
        }
        if(kSS!=-1){
            A.insert(pm.nkx+k, pm.nkx+pm.nkr+kSS) = -aSSr[i][j];
        }
    }

    //Ep-足し込み
    for (int k=0;k<pm.nkp;k++){
        int i = gk.ikp[k];
        int j = gk.jkp[k];

        double kE  = gc.kp[i+1][j];
        double kW  = gc.kp[i-1][j];
        double kN  = gc.kp[i][j+1];
        double kS  = gc.kp[i][j-1];

        double kEE  = gx.kx[i+1][j];
        double kWW  = gx.kx[i][j];

        double kNN  = gr.kr[i][j+1];
        double kSS  = gr.kr[i][j];

        //Ep用
        A.insert(pm.nkx+pm.nkr+k, pm.nkx+pm.nkr+k)      =   aPp[i][j];
        b[pm.nkx+pm.nkr+k]                    =   bp[i][j];
        if(kN!=-1){
            A.insert(pm.nkx+pm.nkr+k, pm.nkx+pm.nkr+kN) = -aNp[i][j];
        }
        if(kS!=-1){
            A.insert(pm.nkx+pm.nkr+k, pm.nkx+pm.nkr+kS) = -aSp[i][j];
        }
        if(kE!=-1){
            A.insert(pm.nkx+pm.nkr+k, pm.nkx+pm.nkr+kE) = -aEp[i][j];
        }
        if(kW!=-1){
            A.insert(pm.nkx+pm.nkr+k, pm.nkx+pm.nkr+kW) = -aWp[i][j];
        }
        //Ex用
        if(kEE!=-1){
            A.insert(pm.nkx+pm.nkr+k, kEE) = -aEEp[i][j];
        }
        if(kWW!=-1){
            A.insert(pm.nkx+pm.nkr+k, kWW) = -aWWp[i][j];
        }
        //Er用
        if(kNN!=-1){
            A.insert(pm.nkx+pm.nkr+k, pm.nkx+kNN) = -aNNp[i][j];
        }
        if(kSS!=-1){
            A.insert(pm.nkx+pm.nkr+k, pm.nkx+kSS) = -aSSp[i][j];
        }
    }
 

    // LU分解でAx = bを解く
    Eigen::SparseLU<Eigen::SparseMatrix<complex<double> > > solver;
    //cout << "solver start!" << endl;
    solver.compute(A);
    xv = solver.solve(b);


    // 実際の誤差を計算
    //Eigen::VectorXcd residual = A*xv- b;
    //cout << residual<< endl;
    //double actual_error = residual.norm()/(b.norm()+1e-100);
    //cout << "Actual error: " << actual_error << endl;

    //Ex-結果を戻す
    for (int k=0;k<pm.nkx;k++){
        int i = gk.ikx[k];
        int j = gk.jkx[k];
        gx.E1x[i][j] = xv[k];
    }
    //Er-結果を戻す
    for (int k=0;k<pm.nkr;k++){
        int i = gk.ikr[k];
        int j = gk.jkr[k];
        gr.E1r[i][j] = xv[pm.nkx+k];
    }
    //Ep-結果を戻す
    for (int k=0;k<pm.nkp;k++){
        int i = gk.ikp[k];
        int j = gk.jkp[k];
        gc.E1p[i][j] = xv[pm.nkx+pm.nkr+k];
    }

    /************************境界条件後処理**************************/
    //Ex-境界条件
    //********* Ex境界条件 (左) *********
    for (int k=0;k<mb.iBndWx.size();k++){
        int i = mb.iBndWx[k];
        int j = mb.jBndWx[k];

        if(mb.sBndWx[k] == 0){ //ディリクレ
            gx.E1x[i-1][j] = 0.0;
        }else if(mb.sBndWx[k] == 2){ //ガウス
            double qL = (gc.r[j-1] + gc.r[j])/2.0/(gc.r[j]+1e-100);
            double qR = (gc.r[j+1] + gc.r[j])/2.0/(gc.r[j]+1e-100);

            gx.E1x[i-1][j] = gx.E1x[i][j] + pm.dx/pm.dr*(qR*gr.E1r[i-1][j+1] - qL*gr.E1r[i-1][j]);
        }else if(mb.sBndWx[k] == 3){ //凹角 下側壁
            double qL = (gc.r[j-1] + gc.r[j])/2.0/(gc.r[j]+1e-100);
            double qR = (gc.r[j+1] + gc.r[j])/2.0/(gc.r[j]+1e-100);
            double deno = 1.0 + pm.dx/pm.dr*qL;

            gx.E1x[i-1][j] = gx.E1x[i][j]/deno + pm.dx/pm.dr*(qR*gr.E1r[i-1][j+1])/deno;
        }else if(mb.sBndWx[k] == 4){ //凹角 上側壁
            double qL = (gc.r[j-1] + gc.r[j])/2.0/(gc.r[j]+1e-100);
            double qR = (gc.r[j+1] + gc.r[j])/2.0/(gc.r[j]+1e-100);
            double deno = 1.0 + pm.dx/pm.dr*qR;

            gx.E1x[i-1][j] = gx.E1x[i][j]/deno - pm.dx/pm.dr*(qL*gr.E1r[i-1][j])/deno;
        }
    }
    //********* Ex境界条件 (右) *********
    for (int k=0;k<mb.iBndEx.size();k++){
        int i = mb.iBndEx[k];
        int j = mb.jBndEx[k];

        if(mb.sBndEx[k] == 2){ //ガウス
            double qL = (gc.r[j-1] + gc.r[j])/2.0/(gc.r[j]+1e-100);
            double qR = (gc.r[j+1] + gc.r[j])/2.0/(gc.r[j]+1e-100);
            
            gx.E1x[i+1][j] = gx.E1x[i][j] - pm.dx/pm.dr*(qR*gr.E1r[i][j+1] - qL*gr.E1r[i][j]);
        }else if(mb.sBndEx[k] == 4){ //凹角 上側壁
            double qL = (gc.r[j-1] + gc.r[j])/2.0/(gc.r[j]+1e-100);
            double qR = (gc.r[j+1] + gc.r[j])/2.0/(gc.r[j]+1e-100);
            double deno = 1.0 + pm.dx/pm.dr*qR;

            gx.E1x[i+1][j] = gx.E1x[i][j]/deno + pm.dx/pm.dr*qL*gr.E1r[i][j]/deno;
        }
    }
    //********* Ex境界条件 (下) *********
    for (int k=0;k<mb.iBndSx.size();k++){
        int i = mb.iBndSx[k];
        int j = mb.jBndSx[k];
        if(mb.sBndSx[k] == 0){ //ディリクレ
            gx.E1x[i][j-1] = -gc.r[j]/gc.r[j-1]*gx.E1x[i][j];
        }else if(mb.sBndSx[k] == 6){ //凸角 右側Open
            //処理なし
        }
    }
    //********* Ex境界条件 (上) *********
    for (int k=0;k<mb.iBndNx.size();k++){
        int i = mb.iBndNx[k];
        int j = mb.jBndNx[k];
        if(mb.sBndNx[k] == 0){ //ディリクレ
            gx.E1x[i][j+1] = -gc.r[j]/gc.r[j+1]*gx.E1x[i][j];
        }else if(mb.sBndNx[k] == 1){ //開放
            double c0 = sqrt(1.0/(ph::eps0*gc.epsr[i][j]*ph::mu0));
            double r_tmp = (gc.r[j+1] + gc.r[j])/2.0;
            double qL = gc.r[j]/(r_tmp+1e-100);
            double qR = gc.r[j+1]  /(r_tmp+1e-100);
            double deno = pow(2.0*c0*sqrt(qR),2) + pow(pm.omegam*pm.dr*qR,2);
            complex<double> R = complex<double>
                (pow(2.0*c0,2)*sqrt(qL*qR) - pow(pm.omegam*pm.dr,2)*qL*qR
                ,-2.0*c0*pm.omegam*pm.dr*(qL*sqrt(qR)+qR*sqrt(qL)));
            R = R/(deno+1e-100);
            gx.E1x[i][j+1] = R*gx.E1x[i][j];
        }else if(mb.sBndNx[k] == 5){ //凸角 左側Open
            //処理なし
        }else if(mb.sBndNx[k] == 6){ //凸角 右側Open
            //処理なし
        }
    }

    //********* Er境界条件 (左) *********
    for (int k=0;k<mb.iBndWr.size();k++){
        int i = mb.iBndWr[k];
        int j = mb.jBndWr[k];

        if(mb.sBndWr[k] == 0){ //ディリクレ
            gr.E1r[i-1][j] = -gr.E1r[i][j];
        }else if(mb.sBndWr[k] == 1){ //開放
            double c0 = sqrt(1.0/(ph::eps0*gc.epsr[i][j]*ph::mu0));
            double deno = pow(2.0*c0,2) + pow(pm.omegam*pm.dx,2);
            complex<double> R = complex<double>
                (pow(2.0*c0,2) - pow(pm.omegam*pm.dx,2),-4.0*c0*pm.omegam*pm.dx);
            R = R/(deno+1e-100);

            gr.E1r[i-1][j] = R*gr.E1r[i][j];
        }else if(mb.sBndWr[k] == 5){ //凸角 下側Open
            //処理なし
        }else if(mb.sBndWr[k] == 6){ //凸角 上側Open
            //処理なし
        }
    }
    //********* Er境界条件 (右) *********
    for (int k=0;k<mb.iBndEr.size();k++){
        int i = mb.iBndEr[k];
        int j = mb.jBndEr[k];

        if(mb.sBndEr[k] == 0){ //ディリクレ
            gr.E1r[i+1][j] = -gr.E1r[i][j];
        }else if(mb.sBndEr[k] == 5){ //凸角 下側Open
            //処理なし
        }
    }
    //********* Er境界条件 (下) *********
    for (int k=0;k<mb.iBndSr.size();k++){
        int i = mb.iBndSr[k];
        int j = mb.jBndSr[k];

        if(mb.sBndSr[k] == 2){ //ガウス
            double rR = (gc.r[j] + gc.r[j-1])/2.0;
            double rL = (gc.r[j-1] + gc.r[j-2])/2.0;
            
            gr.E1r[i][j-1] = rR/rL*gr.E1r[i][j] + pm.dr/pm.dx*gc.r[j-1]/rL*(gx.E1x[i+1][j-1]-gx.E1x[i][j-1]);
        }else if(mb.sBndSr[k] == 3){ //凹角 左側壁
            double rR = (gc.r[j] + gc.r[j-1])/2.0;
            double rL = (gc.r[j-1] + gc.r[j-2])/2.0;
            double deno = pm.dr/pm.dx*gc.r[j-1]/rL;
            if(gc.r[j] < gc.j_flc_bl[0][0]){ //マイクロ波の入り口部の特別処理
                deno = 1.0;
            }

            gr.E1r[i][j-1] = rR/rL*gr.E1r[i][j]/deno + pm.dr/pm.dx*gc.r[j-1]/rL*(gx.E1x[i+1][j-1])/deno;
        }
    }
    //********* Er境界条件 (上) *********
    for (int k=0;k<mb.iBndNr.size();k++){
        int i = mb.iBndNr[k];
        int j = mb.jBndNr[k];

        if(mb.sBndNr[k] == 2){ //ガウス
            double rR = (gc.r[j+1] + gc.r[j])/2.0;
            double rL = (gc.r[j] + gc.r[j-1])/2.0;

            gr.E1r[i][j+1] = rL/rR*gr.E1r[i][j] - pm.dr/pm.dx*gc.r[j]/rR*(gx.E1x[i+1][j]-gx.E1x[i][j]);
        }if(mb.sBndNr[k] == 3){ //凹角 左側壁
            double rR = (gc.r[j+1] + gc.r[j])/2.0;
            double rL = (gc.r[j] + gc.r[j-1])/2.0;
            double deno = pm.dr/pm.dx*gc.r[j]/rR;

            gr.E1r[i][j+1] = rL/rR*gr.E1r[i][j]/deno - pm.dr/pm.dx*gc.r[j]/rR*(gx.E1x[i+1][j])/deno;
        }if(mb.sBndNr[k] == 4){ //凹角 右側壁
            double rR = (gc.r[j+1] + gc.r[j])/2.0;
            double rL = (gc.r[j] + gc.r[j-1])/2.0;
            double deno = pm.dr/pm.dx*gc.r[j]/rR;

            gr.E1r[i][j+1] = rL/rR*gr.E1r[i][j]/deno + pm.dr/pm.dx*gc.r[j]/rR*(gx.E1x[i][j])/deno;
        }
    }

    //********* Ep境界条件 (左) *********
    for (int k=0;k<mb.iBndWp.size();k++){
        int i = mb.iBndWp[k];
        int j = mb.jBndWp[k];

        if(mb.sBndWp[k] == 0){ //ディリクレ
            gc.E1p[i-1][j] = -gc.E1p[i][j];
        }else if(mb.sBndWp[k] == 1){ //開放
            double c0 = sqrt(1.0/(ph::eps0*gc.epsr[i][j]*ph::mu0));
            double deno = pow(2.0*c0,2) + pow(pm.omegam*pm.dx,2);
            complex<double> R = complex<double>
                (pow(2.0*c0,2) - pow(pm.omegam*pm.dx,2),-4.0*c0*pm.omegam*pm.dx);
            R = R/(deno+1e-100);

            gc.E1p[i-1][j] = R*gc.E1p[i][j];
        }
    }
    //********* Ep境界条件 (右) *********
    for (int k=0;k<mb.iBndEp.size();k++){
        int i = mb.iBndEp[k];
        int j = mb.jBndEp[k];

        gc.E1p[i+1][j] = -gc.E1p[i][j];
    }
    //********* Ep境界条件 (下) *********
    for (int k=0;k<mb.iBndSp.size();k++){
        int i = mb.iBndSp[k];
        int j = mb.jBndSp[k];

        gc.E1p[i][j-1] = -gc.r[j]/gc.r[j-1]*gc.E1p[i][j];
    }
    //********* Ep境界条件 (上) *********
    for (int k=0;k<mb.iBndNp.size();k++){
        int i = mb.iBndNp[k];
        int j = mb.jBndNp[k];

        if(mb.sBndNp[k] == 0){ //ディリクレ
            gc.E1p[i][j+1] = -gc.r[j]/gc.r[j+1]*gc.E1p[i][j];
        }else if(mb.sBndNp[k] == 1){ //開放
            double c0 = sqrt(1.0/(ph::eps0*gc.epsr[i][j]*ph::mu0));
            double r_tmp = (gc.r[j+1] + gc.r[j])/2.0;
            double qL = gc.r[j]/(r_tmp+1e-100);
            double qR = gc.r[j+1]  /(r_tmp+1e-100);
            double deno = pow(2.0*c0*sqrt(qR),2) + pow(pm.omegam*pm.dr*qR,2);
            complex<double> R = complex<double>
                (pow(2.0*c0,2)*sqrt(qL*qR) - pow(pm.omegam*pm.dr,2)*qL*qR
                ,-2.0*c0*pm.omegam*pm.dr*(qL*sqrt(qR)+qR*sqrt(qL)));
            R = R/(deno+1e-100);
            
            gc.E1p[i][j+1] = R*gc.E1p[i][j];
        }
    }


    //check coefficients
    if(0==1){
        std::ofstream outputfile1("results/mw_coef0.csv");
        //std::ofstream outputfile1_1(char1+"_tmp"+char2+char_csv);
        //outputfile1<<"x,rho,Ui,rhoUix,Ue,rhoUex,rhoUex_E,rhoUex_D,E,phi,psi,rho_th,U_th,E_th,phi_th,psi_th,rate_ionize,divi,dive1,dive2,dive3, div_poisson, LHS, RHS,zero" << std::endl;
        //outputfile1<<"x,rhon,rho,rhoe(test),Uix,rhoUix,Uex,Uey,rhoUex,rhoUey,Te,heat_flux,Bz,E,phi,nu_m,nu_en,nu_ei,nu_ionz,nu_exc,nu_wall,nu_ano,mue_parae,mue_perp,Halle,jd,Id,rate_eloss,rate_ionize,divi,dive,divn,engy_LHS,engy_LHS1,engy_LHS2,engy_LHS3,engy_LHS4,engy_RHS,engy_RHS1,engy_RHS2,engy_RHS3,zero" << std::endl;
        outputfile1<<"i,j,x,r,aPx(R),aPx(I),aWx(R),aWx(I),aEx(R),aEx(I),aNx(R),aNx(I),aSx(R),aSx(I),aNWx(R),aNWx(I),aNEx(R),aNEx(I),aSWx(R),aSWx(I),aSEx(R),aSEx(I),aEEx(R),aEEx(I),aWWx(R),aWWx(I),bx(R),bx(I),aPr(R),aPr(I),aWr(R),aWr(I),aEr(R),aEr(I),aNr(R),aNr(I),aSr(R),aSr(I),aNWr(R),aNWr(I),aNEr(R),aNEr(I),aSWr(R),aSWr(I),aSEr(R),aSEr(I),aNNr(R),aNNr(I),aSSr(R),aSSr(I),br(R),br(I),aPp(R),aPp(I),aWp(R),aWp(I),aEp(R),aEp(I),aNp(R),aNp(I),aSp(R),aSp(I),aEEp(R),aEEp(I),aWWp(R),aWWp(I),aNNp(R),aNNp(I),aSSp(R),aSSp(I),bp(R),bp(I),zero" << std::endl;
        //outputfile1_1 << "i,j,x,r,rhom,rhoUmx,rhoUmr,nabla_rhoUm,rhon,rhoUnx,rhoUnr,nabla_rhoUn,zero" << std::endl;
        
        for(int i=1;i<=pm.ni;i++){
            for(int j=1;j<=pm.nj;j++){
                outputfile1 << i << ","<< j << "," << gc.x[i]<< ","<< gc.r[j]
                    << "," << real(aPx[i][j] )<< "," << imag(aPx[i][j] )
                    << "," << real(aWx[i][j] )<< "," << imag(aWx[i][j] )
                    << "," << real(aEx[i][j] )<< "," << imag(aEx[i][j] )
                    << "," << real(aNx[i][j] )<< "," << imag(aNx[i][j] )
                    << "," << real(aSx[i][j] )<< "," << imag(aSx[i][j] )
                    << "," << real(aNWx[i][j])<< "," << imag(aNWx[i][j])
                    << "," << real(aNEx[i][j])<< "," << imag(aNEx[i][j])
                    << "," << real(aSWx[i][j])<< "," << imag(aSWx[i][j])
                    << "," << real(aSEx[i][j])<< "," << imag(aSEx[i][j])
                    << "," << real(aEEx[i][j])<< "," << imag(aEEx[i][j])
                    << "," << real(aWWx[i][j])<< "," << imag(aWWx[i][j])
                    << "," << real(bx[i][j]  )<< "," << imag(bx[i][j]  )

                    << "," << real(aPr[i][j] )<< "," << imag(aPr[i][j] )
                    << "," << real(aWr[i][j] )<< "," << imag(aWr[i][j] )
                    << "," << real(aEr[i][j] )<< "," << imag(aEr[i][j] )
                    << "," << real(aNr[i][j] )<< "," << imag(aNr[i][j] )
                    << "," << real(aSr[i][j] )<< "," << imag(aSr[i][j] )
                    << "," << real(aNWr[i][j])<< "," << imag(aNWr[i][j])
                    << "," << real(aNEr[i][j])<< "," << imag(aNEr[i][j])
                    << "," << real(aSWr[i][j])<< "," << imag(aSWr[i][j])
                    << "," << real(aSEr[i][j])<< "," << imag(aSEr[i][j])
                    << "," << real(aNNr[i][j])<< "," << imag(aNNr[i][j])
                    << "," << real(aSSr[i][j])<< "," << imag(aSSr[i][j])
                    << "," << real(br[i][j]  )<< "," << imag(br[i][j]  )

                    << "," << real(aPp[i][j] )<< "," << imag(aPp[i][j] )
                    << "," << real(aWp[i][j] )<< "," << imag(aWp[i][j] )
                    << "," << real(aEp[i][j] )<< "," << imag(aEp[i][j] )
                    << "," << real(aNp[i][j] )<< "," << imag(aNp[i][j] )
                    << "," << real(aSp[i][j] )<< "," << imag(aSp[i][j] )
                    << "," << real(aEEp[i][j])<< "," << imag(aEEp[i][j])
                    << "," << real(aWWp[i][j])<< "," << imag(aWWp[i][j])
                    << "," << real(aNNp[i][j])<< "," << imag(aNNp[i][j])
                    << "," << real(aSSp[i][j])<< "," << imag(aSSp[i][j])
                    << "," << real(bp[i][j]  )<< "," << imag(bp[i][j]  )
                    << "," << 0.0
                    <<std::endl;
            }
        }
    }


}

//*****************************************************************
//**                                                             **
//**           void update_energy_profile()                      **
//**                                                             **
//*****************************************************************
void EmfieldModule::update_energy_profile(Params &pm, GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr){

    double power_sum = 0.0;
    for (int iblock=0;iblock<5;iblock++){ 
        for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
            for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){

                double rL = (gc.r[j]+gc.r[j-1])/2.0;
                double rR = (gc.r[j]+gc.r[j+1])/2.0;

                complex<double> nu_cmp(gc.nu_m1[i][j],pm.omegam);
                complex<double> Hallcmpx = ph::e0*gc.Bx[i][j]/(pm.masse*nu_cmp);
                complex<double> Hallcmpr = ph::e0*gc.Br[i][j]/(pm.masse*nu_cmp);
                complex<double> HallcmpMag2 =  Hallcmpx*Hallcmpx + Hallcmpr*Hallcmpr;
                complex<double> coef = gc.rhoe[i][j]*ph::e0*ph::e0/(pm.masse*nu_cmp+1e-100)/(1.0 + HallcmpMag2);
                complex<double> sigmaexr = coef*(Hallcmpx*Hallcmpr);
                complex<double> sigmaexp = coef*(Hallcmpr);
                complex<double> sigmaexx = coef*(1.0 + Hallcmpx*Hallcmpx);
                complex<double> sigmaerr = coef*(1.0 + Hallcmpr*Hallcmpr);
                complex<double> sigmaerp = coef*(-Hallcmpx);
                complex<double> sigmaerx = coef*(Hallcmpx*Hallcmpr);
                complex<double> sigmaepr = coef*(Hallcmpx);
                complex<double> sigmaepp = coef*(1.0);
                complex<double> sigmaepx = coef*(-Hallcmpr);

                complex<double> iomegaMu0(0.0,pm.omegam*ph::mu0);

                complex<double> Ex_tmp = (gx.E1x[i][j] + gx.E1x[i+1][j])/2.0;
                complex<double> Er_tmp = (rL*gr.E1r[i][j] + rR*gr.E1r[i][j+1])/(2.0*gc.r[j]);
                complex<double> Ep_tmp = gc.E1p[i][j];
                complex<double> Jx_tmp = sigmaexr*Er_tmp + sigmaexp*Ep_tmp + sigmaexx*Ex_tmp;
                complex<double> Jr_tmp = sigmaerr*Er_tmp + sigmaerp*Ep_tmp + sigmaerx*Ex_tmp;
                complex<double> Jp_tmp = sigmaepr*Er_tmp + sigmaepp*Ep_tmp + sigmaepx*Ex_tmp;

                gc.Pabs[i][j] = 0.5*real(Jx_tmp*conj(Ex_tmp) + Jr_tmp*conj(Er_tmp) + Jp_tmp*conj(Ep_tmp));
                power_sum = power_sum + gc.Pabs[i][j]*(pm.dx*pm.dr*2.0*M_PI*gc.r[j]);
            }
        }
    }

    cout  << endl;
    cout << "---------------------------------" << endl;
    cout << "Total power = " << power_sum << " W" << endl;

    //calculate reflection
    if(pm.icon_mwRef == 1){

        //set refelence plane
        //---------------------------------
        double x_ref = 0.004;
        int i_ref = int(x_ref/pm.dx + 1.5);
        //---------------------------------

        int i = i_ref;
        double x_tmp = (gc.x[i] + gc.x[i-1])/2.0;
        //calculate current I (use averaged value in r-direction for robustness)
        //---------------------------------
        int ncount = 0;
        complex<double>I_ref_tmp(0.0,0.0);
        for (int j=gr.j_flr_bl[5][0];j<=gr.j_flr_bl[5][1];j++){
            double r_tmp = (gc.r[j] + gc.r[j-1])/2.0;
            //double qR = (gc.r[j] + gc.r[j+1])/(2.0*gc.r[j]);

            complex<double> H_ref = (gr.E1r[i][j] - gr.E1r[i-1][j])/pm.dx/(pm.omegam*ph::mu0)*complex<double>(0,1.0);
            I_ref_tmp = I_ref_tmp + 2.0*M_PI*r_tmp*H_ref;

            //cout << "j = " << j << " - 1/2, r = "<< r_tmp << " , H_ref = "<<H_ref<<", I_ref = "<<2.0*M_PI*r_tmp*H_ref<< endl;
            ncount = ncount + 1;
        }
        I_ref_tmp = I_ref_tmp/double(ncount);
        //---------------------------------

        //calculate current V
        //---------------------------------
        complex<double> V_ref_tmp(0.0,0.0);
        for (int j=gc.j_flc_bl[5][0];j<=gc.j_flc_bl[5][1];j++){
            //double qL = (gc.r[j] + gc.r[j-1])/(2.0*gc.r[j]);
            //double qR = (gc.r[j] + gc.r[j+1])/(2.0*gc.r[j]);
            //double r_tmp = (gc.r[j] + gc.r[j-1])/2.0;
            double rL = (gc.r[j] + gc.r[j-1])/2.0;
            double rR = (gc.r[j] + gc.r[j+1])/2.0;

            //complex<double> Er_tmp = (gr.E1r[i][j] + gr.E1r[i-1][j])/2.0;
            complex<double> Er_tmp = (rR*(gr.E1r[i][j+1] + gr.E1r[i-1][j+1]) + rL*(gr.E1r[i][j] + gr.E1r[i-1][j]))/(4.0*gc.r[j]);
            
            double c0 = sqrt(1.0/(ph::eps0*gc.epsr[i][j]*ph::mu0));
            
            //cout << "Check ,i = "<<i <<", j = "<< j << ", "<< abs(Er_tmp)*abs(Er_tmp)
            //    //<< ","<<abs(Er_Fw[j+1])<< ","<<abs(Er_Fw[j])
            //    << ","<<(Er_tmp)
            //    << ","<<pow(abs((Er_tmp - Er_anly)/Er_anly),2)*100<< " %"
            //    << endl;

            V_ref_tmp = V_ref_tmp + Er_tmp*pm.dr;
            //cout << "j = " << j <<", Er = "<< Er_tmp << ", rEr = "<< r_tmp*Er_tmp << ", V_ref = "<<V_ref << ", V_ref_th = "<<r_tmp*Er_tmp*log(r2/r1)<< endl;
        }
        //---------------------------------

        complex<double> Z_ref_tmp = V_ref_tmp/(I_ref_tmp + 1e-100);
        
        cout << "Impedance at reference plane x = " << x_tmp << endl;
        cout << "Z_ref = " << Z_ref_tmp << ", V_ref = " << V_ref_tmp << ", I_ref = " << I_ref_tmp 
            << ", |Z_ref| = " << abs(Z_ref_tmp) 
            << ", arg(Z_ref) = " << arg(Z_ref_tmp)*180/M_PI << " deg"<< endl;

        complex<double> a2 = (V_ref_tmp - pm.Z0_base*I_ref_tmp)/(2.0*sqrt(pm.Z0_base));
        complex<double> b2 = (V_ref_tmp + pm.Z0_base*I_ref_tmp)/(2.0*sqrt(pm.Z0_base));

        complex<double> S11(pm.S11_mag * cos(pm.S11_arg_deg/180.0*M_PI), pm.S11_mag * sin(pm.S11_arg_deg/180.0*M_PI));
        complex<double> S21(pm.S21_mag * cos(pm.S21_arg_deg/180.0*M_PI), pm.S21_mag * sin(pm.S21_arg_deg/180.0*M_PI));
        complex<double> S12(pm.S12_mag * cos(pm.S12_arg_deg/180.0*M_PI), pm.S12_mag * sin(pm.S12_arg_deg/180.0*M_PI));
        complex<double> S22(pm.S22_mag * cos(pm.S22_arg_deg/180.0*M_PI), pm.S22_mag * sin(pm.S22_arg_deg/180.0*M_PI));

        complex<double> a1 = (b2 - S22*a2)/S21;
        complex<double> b1 = S11*a1 + S12*a2;
        double P_fwd = norm(a1);
        double P_ref = norm(b1);

        cout << endl;
        cout << "Fwd power = " << P_fwd << " W" << endl;
        cout << "Ref power = " << P_ref << " W" << endl;
        cout << "Transmitted power = " << P_fwd - P_ref << " W" << endl;
        
        double ratio = pm.Pmw/P_fwd;

        cout << "power ratio = " << ratio << endl;
        cout << "J1r_exc ="  << pm.J1r_exc;
        pm.J1r_exc = pm.J1r_exc*sqrt(ratio);
        cout << " -> " << pm.J1r_exc << endl;
        
        //電力調整
        for (int i=0;i<pm.ni+1;i++){
            for (int j=0;j<pm.nj+1;j++){
                gc.Pabs[i][j] = gc.Pabs[i][j]*ratio;
                gr.J1r[i][j] = gr.J1r[i][j]*sqrt(ratio);
            }
        }
       
    }else{

        double ratio = pm.Pmw/power_sum;
        
        //cout << "J1r_exc1 =" << J1r_exc <<","<<real(J1r_exc) <<","<<sqrt(ratio)<<","<<J1r_exc*sqrt(ratio)<< endl;
        //電力調整
        cout << "power ratio = " << ratio << endl;
        cout << "J1r_exc = " << pm.J1r_exc;
        pm.J1r_exc = pm.J1r_exc*sqrt(ratio);
        cout << " -> " << pm.J1r_exc << endl;

        for (int i=0;i<pm.ni+1;i++){
            for (int j=0;j<pm.nj+1;j++){
                gc.Pabs[i][j] = gc.Pabs[i][j]*ratio;
                gr.J1r[i][j] = gr.J1r[i][j]*sqrt(ratio);
            }
        }
        //cout << "J1r_exc2 =" << J1r_exc <<","<<real(J1r_exc) <<","<<sqrt(ratio)<<","<<J1r_exc*sqrt(ratio)<< endl;
        

    }

    cout << "---------------------------------" << endl;
    cout  << endl;

}

