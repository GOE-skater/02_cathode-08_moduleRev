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

//*****************************************************************
//**                                                             **
//**           void solve_Efield                                 **
//**                                                             **
//*****************************************************************
void solve_Microwave(){

    //テスト用 ※ 通常はコメントアウトすること
    for (int i=0;i<ni+2;i++){
        for (int j=0;j<nj+2;j++){
            rhoe[i][j] = 0.0;
            nu_m1[i][j] = nu_eff;
        }
    }

    //Ex用 係数
    std::vector<std::vector<std::complex<double> > > aPx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aWx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aEx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aNx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aSx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aNEx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aNWx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aSEx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aSWx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aEEx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aWWx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > bx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数

    //Er用 係数
    std::vector<std::vector<std::complex<double> > > aPr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aWr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aEr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aNr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aSr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aNEr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aNWr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aSEr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aSWr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aNNr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aSSr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > br(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数

    //Ephi用 係数
    std::vector<std::vector<std::complex<double> > > aPp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aWp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aEp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aNp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aSp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aNNp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aSSp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aEEp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aWWp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > bp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数

    /******************** 係数作成 いらないところも含めて普通に全部やる方が速い ********************/
    std::complex<double> EPS(1e-100,1e-100); //微小複素数
    //for Ex
    for (int iblock=0;iblock<6;iblock++){ 
        for (int i=i_flx_bl[iblock][0];i<=i_flx_bl[iblock][1];i++){ 
            for (int j=j_flx_bl[iblock][0];j<=j_flx_bl[iblock][1];j++){
                double qL = (r[j-1] + r[j])/2.0/(r[j]+1e-100);
                double qR = (r[j+1] + r[j])/2.0/(r[j]+1e-100);
                
                double nu_m_tmp = (nu_m1[i][j] + nu_m1[i-1][j])/2.0;
                double rho_tmp  = (rhoe[i][j]  + rhoe[i-1][j] )/2.0;
                double Bx_tmp   = (Bx[i][j]   + Bx[i-1][j]  )/2.0;
                double Br_tmp   = (Br[i][j]   + Br[i-1][j]  )/2.0;

                std::complex<double> nu_cmp(nu_m_tmp,omegam);
                std::complex<double> Hallcmpx = e0*Bx_tmp/(masse*nu_cmp + EPS);
                std::complex<double> Hallcmpr = e0*Br_tmp/(masse*nu_cmp + EPS);
                std::complex<double> HallcmpMag2 =  Hallcmpx*Hallcmpx + Hallcmpr*Hallcmpr;
                std::complex<double> coef = rho_tmp*e0*e0/(masse*nu_cmp + EPS)/(1.0 + HallcmpMag2 + EPS);
                std::complex<double> sigmaexr = coef*(Hallcmpx*Hallcmpr);
                std::complex<double> sigmaexp = coef*(Hallcmpr);
                std::complex<double> sigmaexx = coef*(1.0 + Hallcmpx*Hallcmpx);
                std::complex<double> iomegaMu0(0.0,omegam*mu0);

                aPx[i][j] = 2.0*dx*dx + 2.0*dr*dr
                    - eps0*epsr[i][j]*mu0*omegam*omegam*dx*dx*dr*dr + iomegaMu0*sigmaexx*dx*dx*dr*dr;
                aEx[i][j] = dr*dr;
                aWx[i][j] = dr*dr;
                aNx[i][j] = qR*dx*dx;
                aSx[i][j] = qL*dx*dx;
                bx[i][j] = std::complex<double>(0,-omegam*mu0*dx*dx*dr*dr)*J1x[i][j];

                //Epの要素
                aEEx[i][j] = -iomegaMu0*sigmaexp*0.5*dx*dx*dr*dr;
                aWWx[i][j] = -iomegaMu0*sigmaexp*0.5*dx*dx*dr*dr;

                //Erの要素
                aNEx[i][j] = -iomegaMu0*sigmaexr*0.25*qR*dx*dx*dr*dr;
                aNWx[i][j] = -iomegaMu0*sigmaexr*0.25*qR*dx*dx*dr*dr;
                aSEx[i][j] = -iomegaMu0*sigmaexr*0.25*qL*dx*dx*dr*dr;
                aSWx[i][j] = -iomegaMu0*sigmaexr*0.25*qL*dx*dx*dr*dr;
            }
        }
    }

    //for Er
    for (int iblock=0;iblock<6;iblock++){ 
        for (int i=i_flr_bl[iblock][0];i<=i_flr_bl[iblock][1];i++){ 
            for (int j=j_flr_bl[iblock][0];j<=j_flr_bl[iblock][1];j++){
                double r_tmp = (r[j-1] + r[j])/2.0;
                double qL = r[j-1]/(r_tmp+1e-100);
                double qR = r[j]  /(r_tmp+1e-100);

                double nu_m_tmp = (nu_m1[i][j] + nu_m1[i][j-1])/2.0;
                double rho_tmp  = (rhoe[i][j]  + rhoe[i][j-1] )/2.0;
                double Bx_tmp   = (r[j]*Bx[i][j]   + r[j-1]*Bx[i][j-1]  )/(2.0*r_tmp);
                double Br_tmp   = (r[j]*Br[i][j]   + r[j-1]*Br[i][j-1]  )/(2.0*r_tmp);

                std::complex<double> nu_cmp(nu_m_tmp,omegam);
                std::complex<double> Hallcmpx = e0*Bx_tmp/(masse*nu_cmp + EPS);
                std::complex<double> Hallcmpr = e0*Br_tmp/(masse*nu_cmp + EPS);
                std::complex<double> HallcmpMag2 =  Hallcmpx*Hallcmpx + Hallcmpr*Hallcmpr;
                std::complex<double> coef = rho_tmp*e0*e0/(masse*nu_cmp + EPS)/(1.0 + HallcmpMag2 + EPS);
                std::complex<double> sigmaerr = coef*(1.0+Hallcmpr*Hallcmpr);
                std::complex<double> sigmaerp = coef*(-Hallcmpx);
                std::complex<double> sigmaerx = coef*(Hallcmpx*Hallcmpr);
                std::complex<double> iomegaMu0(0.0,omegam*mu0);

                aPr[i][j] = 2.0*dx*dx + 2.0*dr*dr + dx*dx*dr*dr/(r_tmp*r_tmp)
                    - eps0*epsr[i][j]*mu0*omegam*omegam*dx*dx*dr*dr + iomegaMu0*sigmaerr*dx*dx*dr*dr;
                aEr[i][j] = dr*dr;
                aWr[i][j] = dr*dr;
                aNr[i][j] = qR*dx*dx;
                aSr[i][j] = qL*dx*dx;
                br[i][j] = std::complex<double>(0,-omegam*mu0*dx*dx*dr*dr)*J1r[i][j];

                //Epの要素
                aNNr[i][j] = -iomegaMu0*sigmaerp*0.5*qR*dx*dx*dr*dr;
                aSSr[i][j] = -iomegaMu0*sigmaerp*0.5*qL*dx*dx*dr*dr;

                //Exの要素
                aNEr[i][j] = -iomegaMu0*sigmaerx*0.25*qR*dx*dx*dr*dr;
                aNWr[i][j] = -iomegaMu0*sigmaerx*0.25*qR*dx*dx*dr*dr;
                aSEr[i][j] = -iomegaMu0*sigmaerx*0.25*qL*dx*dx*dr*dr;
                aSWr[i][j] = -iomegaMu0*sigmaerx*0.25*qL*dx*dx*dr*dr;
            }
        }
    }

    //for Ephi
    for (int iblock=0;iblock<6;iblock++){ 
        for (int i=i_flc_bl[iblock][0];i<=i_flc_bl[iblock][1];i++){ 
            for (int j=j_flc_bl[iblock][0];j<=j_flc_bl[iblock][1];j++){
                double qL = (r[j-1] + r[j])/2.0/(r[j]+1e-100);
                double qR = (r[j+1] + r[j])/2.0/(r[j]+1e-100);

                double nu_m_tmp = nu_m1[i][j];
                double rho_tmp  = rhoe[i][j];
                double Bx_tmp   = Bx[i][j];
                double Br_tmp   = Br[i][j];

                std::complex<double> nu_cmp(nu_m_tmp,omegam);
                std::complex<double> Hallcmpx = e0*Bx_tmp/(masse*nu_cmp + EPS);
                std::complex<double> Hallcmpr = e0*Br_tmp/(masse*nu_cmp + EPS);
                std::complex<double> HallcmpMag2 =  Hallcmpx*Hallcmpx + Hallcmpr*Hallcmpr;
                std::complex<double> coef = rho_tmp*e0*e0/(masse*nu_cmp + EPS)/(1.0 + HallcmpMag2 + EPS);
                std::complex<double> sigmaepr = coef*(Hallcmpx);
                std::complex<double> sigmaepp = coef*(1.0);
                std::complex<double> sigmaepx = coef*(-Hallcmpr);
                std::complex<double> iomegaMu0(0.0,omegam*mu0);

                aPp[i][j] = 2.0*dx*dx + 2.0*dr*dr + dx*dx*dr*dr/(r[j]*r[j])
                    - eps0*epsr[i][j]*mu0*omegam*omegam*dx*dx*dr*dr + iomegaMu0*sigmaepp*dx*dx*dr*dr;
                aEp[i][j] = dr*dr;
                aWp[i][j] = dr*dr;
                aNp[i][j] = qR*dx*dx;
                aSp[i][j] = qL*dx*dx;
                bp[i][j] = std::complex<double>(0,-omegam*mu0*dx*dx*dr*dr)*J1p[i][j];

                //Exの要素
                aEEp[i][j] = -iomegaMu0*sigmaepx*0.5*dx*dx*dr*dr;
                aWWp[i][j] = -iomegaMu0*sigmaepx*0.5*dx*dx*dr*dr;

                //Erの要素
                aNNp[i][j] = -iomegaMu0*sigmaepr*0.5*qR*dx*dx*dr*dr;
                aSSp[i][j] = -iomegaMu0*sigmaepr*0.5*qL*dx*dx*dr*dr;
            }
        }
    }

    /******************** 境界条件設定 (係数修正) ********************/
    //********* Ex境界条件 (左) *********
    for (int k=0;k<iBndWx.size();k++){
        int i = iBndWx[k];
        int j = jBndWx[k];

        if(sBndWx[k] == 0){ //ディリクレ
            aWx[i][j] = 0.0;
        }else if(sBndWx[k] == 2){ //ガウス
            double qL = (r[j-1] + r[j])/2.0/(r[j]+1e-100);
            double qR = (r[j+1] + r[j])/2.0/(r[j]+1e-100);

            aPx[i][j] = aPx[i][j] - aWx[i][j];
            aNWx[i][j] = aNWx[i][j] + aWx[i][j]*dx/dr*qR;
            aSWx[i][j] = aSWx[i][j] - aWx[i][j]*dx/dr*qL;
            aWx[i][j] = 0.0;
        }else if(sBndWx[k] == 3){ //凹角 下側壁
            double qL = (r[j-1] + r[j])/2.0/(r[j]+1e-100);
            double qR = (r[j+1] + r[j])/2.0/(r[j]+1e-100);
            double deno = 1.0 + dx/dr*qL;

            aPx[i][j] = aPx[i][j] - aWx[i][j]/deno;
            aNWx[i][j] = aNWx[i][j] + aWx[i][j]*(dx/dr*qR)/deno;
            aWx[i][j] = 0.0;
        }else if(sBndWx[k] == 4){ //凹角 上側壁
            double qL = (r[j-1] + r[j])/2.0/(r[j]+1e-100);
            double qR = (r[j+1] + r[j])/2.0/(r[j]+1e-100);
            double deno = 1.0 + dx/dr*qR;

            aPx[i][j] = aPx[i][j] - aWx[i][j]/deno;
            aSWx[i][j] = aSWx[i][j] - aWx[i][j]*(dx/dr*qL)/deno;
            aWx[i][j] = 0.0;
        }
    }
    //********* Ex境界条件 (右) *********
    for (int k=0;k<iBndEx.size();k++){
        int i = iBndEx[k];
        int j = jBndEx[k];

        if(sBndEx[k] == 2){ //ガウス
            double qL = (r[j-1] + r[j])/2.0/(r[j]+1e-100);
            double qR = (r[j+1] + r[j])/2.0/(r[j]+1e-100);

            aPx[i][j] = aPx[i][j] - aEx[i][j];
            aNEx[i][j] = aNEx[i][j] - aEx[i][j]*dx/dr*qR;
            aSEx[i][j] = aSEx[i][j] + aEx[i][j]*dx/dr*qL;
            aEx[i][j] = 0.0;
        }else if(sBndEx[k] == 4){ //凹角 上側壁
            double qL = (r[j-1] + r[j])/2.0/(r[j]+1e-100);
            double qR = (r[j+1] + r[j])/2.0/(r[j]+1e-100);
            double deno = 1.0 + dx/dr*qR;
            
            aPx[i][j] = aPx[i][j] - aEx[i][j]/deno;
            aSEx[i][j] = aSEx[i][j] + aEx[i][j]*(dx/dr*qL)/deno;
            aEx[i][j] = 0.0;
        }
    }
    //********* Ex境界条件 (下) *********
    for (int k=0;k<iBndSx.size();k++){
        int i = iBndSx[k];
        int j = jBndSx[k];
        if(sBndSx[k] == 0){ //ディリクレ
            aPx[i][j] = aPx[i][j] + r[j]/r[j-1]*aSx[i][j];
            aSx[i][j] = 0.0;
        }else if(sBndSx[k] == 6){ //凸角 右側Open
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;

            aPx[i][j] = aPx[i][j] - aSx[i][j]*( dr/dx*r[j]/rL);
            aWx[i][j] = aWx[i][j] + aSx[i][j]*(-dr/dx*r[j]/rL);
            aNWx[i][j] = aNWx[i][j] + aSx[i][j]*(rR/rL);
            aSx[i][j] = 0.0;
        }
    }
    //********* Ex境界条件 (上) *********
    for (int k=0;k<iBndNx.size();k++){
        int i = iBndNx[k];
        int j = jBndNx[k];
        if(sBndNx[k] == 0){ //ディリクレ
            aPx[i][j] = aPx[i][j] + r[j]/r[j+1]*aNx[i][j];
            aNx[i][j] = 0.0;
        }if(sBndNx[k] == 1){ //開放
            double c0 = sqrt(1.0/(eps0*epsr[i][j]*mu0));
            double r_tmp = (r[j+1] + r[j])/2.0;
            double qL = r[j]/(r_tmp+1e-100);
            double qR = r[j+1]  /(r_tmp+1e-100);
            double deno = pow(2.0*c0*sqrt(qR),2) + pow(omegam*dr*qR,2);
            std::complex<double> R = std::complex<double>
                (pow(2.0*c0,2)*sqrt(qL*qR) - pow(omegam*dr,2)*qL*qR
                ,-2.0*c0*omegam*dr*(qL*sqrt(qR)+qR*sqrt(qL)));
            R = R/(deno+1e-100);

            aPx[i][j] = aPx[i][j] - aNx[i][j]*R;
            aNx[i][j] = 0.0;
        }else if(sBndNx[k] == 5){ //凸角 左側Open (*)
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;

            aPx[i][j] = aPx[i][j] - aNx[i][j]*( dr/dx*r[j]/rR);
            aEx[i][j] = aEx[i][j] + aNx[i][j]*(-dr/dx*r[j]/rR);
            aSEx[i][j] = aSEx[i][j] + aNx[i][j]*( rL/rR);
            aNx[i][j] = 0.0;
        }else if(sBndNx[k] == 6){ //凸角 右側Open
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;

            aPx[i][j] = aPx[i][j] - aNx[i][j]*( dr/dx*r[j]/rR);
            aWx[i][j] = aWx[i][j] + aNx[i][j]*(-dr/dx*r[j]/rR);
            aSWx[i][j] = aSWx[i][j] + aNx[i][j]*(-rL/rR);
            aNx[i][j] = 0.0;
        }
    }

    //********* Er境界条件 (左) *********
    for (int k=0;k<iBndWr.size();k++){
        int i = iBndWr[k];
        int j = jBndWr[k];

        if(sBndWr[k] == 0){ //ディリクレ
            aPr[i][j] = aPr[i][j] + aWr[i][j];
            aWr[i][j] = 0.0;
        }else if(sBndWr[k] == 1){ //開放
            double c0 = sqrt(1.0/(eps0*epsr[i][j]*mu0));
            double deno = pow(2.0*c0,2) + pow(omegam*dx,2);
            std::complex<double> R = std::complex<double>(pow(2.0*c0,2) - pow(omegam*dx,2),-4.0*c0*omegam*dx);
            R = R/(deno+1e-100);
            
            aPr[i][j] = aPr[i][j] - aWr[i][j]*R;
            aWr[i][j] = 0.0;
        }else if(sBndWr[k] == 5){ //凸角 下側Open
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;
            
            aPr[i][j] = aPr[i][j] - aWr[i][j]*( dx/dr*rL/r[j]);
            aNr[i][j] = aNr[i][j] + aWr[i][j]*(-dx/dr*rR/r[j]);
            aNEr[i][j] = aNEr[i][j] + aWr[i][j]*(-1.0);
            aWr[i][j] = 0.0;
        }else if(sBndWr[k] == 6){ //凸角 上側Open
            double rR = (r[j] + r[j-1])/2.0;
            double rL = (r[j-1] + r[j-2])/2.0;

            aPr[i][j] = aPr[i][j] - aWr[i][j]*( dx/dr*rR/r[j-1]);
            aSr[i][j] = aSr[i][j] + aWr[i][j]*(-dx/dr*rL/r[j-1]);
            aSEr[i][j] = aSEr[i][j] + aWr[i][j]*( 1.0);
            aWr[i][j] = 0.0;
        }
    }
    //********* Er境界条件 (右) *********
    for (int k=0;k<iBndEr.size();k++){
        int i = iBndEr[k];
        int j = jBndEr[k];

        if(sBndEr[k] == 0){ //ディリクレ
            aPr[i][j] = aPr[i][j] + aEr[i][j];
            aEr[i][j] = 0.0;
        }else if(sBndEr[k] == 5){ //凸角 下側Open (*)
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;

            aPr[i][j] = aPr[i][j] - aEr[i][j]*( dx/dr*rL/r[j]);
            aNr[i][j] = aNr[i][j] + aEr[i][j]*(-dx/dr*rR/r[j]);
            aNWr[i][j] = aNWr[i][j] + aEr[i][j]*( 1.0);
            aEr[i][j] = 0.0;
        }
    }
    //********* Er境界条件 (下) *********
    for (int k=0;k<iBndSr.size();k++){
        int i = iBndSr[k];
        int j = jBndSr[k];

        if(sBndSr[k] == 2){ //ガウス
            double rR = (r[j] + r[j-1])/2.0;
            double rL = (r[j-1] + r[j-2])/2.0;
            aPr[i][j] = aPr[i][j] - aSr[i][j]*rR/rL;
            aSEr[i][j] = aSEr[i][j] + aSr[i][j]*dr/dx*r[j-1]/rL;
            aSWr[i][j] = aSWr[i][j] - aSr[i][j]*dr/dx*r[j-1]/rL;
            aSr[i][j] = 0.0;
        }else if(sBndSr[k] == 3){ //凹角 左側壁
            double rR = (r[j] + r[j-1])/2.0;
            double rL = (r[j-1] + r[j-2])/2.0;
            double deno = dr/dx*r[j-1]/rL;
            if(r[j] < j_flc_bl[0][0]){ //マイクロ波の入り口部の特別処理
                deno = 1.0;
            }

            aPr[i][j] = aPr[i][j] - aSr[i][j]*rR/rL/deno;
            aNEr[i][j] = aNEr[i][j] + aSr[i][j]*dr/dx*r[j-1]/rL/deno;
            aSr[i][j] = 0.0;
        }
    }
    //********* Er境界条件 (上) *********
    for (int k=0;k<iBndNr.size();k++){
        int i = iBndNr[k];
        int j = jBndNr[k];

        if(sBndNr[k] == 2){ //ガウス
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;
            aPr[i][j] = aPr[i][j] - aNr[i][j]*rL/rR;
            aNEr[i][j] = aNEr[i][j] - aNr[i][j]*dr/dx*r[j]/rR;
            aNWr[i][j] = aNWr[i][j] + aNr[i][j]*dr/dx*r[j]/rR;
            aNr[i][j] = 0.0;
        }if(sBndNr[k] == 3){ //凹角 左側壁
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;
            double deno = dr/dx*r[j]/rR;

            aPr[i][j] = aPr[i][j] - aNr[i][j]*rL/rR/deno;
            aSEr[i][j] = aSEr[i][j] - aNr[i][j]*dr/dx*r[j]/rR/deno;
            aNr[i][j] = 0.0;
        }if(sBndNr[k] == 4){ //凹角 右側壁
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;
            double deno = dr/dx*r[j]/rR;

            aPr[i][j] = aPr[i][j] - aNr[i][j]*rL/rR/deno;
            aSWr[i][j] = aSWr[i][j] + aNr[i][j]*dr/dx*r[j]/rR/deno;
            aNr[i][j] = 0.0;
        }
    }

    //********* Ep境界条件 (左) *********
    for (int k=0;k<iBndWp.size();k++){
        int i = iBndWp[k];
        int j = jBndWp[k];

        if(sBndWp[k] == 0){ //ディリクレ
            aPp[i][j] = aPp[i][j] + aWp[i][j];
            aWp[i][j] = 0.0;
        }else if(sBndWp[k] == 1){ //開放
            double c0 = sqrt(1.0/(eps0*epsr[i][j]*mu0));
            double deno = pow(2.0*c0,2) + pow(omegam*dx,2);
            std::complex<double> R = std::complex<double>
                (pow(2.0*c0,2) - pow(omegam*dx,2),-4.0*c0*omegam*dx);
            R = R/(deno+1e-100);
            aPp[i][j] = aPp[i][j] - aWp[i][j]*R;
            aWp[i][j] = 0.0;
        }
    }
    //********* Ep境界条件 (右) *********
    for (int k=0;k<iBndEp.size();k++){
        int i = iBndEp[k];
        int j = jBndEp[k];

        if(sBndEp[k] == 0){ //ディリクレ
            aPp[i][j] = aPp[i][j] + aEp[i][j];
            aEp[i][j] = 0.0;
        }
    }
    //********* Ep境界条件 (下) *********
    for (int k=0;k<iBndSp.size();k++){
        int i = iBndSp[k];
        int j = jBndSp[k];

        if(sBndSp[k] == 0){ //ディリクレ
            aPp[i][j] = aPp[i][j] + r[j]/r[j-1]*aSp[i][j];
            aSp[i][j] = 0.0;
        }
    }
    //********* Ep境界条件 (上) *********
    for (int k=0;k<iBndNp.size();k++){
        int i = iBndNp[k];
        int j = jBndNp[k];

        if(sBndNp[k] == 0){ //ディリクレ
            aPp[i][j] = aPp[i][j] + r[j]/r[j+1]*aNp[i][j];
            aNp[i][j] = 0.0;
        }else if(sBndNp[k] == 1){ //開放
            double c0 = sqrt(1.0/(eps0*epsr[i][j]*mu0));
            double r_tmp = (r[j+1] + r[j])/2.0;
            double qL = r[j]/(r_tmp+1e-100);
            double qR = r[j+1]  /(r_tmp+1e-100);
            double deno = pow(2.0*c0*sqrt(qR),2) + pow(omegam*dr*qR,2);
            std::complex<double> R = std::complex<double>
                (pow(2.0*c0,2)*sqrt(qL*qR) - pow(omegam*dr,2)*qL*qR
                ,-2.0*c0*omegam*dr*(qL*sqrt(qR)+qR*sqrt(qL)));
            R = R/(deno+1e-100);

            aPp[i][j] = aPp[i][j] - aNp[i][j]*R;
            aNp[i][j] = 0.0;
        }
    }

    /******************** LUソルバーの準備 ********************/
    Eigen::SparseMatrix<std::complex<double> > A(nk, nk);
    A.reserve(Eigen::VectorXi::Constant(nk,11)); //ここの数字を変えて帯域幅を確保する
    Eigen::VectorXcd b(nk);
    Eigen::VectorXcd xv(nk);

    //Ex-足し込み
    for (int k=0;k<nkx;k++){

        int i = ikx[k];
        int j = jkx[k];

        double kE = kx[i+1][j];
        double kW = kx[i-1][j];
        double kN = kx[i][j+1];
        double kS = kx[i][j-1];

        double kNE = kr[i][j+1];
        double kSE = kr[i][j];
        double kNW = kr[i-1][j+1];
        double kSW = kr[i-1][j];

        double kEE = kp[i][j];
        double kWW = kp[i-1][j];
        
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
            A.insert(k, nkx+kNE) = -aNEx[i][j];
        }
        if(kSE!=-1){
            A.insert(k, nkx+kSE) = -aSEx[i][j];
        }
        if(kNW!=-1){
            A.insert(k, nkx+kNW) = -aNWx[i][j];
        }
        if(kSW!=-1){
            A.insert(k, nkx+kSW) = -aSWx[i][j];
        }
        //Ep用
        if(kEE!=-1){
            A.insert(k, nkx+nkr+kEE) = -aEEx[i][j];
        }
        if(kWW!=-1){
            A.insert(k, nkx+nkr+kWW) = -aWWx[i][j];
        }
    }
    
    //Er-足し込み
    for (int k=0;k<nkr;k++){
        int i = ikr[k];
        int j = jkr[k];

        double kE  = kr[i+1][j];
        double kW  = kr[i-1][j];
        double kN  = kr[i][j+1];
        double kS  = kr[i][j-1];

        double kNE = kx[i+1][j];
        double kSE = kx[i+1][j-1];
        double kNW = kx[i][j];
        double kSW = kx[i][j-1];

        double kNN = kp[i][j];
        double kSS = kp[i][j-1];

        //Er用
        A.insert(nkx+k, nkx+k)      =   aPr[i][j];
        b[nkx+k]                    =   br[i][j];
        if(kN!=-1){
            A.insert(nkx+k, nkx+kN) = -aNr[i][j];
        }
        if(kS!=-1){
            A.insert(nkx+k, nkx+kS) = -aSr[i][j];
        }
        if(kE!=-1){
            A.insert(nkx+k, nkx+kE) = -aEr[i][j];
        }
        if(kW!=-1){
            A.insert(nkx+k, nkx+kW) = -aWr[i][j];
        }
        //Ex用
        if(kNE!=-1){
            A.insert(nkx+k, kNE) = -aNEr[i][j];
        }
        if(kSE!=-1){
            A.insert(nkx+k, kSE) = -aSEr[i][j];
        }
        if(kNW!=-1){
            A.insert(nkx+k, kNW) = -aNWr[i][j];
        }
        if(kSW!=-1){
            A.insert(nkx+k, kSW) = -aSWr[i][j];
        }
        //Ep用
        if(kNN!=-1){
            A.insert(nkx+k, nkx+nkr+kNN) = -aNNr[i][j];
        }
        if(kSS!=-1){
            A.insert(nkx+k, nkx+nkr+kSS) = -aSSr[i][j];
        }
    }

    //Ep-足し込み
    for (int k=0;k<nkp;k++){
        int i = ikp[k];
        int j = jkp[k];

        double kE  = kp[i+1][j];
        double kW  = kp[i-1][j];
        double kN  = kp[i][j+1];
        double kS  = kp[i][j-1];

        double kEE  = kx[i+1][j];
        double kWW  = kx[i][j];

        double kNN  = kr[i][j+1];
        double kSS  = kr[i][j];

        //Ep用
        A.insert(nkx+nkr+k, nkx+nkr+k)      =   aPp[i][j];
        b[nkx+nkr+k]                    =   bp[i][j];
        if(kN!=-1){
            A.insert(nkx+nkr+k, nkx+nkr+kN) = -aNp[i][j];
        }
        if(kS!=-1){
            A.insert(nkx+nkr+k, nkx+nkr+kS) = -aSp[i][j];
        }
        if(kE!=-1){
            A.insert(nkx+nkr+k, nkx+nkr+kE) = -aEp[i][j];
        }
        if(kW!=-1){
            A.insert(nkx+nkr+k, nkx+nkr+kW) = -aWp[i][j];
        }
        //Ex用
        if(kEE!=-1){
            A.insert(nkx+nkr+k, kEE) = -aEEp[i][j];
        }
        if(kWW!=-1){
            A.insert(nkx+nkr+k, kWW) = -aWWp[i][j];
        }
        //Er用
        if(kNN!=-1){
            A.insert(nkx+nkr+k, nkx+kNN) = -aNNp[i][j];
        }
        if(kSS!=-1){
            A.insert(nkx+nkr+k, nkx+kSS) = -aSSp[i][j];
        }
    }
 

    // LU分解でAx = bを解く
    Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double> > > solver;
    //std::cout << "solver start!" << std::endl;
    solver.compute(A);
    xv = solver.solve(b);

    /*
    // 実際の誤差を計算
    Eigen::VectorXcd residual = A*xv- b;
    std::cout << residual<< std::endl;
    double actual_error = residual.norm()/(b.norm()+1e-100);
    std::cout << "Actual error: " << actual_error << std::endl;
    */

    //Ex-結果を戻す
    for (int k=0;k<nkx;k++){
        int i = ikx[k];
        int j = jkx[k];
        E1x[i][j] = xv[k];
    }
    //Er-結果を戻す
    for (int k=0;k<nkr;k++){
        int i = ikr[k];
        int j = jkr[k];
        E1r[i][j] = xv[nkx+k];
    }
    //Ep-結果を戻す
    for (int k=0;k<nkp;k++){
        int i = ikp[k];
        int j = jkp[k];
        E1p[i][j] = xv[nkx+nkr+k];
    }

    /************************境界条件後処理**************************/
    //Ex-境界条件
    //********* Ex境界条件 (左) *********
    for (int k=0;k<iBndWx.size();k++){
        int i = iBndWx[k];
        int j = jBndWx[k];

        if(sBndWx[k] == 0){ //ディリクレ
            E1x[i-1][j] = 0.0;
        }else if(sBndWx[k] == 2){ //ガウス
            double qL = (r[j-1] + r[j])/2.0/(r[j]+1e-100);
            double qR = (r[j+1] + r[j])/2.0/(r[j]+1e-100);

            E1x[i-1][j] = E1x[i][j] + dx/dr*(qR*E1r[i-1][j+1] - qL*E1r[i-1][j]);
        }else if(sBndWx[k] == 3){ //凹角 下側壁
            double qL = (r[j-1] + r[j])/2.0/(r[j]+1e-100);
            double qR = (r[j+1] + r[j])/2.0/(r[j]+1e-100);
            double deno = 1.0 + dx/dr*qL;

            E1x[i-1][j] = E1x[i][j]/deno + dx/dr*(qR*E1r[i-1][j+1])/deno;
        }else if(sBndWx[k] == 4){ //凹角 上側壁
            double qL = (r[j-1] + r[j])/2.0/(r[j]+1e-100);
            double qR = (r[j+1] + r[j])/2.0/(r[j]+1e-100);
            double deno = 1.0 + dx/dr*qR;

            E1x[i-1][j] = E1x[i][j]/deno - dx/dr*(qL*E1r[i-1][j])/deno;
        }
    }
    //********* Ex境界条件 (右) *********
    for (int k=0;k<iBndEx.size();k++){
        int i = iBndEx[k];
        int j = jBndEx[k];

        if(sBndEx[k] == 2){ //ガウス
            double qL = (r[j-1] + r[j])/2.0/(r[j]+1e-100);
            double qR = (r[j+1] + r[j])/2.0/(r[j]+1e-100);
            
            E1x[i+1][j] = E1x[i][j] - dx/dr*(qR*E1r[i][j+1] - qL*E1r[i][j]);
        }else if(sBndEx[k] == 4){ //凹角 上側壁
            double qL = (r[j-1] + r[j])/2.0/(r[j]+1e-100);
            double qR = (r[j+1] + r[j])/2.0/(r[j]+1e-100);
            double deno = 1.0 + dx/dr*qR;

            E1x[i+1][j] = E1x[i][j]/deno + dx/dr*qL*E1r[i][j]/deno;
        }
    }
    //********* Ex境界条件 (下) *********
    for (int k=0;k<iBndSx.size();k++){
        int i = iBndSx[k];
        int j = jBndSx[k];
        if(sBndSx[k] == 0){ //ディリクレ
            E1x[i][j-1] = -r[j]/r[j-1]*E1x[i][j];
        }else if(sBndSx[k] == 6){ //凸角 右側Open
            //処理なし
        }
    }
    //********* Ex境界条件 (上) *********
    for (int k=0;k<iBndNx.size();k++){
        int i = iBndNx[k];
        int j = jBndNx[k];
        if(sBndNx[k] == 0){ //ディリクレ
            E1x[i][j+1] = -r[j]/r[j+1]*E1x[i][j];
        }else if(sBndNx[k] == 1){ //開放
            double c0 = sqrt(1.0/(eps0*epsr[i][j]*mu0));
            double r_tmp = (r[j+1] + r[j])/2.0;
            double qL = r[j]/(r_tmp+1e-100);
            double qR = r[j+1]  /(r_tmp+1e-100);
            double deno = pow(2.0*c0*sqrt(qR),2) + pow(omegam*dr*qR,2);
            std::complex<double> R = std::complex<double>
                (pow(2.0*c0,2)*sqrt(qL*qR) - pow(omegam*dr,2)*qL*qR
                ,-2.0*c0*omegam*dr*(qL*sqrt(qR)+qR*sqrt(qL)));
            R = R/(deno+1e-100);
            E1x[i][j+1] = R*E1x[i][j];
        }else if(sBndNx[k] == 5){ //凸角 左側Open
            //処理なし
        }else if(sBndNx[k] == 6){ //凸角 右側Open
            //処理なし
        }
    }

    //********* Er境界条件 (左) *********
    for (int k=0;k<iBndWr.size();k++){
        int i = iBndWr[k];
        int j = jBndWr[k];

        if(sBndWr[k] == 0){ //ディリクレ
            E1r[i-1][j] = -E1r[i][j];
        }else if(sBndWr[k] == 1){ //開放
            double c0 = sqrt(1.0/(eps0*epsr[i][j]*mu0));
            double deno = pow(2.0*c0,2) + pow(omegam*dx,2);
            std::complex<double> R = std::complex<double>
                (pow(2.0*c0,2) - pow(omegam*dx,2),-4.0*c0*omegam*dx);
            R = R/(deno+1e-100);

            E1r[i-1][j] = R*E1r[i][j];
        }else if(sBndWr[k] == 5){ //凸角 下側Open
            //処理なし
        }else if(sBndWr[k] == 6){ //凸角 上側Open
            //処理なし
        }
    }
    //********* Er境界条件 (右) *********
    for (int k=0;k<iBndEr.size();k++){
        int i = iBndEr[k];
        int j = jBndEr[k];

        if(sBndEr[k] == 0){ //ディリクレ
            E1r[i+1][j] = -E1r[i][j];
        }else if(sBndEr[k] == 5){ //凸角 下側Open
            //処理なし
        }
    }
    //********* Er境界条件 (下) *********
    for (int k=0;k<iBndSr.size();k++){
        int i = iBndSr[k];
        int j = jBndSr[k];

        if(sBndSr[k] == 2){ //ガウス
            double rR = (r[j] + r[j-1])/2.0;
            double rL = (r[j-1] + r[j-2])/2.0;
            
            E1r[i][j-1] = rR/rL*E1r[i][j] + dr/dx*r[j-1]/rL*(E1x[i+1][j-1]-E1x[i][j-1]);
        }else if(sBndSr[k] == 3){ //凹角 左側壁
            double rR = (r[j] + r[j-1])/2.0;
            double rL = (r[j-1] + r[j-2])/2.0;
            double deno = dr/dx*r[j-1]/rL;
            if(r[j] < j_flc_bl[0][0]){ //マイクロ波の入り口部の特別処理
                deno = 1.0;
            }

            E1r[i][j-1] = rR/rL*E1r[i][j]/deno + dr/dx*r[j-1]/rL*(E1x[i+1][j-1])/deno;
        }
    }
    //********* Er境界条件 (上) *********
    for (int k=0;k<iBndNr.size();k++){
        int i = iBndNr[k];
        int j = jBndNr[k];

        if(sBndNr[k] == 2){ //ガウス
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;

            E1r[i][j+1] = rL/rR*E1r[i][j] - dr/dx*r[j]/rR*(E1x[i+1][j]-E1x[i][j]);
        }if(sBndNr[k] == 3){ //凹角 左側壁
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;
            double deno = dr/dx*r[j]/rR;

            E1r[i][j+1] = rL/rR*E1r[i][j]/deno - dr/dx*r[j]/rR*(E1x[i+1][j])/deno;
        }if(sBndNr[k] == 4){ //凹角 右側壁
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;
            double deno = dr/dx*r[j]/rR;

            E1r[i][j+1] = rL/rR*E1r[i][j]/deno + dr/dx*r[j]/rR*(E1x[i][j])/deno;
        }
    }

    //********* Ep境界条件 (左) *********
    for (int k=0;k<iBndWp.size();k++){
        int i = iBndWp[k];
        int j = jBndWp[k];

        if(sBndWp[k] == 0){ //ディリクレ
            E1p[i-1][j] = -E1p[i][j];
        }else if(sBndWp[k] == 1){ //開放
            double c0 = sqrt(1.0/(eps0*epsr[i][j]*mu0));
            double deno = pow(2.0*c0,2) + pow(omegam*dx,2);
            std::complex<double> R = std::complex<double>
                (pow(2.0*c0,2) - pow(omegam*dx,2),-4.0*c0*omegam*dx);
            R = R/(deno+1e-100);

            E1p[i-1][j] = R*E1p[i][j];
        }
    }
    //********* Ep境界条件 (右) *********
    for (int k=0;k<iBndEp.size();k++){
        int i = iBndEp[k];
        int j = jBndEp[k];

        E1p[i+1][j] = -E1p[i][j];
    }
    //********* Ep境界条件 (下) *********
    for (int k=0;k<iBndSp.size();k++){
        int i = iBndSp[k];
        int j = jBndSp[k];

        E1p[i][j-1] = -r[j]/r[j-1]*E1p[i][j];
    }
    //********* Ep境界条件 (上) *********
    for (int k=0;k<iBndNp.size();k++){
        int i = iBndNp[k];
        int j = jBndNp[k];

        if(sBndNp[k] == 0){ //ディリクレ
            E1p[i][j+1] = -r[j]/r[j+1]*E1p[i][j];
        }else if(sBndNp[k] == 1){ //開放
            double c0 = sqrt(1.0/(eps0*epsr[i][j]*mu0));
            double r_tmp = (r[j+1] + r[j])/2.0;
            double qL = r[j]/(r_tmp+1e-100);
            double qR = r[j+1]  /(r_tmp+1e-100);
            double deno = pow(2.0*c0*sqrt(qR),2) + pow(omegam*dr*qR,2);
            std::complex<double> R = std::complex<double>
                (pow(2.0*c0,2)*sqrt(qL*qR) - pow(omegam*dr,2)*qL*qR
                ,-2.0*c0*omegam*dr*(qL*sqrt(qR)+qR*sqrt(qL)));
            R = R/(deno+1e-100);
            
            E1p[i][j+1] = R*E1p[i][j];
        }
    }
}

//*****************************************************************
//**                                                             **
//**           void solve_Efield                                 **
//**                                                             **
//*****************************************************************
void solve_Microwave_impedanceTest_old(){

    //テスト用
    for (int i=0;i<ni+2;i++){
        for (int j=0;j<nj+2;j++){
            //rhoe[i][j] = 0.0;
            nu_m1[i][j] = nu_eff;
        }
    }

    //Ex用 係数
    std::vector<std::vector<std::complex<double> > > aPx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aWx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aEx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aNx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aSx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aNEx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aNWx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aSEx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aSWx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aEEx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aWWx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > bx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数

    //Er用 係数
    std::vector<std::vector<std::complex<double> > > aPr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aWr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aEr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aNr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aSr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aNEr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aNWr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aSEr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aSWr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aNNr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aSSr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > br(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数

    //Ephi用 係数
    std::vector<std::vector<std::complex<double> > > aPp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aWp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aEp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aNp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aSp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aNNp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aSSp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aEEp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aWWp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > bp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数

    /******************** 係数作成 いらないところも含めて普通に全部やる方が速い ********************/
    std::complex<double> EPS(1e-100,1e-100); //微小複素数
    //for Ex
    for (int iblock=5;iblock<=5;iblock++){ 
        for (int i=i_flx_bl[iblock][0];i<=i_flx_bl[iblock][1];i++){ 
            for (int j=j_flx_bl[iblock][0];j<=j_flx_bl[iblock][1];j++){
                double qL = (r[j-1] + r[j])/2.0/(r[j]+1e-100);
                double qR = (r[j+1] + r[j])/2.0/(r[j]+1e-100);
                
                double nu_m_tmp = (nu_m1[i][j] + nu_m1[i-1][j])/2.0;
                double rho_tmp  = (rhoe[i][j]  + rhoe[i-1][j] )/2.0;
                double Bx_tmp   = (Bx[i][j]   + Bx[i-1][j]  )/2.0;
                double Br_tmp   = (Br[i][j]   + Br[i-1][j]  )/2.0;

                std::complex<double> nu_cmp(nu_m_tmp,omegam);
                std::complex<double> Hallcmpx = e0*Bx_tmp/(masse*nu_cmp + EPS);
                std::complex<double> Hallcmpr = e0*Br_tmp/(masse*nu_cmp + EPS);
                std::complex<double> HallcmpMag2 =  Hallcmpx*Hallcmpx + Hallcmpr*Hallcmpr;
                std::complex<double> coef = rho_tmp*e0*e0/(masse*nu_cmp + EPS)/(1.0 + HallcmpMag2 + EPS);
                std::complex<double> sigmaexr = coef*(Hallcmpx*Hallcmpr);
                std::complex<double> sigmaexp = coef*(Hallcmpr);
                std::complex<double> sigmaexx = coef*(1.0 + Hallcmpx*Hallcmpx);
                std::complex<double> iomegaMu0(0.0,omegam*mu0);

                aPx[i][j] = 2.0*dx*dx + 2.0*dr*dr
                    - eps0*epsr[i][j]*mu0*omegam*omegam*dx*dx*dr*dr + iomegaMu0*sigmaexx*dx*dx*dr*dr;
                aEx[i][j] = dr*dr;
                aWx[i][j] = dr*dr;
                aNx[i][j] = qR*dx*dx;
                aSx[i][j] = qL*dx*dx;
                bx[i][j] = std::complex<double>(0,-omegam*mu0*dx*dx*dr*dr)*J1x[i][j];

                //Epの要素
                aEEx[i][j] = -iomegaMu0*sigmaexp*0.5*dx*dx*dr*dr;
                aWWx[i][j] = -iomegaMu0*sigmaexp*0.5*dx*dx*dr*dr;

                //Erの要素
                aNEx[i][j] = -iomegaMu0*sigmaexr*0.25*qR*dx*dx*dr*dr;
                aNWx[i][j] = -iomegaMu0*sigmaexr*0.25*qR*dx*dx*dr*dr;
                aSEx[i][j] = -iomegaMu0*sigmaexr*0.25*qL*dx*dx*dr*dr;
                aSWx[i][j] = -iomegaMu0*sigmaexr*0.25*qL*dx*dx*dr*dr;
            }
        }
    }

    //for Er
    for (int iblock=5;iblock<=5;iblock++){ 
        for (int i=i_flr_bl[iblock][0];i<=i_flr_bl[iblock][1];i++){ 
            for (int j=j_flr_bl[iblock][0];j<=j_flr_bl[iblock][1];j++){
                double r_tmp = (r[j-1] + r[j])/2.0;
                double qL = r[j-1]/(r_tmp+1e-100);
                double qR = r[j]  /(r_tmp+1e-100);

                double nu_m_tmp = (nu_m1[i][j] + nu_m1[i][j-1])/2.0;
                double rho_tmp  = (rhoe[i][j]  + rhoe[i][j-1] )/2.0;
                double Bx_tmp   = (r[j]*Bx[i][j]   + r[j-1]*Bx[i][j-1]  )/(2.0*r_tmp);
                double Br_tmp   = (r[j]*Br[i][j]   + r[j-1]*Br[i][j-1]  )/(2.0*r_tmp);

                std::complex<double> nu_cmp(nu_m_tmp,omegam);
                std::complex<double> Hallcmpx = e0*Bx_tmp/(masse*nu_cmp + EPS);
                std::complex<double> Hallcmpr = e0*Br_tmp/(masse*nu_cmp + EPS);
                std::complex<double> HallcmpMag2 =  Hallcmpx*Hallcmpx + Hallcmpr*Hallcmpr;
                std::complex<double> coef = rho_tmp*e0*e0/(masse*nu_cmp + EPS)/(1.0 + HallcmpMag2 + EPS);
                std::complex<double> sigmaerr = coef*(1.0+Hallcmpr*Hallcmpr);
                std::complex<double> sigmaerp = coef*(-Hallcmpx);
                std::complex<double> sigmaerx = coef*(Hallcmpx*Hallcmpr);
                std::complex<double> iomegaMu0(0.0,omegam*mu0);

                aPr[i][j] = 2.0*dx*dx + 2.0*dr*dr + dx*dx*dr*dr/(r_tmp*r_tmp)
                    - eps0*epsr[i][j]*mu0*omegam*omegam*dx*dx*dr*dr + iomegaMu0*sigmaerr*dx*dx*dr*dr;
                aEr[i][j] = dr*dr;
                aWr[i][j] = dr*dr;
                aNr[i][j] = qR*dx*dx;
                aSr[i][j] = qL*dx*dx;
                br[i][j] = std::complex<double>(0,-omegam*mu0*dx*dx*dr*dr)*J1r[i][j];

                //Epの要素
                aNNr[i][j] = -iomegaMu0*sigmaerp*0.5*qR*dx*dx*dr*dr;
                aSSr[i][j] = -iomegaMu0*sigmaerp*0.5*qL*dx*dx*dr*dr;

                //Exの要素
                aNEr[i][j] = -iomegaMu0*sigmaerx*0.25*qR*dx*dx*dr*dr;
                aNWr[i][j] = -iomegaMu0*sigmaerx*0.25*qR*dx*dx*dr*dr;
                aSEr[i][j] = -iomegaMu0*sigmaerx*0.25*qL*dx*dx*dr*dr;
                aSWr[i][j] = -iomegaMu0*sigmaerx*0.25*qL*dx*dx*dr*dr;
            }
        }
    }

    //for Ephi
    for (int iblock=5;iblock<=5;iblock++){ 
        for (int i=i_flc_bl[iblock][0];i<=i_flc_bl[iblock][1];i++){ 
            for (int j=j_flc_bl[iblock][0];j<=j_flc_bl[iblock][1];j++){
                double qL = (r[j-1] + r[j])/2.0/(r[j]+1e-100);
                double qR = (r[j+1] + r[j])/2.0/(r[j]+1e-100);

                double nu_m_tmp = nu_m1[i][j];
                double rho_tmp  = rhoe[i][j];
                double Bx_tmp   = Bx[i][j];
                double Br_tmp   = Br[i][j];

                std::complex<double> nu_cmp(nu_m_tmp,omegam);
                std::complex<double> Hallcmpx = e0*Bx_tmp/(masse*nu_cmp + EPS);
                std::complex<double> Hallcmpr = e0*Br_tmp/(masse*nu_cmp + EPS);
                std::complex<double> HallcmpMag2 =  Hallcmpx*Hallcmpx + Hallcmpr*Hallcmpr;
                std::complex<double> coef = rho_tmp*e0*e0/(masse*nu_cmp + EPS)/(1.0 + HallcmpMag2 + EPS);
                std::complex<double> sigmaepr = coef*(Hallcmpx);
                std::complex<double> sigmaepp = coef*(1.0);
                std::complex<double> sigmaepx = coef*(-Hallcmpr);
                std::complex<double> iomegaMu0(0.0,omegam*mu0);

                aPp[i][j] = 2.0*dx*dx + 2.0*dr*dr + dx*dx*dr*dr/(r[j]*r[j])
                    - eps0*epsr[i][j]*mu0*omegam*omegam*dx*dx*dr*dr + iomegaMu0*sigmaepp*dx*dx*dr*dr;
                aEp[i][j] = dr*dr;
                aWp[i][j] = dr*dr;
                aNp[i][j] = qR*dx*dx;
                aSp[i][j] = qL*dx*dx;
                bp[i][j] = std::complex<double>(0,-omegam*mu0*dx*dx*dr*dr)*J1p[i][j];

                //Exの要素
                aEEp[i][j] = -iomegaMu0*sigmaepx*0.5*dx*dx*dr*dr;
                aWWp[i][j] = -iomegaMu0*sigmaepx*0.5*dx*dx*dr*dr;

                //Erの要素
                aNNp[i][j] = -iomegaMu0*sigmaepr*0.5*qR*dx*dx*dr*dr;
                aSSp[i][j] = -iomegaMu0*sigmaepr*0.5*qL*dx*dx*dr*dr;
            }
        }
    }

    /******************** 境界条件設定 (係数修正) ********************/
    //********* Ex境界条件 (左) *********
    for (int k=0;k<iBndWx.size();k++){
        int i = iBndWx[k];
        int j = jBndWx[k];

        if(sBndWx[k] == 0){ //ディリクレ
            aWx[i][j] = 0.0;
        }else if(sBndWx[k] == 2){ //ガウス
            double qL = (r[j-1] + r[j])/2.0/(r[j]+1e-100);
            double qR = (r[j+1] + r[j])/2.0/(r[j]+1e-100);

            aPx[i][j] = aPx[i][j] - aWx[i][j];
            aNWx[i][j] = aNWx[i][j] + aWx[i][j]*dx/dr*qR;
            aSWx[i][j] = aSWx[i][j] - aWx[i][j]*dx/dr*qL;
            aWx[i][j] = 0.0;
        }else if(sBndWx[k] == 3){ //凹角 下側壁
            double qL = (r[j-1] + r[j])/2.0/(r[j]+1e-100);
            double qR = (r[j+1] + r[j])/2.0/(r[j]+1e-100);
            double deno = 1.0 + dx/dr*qL;

            aPx[i][j] = aPx[i][j] - aWx[i][j]/deno;
            aNWx[i][j] = aNWx[i][j] + aWx[i][j]*(dx/dr*qR)/deno;
            aWx[i][j] = 0.0;
        }else if(sBndWx[k] == 4){ //凹角 上側壁
            double qL = (r[j-1] + r[j])/2.0/(r[j]+1e-100);
            double qR = (r[j+1] + r[j])/2.0/(r[j]+1e-100);
            double deno = 1.0 + dx/dr*qR;

            aPx[i][j] = aPx[i][j] - aWx[i][j]/deno;
            aSWx[i][j] = aSWx[i][j] - aWx[i][j]*(dx/dr*qL)/deno;
            aWx[i][j] = 0.0;
        }
    }
    //********* Ex境界条件 (右) *********
    for (int k=0;k<iBndEx.size();k++){
        int i = iBndEx[k];
        int j = jBndEx[k];

        if(sBndEx[k] == 0){ //ディリクレ
            aEx[i][j] = 0.0;
        }else if(sBndEx[k] == 2){ //ガウス
            double qL = (r[j-1] + r[j])/2.0/(r[j]+1e-100);
            double qR = (r[j+1] + r[j])/2.0/(r[j]+1e-100);

            aPx[i][j] = aPx[i][j] - aEx[i][j];
            aNEx[i][j] = aNEx[i][j] - aEx[i][j]*dx/dr*qR;
            aSEx[i][j] = aSEx[i][j] + aEx[i][j]*dx/dr*qL;
            aEx[i][j] = 0.0;
        }if(sBndEx[k] == 3){ //凹角 下側壁
            double qL = (r[j-1] + r[j])/2.0/(r[j]+1e-100);
            double qR = (r[j+1] + r[j])/2.0/(r[j]+1e-100);
            double deno = 1.0 + dx/dr*qL;

            aPx[i][j] = aPx[i][j] - aEx[i][j]/deno;
            aNEx[i][j] = aNEx[i][j] - aEx[i][j]*(dx/dr*qR)/deno;
            aEx[i][j] = 0.0;
        }else if(sBndEx[k] == 4){ //凹角 上側壁
            double qL = (r[j-1] + r[j])/2.0/(r[j]+1e-100);
            double qR = (r[j+1] + r[j])/2.0/(r[j]+1e-100);
            double deno = 1.0 + dx/dr*qR;
            
            aPx[i][j] = aPx[i][j] - aEx[i][j]/deno;
            aSEx[i][j] = aSEx[i][j] + aEx[i][j]*(dx/dr*qL)/deno;
            aEx[i][j] = 0.0;
        }
    }
    //********* Ex境界条件 (下) *********
    for (int k=0;k<iBndSx.size();k++){
        int i = iBndSx[k];
        int j = jBndSx[k];
        if(sBndSx[k] == 0){ //ディリクレ
            aPx[i][j] = aPx[i][j] + r[j]/r[j-1]*aSx[i][j];
            aSx[i][j] = 0.0;
        }else if(sBndSx[k] == 6){ //凸角 右側Open
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;

            aPx[i][j] = aPx[i][j] - aSx[i][j]*( dr/dx*r[j]/rL);
            aWx[i][j] = aWx[i][j] + aSx[i][j]*(-dr/dx*r[j]/rL);
            aNWx[i][j] = aNWx[i][j] + aSx[i][j]*(rR/rL);
            aSx[i][j] = 0.0;
        }
    }
    //********* Ex境界条件 (上) *********
    for (int k=0;k<iBndNx.size();k++){
        int i = iBndNx[k];
        int j = jBndNx[k];
        if(sBndNx[k] == 0){ //ディリクレ
            aPx[i][j] = aPx[i][j] + r[j]/r[j+1]*aNx[i][j];
            aNx[i][j] = 0.0;
        }if(sBndNx[k] == 1){ //開放
            double c0 = sqrt(1.0/(eps0*epsr[i][j]*mu0));
            double r_tmp = (r[j+1] + r[j])/2.0;
            double qL = r[j]/(r_tmp+1e-100);
            double qR = r[j+1]  /(r_tmp+1e-100);
            double deno = pow(2.0*c0*sqrt(qR),2) + pow(omegam*dr*qR,2);
            std::complex<double> R = std::complex<double>
                (pow(2.0*c0,2)*sqrt(qL*qR) - pow(omegam*dr,2)*qL*qR
                ,-2.0*c0*omegam*dr*(qL*sqrt(qR)+qR*sqrt(qL)));
            R = R/(deno+1e-100);

            aPx[i][j] = aPx[i][j] - aNx[i][j]*R;
            aNx[i][j] = 0.0;
        }else if(sBndNx[k] == 5){ //凸角 左側Open (*)
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;

            aPx[i][j] = aPx[i][j] - aNx[i][j]*( dr/dx*r[j]/rR);
            aEx[i][j] = aEx[i][j] + aNx[i][j]*(-dr/dx*r[j]/rR);
            aSEx[i][j] = aSEx[i][j] + aNx[i][j]*( rL/rR);
            aNx[i][j] = 0.0;
        }else if(sBndNx[k] == 6){ //凸角 右側Open
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;

            aPx[i][j] = aPx[i][j] - aNx[i][j]*( dr/dx*r[j]/rR);
            aWx[i][j] = aWx[i][j] + aNx[i][j]*(-dr/dx*r[j]/rR);
            aSWx[i][j] = aSWx[i][j] + aNx[i][j]*(-rL/rR);
            aNx[i][j] = 0.0;
        }
    }

    //********* Er境界条件 (左) *********
    for (int k=0;k<iBndWr.size();k++){
        int i = iBndWr[k];
        int j = jBndWr[k];

        if(sBndWr[k] == 0){ //ディリクレ
            aPr[i][j] = aPr[i][j] + aWr[i][j];
            aWr[i][j] = 0.0;
        }else if(sBndWr[k] == 1){ //開放
            double c0 = sqrt(1.0/(eps0*epsr[i][j]*mu0));
            double deno = pow(2.0*c0,2) + pow(omegam*dx,2);
            std::complex<double> R = std::complex<double>(pow(2.0*c0,2) - pow(omegam*dx,2),-4.0*c0*omegam*dx);
            R = R/(deno+1e-100);
            
            aPr[i][j] = aPr[i][j] - aWr[i][j]*R;
            aWr[i][j] = 0.0;
        }else if(sBndWr[k] == 5){ //凸角 下側Open
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;
            
            aPr[i][j] = aPr[i][j] - aWr[i][j]*( dx/dr*rL/r[j]);
            aNr[i][j] = aNr[i][j] + aWr[i][j]*(-dx/dr*rR/r[j]);
            aNEr[i][j] = aNEr[i][j] + aWr[i][j]*(-1.0);
            aWr[i][j] = 0.0;
        }else if(sBndWr[k] == 6){ //凸角 上側Open
            double rR = (r[j] + r[j-1])/2.0;
            double rL = (r[j-1] + r[j-2])/2.0;

            aPr[i][j] = aPr[i][j] - aWr[i][j]*( dx/dr*rR/r[j-1]);
            aSr[i][j] = aSr[i][j] + aWr[i][j]*(-dx/dr*rL/r[j-1]);
            aSEr[i][j] = aSEr[i][j] + aWr[i][j]*( 1.0);
            aWr[i][j] = 0.0;
        }
    }
    //********* Er境界条件 (右) *********
    for (int k=0;k<iBndEr.size();k++){
        int i = iBndEr[k];
        int j = jBndEr[k];

        if(sBndEr[k] == 0){ //ディリクレ
            aPr[i][j] = aPr[i][j] + aEr[i][j];
            aEr[i][j] = 0.0;
        }if(sBndEr[k] == 1){ //開放
            double c0 = sqrt(1.0/(eps0*mu0));
            
            //std::complex<double> deno(2.0*c0,omegam*dx);
            //std::complex<double> R(2.0*c0,-omegam*dx);
            //aPr[i][j] = aPr[i][j] - aEr[i][j]*R/(deno+1e-100);
            //aEr[i][j] = 0.0;

            double deno = pow(2.0*c0,2) + pow(omegam*dx,2);
            std::complex<double> R = std::complex<double>
                (pow(2.0*c0,2) - pow(omegam*dx,2),-4.0*c0*omegam*dx);
            aPr[i][j] = aPr[i][j] - aEr[i][j]*R/(deno+1e-100);
            aEr[i][j] = 0.0;


        }else if(sBndEr[k] == 5){ //凸角 下側Open (*)
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;

            aPr[i][j] = aPr[i][j] - aEr[i][j]*( dx/dr*rL/r[j]);
            aNr[i][j] = aNr[i][j] + aEr[i][j]*(-dx/dr*rR/r[j]);
            aNWr[i][j] = aNWr[i][j] + aEr[i][j]*( 1.0);
            aEr[i][j] = 0.0;
        }
    }
    //********* Er境界条件 (下) *********
    for (int k=0;k<iBndSr.size();k++){
        int i = iBndSr[k];
        int j = jBndSr[k];

        if(sBndSr[k] == 2){ //ガウス
            double rR = (r[j] + r[j-1])/2.0;
            double rL = (r[j-1] + r[j-2])/2.0;
            aPr[i][j] = aPr[i][j] - aSr[i][j]*rR/rL;
            aSEr[i][j] = aSEr[i][j] + aSr[i][j]*dr/dx*r[j-1]/rL;
            aSWr[i][j] = aSWr[i][j] - aSr[i][j]*dr/dx*r[j-1]/rL;
            aSr[i][j] = 0.0;
        }else if(sBndSr[k] == 3){ //凹角 左側壁
            double rR = (r[j] + r[j-1])/2.0;
            double rL = (r[j-1] + r[j-2])/2.0;
            double deno = dr/dx*r[j-1]/rL;
            if(r[j] < j_flc_bl[0][0]){ //マイクロ波の入り口部の特別処理
                deno = 1.0;
            }

            aPr[i][j] = aPr[i][j] - aSr[i][j]*rR/rL/deno;
            aNEr[i][j] = aNEr[i][j] + aSr[i][j]*dr/dx*r[j-1]/rL/deno;
            aSr[i][j] = 0.0;
        }
    }
    //********* Er境界条件 (上) *********
    for (int k=0;k<iBndNr.size();k++){
        int i = iBndNr[k];
        int j = jBndNr[k];

        if(sBndNr[k] == 2){ //ガウス
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;
            aPr[i][j] = aPr[i][j] - aNr[i][j]*rL/rR;
            aNEr[i][j] = aNEr[i][j] - aNr[i][j]*dr/dx*r[j]/rR;
            aNWr[i][j] = aNWr[i][j] + aNr[i][j]*dr/dx*r[j]/rR;
            aNr[i][j] = 0.0;
        }if(sBndNr[k] == 3){ //凹角 左側壁
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;
            double deno = dr/dx*r[j]/rR;

            aPr[i][j] = aPr[i][j] - aNr[i][j]*rL/rR/deno;
            aSEr[i][j] = aSEr[i][j] - aNr[i][j]*dr/dx*r[j]/rR/deno;
            aNr[i][j] = 0.0;
        }if(sBndNr[k] == 4){ //凹角 右側壁
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;
            double deno = dr/dx*r[j]/rR;

            aPr[i][j] = aPr[i][j] - aNr[i][j]*rL/rR/deno;
            aSWr[i][j] = aSWr[i][j] + aNr[i][j]*dr/dx*r[j]/rR/deno;
            aNr[i][j] = 0.0;
        }
    }

    //********* Ep境界条件 (左) *********
    for (int k=0;k<iBndWp.size();k++){
        int i = iBndWp[k];
        int j = jBndWp[k];

        if(sBndWp[k] == 0){ //ディリクレ
            aPp[i][j] = aPp[i][j] + aWp[i][j];
            aWp[i][j] = 0.0;
        }else if(sBndWp[k] == 1){ //開放
            double c0 = sqrt(1.0/(eps0*epsr[i][j]*mu0));
            double deno = pow(2.0*c0,2) + pow(omegam*dx,2);
            std::complex<double> R = std::complex<double>
                (pow(2.0*c0,2) - pow(omegam*dx,2),-4.0*c0*omegam*dx);
            R = R/(deno+1e-100);
            aPp[i][j] = aPp[i][j] - aWp[i][j]*R;
            aWp[i][j] = 0.0;
        }
    }
    //********* Ep境界条件 (右) *********
    for (int k=0;k<iBndEp.size();k++){
        int i = iBndEp[k];
        int j = jBndEp[k];

        if(sBndEp[k] == 0){ //ディリクレ
            aPp[i][j] = aPp[i][j] + aEp[i][j];
            aEp[i][j] = 0.0;
        }else if(sBndEp[k] == 1){ //開放
            double c0 = sqrt(1.0/(eps0*mu0));
            double deno = pow(2.0*c0,2) + pow(omegam*dx,2);
            std::complex<double> R = std::complex<double>
                (pow(2.0*c0,2) - pow(omegam*dx,2),-4.0*c0*omegam*dx);
            R = R/(deno+1e-100);
            
            aPp[i][j] = aPp[i][j] - aEp[i][j]*R;
            aEp[i][j] = 0.0;
        }
    }
    //********* Ep境界条件 (下) *********
    for (int k=0;k<iBndSp.size();k++){
        int i = iBndSp[k];
        int j = jBndSp[k];

        if(sBndSp[k] == 0){ //ディリクレ
            aPp[i][j] = aPp[i][j] + r[j]/r[j-1]*aSp[i][j];
            aSp[i][j] = 0.0;
        }
    }
    //********* Ep境界条件 (上) *********
    for (int k=0;k<iBndNp.size();k++){
        int i = iBndNp[k];
        int j = jBndNp[k];

        if(sBndNp[k] == 0){ //ディリクレ
            aPp[i][j] = aPp[i][j] + r[j]/r[j+1]*aNp[i][j];
            aNp[i][j] = 0.0;
        }else if(sBndNp[k] == 1){ //開放
            double c0 = sqrt(1.0/(eps0*epsr[i][j]*mu0));
            double r_tmp = (r[j+1] + r[j])/2.0;
            double qL = r[j]/(r_tmp+1e-100);
            double qR = r[j+1]  /(r_tmp+1e-100);
            double deno = pow(2.0*c0*sqrt(qR),2) + pow(omegam*dr*qR,2);
            std::complex<double> R = std::complex<double>
                (pow(2.0*c0,2)*sqrt(qL*qR) - pow(omegam*dr,2)*qL*qR
                ,-2.0*c0*omegam*dr*(qL*sqrt(qR)+qR*sqrt(qL)));
            R = R/(deno+1e-100);

            aPp[i][j] = aPp[i][j] - aNp[i][j]*R;
            aNp[i][j] = 0.0;
        }
    }

    /******************** LUソルバーの準備 ********************/
    Eigen::SparseMatrix<std::complex<double> > A(nk, nk);
    A.reserve(Eigen::VectorXi::Constant(nk,11)); //ここの数字を変えて帯域幅を確保する
    Eigen::VectorXcd b(nk);
    Eigen::VectorXcd xv(nk);

    //Ex-足し込み
    for (int k=0;k<nkx;k++){

        int i = ikx[k];
        int j = jkx[k];

        double kE = kx[i+1][j];
        double kW = kx[i-1][j];
        double kN = kx[i][j+1];
        double kS = kx[i][j-1];

        double kNE = kr[i][j+1];
        double kSE = kr[i][j];
        double kNW = kr[i-1][j+1];
        double kSW = kr[i-1][j];

        double kEE = kp[i][j];
        double kWW = kp[i-1][j];
        
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
            A.insert(k, nkx+kNE) = -aNEx[i][j];
        }
        if(kSE!=-1){
            A.insert(k, nkx+kSE) = -aSEx[i][j];
        }
        if(kNW!=-1){
            A.insert(k, nkx+kNW) = -aNWx[i][j];
        }
        if(kSW!=-1){
            A.insert(k, nkx+kSW) = -aSWx[i][j];
        }
        //Ep用
        if(kEE!=-1){
            A.insert(k, nkx+nkr+kEE) = -aEEx[i][j];
        }
        if(kWW!=-1){
            A.insert(k, nkx+nkr+kWW) = -aWWx[i][j];
        }
    }
    
    //Er-足し込み
    for (int k=0;k<nkr;k++){
        int i = ikr[k];
        int j = jkr[k];

        double kE  = kr[i+1][j];
        double kW  = kr[i-1][j];
        double kN  = kr[i][j+1];
        double kS  = kr[i][j-1];

        double kNE = kx[i+1][j];
        double kSE = kx[i+1][j-1];
        double kNW = kx[i][j];
        double kSW = kx[i][j-1];

        double kNN = kp[i][j];
        double kSS = kp[i][j-1];

        //Er用
        A.insert(nkx+k, nkx+k)      =   aPr[i][j];
        b[nkx+k]                    =   br[i][j];
        if(kN!=-1){
            A.insert(nkx+k, nkx+kN) = -aNr[i][j];
        }
        if(kS!=-1){
            A.insert(nkx+k, nkx+kS) = -aSr[i][j];
        }
        if(kE!=-1){
            A.insert(nkx+k, nkx+kE) = -aEr[i][j];
        }
        if(kW!=-1){
            A.insert(nkx+k, nkx+kW) = -aWr[i][j];
        }
        //Ex用
        if(kNE!=-1){
            A.insert(nkx+k, kNE) = -aNEr[i][j];
        }
        if(kSE!=-1){
            A.insert(nkx+k, kSE) = -aSEr[i][j];
        }
        if(kNW!=-1){
            A.insert(nkx+k, kNW) = -aNWr[i][j];
        }
        if(kSW!=-1){
            A.insert(nkx+k, kSW) = -aSWr[i][j];
        }
        //Ep用
        if(kNN!=-1){
            A.insert(nkx+k, nkx+nkr+kNN) = -aNNr[i][j];
        }
        if(kSS!=-1){
            A.insert(nkx+k, nkx+nkr+kSS) = -aSSr[i][j];
        }
    }

    //Ep-足し込み
    for (int k=0;k<nkp;k++){
        int i = ikp[k];
        int j = jkp[k];

        double kE  = kp[i+1][j];
        double kW  = kp[i-1][j];
        double kN  = kp[i][j+1];
        double kS  = kp[i][j-1];

        double kEE  = kx[i+1][j];
        double kWW  = kx[i][j];

        double kNN  = kr[i][j+1];
        double kSS  = kr[i][j];

        //Ep用
        A.insert(nkx+nkr+k, nkx+nkr+k)      =   aPp[i][j];
        b[nkx+nkr+k]                    =   bp[i][j];
        if(kN!=-1){
            A.insert(nkx+nkr+k, nkx+nkr+kN) = -aNp[i][j];
        }
        if(kS!=-1){
            A.insert(nkx+nkr+k, nkx+nkr+kS) = -aSp[i][j];
        }
        if(kE!=-1){
            A.insert(nkx+nkr+k, nkx+nkr+kE) = -aEp[i][j];
        }
        if(kW!=-1){
            A.insert(nkx+nkr+k, nkx+nkr+kW) = -aWp[i][j];
        }
        //Ex用
        if(kEE!=-1){
            A.insert(nkx+nkr+k, kEE) = -aEEp[i][j];
        }
        if(kWW!=-1){
            A.insert(nkx+nkr+k, kWW) = -aWWp[i][j];
        }
        //Er用
        if(kNN!=-1){
            A.insert(nkx+nkr+k, nkx+kNN) = -aNNp[i][j];
        }
        if(kSS!=-1){
            A.insert(nkx+nkr+k, nkx+kSS) = -aSSp[i][j];
        }
    }
 

    // LU分解でAx = bを解く
    Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double> > > solver;
    //std::cout << "solver start!" << std::endl;
    solver.compute(A);
    xv = solver.solve(b);

    /*
    // 実際の誤差を計算
    Eigen::VectorXcd residual = A*xv- b;
    std::cout << residual<< std::endl;
    double actual_error = residual.norm()/(b.norm()+1e-100);
    std::cout << "Actual error: " << actual_error << std::endl;
    */

    //Ex-結果を戻す
    for (int k=0;k<nkx;k++){
        int i = ikx[k];
        int j = jkx[k];
        E1x[i][j] = xv[k];
    }
    //Er-結果を戻す
    for (int k=0;k<nkr;k++){
        int i = ikr[k];
        int j = jkr[k];
        E1r[i][j] = xv[nkx+k];
    }
    //Ep-結果を戻す
    for (int k=0;k<nkp;k++){
        int i = ikp[k];
        int j = jkp[k];
        E1p[i][j] = xv[nkx+nkr+k];
    }

    /************************境界条件後処理**************************/
    //Ex-境界条件
    //********* Ex境界条件 (左) *********
    for (int k=0;k<iBndWx.size();k++){
        int i = iBndWx[k];
        int j = jBndWx[k];

        if(sBndWx[k] == 0){ //ディリクレ
            E1x[i-1][j] = 0.0;
        }else if(sBndWx[k] == 2){ //ガウス
            double qL = (r[j-1] + r[j])/2.0/(r[j]+1e-100);
            double qR = (r[j+1] + r[j])/2.0/(r[j]+1e-100);

            E1x[i-1][j] = E1x[i][j] + dx/dr*(qR*E1r[i-1][j+1] - qL*E1r[i-1][j]);
        }else if(sBndWx[k] == 3){ //凹角 下側壁
            double qL = (r[j-1] + r[j])/2.0/(r[j]+1e-100);
            double qR = (r[j+1] + r[j])/2.0/(r[j]+1e-100);
            double deno = 1.0 + dx/dr*qL;

            E1x[i-1][j] = E1x[i][j]/deno + dx/dr*(qR*E1r[i-1][j+1])/deno;
        }else if(sBndWx[k] == 4){ //凹角 上側壁
            double qL = (r[j-1] + r[j])/2.0/(r[j]+1e-100);
            double qR = (r[j+1] + r[j])/2.0/(r[j]+1e-100);
            double deno = 1.0 + dx/dr*qR;

            E1x[i-1][j] = E1x[i][j]/deno - dx/dr*(qL*E1r[i-1][j])/deno;
        }
    }
    //********* Ex境界条件 (右) *********
    for (int k=0;k<iBndEx.size();k++){
        int i = iBndEx[k];
        int j = jBndEx[k];

        if(sBndEx[k] == 2){ //ガウス
            double qL = (r[j-1] + r[j])/2.0/(r[j]+1e-100);
            double qR = (r[j+1] + r[j])/2.0/(r[j]+1e-100);
            
            E1x[i+1][j] = E1x[i][j] - dx/dr*(qR*E1r[i][j+1] - qL*E1r[i][j]);
        }else if(sBndEx[k] == 3){ //凹角 下側壁
            double qL = (r[j-1] + r[j])/2.0/(r[j]+1e-100);
            double qR = (r[j+1] + r[j])/2.0/(r[j]+1e-100);
            double deno = 1.0 + dx/dr*qL;

            E1x[i+1][j] = E1x[i][j]/deno - dx/dr*qR*E1r[i][j+1]/deno;
        }else if(sBndEx[k] == 4){ //凹角 上側壁
            double qL = (r[j-1] + r[j])/2.0/(r[j]+1e-100);
            double qR = (r[j+1] + r[j])/2.0/(r[j]+1e-100);
            double deno = 1.0 + dx/dr*qR;

            E1x[i+1][j] = E1x[i][j]/deno + dx/dr*qL*E1r[i][j]/deno;
        }
    }
    //********* Ex境界条件 (下) *********
    for (int k=0;k<iBndSx.size();k++){
        int i = iBndSx[k];
        int j = jBndSx[k];
        if(sBndSx[k] == 0){ //ディリクレ
            E1x[i][j-1] = -r[j]/r[j-1]*E1x[i][j];
        }else if(sBndSx[k] == 6){ //凸角 右側Open
            //処理なし
        }
    }
    //********* Ex境界条件 (上) *********
    for (int k=0;k<iBndNx.size();k++){
        int i = iBndNx[k];
        int j = jBndNx[k];
        if(sBndNx[k] == 0){ //ディリクレ
            E1x[i][j+1] = -r[j]/r[j+1]*E1x[i][j];
        }else if(sBndNx[k] == 1){ //開放
            double c0 = sqrt(1.0/(eps0*epsr[i][j]*mu0));
            double r_tmp = (r[j+1] + r[j])/2.0;
            double qL = r[j]/(r_tmp+1e-100);
            double qR = r[j+1]  /(r_tmp+1e-100);
            double deno = pow(2.0*c0*sqrt(qR),2) + pow(omegam*dr*qR,2);
            std::complex<double> R = std::complex<double>
                (pow(2.0*c0,2)*sqrt(qL*qR) - pow(omegam*dr,2)*qL*qR
                ,-2.0*c0*omegam*dr*(qL*sqrt(qR)+qR*sqrt(qL)));
            R = R/(deno+1e-100);
            E1x[i][j+1] = R*E1x[i][j];
        }else if(sBndNx[k] == 5){ //凸角 左側Open
            //処理なし
        }else if(sBndNx[k] == 6){ //凸角 右側Open
            //処理なし
        }
    }

    //********* Er境界条件 (左) *********
    for (int k=0;k<iBndWr.size();k++){
        int i = iBndWr[k];
        int j = jBndWr[k];

        if(sBndWr[k] == 0){ //ディリクレ
            E1r[i-1][j] = -E1r[i][j];
        }else if(sBndWr[k] == 1){ //開放
            double c0 = sqrt(1.0/(eps0*epsr[i][j]*mu0));
            double deno = pow(2.0*c0,2) + pow(omegam*dx,2);
            std::complex<double> R = std::complex<double>
                (pow(2.0*c0,2) - pow(omegam*dx,2),-4.0*c0*omegam*dx);
            R = R/(deno+1e-100);

            E1r[i-1][j] = R*E1r[i][j];
        }else if(sBndWr[k] == 5){ //凸角 下側Open
            //処理なし
        }else if(sBndWr[k] == 6){ //凸角 上側Open
            //処理なし
        }
    }
    //********* Er境界条件 (右) *********
    for (int k=0;k<iBndEr.size();k++){
        int i = iBndEr[k];
        int j = jBndEr[k];

        if(sBndEr[k] == 0){ //ディリクレ
            E1r[i+1][j] = -E1r[i][j];
        }else if(sBndEr[k] == 1){ //開放
            double c0 = sqrt(1.0/(eps0*mu0));
            //std::complex<double> deno(2.0*c0,omegam*dx);
            //std::complex<double> R(2.0*c0,-omegam*dx);

            double deno = pow(2.0*c0,2) + pow(omegam*dx,2);
            std::complex<double> R = std::complex<double>
                (pow(2.0*c0,2) - pow(omegam*dx,2),-4.0*c0*omegam*dx);
            R = R/(deno+1e-100);

            E1r[i+1][j] = R/(deno+1e-100)*E1r[i][j];
        }else if(sBndEr[k] == 5){ //凸角 下側Open
            //処理なし
        }
    }
    //********* Er境界条件 (下) *********
    for (int k=0;k<iBndSr.size();k++){
        int i = iBndSr[k];
        int j = jBndSr[k];

        if(sBndSr[k] == 2){ //ガウス
            double rR = (r[j] + r[j-1])/2.0;
            double rL = (r[j-1] + r[j-2])/2.0;
            
            E1r[i][j-1] = rR/rL*E1r[i][j] + dr/dx*r[j-1]/rL*(E1x[i+1][j-1]-E1x[i][j-1]);
        }else if(sBndSr[k] == 3){ //凹角 左側壁
            double rR = (r[j] + r[j-1])/2.0;
            double rL = (r[j-1] + r[j-2])/2.0;
            double deno = dr/dx*r[j-1]/rL;
            if(r[j] < j_flc_bl[0][0]){ //マイクロ波の入り口部の特別処理
                deno = 1.0;
            }

            E1r[i][j-1] = rR/rL*E1r[i][j]/deno + dr/dx*r[j-1]/rL*(E1x[i+1][j-1])/deno;
        }
    }
    //********* Er境界条件 (上) *********
    for (int k=0;k<iBndNr.size();k++){
        int i = iBndNr[k];
        int j = jBndNr[k];

        if(sBndNr[k] == 2){ //ガウス
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;

            E1r[i][j+1] = rL/rR*E1r[i][j] - dr/dx*r[j]/rR*(E1x[i+1][j]-E1x[i][j]);
        }if(sBndNr[k] == 3){ //凹角 左側壁
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;
            double deno = dr/dx*r[j]/rR;

            E1r[i][j+1] = rL/rR*E1r[i][j]/deno - dr/dx*r[j]/rR*(E1x[i+1][j])/deno;
        }if(sBndNr[k] == 4){ //凹角 右側壁
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;
            double deno = dr/dx*r[j]/rR;

            E1r[i][j+1] = rL/rR*E1r[i][j]/deno + dr/dx*r[j]/rR*(E1x[i][j])/deno;
        }
    }

    //********* Ep境界条件 (左) *********
    for (int k=0;k<iBndWp.size();k++){
        int i = iBndWp[k];
        int j = jBndWp[k];

        if(sBndWp[k] == 0){ //ディリクレ
            E1p[i-1][j] = -E1p[i][j];
        }else if(sBndWp[k] == 1){ //開放
            double c0 = sqrt(1.0/(eps0*epsr[i][j]*mu0));
            double deno = pow(2.0*c0,2) + pow(omegam*dx,2);
            std::complex<double> R = std::complex<double>
                (pow(2.0*c0,2) - pow(omegam*dx,2),-4.0*c0*omegam*dx);
            R = R/(deno+1e-100);

            E1p[i-1][j] = R*E1p[i][j];
        }
    }
    //********* Ep境界条件 (右) *********
    for (int k=0;k<iBndEp.size();k++){
        int i = iBndEp[k];
        int j = jBndEp[k];

        if(sBndEp[k] == 0){ //ディリクレ
            E1p[i+1][j] = -E1p[i][j];
        }else if(sBndEp[k] == 1){ //開放
            double c0 = sqrt(1.0/(eps0*mu0));
            //double deno = pow(2.0*c0,2) + pow(omegam*dx,2);
            std::complex<double> deno(2.0*c0,omegam*dx);
            std::complex<double> R(2.0*c0,-omegam*dx);
            R = R/(deno+1e-100);

            E1p[i+1][j] = R*E1p[i][j];
        }
    }
    //********* Ep境界条件 (下) *********
    for (int k=0;k<iBndSp.size();k++){
        int i = iBndSp[k];
        int j = jBndSp[k];

        E1p[i][j-1] = -r[j]/r[j-1]*E1p[i][j];
    }
    //********* Ep境界条件 (上) *********
    for (int k=0;k<iBndNp.size();k++){
        int i = iBndNp[k];
        int j = jBndNp[k];

        if(sBndNp[k] == 0){ //ディリクレ
            E1p[i][j+1] = -r[j]/r[j+1]*E1p[i][j];
        }else if(sBndNp[k] == 1){ //開放
            double c0 = sqrt(1.0/(eps0*epsr[i][j]*mu0));
            double r_tmp = (r[j+1] + r[j])/2.0;
            double qL = r[j]/(r_tmp+1e-100);
            double qR = r[j+1]  /(r_tmp+1e-100);
            double deno = pow(2.0*c0*sqrt(qR),2) + pow(omegam*dr*qR,2);
            std::complex<double> R = std::complex<double>
                (pow(2.0*c0,2)*sqrt(qL*qR) - pow(omegam*dr,2)*qL*qR
                ,-2.0*c0*omegam*dr*(qL*sqrt(qR)+qR*sqrt(qL)));
            R = R/(deno+1e-100);
            
            E1p[i][j+1] = R*E1p[i][j];
        }
    }


    {

        //set refelence plane
        //---------------------------------
        double x_ref = 0.004;
        int i_ref = int(x_ref/dx + 1.5);
        //---------------------------------

        int i = i_ref;
        double x_tmp = (x[i] + x[i-1])/2.0;
        //calculate current I (use averaged value in r-direction for robustness)
        //---------------------------------
        int ncount = 0;
        std::complex<double>I_ref_tmp(0.0,0.0);
        for (int j=j_flr_bl[5][0];j<=j_flr_bl[5][1];j++){
            double r_tmp = (r[j] + r[j-1])/2.0;
            //double qR = (r[j] + r[j+1])/(2.0*r[j]);

            std::complex<double> H_ref = (E1r[i][j] - E1r[i-1][j])/dx/(omegam*mu0)*std::complex<double>(0,1.0);
            I_ref_tmp = I_ref_tmp + 2.0*M_PI*r_tmp*H_ref;

            //std::cout << "j = " << j << " - 1/2, r = "<< r_tmp << " , H_ref = "<<H_ref<<", I_ref = "<<2.0*M_PI*r_tmp*H_ref<< std::endl;
            ncount = ncount + 1;
        }
        I_ref_tmp = I_ref_tmp/double(ncount);
        //---------------------------------

        //calculate current V
        //---------------------------------
        std::complex<double> V_ref_tmp(0.0,0.0);
        for (int j=j_flc_bl[5][0];j<=j_flc_bl[5][1];j++){
            //double qL = (r[j] + r[j-1])/(2.0*r[j]);
            //double qR = (r[j] + r[j+1])/(2.0*r[j]);
            //double r_tmp = (r[j] + r[j-1])/2.0;
            double rL = (r[j] + r[j-1])/2.0;
            double rR = (r[j] + r[j+1])/2.0;

            //std::complex<double> Er_tmp = (E1r[i][j] + E1r[i-1][j])/2.0;
            std::complex<double> Er_tmp = (rR*(E1r[i][j+1] + E1r[i-1][j+1]) + rL*(E1r[i][j] + E1r[i-1][j]))/(4.0*r[j]);
            
            double c0 = sqrt(1.0/(eps0*epsr[i][j]*mu0));
            
            //std::cout << "Check ,i = "<<i <<", j = "<< j << ", "<< std::abs(Er_tmp)*std::abs(Er_tmp)
            //    //<< ","<<std::abs(Er_Fw[j+1])<< ","<<std::abs(Er_Fw[j])
            //    << ","<<(Er_tmp)
            //    << ","<<pow(std::abs((Er_tmp - Er_anly)/Er_anly),2)*100<< " %"
            //    << std::endl;

            V_ref_tmp = V_ref_tmp + Er_tmp*dr;
            //std::cout << "j = " << j <<", Er = "<< Er_tmp << ", rEr = "<< r_tmp*Er_tmp << ", V_ref = "<<V_ref << ", V_ref_th = "<<r_tmp*Er_tmp*log(r2/r1)<< std::endl;
        }
        //---------------------------------

        std::complex<double> Z_ref_tmp = V_ref_tmp/(I_ref_tmp + 1e-100);
        
        std::cout << "Impedance at reference plane x = " << x_tmp << std::endl;
        std::cout << "Z_ref = " << Z_ref_tmp << ", V_ref = " << V_ref_tmp << ", I_ref = " << I_ref_tmp 
            << ", |Z_ref| = " << std::abs(Z_ref_tmp) 
            << ", arg(Z_ref) = " << std::arg(Z_ref_tmp)*180/M_PI << " deg"<< std::endl;

        std::complex<double> a2 = (V_ref_tmp - Z0_base*I_ref_tmp)/(2.0*sqrt(Z0_base));
        std::complex<double> b2 = (V_ref_tmp + Z0_base*I_ref_tmp)/(2.0*sqrt(Z0_base));

        std::complex<double> S11(S11_mag * std::cos(S11_arg_deg/180.0*M_PI), S11_mag * std::sin(S11_arg_deg/180.0*M_PI));
        std::complex<double> S21(S21_mag * std::cos(S21_arg_deg/180.0*M_PI), S21_mag * std::sin(S21_arg_deg/180.0*M_PI));
        std::complex<double> S12(S12_mag * std::cos(S12_arg_deg/180.0*M_PI), S12_mag * std::sin(S12_arg_deg/180.0*M_PI));
        std::complex<double> S22(S22_mag * std::cos(S22_arg_deg/180.0*M_PI), S22_mag * std::sin(S22_arg_deg/180.0*M_PI));

        //std::complex<double> a1 = (b2 - S22*a2)/S21;
        //std::complex<double> b1 = S11*a1 + S12*a2;
        double P_fwd = std::norm(b2);
        double P_ref = std::norm(a2);

        std::cout << std::endl;
        std::cout << "Fwd power = " << P_fwd << " W" << std::endl;
        std::cout << "Ref power = " << P_ref << " W, " <<P_ref/P_fwd*100<<" %"<< std::endl;
        
    }

}


//*****************************************************************
//**                                                             **
//**           void solve_Efield                                 **
//**                                                             **
//*****************************************************************
void solve_Microwave_impedanceTest(){

    //テスト用
    for (int i=0;i<ni+2;i++){
        for (int j=0;j<nj+2;j++){
            //rhoe[i][j] = 0.0;
            nu_m1[i][j] = nu_eff;
        }
    }

    //Ex用 係数
    std::vector<std::vector<std::complex<double> > > aPx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aWx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aEx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aNx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aSx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aNEx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aNWx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aSEx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aSWx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aEEx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aWWx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > bx(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数

    //Er用 係数
    std::vector<std::vector<std::complex<double> > > aPr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aWr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aEr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aNr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aSr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aNEr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aNWr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aSEr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aSWr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aNNr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aSSr(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > br(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数

    //Ephi用 係数
    std::vector<std::vector<std::complex<double> > > aPp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aWp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aEp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aNp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aSp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aNNp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aSSp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aEEp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > aWWp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数
    std::vector<std::vector<std::complex<double> > > bp(ni+2,std::vector<std::complex<double> >(nj+2,0.0)); //係数

    /******************** 係数作成 いらないところも含めて普通に全部やる方が速い ********************/
    std::complex<double> EPS(1e-100,1e-100); //微小複素数
    
    //for Er
    for (int iblock=5;iblock<=5;iblock++){ 
        for (int i=i_flr_bl[iblock][0];i<=i_flr_bl[iblock][1];i++){ 
            for (int j=j_flr_bl[iblock][0];j<=j_flr_bl[iblock][1];j++){
                double r_tmp = (r[j-1] + r[j])/2.0;
                double qL = r[j-1]/(r_tmp+1e-100);
                double qR = r[j]  /(r_tmp+1e-100);

                double nu_m_tmp = (nu_m1[i][j] + nu_m1[i][j-1])/2.0;
                double rho_tmp  = (rhoe[i][j]  + rhoe[i][j-1] )/2.0;
                double Bx_tmp   = (r[j]*Bx[i][j]   + r[j-1]*Bx[i][j-1]  )/(2.0*r_tmp);
                double Br_tmp   = (r[j]*Br[i][j]   + r[j-1]*Br[i][j-1]  )/(2.0*r_tmp);

                std::complex<double> nu_cmp(nu_m_tmp,omegam);
                std::complex<double> Hallcmpx = e0*Bx_tmp/(masse*nu_cmp + EPS);
                std::complex<double> Hallcmpr = e0*Br_tmp/(masse*nu_cmp + EPS);
                std::complex<double> HallcmpMag2 =  Hallcmpx*Hallcmpx + Hallcmpr*Hallcmpr;
                std::complex<double> coef = rho_tmp*e0*e0/(masse*nu_cmp + EPS)/(1.0 + HallcmpMag2 + EPS);
                
                std::complex<double> iomegaMu0(0.0,omegam*mu0);

                aPr[i][j] = 2.0*dx*dx + 2.0*dr*dr + dx*dx*dr*dr/(r_tmp*r_tmp)
                    - eps0*epsr[i][j]*mu0*omegam*omegam*dx*dx*dr*dr;
                aEr[i][j] = dr*dr;
                aWr[i][j] = dr*dr;
                aNr[i][j] = qR*dx*dx;
                aSr[i][j] = qL*dx*dx;
                br[i][j] = std::complex<double>(0,-omegam*mu0*dx*dx*dr*dr)*J1r[i][j];
            }
        }
    }

    /******************** 境界条件設定 (係数修正) ********************/
    

    //********* Er境界条件 (左) *********
    for (int k=0;k<iBndWr.size();k++){
        int i = iBndWr[k];
        int j = jBndWr[k];

        if(sBndWr[k] == 0){ //ディリクレ
            aPr[i][j] = aPr[i][j] + aWr[i][j];
            aWr[i][j] = 0.0;
        }else if(sBndWr[k] == 1){ //開放
            double c0 = sqrt(1.0/(eps0*epsr[i][j]*mu0));
            double deno = pow(2.0*c0,2) + pow(omegam*dx,2);
            std::complex<double> R = std::complex<double>(pow(2.0*c0,2) - pow(omegam*dx,2),-4.0*c0*omegam*dx);
            R = R/(deno+1e-100);
            
            aPr[i][j] = aPr[i][j] - aWr[i][j]*R;
            aWr[i][j] = 0.0;
        }else if(sBndWr[k] == 5){ //凸角 下側Open
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;
            
            aPr[i][j] = aPr[i][j] - aWr[i][j]*( dx/dr*rL/r[j]);
            aNr[i][j] = aNr[i][j] + aWr[i][j]*(-dx/dr*rR/r[j]);
            aNEr[i][j] = aNEr[i][j] + aWr[i][j]*(-1.0);
            aWr[i][j] = 0.0;
        }else if(sBndWr[k] == 6){ //凸角 上側Open
            double rR = (r[j] + r[j-1])/2.0;
            double rL = (r[j-1] + r[j-2])/2.0;

            aPr[i][j] = aPr[i][j] - aWr[i][j]*( dx/dr*rR/r[j-1]);
            aSr[i][j] = aSr[i][j] + aWr[i][j]*(-dx/dr*rL/r[j-1]);
            aSEr[i][j] = aSEr[i][j] + aWr[i][j]*( 1.0);
            aWr[i][j] = 0.0;
        }
    }
    //********* Er境界条件 (右) *********
    for (int k=0;k<iBndEr.size();k++){
        int i = iBndEr[k];
        int j = jBndEr[k];

        if(sBndEr[k] == 0){ //ディリクレ
            aPr[i][j] = aPr[i][j] + aEr[i][j];
            aEr[i][j] = 0.0;
        }if(sBndEr[k] == 1){ //開放
            double c0 = sqrt(1.0/(eps0*mu0));
            
            //std::complex<double> deno(2.0*c0,omegam*dx);
            //std::complex<double> R(2.0*c0,-omegam*dx);
            //aPr[i][j] = aPr[i][j] - aEr[i][j]*R/(deno+1e-100);
            //aEr[i][j] = 0.0;

            double deno = pow(2.0*c0,2) + pow(omegam*dx,2);
            std::complex<double> R = std::complex<double>
                (pow(2.0*c0,2) - pow(omegam*dx,2),-4.0*c0*omegam*dx);
            aPr[i][j] = aPr[i][j] - aEr[i][j]*R/(deno+1e-100);
            aEr[i][j] = 0.0;


        }else if(sBndEr[k] == 5){ //凸角 下側Open (*)
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;

            aPr[i][j] = aPr[i][j] - aEr[i][j]*( dx/dr*rL/r[j]);
            aNr[i][j] = aNr[i][j] + aEr[i][j]*(-dx/dr*rR/r[j]);
            aNWr[i][j] = aNWr[i][j] + aEr[i][j]*( 1.0);
            aEr[i][j] = 0.0;
        }
    }
    //********* Er境界条件 (下) *********
    for (int k=0;k<iBndSr.size();k++){
        int i = iBndSr[k];
        int j = jBndSr[k];

        if(sBndSr[k] == 2){ //ガウス
            double rR = (r[j] + r[j-1])/2.0;
            double rL = (r[j-1] + r[j-2])/2.0;
            aPr[i][j] = aPr[i][j] - aSr[i][j]*rR/rL;
            aSEr[i][j] = aSEr[i][j] + aSr[i][j]*dr/dx*r[j-1]/rL;
            aSWr[i][j] = aSWr[i][j] - aSr[i][j]*dr/dx*r[j-1]/rL;
            aSr[i][j] = 0.0;
        }else if(sBndSr[k] == 3){ //凹角 左側壁
            double rR = (r[j] + r[j-1])/2.0;
            double rL = (r[j-1] + r[j-2])/2.0;
            double deno = dr/dx*r[j-1]/rL;
            if(r[j] < j_flc_bl[0][0]){ //マイクロ波の入り口部の特別処理
                deno = 1.0;
            }

            aPr[i][j] = aPr[i][j] - aSr[i][j]*rR/rL/deno;
            aNEr[i][j] = aNEr[i][j] + aSr[i][j]*dr/dx*r[j-1]/rL/deno;
            aSr[i][j] = 0.0;
        }
    }
    //********* Er境界条件 (上) *********
    for (int k=0;k<iBndNr.size();k++){
        int i = iBndNr[k];
        int j = jBndNr[k];

        if(sBndNr[k] == 2){ //ガウス
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;
            aPr[i][j] = aPr[i][j] - aNr[i][j]*rL/rR;
            aNEr[i][j] = aNEr[i][j] - aNr[i][j]*dr/dx*r[j]/rR;
            aNWr[i][j] = aNWr[i][j] + aNr[i][j]*dr/dx*r[j]/rR;
            aNr[i][j] = 0.0;
        }if(sBndNr[k] == 3){ //凹角 左側壁
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;
            double deno = dr/dx*r[j]/rR;

            aPr[i][j] = aPr[i][j] - aNr[i][j]*rL/rR/deno;
            aSEr[i][j] = aSEr[i][j] - aNr[i][j]*dr/dx*r[j]/rR/deno;
            aNr[i][j] = 0.0;
        }if(sBndNr[k] == 4){ //凹角 右側壁
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;
            double deno = dr/dx*r[j]/rR;

            aPr[i][j] = aPr[i][j] - aNr[i][j]*rL/rR/deno;
            aSWr[i][j] = aSWr[i][j] + aNr[i][j]*dr/dx*r[j]/rR/deno;
            aNr[i][j] = 0.0;
        }
    }

    //********* Ep境界条件 (左) *********
    for (int k=0;k<iBndWp.size();k++){
        int i = iBndWp[k];
        int j = jBndWp[k];

        if(sBndWp[k] == 0){ //ディリクレ
            aPp[i][j] = aPp[i][j] + aWp[i][j];
            aWp[i][j] = 0.0;
        }else if(sBndWp[k] == 1){ //開放
            double c0 = sqrt(1.0/(eps0*epsr[i][j]*mu0));
            double deno = pow(2.0*c0,2) + pow(omegam*dx,2);
            std::complex<double> R = std::complex<double>
                (pow(2.0*c0,2) - pow(omegam*dx,2),-4.0*c0*omegam*dx);
            R = R/(deno+1e-100);
            aPp[i][j] = aPp[i][j] - aWp[i][j]*R;
            aWp[i][j] = 0.0;
        }
    }

    /******************** LUソルバーの準備 ********************/
    Eigen::SparseMatrix<std::complex<double> > A(nk, nk);
    A.reserve(Eigen::VectorXi::Constant(nk,11)); //ここの数字を変えて帯域幅を確保する
    Eigen::VectorXcd b(nk);
    Eigen::VectorXcd xv(nk);
    
    //Er-足し込み
    for (int k=0;k<nkr;k++){
        int i = ikr[k];
        int j = jkr[k];

        double kE  = kr[i+1][j];
        double kW  = kr[i-1][j];
        double kN  = kr[i][j+1];
        double kS  = kr[i][j-1];

        double kNE = kx[i+1][j];
        double kSE = kx[i+1][j-1];
        double kNW = kx[i][j];
        double kSW = kx[i][j-1];

        double kNN = kp[i][j];
        double kSS = kp[i][j-1];

        //Er用
        A.insert(nkx+k, nkx+k)      =   aPr[i][j];
        b[nkx+k]                    =   br[i][j];
        if(kN!=-1){
            A.insert(nkx+k, nkx+kN) = -aNr[i][j];
        }
        if(kS!=-1){
            A.insert(nkx+k, nkx+kS) = -aSr[i][j];
        }
        if(kE!=-1){
            A.insert(nkx+k, nkx+kE) = -aEr[i][j];
        }
        if(kW!=-1){
            A.insert(nkx+k, nkx+kW) = -aWr[i][j];
        }
        //Ex用
        if(kNE!=-1){
            A.insert(nkx+k, kNE) = -aNEr[i][j];
        }
        if(kSE!=-1){
            A.insert(nkx+k, kSE) = -aSEr[i][j];
        }
        if(kNW!=-1){
            A.insert(nkx+k, kNW) = -aNWr[i][j];
        }
        if(kSW!=-1){
            A.insert(nkx+k, kSW) = -aSWr[i][j];
        }
        //Ep用
        if(kNN!=-1){
            A.insert(nkx+k, nkx+nkr+kNN) = -aNNr[i][j];
        }
        if(kSS!=-1){
            A.insert(nkx+k, nkx+nkr+kSS) = -aSSr[i][j];
        }
    }


    // LU分解でAx = bを解く
    Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double> > > solver;
    //std::cout << "solver start!" << std::endl;
    solver.compute(A);
    xv = solver.solve(b);

    /*
    // 実際の誤差を計算
    Eigen::VectorXcd residual = A*xv- b;
    std::cout << residual<< std::endl;
    double actual_error = residual.norm()/(b.norm()+1e-100);
    std::cout << "Actual error: " << actual_error << std::endl;
    */

    //Er-結果を戻す
    for (int k=0;k<nkr;k++){
        int i = ikr[k];
        int j = jkr[k];
        E1r[i][j] = xv[nkx+k];
    }

    /************************境界条件後処理**************************/
   
    //********* Er境界条件 (左) *********
    for (int k=0;k<iBndWr.size();k++){
        int i = iBndWr[k];
        int j = jBndWr[k];

        if(sBndWr[k] == 0){ //ディリクレ
            E1r[i-1][j] = -E1r[i][j];
        }else if(sBndWr[k] == 1){ //開放
            double c0 = sqrt(1.0/(eps0*epsr[i][j]*mu0));
            double deno = pow(2.0*c0,2) + pow(omegam*dx,2);
            std::complex<double> R = std::complex<double>
                (pow(2.0*c0,2) - pow(omegam*dx,2),-4.0*c0*omegam*dx);
            R = R/(deno+1e-100);

            E1r[i-1][j] = R*E1r[i][j];
        }else if(sBndWr[k] == 5){ //凸角 下側Open
            //処理なし
        }else if(sBndWr[k] == 6){ //凸角 上側Open
            //処理なし
        }
    }
    //********* Er境界条件 (右) *********
    for (int k=0;k<iBndEr.size();k++){
        int i = iBndEr[k];
        int j = jBndEr[k];

        if(sBndEr[k] == 0){ //ディリクレ
            E1r[i+1][j] = -E1r[i][j];
        }else if(sBndEr[k] == 1){ //開放
            double c0 = sqrt(1.0/(eps0*mu0));
            //std::complex<double> deno(2.0*c0,omegam*dx);
            //std::complex<double> R(2.0*c0,-omegam*dx);

            double deno = pow(2.0*c0,2) + pow(omegam*dx,2);
            std::complex<double> R = std::complex<double>
                (pow(2.0*c0,2) - pow(omegam*dx,2),-4.0*c0*omegam*dx);
            R = R/(deno+1e-100);

            E1r[i+1][j] = R/(deno+1e-100)*E1r[i][j];
        }else if(sBndEr[k] == 5){ //凸角 下側Open
            //処理なし
        }
    }
    //********* Er境界条件 (下) *********
    for (int k=0;k<iBndSr.size();k++){
        int i = iBndSr[k];
        int j = jBndSr[k];

        if(sBndSr[k] == 2){ //ガウス
            double rR = (r[j] + r[j-1])/2.0;
            double rL = (r[j-1] + r[j-2])/2.0;
            
            E1r[i][j-1] = rR/rL*E1r[i][j] + dr/dx*r[j-1]/rL*(E1x[i+1][j-1]-E1x[i][j-1]);
        }else if(sBndSr[k] == 3){ //凹角 左側壁
            double rR = (r[j] + r[j-1])/2.0;
            double rL = (r[j-1] + r[j-2])/2.0;
            double deno = dr/dx*r[j-1]/rL;
            if(r[j] < j_flc_bl[0][0]){ //マイクロ波の入り口部の特別処理
                deno = 1.0;
            }

            E1r[i][j-1] = rR/rL*E1r[i][j]/deno + dr/dx*r[j-1]/rL*(E1x[i+1][j-1])/deno;
        }
    }
    //********* Er境界条件 (上) *********
    for (int k=0;k<iBndNr.size();k++){
        int i = iBndNr[k];
        int j = jBndNr[k];

        if(sBndNr[k] == 2){ //ガウス
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;

            E1r[i][j+1] = rL/rR*E1r[i][j] - dr/dx*r[j]/rR*(E1x[i+1][j]-E1x[i][j]);
        }if(sBndNr[k] == 3){ //凹角 左側壁
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;
            double deno = dr/dx*r[j]/rR;

            E1r[i][j+1] = rL/rR*E1r[i][j]/deno - dr/dx*r[j]/rR*(E1x[i+1][j])/deno;
        }if(sBndNr[k] == 4){ //凹角 右側壁
            double rR = (r[j+1] + r[j])/2.0;
            double rL = (r[j] + r[j-1])/2.0;
            double deno = dr/dx*r[j]/rR;

            E1r[i][j+1] = rL/rR*E1r[i][j]/deno + dr/dx*r[j]/rR*(E1x[i][j])/deno;
        }
    }


    {

        //set refelence plane
        //---------------------------------
        double x_ref = 0.004;
        int i_ref = int(x_ref/dx + 1.5);
        //---------------------------------

        int i = i_ref;
        double x_tmp = (x[i] + x[i-1])/2.0;
        //calculate current I (use averaged value in r-direction for robustness)
        //---------------------------------
        int ncount = 0;
        std::complex<double>I_ref_tmp(0.0,0.0);
        for (int j=j_flr_bl[5][0];j<=j_flr_bl[5][1];j++){
            double r_tmp = (r[j] + r[j-1])/2.0;
            //double qR = (r[j] + r[j+1])/(2.0*r[j]);

            std::complex<double> H_ref = (E1r[i][j] - E1r[i-1][j])/dx/(omegam*mu0)*std::complex<double>(0,1.0);
            I_ref_tmp = I_ref_tmp + 2.0*M_PI*r_tmp*H_ref;

            //std::cout << "j = " << j << " - 1/2, r = "<< r_tmp << " , H_ref = "<<H_ref<<", I_ref = "<<2.0*M_PI*r_tmp*H_ref<< std::endl;
            ncount = ncount + 1;
        }
        I_ref_tmp = I_ref_tmp/double(ncount);
        //---------------------------------

        //calculate current V
        //---------------------------------
        std::complex<double> V_ref_tmp(0.0,0.0);
        for (int j=j_flc_bl[5][0];j<=j_flc_bl[5][1];j++){
            //double qL = (r[j] + r[j-1])/(2.0*r[j]);
            //double qR = (r[j] + r[j+1])/(2.0*r[j]);
            //double r_tmp = (r[j] + r[j-1])/2.0;
            double rL = (r[j] + r[j-1])/2.0;
            double rR = (r[j] + r[j+1])/2.0;

            //std::complex<double> Er_tmp = (E1r[i][j] + E1r[i-1][j])/2.0;
            std::complex<double> Er_tmp = (rR*(E1r[i][j+1] + E1r[i-1][j+1]) + rL*(E1r[i][j] + E1r[i-1][j]))/(4.0*r[j]);
            
            double c0 = sqrt(1.0/(eps0*epsr[i][j]*mu0));
            
            //std::cout << "Check ,i = "<<i <<", j = "<< j << ", "<< std::abs(Er_tmp)*std::abs(Er_tmp)
            //    //<< ","<<std::abs(Er_Fw[j+1])<< ","<<std::abs(Er_Fw[j])
            //    << ","<<(Er_tmp)
            //    << ","<<pow(std::abs((Er_tmp - Er_anly)/Er_anly),2)*100<< " %"
            //    << std::endl;

            V_ref_tmp = V_ref_tmp + Er_tmp*dr;
            //std::cout << "j = " << j <<", Er = "<< Er_tmp << ", rEr = "<< r_tmp*Er_tmp << ", V_ref = "<<V_ref << ", V_ref_th = "<<r_tmp*Er_tmp*log(r2/r1)<< std::endl;
        }
        //---------------------------------

        std::complex<double> Z_ref_tmp = V_ref_tmp/(I_ref_tmp + 1e-100);
        
        std::cout << "Impedance at reference plane x = " << x_tmp << std::endl;
        std::cout << "Z_ref = " << Z_ref_tmp << ", V_ref = " << V_ref_tmp << ", I_ref = " << I_ref_tmp 
            << ", |Z_ref| = " << std::abs(Z_ref_tmp) 
            << ", arg(Z_ref) = " << std::arg(Z_ref_tmp)*180/M_PI << " deg"<< std::endl;

        std::complex<double> a2 = (V_ref_tmp - Z0_base*I_ref_tmp)/(2.0*sqrt(Z0_base));
        std::complex<double> b2 = (V_ref_tmp + Z0_base*I_ref_tmp)/(2.0*sqrt(Z0_base));

        std::complex<double> S11(S11_mag * std::cos(S11_arg_deg/180.0*M_PI), S11_mag * std::sin(S11_arg_deg/180.0*M_PI));
        std::complex<double> S21(S21_mag * std::cos(S21_arg_deg/180.0*M_PI), S21_mag * std::sin(S21_arg_deg/180.0*M_PI));
        std::complex<double> S12(S12_mag * std::cos(S12_arg_deg/180.0*M_PI), S12_mag * std::sin(S12_arg_deg/180.0*M_PI));
        std::complex<double> S22(S22_mag * std::cos(S22_arg_deg/180.0*M_PI), S22_mag * std::sin(S22_arg_deg/180.0*M_PI));

        //std::complex<double> a1 = (b2 - S22*a2)/S21;
        //std::complex<double> b1 = S11*a1 + S12*a2;
        double P_fwd = std::norm(b2);
        double P_ref = std::norm(a2);

        std::cout << std::endl;
        std::cout << "Fwd power = " << P_fwd << " W" << std::endl;
        std::cout << "Ref power = " << P_ref << " W, " <<P_ref/P_fwd*100<<" %"<< std::endl;
        
    }

}

//*****************************************************************
//**                                                             **
//**           void update_energy_profile()                      **
//**                                                             **
//*****************************************************************
void update_energy_profile(){

    double power_sum = 0.0;
    for (int iblock=0;iblock<5;iblock++){ 
        for (int i=i_flc_bl[iblock][0];i<=i_flc_bl[iblock][1];i++){ 
            for (int j=j_flc_bl[iblock][0];j<=j_flc_bl[iblock][1];j++){

                double rL = (r[j]+r[j-1])/2.0;
                double rR = (r[j]+r[j+1])/2.0;

                std::complex<double> nu_cmp(nu_m1[i][j],omegam);
                std::complex<double> Hallcmpx = e0*Bx[i][j]/(masse*nu_cmp);
                std::complex<double> Hallcmpr = e0*Br[i][j]/(masse*nu_cmp);
                std::complex<double> HallcmpMag2 =  Hallcmpx*Hallcmpx + Hallcmpr*Hallcmpr;
                std::complex<double> coef = rhoe[i][j]*e0*e0/(masse*nu_cmp+1e-100)/(1.0 + HallcmpMag2);
                std::complex<double> sigmaexr = coef*(Hallcmpx*Hallcmpr);
                std::complex<double> sigmaexp = coef*(Hallcmpr);
                std::complex<double> sigmaexx = coef*(1.0 + Hallcmpx*Hallcmpx);
                std::complex<double> sigmaerr = coef*(1.0 + Hallcmpr*Hallcmpr);
                std::complex<double> sigmaerp = coef*(-Hallcmpx);
                std::complex<double> sigmaerx = coef*(Hallcmpx*Hallcmpr);
                std::complex<double> sigmaepr = coef*(Hallcmpx);
                std::complex<double> sigmaepp = coef*(1.0);
                std::complex<double> sigmaepx = coef*(-Hallcmpr);

                std::complex<double> iomegaMu0(0.0,omegam*mu0);

                std::complex<double> Ex_tmp = (E1x[i][j] + E1x[i+1][j])/2.0;
                std::complex<double> Er_tmp = (rL*E1r[i][j] + rR*E1r[i][j+1])/(2.0*r[j]);
                std::complex<double> Ep_tmp = E1p[i][j];
                std::complex<double> Jx_tmp = sigmaexr*Er_tmp + sigmaexp*Ep_tmp + sigmaexx*Ex_tmp;
                std::complex<double> Jr_tmp = sigmaerr*Er_tmp + sigmaerp*Ep_tmp + sigmaerx*Ex_tmp;
                std::complex<double> Jp_tmp = sigmaepr*Er_tmp + sigmaepp*Ep_tmp + sigmaepx*Ex_tmp;

                Pabs[i][j] = 0.5*real(Jx_tmp*conj(Ex_tmp) + Jr_tmp*conj(Er_tmp) + Jp_tmp*conj(Ep_tmp));
                power_sum = power_sum + Pabs[i][j]*(dx*dr*2.0*M_PI*r[j]);
            }
        }
    }

    std::cout  << std::endl;
    std::cout << "---------------------------------" << std::endl;
    std::cout << "Total power = " << power_sum << " W" << std::endl;

    //calculate reflection
    if(icon_mwRef == 1){

        //set refelence plane
        //---------------------------------
        double x_ref = 0.004;
        int i_ref = int(x_ref/dx + 1.5);
        //---------------------------------

        int i = i_ref;
        double x_tmp = (x[i] + x[i-1])/2.0;
        //calculate current I (use averaged value in r-direction for robustness)
        //---------------------------------
        int ncount = 0;
        std::complex<double>I_ref_tmp(0.0,0.0);
        for (int j=j_flr_bl[5][0];j<=j_flr_bl[5][1];j++){
            double r_tmp = (r[j] + r[j-1])/2.0;
            //double qR = (r[j] + r[j+1])/(2.0*r[j]);

            std::complex<double> H_ref = (E1r[i][j] - E1r[i-1][j])/dx/(omegam*mu0)*std::complex<double>(0,1.0);
            I_ref_tmp = I_ref_tmp + 2.0*M_PI*r_tmp*H_ref;

            //std::cout << "j = " << j << " - 1/2, r = "<< r_tmp << " , H_ref = "<<H_ref<<", I_ref = "<<2.0*M_PI*r_tmp*H_ref<< std::endl;
            ncount = ncount + 1;
        }
        I_ref_tmp = I_ref_tmp/double(ncount);
        //---------------------------------

        //calculate current V
        //---------------------------------
        std::complex<double> V_ref_tmp(0.0,0.0);
        for (int j=j_flc_bl[5][0];j<=j_flc_bl[5][1];j++){
            //double qL = (r[j] + r[j-1])/(2.0*r[j]);
            //double qR = (r[j] + r[j+1])/(2.0*r[j]);
            //double r_tmp = (r[j] + r[j-1])/2.0;
            double rL = (r[j] + r[j-1])/2.0;
            double rR = (r[j] + r[j+1])/2.0;

            //std::complex<double> Er_tmp = (E1r[i][j] + E1r[i-1][j])/2.0;
            std::complex<double> Er_tmp = (rR*(E1r[i][j+1] + E1r[i-1][j+1]) + rL*(E1r[i][j] + E1r[i-1][j]))/(4.0*r[j]);
            
            double c0 = sqrt(1.0/(eps0*epsr[i][j]*mu0));
            
            //std::cout << "Check ,i = "<<i <<", j = "<< j << ", "<< std::abs(Er_tmp)*std::abs(Er_tmp)
            //    //<< ","<<std::abs(Er_Fw[j+1])<< ","<<std::abs(Er_Fw[j])
            //    << ","<<(Er_tmp)
            //    << ","<<pow(std::abs((Er_tmp - Er_anly)/Er_anly),2)*100<< " %"
            //    << std::endl;

            V_ref_tmp = V_ref_tmp + Er_tmp*dr;
            //std::cout << "j = " << j <<", Er = "<< Er_tmp << ", rEr = "<< r_tmp*Er_tmp << ", V_ref = "<<V_ref << ", V_ref_th = "<<r_tmp*Er_tmp*log(r2/r1)<< std::endl;
        }
        //---------------------------------

        std::complex<double> Z_ref_tmp = V_ref_tmp/(I_ref_tmp + 1e-100);
        
        std::cout << "Impedance at reference plane x = " << x_tmp << std::endl;
        std::cout << "Z_ref = " << Z_ref_tmp << ", V_ref = " << V_ref_tmp << ", I_ref = " << I_ref_tmp 
            << ", |Z_ref| = " << std::abs(Z_ref_tmp) 
            << ", arg(Z_ref) = " << std::arg(Z_ref_tmp)*180/M_PI << " deg"<< std::endl;

        std::complex<double> a2 = (V_ref_tmp - Z0_base*I_ref_tmp)/(2.0*sqrt(Z0_base));
        std::complex<double> b2 = (V_ref_tmp + Z0_base*I_ref_tmp)/(2.0*sqrt(Z0_base));

        std::complex<double> S11(S11_mag * std::cos(S11_arg_deg/180.0*M_PI), S11_mag * std::sin(S11_arg_deg/180.0*M_PI));
        std::complex<double> S21(S21_mag * std::cos(S21_arg_deg/180.0*M_PI), S21_mag * std::sin(S21_arg_deg/180.0*M_PI));
        std::complex<double> S12(S12_mag * std::cos(S12_arg_deg/180.0*M_PI), S12_mag * std::sin(S12_arg_deg/180.0*M_PI));
        std::complex<double> S22(S22_mag * std::cos(S22_arg_deg/180.0*M_PI), S22_mag * std::sin(S22_arg_deg/180.0*M_PI));

        std::complex<double> a1 = (b2 - S22*a2)/S21;
        std::complex<double> b1 = S11*a1 + S12*a2;
        double P_fwd = std::norm(a1);
        double P_ref = std::norm(b1);

        std::cout << std::endl;
        std::cout << "Fwd power = " << P_fwd << " W" << std::endl;
        std::cout << "Ref power = " << P_ref << " W" << std::endl;
        std::cout << "Transmitted power = " << P_fwd - P_ref << " W" << std::endl;
        
        
        double ratio = Pmw/P_fwd;

        std::cout << "power ratio = " << ratio << std::endl;
        std::cout << "J1r_exc ="  << J1r_exc;
        J1r_exc = J1r_exc*sqrt(ratio);
        std::cout << " -> " << J1r_exc << std::endl;
        
        //電力調整
        for (int i=0;i<ni+1;i++){
            for (int j=0;j<nj+1;j++){
                Pabs[i][j] = Pabs[i][j]*ratio;
                J1r[i][j] = J1r[i][j]*sqrt(ratio);
            }
        }
       
    }else{

        double ratio = Pmw/power_sum;
        
        //std::cout << "J1r_exc1 =" << J1r_exc <<","<<std::real(J1r_exc) <<","<<sqrt(ratio)<<","<<J1r_exc*sqrt(ratio)<< std::endl;
        //電力調整
        std::cout << "power ratio = " << ratio << std::endl;
        std::cout << "J1r_exc = " << J1r_exc;
        J1r_exc = J1r_exc*sqrt(ratio);
        std::cout << " -> " << J1r_exc << std::endl;

        for (int i=0;i<ni+1;i++){
            for (int j=0;j<nj+1;j++){
                Pabs[i][j] = Pabs[i][j]*ratio;
                J1r[i][j] = J1r[i][j]*sqrt(ratio);
            }
        }
        //std::cout << "J1r_exc2 =" << J1r_exc <<","<<std::real(J1r_exc) <<","<<sqrt(ratio)<<","<<J1r_exc*sqrt(ratio)<< std::endl;
        

    }

    std::cout << "---------------------------------" << std::endl;
    std::cout  << std::endl;

}

