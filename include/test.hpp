#pragma once
#include <string>
#include <iostream>
#include <vector>
#include <sys/stat.h>

#include "constants.hpp"
#include "params.hpp"
#include "miniFuncs.hpp"
#include "arrays.hpp"

using namespace std;

//*****************************************************************
//**                                                             **
//**           class InputModule                                 **
//**                                                             **
//*****************************************************************
class InitialFuncs
{
    private:

    public:
        void iniParam(Params &pm,GridCenter &gc,GridInterfaceX &gx,GridInterfaceR &gr);
        void makeBoundary(Params &pm, GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr, GridK &gk, MicrowaveBC &mb);
        void makeBoundary_impedanceTest(Params &pm, GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr, GridK &gk, MicrowaveBC &mb);
};

//*****************************************************************
//**                                                             **
//**           void makeBoundary                                 **
//**                                                             **
//*****************************************************************
void InitialFuncs::makeBoundary_impedanceTest(Params &pm, GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr, GridK &gk, MicrowaveBC &mb)
{
    double x_tmp = 0.0;
    double r_tmp = 0.0;
    

    pm.x1 = 0.0e-3;
    //x2 =  5.6e-3 + 7.8e-3; //5.5 mm
    pm.x2 = pm.xR;
    pm.x3 = 10.0e-3 + 7.8e-3;
    pm.x4 = 15.0e-3 + 7.8e-3;
    pm.x5 = 20.0e-3 + 7.8e-3;
    pm.x6 =  0.0e-3 + 7.8e-3;
    
    pm.r1 =  0.8e-3;
    pm.r2 =  2.0e-3;
    pm.r3 =  4.0e-3;
    pm.r4 =  9.0e-3;
    pm.r5 =  2.0e-3;
    pm.r6 = 14.0e-3;

    gc.i_flc_bl[5][0] = int(pm.x1/pm.dx + 0.5) + 1; //1
    gc.j_flc_bl[5][0] = int(pm.r1/pm.dr + 0.5) + 1; //5
    gc.i_flc_bl[5][1] = int(pm.x2/pm.dx - 0.5) + 1; //20
    gc.j_flc_bl[5][1] = int(pm.r2/pm.dr - 0.5) + 1; //10

    gx.i_flx_bl[5][0] = int(pm.x1/pm.dx + 0.5) + 2; //2
    gx.j_flx_bl[5][0] = int(pm.r1/pm.dr + 0.5) + 1; //5
    gx.i_flx_bl[5][1] = int(pm.x2/pm.dx - 0.5) + 1; //21
    gx.j_flx_bl[5][1] = int(pm.r2/pm.dr - 0.5) + 1; //10

    //std::cout << "flx" << std::endl;
    //std::cout << gx.i_flx_bl[0][0] << " , "<<gx.j_flx_bl[0][0] << " , "<< gx.i_flx_bl[0][1] << " , "<<gx.j_flx_bl[0][1] << std::endl;
    //std::cout << gx.i_flx_bl[1][0] << " , "<<gx.j_flx_bl[1][0] << " , "<< gx.i_flx_bl[1][1] << " , "<<gx.j_flx_bl[1][1] << std::endl;
    //std::cout << gx.i_flx_bl[2][0] << " , "<<gx.j_flx_bl[2][0] << " , "<< gx.i_flx_bl[2][1] << " , "<<gx.j_flx_bl[2][1] << std::endl;
    //std::cout << gx.i_flx_bl[3][0] << " , "<<gx.j_flx_bl[3][0] << " , "<< gx.i_flx_bl[3][1] << " , "<<gx.j_flx_bl[3][1] << std::endl;
    //std::cout << gx.i_flx_bl[4][0] << " , "<<gx.j_flx_bl[4][0] << " , "<< gx.i_flx_bl[4][1] << " , "<<gx.j_flx_bl[4][1] << std::endl;
    //std::cout << gx.i_flx_bl[5][0] << " , "<<gx.j_flx_bl[5][0] << " , "<< gx.i_flx_bl[5][1] << " , "<<gx.j_flx_bl[5][1] << std::endl;


    gr.i_flr_bl[5][0] = int(pm.x1/pm.dx + 0.5) + 1; //1
    gr.j_flr_bl[5][0] = int(pm.r1/pm.dr + 0.5) + 2; //6
    gr.i_flr_bl[5][1] = int(pm.x2/pm.dx - 0.5) + 1; //20
    gr.j_flr_bl[5][1] = int(pm.r2/pm.dr - 0.5) + 1; //10

    //std::cout << "flr" << std::endl;
    //std::cout << gr.i_flr_bl[0][0] << " , "<<gr.j_flr_bl[0][0] << " , "<< gr.i_flr_bl[0][1] << " , "<<gr.j_flr_bl[0][1] << std::endl;
    //std::cout << gr.i_flr_bl[1][0] << " , "<<gr.j_flr_bl[1][0] << " , "<< gr.i_flr_bl[1][1] << " , "<<gr.j_flr_bl[1][1] << std::endl;
    //std::cout << gr.i_flr_bl[2][0] << " , "<<gr.j_flr_bl[2][0] << " , "<< gr.i_flr_bl[2][1] << " , "<<gr.j_flr_bl[2][1] << std::endl;
    //std::cout << gr.i_flr_bl[3][0] << " , "<<gr.j_flr_bl[3][0] << " , "<< gr.i_flr_bl[3][1] << " , "<<gr.j_flr_bl[3][1] << std::endl;
    //std::cout << gr.i_flr_bl[4][0] << " , "<<gr.j_flr_bl[4][0] << " , "<< gr.i_flr_bl[4][1] << " , "<<gr.j_flr_bl[4][1] << std::endl;
    //std::cout << gr.i_flr_bl[5][0] << " , "<<gr.j_flr_bl[5][0] << " , "<< gr.i_flr_bl[5][1] << " , "<<gr.j_flr_bl[5][1] << std::endl;


    for (int iblock=5;iblock<6;iblock++){ 
        for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
            for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                gc.jdgBnd_Ep[i][j] = 1;
            }
        }
    }

    for (int iblock=5;iblock<6;iblock++){ 
        for (int i=gx.i_flx_bl[iblock][0];i<=gx.i_flx_bl[iblock][1];i++){ 
            for (int j=gx.j_flx_bl[iblock][0];j<=gx.j_flx_bl[iblock][1];j++){
                gx.jdgBnd_Ex[i][j] = 1;
            }
        }
    }

    for (int iblock=5;iblock<6;iblock++){ 
        for (int i=gr.i_flr_bl[iblock][0];i<=gr.i_flr_bl[iblock][1];i++){ 
            for (int j=gr.j_flr_bl[iblock][0];j<=gr.j_flr_bl[iblock][1];j++){
                gr.jdgBnd_Er[i][j] = 1;
            }
        }
    }

    //変換係数作成
    int kx_tmp = 0;
    int kr_tmp = 0;
    int kp_tmp = 0;
    int kfc_tmp = 0;
    int kfx_tmp = 0;
    int kfr_tmp = 0;

    for (int i=0;i<=pm.ni+1;i++){
        for (int j=0;j<=pm.nj+1;j++){
            

            if(gr.jdgBnd_Er[i][j]==1){
                gr.kr[i][j] = kr_tmp;
                kr_tmp++;
                gk.ikr.push_back(i);
                gk.jkr.push_back(j);
            }

            
        }
    }


    pm.nkx = gk.ikx.size();
    pm.nkr = gk.ikr.size();
    pm.nkp = gk.ikp.size();
    pm.nk = pm.nkx + pm.nkr + pm.nkp;
    std::cout << "nkx = " << pm.nkx << " nkr = " << pm.nkr << " nkp = " << pm.nkp << " nk = " << pm.nk  << std::endl;

    //Output boundary check file
    if(pm.icon_chk == 1){

        std::ofstream outputfile1("results/boundary0.csv");
        outputfile1<<"i,j,x,r,jdgBnd_Ep,jdgBnd_Ex,jdgBnd_Er,jdgBnd_flc,jdgBnd_flx,jdgBnd_flr,kx,kr,kp,zero" << std::endl;

        for(int i=0;i<=pm.ni+1;i++){
            for(int j=0;j<=pm.nj+1;j++){
                outputfile1<< i << ","<< j << "," << gc.x[i]<< ","<< gc.r[j]
                    << "," << gc.jdgBnd_Ep[i][j] << "," << gx.jdgBnd_Ex[i][j]<< "," << gr.jdgBnd_Er[i][j]
                    << "," << gc.jdgBnd_flc[i][j]<< "," << gx.jdgBnd_flx[i][j]<< "," << gr.jdgBnd_flr[i][j]
                    << "," << gx.kx[i][j]<< "," << gr.kr[i][j]<< "," << gc.kp[i][j]
                    << "," << 0.0<< std::endl;
            }
        }

        outputfile1.close();
    }
    

    //**************** 境界条件保持配列作成 ****************
    
    //Er
    //=============================================
    //********* Er境界条件 (左) *********
    for (int j=gc.j_flc_bl[5][0]+1;j<=gc.j_flc_bl[5][1];j++){
        int i=gc.i_flc_bl[5][0];
        mb.iBndWr.push_back(i);
        mb.jBndWr.push_back(j);
        mb.sBndWr.push_back(1); //開放
        //mb.sBndWr.push_back(7); //励振
    }
    //出力
    if(pm.icon_chk == 1){
        std::ofstream outputfileWr("results/BCr0.csv");
        outputfileWr <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<mb.iBndWr.size();k++){
            int i = mb.iBndWr[k];
            int j = mb.jBndWr[k];
            
            outputfileWr << i << ","<< j 
                << ","<< gc.x[i] << ","<< gc.r[j] << ","<<mb.sBndWr[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileWr.close();
    }
    //********* Er境界条件 (右) *********
    for (int j=gc.j_flc_bl[5][0]+1;j<=gc.j_flc_bl[5][1];j++){
        int i=gc.i_flc_bl[5][1];
        mb.iBndEr.push_back(i);
        mb.jBndEr.push_back(j);
        mb.sBndEr.push_back(1); //開放
    }
    //出力
    if(pm.icon_chk == 1){
        std::ofstream outputfileEr("results/BCr1.csv");
        outputfileEr <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<mb.iBndEr.size();k++){
            int i = mb.iBndEr[k];
            int j = mb.jBndEr[k];
            
            outputfileEr << i << ","<< j 
                << ","<< gc.x[i] << ","<< gc.r[j] << ","<<mb.sBndEr[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileEr.close();
    }
    //********* Er境界条件 (下) *********
    for (int i=gc.i_flc_bl[5][0];i<=gc.i_flc_bl[5][1];i++){
        int j=gc.j_flc_bl[5][0]+1;
        mb.iBndSr.push_back(i);
        mb.jBndSr.push_back(j);
        mb.sBndSr.push_back(2); //ガウス
    }
    //出力
    if(pm.icon_chk == 1){
        std::ofstream outputfileSr("results/BCr2.csv");
        outputfileSr <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<mb.iBndSr.size();k++){
            int i = mb.iBndSr[k];
            int j = mb.jBndSr[k];
            
            outputfileSr << i << ","<< j 
                << ","<< gc.x[i] << ","<< gc.r[j] << ","<<mb.sBndSr[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileSr.close();
    }
    //********* Er境界条件 (上) *********
    for (int i=gc.i_flc_bl[5][0];i<=gc.i_flc_bl[5][1];i++){
        int j=gc.j_flc_bl[5][1];
        mb.iBndNr.push_back(i);
        mb.jBndNr.push_back(j);
        mb.sBndNr.push_back(2); //ガウス
    }
    //出力
    if(pm.icon_chk == 1){
        std::ofstream outputfileNr("results/BCr3.csv");
        outputfileNr <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<mb.iBndNr.size();k++){
            int i = mb.iBndNr[k];
            int j = mb.jBndNr[k];
            
            outputfileNr << i << ","<< j 
                << ","<< gc.x[i] << ","<< gc.r[j] << ","<<mb.sBndNr[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileNr.close();
    }
    //=============================================

}
