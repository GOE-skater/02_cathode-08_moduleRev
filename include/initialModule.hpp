#pragma once
#include <string>
#include <iostream>
#include <vector>
#include <sys/stat.h>

#include "constants.hpp"
#include "params.hpp"
#include "miniFuncs.hpp"
#include "arrays.hpp"
#include "cmpMat2.hpp"

using namespace std;

//*****************************************************************
//**                                                             **
//**           class InputModule                                 **
//**                                                             **
//*****************************************************************
class InitialModule
{
    private:
        void sParamToABCD(
            CmpMat2& T,
            const complex<double>& S11,
            const complex<double>& S21,
            const complex<double>& S12,
            const complex<double>& S22,
            const double Z0
        );

        void ABCDToSParam(
            complex<double>& S11,
            complex<double>& S21,
            complex<double>& S12,
            complex<double>& S22,
            const CmpMat2& T,
            const double Z0
        );

    public:
        void iniParam(Params &pm,GridCenter &gc,GridInterfaceX &gx,GridInterfaceR &gr);
        void makeBoundary(Params &pm, GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr, MicrowaveBC &mb);
        void makeBoundary_impedanceTest(Params &pm, GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr, MicrowaveBC &mb);
        void makeProfile(Params &pm, GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr);
        
};

//*****************************************************************
//**                                                             **
//**           void iniParam                                     **
//**                                                             **
//*****************************************************************
void InitialModule::iniParam(Params &pm,GridCenter &gc,GridInterfaceX &gx,GridInterfaceR &gr)
{
    
    //******************* Initialization **************************

    std::cout << "Inizialization" << std::endl;
    cout << "---------------------------------" << endl;

    //domain setting
    //--------------------------------
    pm.Lx = pm.xR - pm.xL;
    pm.Lr = pm.rmax - pm.rmin;
    pm.dx = pm.Lx/double(pm.ni);
    pm.dr = pm.Lr/double(pm.nj);
    std::cout << "dx = " << pm.dx<< " dr = " << pm.dr << std::endl;
    //--------------------------------

    //time step
    //--------------------------------
    pm.dt = pm.dt_ini;
    //--------------------------------

    //microwave angular frequency
    //--------------------------------
    pm.omegam = 2.0*M_PI*pm.fmw;
    //--------------------------------

    //mesh Gegeration
    //--------------------------------
    for (int i = 0; i < pm.ni+2; i++){
        gc.x[i] = pm.xL+(i-0.5)*pm.dx;
    }
    for (int j = 0; j < pm.nj+2; j++){
        gc.r[j] = pm.rmin+(j-0.5)*pm.dr;
    }
    //--------------------------------

    /*
    //Give the r-directin current density Jr at x = 0.01 + dx/2, 0.003 < r < 0.006
    double i_exc = int((0.0002 - pm.xL)/pm.dx + 0.5); //5
    double j_exc_min = int((0.0008 - pm.rmin)/pm.dr + 0.5) + 2; //6
    double j_exc_max = int((0.0020 - pm.rmin)/pm.dr + 0.5); //10
    
    for (int j=j_exc_min;j<=j_exc_max;j++){
        //J1r[i_exc][j] = J1r_exc;

        double r_tmp = (gc.r[j] + gc.r[j-1])/2.0;
        gr.J1r_exc[i_exc][j] = pm.J1r_exc/r_tmp/pm.dx;
        //cout << j << ","<< gc.r[j] << ", gr.J1r[i_exc][j] = " << gr.J1r[i_exc][j] << endl;
        //std::cout << i_exc << ","<< j <<","<<gr.J1r[i_exc][j] << std::endl;
    }
    */

    //initialize S-parameter
    //--------------------------------
    pm.S11 = complex<double>(pm.S11_mag*cos(pm.S11_arg_deg/180.0*M_PI), pm.S11_mag*sin(pm.S11_arg_deg/180.0*M_PI));
    pm.S21 = complex<double>(pm.S21_mag*cos(pm.S21_arg_deg/180.0*M_PI), pm.S21_mag*sin(pm.S21_arg_deg/180.0*M_PI));
    pm.S12 = complex<double>(pm.S12_mag*cos(pm.S12_arg_deg/180.0*M_PI), pm.S12_mag*sin(pm.S12_arg_deg/180.0*M_PI));
    pm.S22 = complex<double>(pm.S22_mag*cos(pm.S22_arg_deg/180.0*M_PI), pm.S22_mag*sin(pm.S22_arg_deg/180.0*M_PI));
    //--------------------------------

    cout << "S-para (original), S11 = " << pm.S11 << ", S21 = " <<pm.S21 << ", S12 = " << pm.S12<< ", S22 = " << pm.S22<<endl;
    double phase_pre = pm.phase_pre_deg/180.0*M_PI;

    CmpMat2 Tpre{
        complex<double>(cos(phase_pre), 0.0),          complex<double>(0.0,1.0)*pm.Z0*sin(phase_pre),
        complex<double>(0.0,1.0)*sin(phase_pre)/pm.Z0, complex<double>(cos(phase_pre), 0.0)
    };

    //Convert measured/base S-parameters to ABCD matrix.
    //---------------------------------
    CmpMat2 Tbase;
    sParamToABCD(Tbase, pm.S11, pm.S21, pm.S12, pm.S22, pm.Z0);
    //---------------------------------

    //Compute the short-stub admittance before and after changing the stub length.
    //---------------------------------
    double c0 = sqrt(1.0/(ph::eps0*pm.epsr_diele*ph::mu0));
    double beta = pm.omegam/c0;
    const double Y0Stub = 1.0/pm.Z0;

    double tan0 = tan(beta*pm.h0_stab);
    double tanNew = tan(beta*(pm.h0_stab + pm.dh_stab));

    if (abs(tan0) < 1.0e-14 || abs(tanNew) < 1.0e-14) {
        throw runtime_error("short-stub admittance: tan(beta*l) is too small.");
    }

    const complex<double> Ys0 =
        -complex<double>(0.0,1.0)*Y0Stub/tan0;

    const complex<double> YsNew =
        -complex<double>(0.0,1.0)*Y0Stub/tanNew;

    const complex<double> dYs = YsNew - Ys0;
    //---------------------------------

    // ABCD matrix of the shunt admittance correction.
    //---------------------------------
    const CmpMat2 Tcorr{
        1.0,  0.0,
        dYs,  1.0
    };
    //---------------------------------

    //Apply the correction at the original stub junction.
    //---------------------------------
    const CmpMat2 Tstub =
        Tpre*Tcorr*Tpre.inv()*Tbase;
    //const CmpMat2 Tstub =
    //    Tcorr*Tpre.inv()*Tbase;
    //---------------------------------

    //Convert corrected ABCD matrix back to S-parameters.
    //---------------------------------
    ABCDToSParam(pm.S11, pm.S21, pm.S12, pm.S22, Tstub, pm.Z0);
    //---------------------------------

    cout << "S-para (w/ stab), S11 = " << pm.S11 << ", S21 = " <<pm.S21 << ", S12 = " << pm.S12<< ", S22 = " << pm.S22<<endl;
    
    //S-parameter change by stab length dH
    //--------------------------------
    //complex<double> A_base = ((1.0 + pm.S11)*(1.0 - pm.S22) + pm.S12*pm.S21)/(2.0*pm.S21 + 1e-100);
    //complex<double> B_base = pm.Z0*((1.0 + pm.S11)*(1.0 + pm.S22) - pm.S12*pm.S21)/(2.0*pm.S21 + 1e-100);
    //complex<double> C_base = 1.0/pm.Z0*((1.0 - pm.S11)*(1.0 - pm.S22) - pm.S12*pm.S21)/(2.0*pm.S21 + 1e-100);
    //complex<double> D_base = ((1.0 - pm.S11)*(1.0 + pm.S22) + pm.S12*pm.S21)/(2.0*pm.S21 + 1e-100);
    //--------------------------------

    //S-parameter change by cable length dL
    //--------------------------------
    complex<double> p = exp(complex<double>(0,-beta)*pm.dl_cable);

    pm.S11 = pm.S11;
    pm.S21 = pm.S21*p;
    pm.S12 = pm.S12*p;
    pm.S22 = pm.S22*p*p;
    cout << "S-para (w/ stab & cable), S11 = " << pm.S11 << ", S21 = " <<pm.S21 << ", S12 = " << pm.S12<< ", S22 = " << pm.S22<<endl;
    

    cout << "Phase Deley (by cable length) = " << beta*pm.dl_cable/M_PI*180 << " deg" << endl;
    //--------------------------------
    cout << "---------------------------------" << endl;
    
    std::string folder_name = ("./results");
    mkdir(folder_name.c_str(), 0777);
    
}

//*****************************************************************
//**                                                             **
//**           void makeBoundary                                 **
//**                                                             **
//*****************************************************************
void InitialModule::makeBoundary(Params &pm, GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr, MicrowaveBC &mb)
{
    
    double x_tmp = 0.0;
    double r_tmp = 0.0;
    
    //double x1 =  0.0e-3;
    //double x2 =  4.0e-3;
    //double x3 = 10.0e-3;
    //double x4 = 15.6e-3;
    //double x5 = 20.6e-3;
    //double x6 =  0.0e-3;

    pm.x1 =  0.0e-3;
    pm.x2 =  5.6e-3 + 7.8e-3; //5.5 mm
    pm.x3 = 10.0e-3 + 7.8e-3;
    pm.x4 = 15.0e-3 + 7.8e-3;
    pm.x5 = 20.0e-3 + 7.8e-3;
    pm.x6 =  0.0e-3 + 7.8e-3;
    
    pm.r1 =  0.6e-3;
    pm.r2 =  2.0e-3;
    pm.r3 =  4.0e-3;
    pm.r4 =  9.0e-3;
    pm.r5 =  2.0e-3;
    pm.r6 = 14.0e-3;

    //中性粒子の流入条件
    double Q_neutIn = pm.Q_neutIn_mgs*1e-6/pm.massi; //mg/s -> kg/s -> 個/s
    double S_in = 2.0*M_PI*pm.r4*pm.width_neutIn;
    pm.fn_In = Q_neutIn/S_in;
    
    gc.i_flc_bl[0][0] = int(pm.x6/pm.dx + 0.5) + 1; //1
    gc.j_flc_bl[0][0] = int(pm.r3/pm.dr + 0.5) + 1; //21
    gc.i_flc_bl[0][1] = int(pm.x2/pm.dx - 0.5) + 1; //20
    gc.j_flc_bl[0][1] = int(pm.r4/pm.dr - 0.5) + 1; //45

    gc.i_flc_bl[1][0] = int(pm.x2/pm.dx + 0.5) + 1; //21
    gc.j_flc_bl[1][0] = int(pm.r5/pm.dr + 0.5) + 1; //11
    gc.i_flc_bl[1][1] = int(pm.x4/pm.dx - 0.5) + 1; //75
    gc.j_flc_bl[1][1] = int(pm.r4/pm.dr - 0.5) + 1; //45

    gc.i_flc_bl[2][0] = int(pm.x2/pm.dx + 0.5) + 1; //21
    gc.j_flc_bl[2][0] = int(pm.r1/pm.dr + 0.5) + 1; //5
    gc.i_flc_bl[2][1] = int(pm.x3/pm.dx - 0.5) + 1; //50
    gc.j_flc_bl[2][1] = int(pm.r5/pm.dr - 0.5) + 1; //10

    gc.i_flc_bl[3][0] = int(pm.x3/pm.dx + 0.5) + 1; //51
    gc.j_flc_bl[3][0] = int(0.0/pm.dr + 0.5) + 1; //1
    gc.i_flc_bl[3][1] = int(pm.x5/pm.dx - 0.5) + 1; //100
    gc.j_flc_bl[3][1] = int(pm.r5/pm.dr - 0.5) + 1; //10

    gc.i_flc_bl[4][0] = int(pm.x5/pm.dx + 0.5) + 1; //101
    gc.j_flc_bl[4][0] = int(0.0/pm.dr + 0.5) + 1; //1
    gc.i_flc_bl[4][1] = int(pm.xR/pm.dx - 0.5) + 1; //158
    gc.j_flc_bl[4][1] = int(pm.rmax/pm.dr - 0.5) + 1; //70

    gc.i_flc_bl[5][0] = int(pm.x1/pm.dx + 0.5) + 1; //1
    gc.j_flc_bl[5][0] = int(pm.r1/pm.dr + 0.5) + 1; //5
    gc.i_flc_bl[5][1] = int(pm.x2/pm.dx - 0.5) + 1; //20
    gc.j_flc_bl[5][1] = int(pm.r2/pm.dr - 0.5) + 1; //10

    //std::cout << "flc" << std::endl;
    //std::cout << gc.i_flc_bl[0][0] << " , "<<gc.j_flc_bl[0][0] << " , "<< gc.i_flc_bl[0][1] << " , "<<gc.j_flc_bl[0][1] << std::endl;
    //std::cout << gc.i_flc_bl[1][0] << " , "<<gc.j_flc_bl[1][0] << " , "<< gc.i_flc_bl[1][1] << " , "<<gc.j_flc_bl[1][1] << std::endl;
    //std::cout << gc.i_flc_bl[2][0] << " , "<<gc.j_flc_bl[2][0] << " , "<< gc.i_flc_bl[2][1] << " , "<<gc.j_flc_bl[2][1] << std::endl;
    //std::cout << gc.i_flc_bl[3][0] << " , "<<gc.j_flc_bl[3][0] << " , "<< gc.i_flc_bl[3][1] << " , "<<gc.j_flc_bl[3][1] << std::endl;
    //std::cout << gc.i_flc_bl[4][0] << " , "<<gc.j_flc_bl[4][0] << " , "<< gc.i_flc_bl[4][1] << " , "<<gc.j_flc_bl[4][1] << std::endl;
    //std::cout << gc.i_flc_bl[5][0] << " , "<<gc.j_flc_bl[5][0] << " , "<< gc.i_flc_bl[5][1] << " , "<<gc.j_flc_bl[5][1] << std::endl;

    gx.i_flx_bl[0][0] = int(pm.x6/pm.dx + 0.5) + 2; //2
    gx.j_flx_bl[0][0] = int(pm.r3/pm.dr + 0.5) + 1; //21
    gx.i_flx_bl[0][1] = int(pm.x2/pm.dx - 0.5) + 2; //21
    gx.j_flx_bl[0][1] = int(pm.r4/pm.dr - 0.5) + 1; //47

    gx.i_flx_bl[1][0] = int(pm.x2/pm.dx + 0.5) + 2; //22
    gx.j_flx_bl[1][0] = int(pm.r5/pm.dr + 0.5) + 1; //11
    gx.i_flx_bl[1][1] = int(pm.x4/pm.dx - 0.5) + 1; //78
    gx.j_flx_bl[1][1] = int(pm.r4/pm.dr - 0.5) + 1; //47

    gx.i_flx_bl[2][0] = int(pm.x2/pm.dx + 0.5) + 2; //22
    gx.j_flx_bl[2][0] = int(pm.r1/pm.dr + 0.5) + 1; //5
    gx.i_flx_bl[2][1] = int(pm.x3/pm.dx - 0.5) + 2; //51
    gx.j_flx_bl[2][1] = int(pm.r5/pm.dr - 0.5) + 1; //10

    gx.i_flx_bl[3][0] = int(pm.x3/pm.dx + 0.5) + 2; //52
    gx.j_flx_bl[3][0] = int(   0.0/pm.dr + 0.5) + 1; //1
    gx.i_flx_bl[3][1] = int(pm.x5/pm.dx - 0.5) + 2; //104
    gx.j_flx_bl[3][1] = int(pm.r5/pm.dr - 0.5) + 1; //10

    gx.i_flx_bl[4][0] = int(pm.x5/pm.dx + 0.5) + 2; //105
    gx.j_flx_bl[4][0] = int(   0.0/pm.dr + 0.5) + 1; //1
    gx.i_flx_bl[4][1] = int(pm.xR/pm.dx - 0.5) + 1; //158
    gx.j_flx_bl[4][1] = int(pm.rmax/pm.dr - 0.5) + 1; //70

    gx.i_flx_bl[5][0] = int(pm.x1/pm.dx + 0.5) + 2; //2
    gx.j_flx_bl[5][0] = int(pm.r1/pm.dr + 0.5) + 1; //5
    gx.i_flx_bl[5][1] = int(pm.x2/pm.dx - 0.5) + 2; //21
    gx.j_flx_bl[5][1] = int(pm.r2/pm.dr - 0.5) + 1; //10

    //std::cout << "flx" << std::endl;
    //std::cout << i_flx_bl[0][0] << " , "<<j_flx_bl[0][0] << " , "<< i_flx_bl[0][1] << " , "<<j_flx_bl[0][1] << std::endl;
    //std::cout << i_flx_bl[1][0] << " , "<<j_flx_bl[1][0] << " , "<< i_flx_bl[1][1] << " , "<<j_flx_bl[1][1] << std::endl;
    //std::cout << i_flx_bl[2][0] << " , "<<j_flx_bl[2][0] << " , "<< i_flx_bl[2][1] << " , "<<j_flx_bl[2][1] << std::endl;
    //std::cout << i_flx_bl[3][0] << " , "<<j_flx_bl[3][0] << " , "<< i_flx_bl[3][1] << " , "<<j_flx_bl[3][1] << std::endl;
    //std::cout << i_flx_bl[4][0] << " , "<<j_flx_bl[4][0] << " , "<< i_flx_bl[4][1] << " , "<<j_flx_bl[4][1] << std::endl;
    //std::cout << i_flx_bl[5][0] << " , "<<j_flx_bl[5][0] << " , "<< i_flx_bl[5][1] << " , "<<j_flx_bl[5][1] << std::endl;

    gr.i_flr_bl[0][0] = int(pm.x6/pm.dx + 0.5) + 1; //1
    gr.j_flr_bl[0][0] = int(pm.r3/pm.dr + 0.5) + 2; //22
    gr.i_flr_bl[0][1] = int(pm.x2/pm.dx - 0.5) + 1; //20
    gr.j_flr_bl[0][1] = int(pm.r4/pm.dr - 0.5) + 1; //47

    gr.i_flr_bl[1][0] = int(pm.x2/pm.dx + 0.5) + 1; //21
    gr.j_flr_bl[1][0] = int(pm.r5/pm.dr + 0.5) + 1; //11
    gr.i_flr_bl[1][1] = int(pm.x4/pm.dx - 0.5) + 1; //78
    gr.j_flr_bl[1][1] = int(pm.r4/pm.dr - 0.5) + 1; //47

    gr.i_flr_bl[2][0] = int(pm.x2/pm.dx + 0.5) + 1; //21
    gr.j_flr_bl[2][0] = int(pm.r1/pm.dr + 0.5) + 2; //6
    gr.i_flr_bl[2][1] = int(pm.x3/pm.dx - 0.5) + 1; //50
    gr.j_flr_bl[2][1] = int(pm.r5/pm.dr - 0.5) + 1; //10

    gr.i_flr_bl[3][0] = int(pm.x3/pm.dx + 0.5) + 1; //51
    gr.j_flr_bl[3][0] = int(   0.0/pm.dr + 0.5) + 2; //2
    gr.i_flr_bl[3][1] = int(pm.x5/pm.dx - 0.5) + 1; //103
    gr.j_flr_bl[3][1] = int(pm.r5/pm.dr - 0.5) + 1; //10

    gr.i_flr_bl[4][0] = int(pm.x5/pm.dx + 0.5) + 1; //104
    gr.j_flr_bl[4][0] = int(   0.0/pm.dr + 0.5) + 2; //2
    gr.i_flr_bl[4][1] = int(pm.xR/pm.dx - 0.5) + 1; //158
    gr.j_flr_bl[4][1] = int(pm.rmax/pm.dr - 0.5) + 1; //70

    gr.i_flr_bl[5][0] = int(pm.x1/pm.dx + 0.5) + 1; //1
    gr.j_flr_bl[5][0] = int(pm.r1/pm.dr + 0.5) + 2; //6
    gr.i_flr_bl[5][1] = int(pm.x2/pm.dx - 0.5) + 1; //20
    gr.j_flr_bl[5][1] = int(pm.r2/pm.dr - 0.5) + 1; //10

    //std::cout << "flr" << std::endl;
    //std::cout << i_flr_bl[0][0] << " , "<<j_flr_bl[0][0] << " , "<< i_flr_bl[0][1] << " , "<<j_flr_bl[0][1] << std::endl;
    //std::cout << i_flr_bl[1][0] << " , "<<j_flr_bl[1][0] << " , "<< i_flr_bl[1][1] << " , "<<j_flr_bl[1][1] << std::endl;
    //std::cout << i_flr_bl[2][0] << " , "<<j_flr_bl[2][0] << " , "<< i_flr_bl[2][1] << " , "<<j_flr_bl[2][1] << std::endl;
    //std::cout << i_flr_bl[3][0] << " , "<<j_flr_bl[3][0] << " , "<< i_flr_bl[3][1] << " , "<<j_flr_bl[3][1] << std::endl;
    //std::cout << i_flr_bl[4][0] << " , "<<j_flr_bl[4][0] << " , "<< i_flr_bl[4][1] << " , "<<j_flr_bl[4][1] << std::endl;
    //std::cout << i_flr_bl[5][0] << " , "<<j_flr_bl[5][0] << " , "<< i_flr_bl[5][1] << " , "<<j_flr_bl[5][1] << std::endl;

    for (int iblock=0;iblock<5;iblock++){ 
        for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
            for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                gc.jdgBnd_flc[i][j] = 1;
            }
        }
    }

    for (int iblock=0;iblock<5;iblock++){ 
        for (int i=gx.i_flx_bl[iblock][0];i<=gx.i_flx_bl[iblock][1];i++){ 
            for (int j=gx.j_flx_bl[iblock][0];j<=gx.j_flx_bl[iblock][1];j++){
                gx.jdgBnd_flx[i][j] = 1;
            }
        }
    }

    for (int iblock=0;iblock<5;iblock++){ 
        for (int i=gr.i_flr_bl[iblock][0];i<=gr.i_flr_bl[iblock][1];i++){ 
            for (int j=gr.j_flr_bl[iblock][0];j<=gr.j_flr_bl[iblock][1];j++){
                gr.jdgBnd_flr[i][j] = 1;
            }
        }
    }

    for (int iblock=0;iblock<6;iblock++){ 
        for (int i=gc.i_flc_bl[iblock][0];i<=gc.i_flc_bl[iblock][1];i++){ 
            for (int j=gc.j_flc_bl[iblock][0];j<=gc.j_flc_bl[iblock][1];j++){
                gc.jdgBnd_Ep[i][j] = 1;
            }
        }
    }

    for (int iblock=0;iblock<6;iblock++){ 
        for (int i=gx.i_flx_bl[iblock][0];i<=gx.i_flx_bl[iblock][1];i++){ 
            for (int j=gx.j_flx_bl[iblock][0];j<=gx.j_flx_bl[iblock][1];j++){
                gx.jdgBnd_Ex[i][j] = 1;
            }
        }
    }

    for (int iblock=0;iblock<6;iblock++){ 
        for (int i=gr.i_flr_bl[iblock][0];i<=gr.i_flr_bl[iblock][1];i++){ 
            for (int j=gr.j_flr_bl[iblock][0];j<=gr.j_flr_bl[iblock][1];j++){
                gr.jdgBnd_Er[i][j] = 1;
            }
        }
    }

    //Output boundary check file
    if(pm.flag_chk == 1){

        std::ofstream outputfile1("results/boundary0.csv");
        //outputfile1<<"i,j,x,r,jdgBnd_Ep,jdgBnd_Ex,jdgBnd_Er,jdgBnd_flc,jdgBnd_flx,jdgBnd_flr,kx,kr,kp,zero" << std::endl;
        outputfile1<<"i,j,x,r,jdgBnd_Ep,jdgBnd_Ex,jdgBnd_Er,jdgBnd_flc,jdgBnd_flx,jdgBnd_flr,zero" << std::endl;

        for(int i=0;i<=pm.ni+1;i++){
            for(int j=0;j<=pm.nj+1;j++){
                outputfile1<< i << ","<< j << "," << gc.x[i]<< ","<< gc.r[j]
                    << "," << gc.jdgBnd_Ep[i][j] << "," << gx.jdgBnd_Ex[i][j]<< "," << gr.jdgBnd_Er[i][j]
                    << "," << gc.jdgBnd_flc[i][j]<< "," << gx.jdgBnd_flx[i][j]<< "," << gr.jdgBnd_flr[i][j]
                    //<< "," << gx.kx[i][j]<< "," << gr.kr[i][j]<< "," << gc.kp[i][j]
                    << "," << 0.0<< std::endl;
            }
        }

        outputfile1.close();
    }
    
    //**************** 境界条件保持配列作成 ****************
    //********* Ex境界条件 (左) *********
    //z0
    for (int j=gc.j_flc_bl[0][0]+1;j<=gc.j_flc_bl[0][1]-1;j++){
        int i=gc.i_flc_bl[0][0]+1;
        mb.iBndWx.push_back(i);
        mb.jBndWx.push_back(j);
        mb.sBndWx.push_back(2); //ガウス
    }
    //z0
    {
        int i=gc.i_flc_bl[0][0]+1;
        int j=gc.j_flc_bl[0][0];
        mb.iBndWx.push_back(i);
        mb.jBndWx.push_back(j);
        mb.sBndWx.push_back(3); //凹角 下側壁
    }
    //z0
    {
        int i=gc.i_flc_bl[0][0]+1;
        int j=gc.j_flc_bl[0][1];
        mb.iBndWx.push_back(i);
        mb.jBndWx.push_back(j);
        mb.sBndWx.push_back(4); //凹角 上側壁
    }
    //z1
    for (int j=gc.j_flc_bl[5][1]+1;j<=gc.j_flc_bl[0][0]-1;j++){
        int i=gc.i_flc_bl[2][0]+1;
        mb.iBndWx.push_back(i);
        mb.jBndWx.push_back(j);
        mb.sBndWx.push_back(2); //ガウス
    }
    //z2
    for (int j=gc.j_flc_bl[3][0];j<=gc.j_flc_bl[2][0]-1;j++){
        int i=gc.i_flc_bl[3][0]+1;
        mb.iBndWx.push_back(i);
        mb.jBndWx.push_back(j);
        mb.sBndWx.push_back(2); //ガウス
    }
    //z4
    for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[4][1]-1;j++){
        int i=gc.i_flc_bl[4][0]+1;
        mb.iBndWx.push_back(i);
        mb.jBndWx.push_back(j);
        mb.sBndWx.push_back(2); //ガウス
    }
    //z4
    {
        int i=gc.i_flc_bl[4][0]+1;
        int j=gc.j_flc_bl[4][1];
        mb.iBndWx.push_back(i);
        mb.jBndWx.push_back(j);
        mb.sBndWx.push_back(4); //凹角 上側壁
    }
    //z6
    for (int j=gc.j_flc_bl[5][0];j<=gc.j_flc_bl[5][1];j++){
        int i=gc.i_flc_bl[5][0]+1;
        mb.iBndWx.push_back(i);
        mb.jBndWx.push_back(j);
        mb.sBndWx.push_back(0); //ディリクレ
    }
    //出力
    if(pm.flag_chk == 1){
        std::ofstream outputfileWx("results/BCx0.csv");
        outputfileWx <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<mb.iBndWx.size();k++){
            int i = mb.iBndWx[k];
            int j = mb.jBndWx[k];
            
            outputfileWx << i << ","<< j 
                << ","<< gc.x[i] << ","<<gc.r[j] << ","<<mb.sBndWx[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileWx.close();
    }
    //********* Ex境界条件 (右) *********
    //z3
    for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[1][1]-1;j++){
        int i=gc.i_flc_bl[1][1];
        mb.iBndEx.push_back(i);
        mb.jBndEx.push_back(j);
        mb.sBndEx.push_back(2); //ガウス
    }
    //z3
    {
        int i=gc.i_flc_bl[1][1];
        int j=gc.j_flc_bl[1][1];
        mb.iBndEx.push_back(i);
        mb.jBndEx.push_back(j);
        mb.sBndEx.push_back(4); //凹角 上側壁
    }
    //z5
    for (int j=gc.j_flc_bl[4][0];j<=gc.j_flc_bl[4][1]-1;j++){
        int i=gc.i_flc_bl[4][1];
        mb.jBndEx.push_back(j);
        mb.iBndEx.push_back(i);
        mb.sBndEx.push_back(2); //ガウス
    }
    //z5
    {
        int i=gc.i_flc_bl[4][1];
        int j=gc.j_flc_bl[4][1];
        mb.iBndEx.push_back(i);
        mb.jBndEx.push_back(j);
        mb.sBndEx.push_back(4); //凹角 上側壁
    }
    //出力
    if(pm.flag_chk == 1){
        std::ofstream outputfileEx("results/BCx1.csv");
        outputfileEx <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<mb.iBndEx.size();k++){
            int i = mb.iBndEx[k];
            int j = mb.jBndEx[k];
            
            outputfileEx << i << ","<< j 
                << ","<< gc.x[i] << ","<<gc.r[j] << ","<<mb.sBndEx[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileEx.close();
    }
    //********* Ex境界条件 (下) *********
    //x0
    for (int i=gc.i_flc_bl[0][0]+1;i<=gc.i_flc_bl[0][1];i++){
        int j=gc.j_flc_bl[0][0];
        mb.iBndSx.push_back(i);
        mb.jBndSx.push_back(j);
        mb.sBndSx.push_back(0); //ディリクレ
    }
    //x0
    {
        int i=gc.i_flc_bl[0][1]+1;
        int j=gc.j_flc_bl[0][0];
        mb.iBndSx.push_back(i);
        mb.jBndSx.push_back(j);
        mb.sBndSx.push_back(6); //凸角 右側Open
    }
    //x2
    for (int i=gc.i_flc_bl[5][0]+1;i<=gc.i_flc_bl[2][1];i++){
        int j=gc.j_flc_bl[5][0];
        mb.iBndSx.push_back(i);
        mb.jBndSx.push_back(j);
        mb.sBndSx.push_back(0); //ディリクレ
    }
    //x2
    {
        int i=gc.i_flc_bl[2][1]+1;
        int j=gc.j_flc_bl[5][0];
        mb.iBndSx.push_back(i);
        mb.jBndSx.push_back(j);
        mb.sBndSx.push_back(6); //凸角 右側Open
    }
    //出力
    if(pm.flag_chk == 1){
        std::ofstream outputfileSx("results/BCx2.csv");
        outputfileSx <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<mb.iBndSx.size();k++){
            int i = mb.iBndSx[k];
            int j = mb.jBndSx[k];
            
            outputfileSx << i << ","<< j 
                << ","<< gc.x[i] << ","<<gc.r[j] << ","<<mb.sBndSx[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileSx.close();
    }
    //********* Ex境界条件 (上) *********
    //x1
    for (int i=gc.i_flc_bl[0][0]+1;i<=gc.i_flc_bl[1][1];i++){
        int j=gc.j_flc_bl[0][1];
        mb.iBndNx.push_back(i);
        mb.jBndNx.push_back(j);
        mb.sBndNx.push_back(0); //ディリクレ
    }
    //x3
    for (int i=gc.i_flc_bl[5][0]+1;i<=gc.i_flc_bl[5][1];i++){
        int j=gc.j_flc_bl[5][1];
        mb.iBndNx.push_back(i);
        mb.jBndNx.push_back(j);
        mb.sBndNx.push_back(0); //ディリクレ
    }
    //x3
    {
        int i=gc.i_flc_bl[5][1]+1;
        int j=gc.j_flc_bl[5][1];
        mb.iBndNx.push_back(i);
        mb.jBndNx.push_back(j);
        mb.sBndNx.push_back(6); //凸角 右側Open
    }
    //x4
    for (int i=gc.i_flc_bl[1][1]+2;i<=gc.i_flc_bl[3][1];i++){
        int j=gc.j_flc_bl[3][1];
        mb.iBndNx.push_back(i);
        mb.jBndNx.push_back(j);
        mb.sBndNx.push_back(0); //ディリクレ
    }
    //x4
    {
        int i=gc.i_flc_bl[1][1]+1;
        int j=gc.j_flc_bl[3][1];
        mb.iBndNx.push_back(i);
        mb.jBndNx.push_back(j);
        mb.sBndNx.push_back(5); //凸角 左側Open
    }
    //x4
    {
        int i=gc.i_flc_bl[3][1]+1;
        int j=gc.j_flc_bl[3][1];
        mb.iBndNx.push_back(i);
        mb.jBndNx.push_back(j);
        mb.sBndNx.push_back(6); //凸角 右側Open
    }
    //x5
    for (int i=gc.i_flc_bl[4][0]+1;i<=gc.i_flc_bl[4][1];i++){
        int j=gc.j_flc_bl[4][1];
        mb.iBndNx.push_back(i);
        mb.jBndNx.push_back(j);
        mb.sBndNx.push_back(1); //開放
    }
    //x5 (仮)
    //for (int i=gc.i_flc_bl[4][0]+1;i<=gc.i_flc_bl[4][1];i++){
    //    int j=gc.j_flc_bl[4][1];
    //    mb.iBndNx.push_back(i);
    //    mb.jBndNx.push_back(j);
    //    mb.sBndNx.push_back(0); //ディリクレ
    //}

    //出力
    if(pm.flag_chk == 1){
        std::ofstream outputfileNx("results/BCx3.csv");
        outputfileNx <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<mb.iBndNx.size();k++){
            int i = mb.iBndNx[k];
            int j = mb.jBndNx[k];
            
            outputfileNx << i << ","<< j 
                << ","<< gc.x[i] << ","<<gc.r[j] << ","<<mb.sBndNx[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileNx.close();
    }

    //********* Er境界条件 (左) *********
    //z0
    for (int j=gc.j_flc_bl[0][0]+1;j<=gc.j_flc_bl[0][1];j++){
        int i=gc.i_flc_bl[0][0];
        mb.iBndWr.push_back(i);
        mb.jBndWr.push_back(j);
        mb.sBndWr.push_back(0); //ディリクレ
    }
    //z1
    for (int j=gc.j_flc_bl[5][1]+2;j<=gc.j_flc_bl[0][0]-1;j++){
        int i=gc.i_flc_bl[1][0];
        mb.iBndWr.push_back(i);
        mb.jBndWr.push_back(j);
        mb.sBndWr.push_back(0); //ディリクレ
    }
    //z1
    {
        int i=gc.i_flc_bl[1][0];
        int j=gc.j_flc_bl[5][1]+1;
        mb.iBndWr.push_back(i);
        mb.jBndWr.push_back(j);
        mb.sBndWr.push_back(5); //凸角 下側Open
    }
    //z1
    {
        int i=gc.i_flc_bl[1][0];
        int j=gc.j_flc_bl[0][0];
        mb.iBndWr.push_back(i);
        mb.jBndWr.push_back(j);
        mb.sBndWr.push_back(6); //凸角 上側Open
    }
    //z2
    for (int j=gc.j_flc_bl[3][0]+1;j<=gc.j_flc_bl[2][0]-1;j++){
        int i=gc.i_flc_bl[3][0];
        mb.iBndWr.push_back(i);
        mb.jBndWr.push_back(j);
        mb.sBndWr.push_back(0); //ディリクレ
    }
    //z2
    {
        int i=gc.i_flc_bl[3][0];
        int j=gc.j_flc_bl[2][0];
        mb.iBndWr.push_back(i);
        mb.jBndWr.push_back(j);
        mb.sBndWr.push_back(6); //凸角 上側Open
    }
    //z4
    for (int j=gc.j_flc_bl[1][0]+1;j<=gc.j_flc_bl[4][1];j++){
        int i=gc.i_flc_bl[4][0];
        mb.iBndWr.push_back(i);
        mb.jBndWr.push_back(j);
        mb.sBndWr.push_back(0); //ディリクレ
    }
    //z4
    {
        int i=gc.i_flc_bl[4][0];
        int j=gc.j_flc_bl[1][0];
        mb.iBndWr.push_back(i);
        mb.jBndWr.push_back(j);
        mb.sBndWr.push_back(5); //凸角 下側Open
    }
    //z6
    for (int j=gc.j_flc_bl[5][0]+1;j<=gc.j_flc_bl[5][1];j++){
        int i=gc.i_flc_bl[5][0];
        mb.iBndWr.push_back(i);
        mb.jBndWr.push_back(j);
        //mb.sBndWr.push_back(1); //開放
        mb.sBndWr.push_back(7); //励振
    }
    //出力
    if(pm.flag_chk == 1){
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
    //z3
    for (int j=gc.j_flc_bl[1][0]+1;j<=gc.j_flc_bl[1][1];j++){
        int i=gc.i_flc_bl[1][1];
        mb.iBndEr.push_back(i);
        mb.jBndEr.push_back(j);
        mb.sBndEr.push_back(0); //ディリクレ
    }
    //z3
    {
        int i=gc.i_flc_bl[1][1];
        int j=gc.j_flc_bl[1][0];
        mb.iBndEr.push_back(i);
        mb.jBndEr.push_back(j);
        mb.sBndEr.push_back(5); //凸角 下側Open
    }
    //z5
    for (int j=gc.j_flc_bl[4][0]+1;j<=gc.j_flc_bl[4][1];j++){
        int i=gc.i_flc_bl[4][1];
        mb.iBndEr.push_back(i);
        mb.jBndEr.push_back(j);
        mb.sBndEr.push_back(0); //ディリクレ
    }
    //出力
    if(pm.flag_chk == 1){
        std::ofstream outputfileEr("results/BCr1.csv");
        outputfileEr <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<mb.iBndEr.size();k++){
            int j = mb.jBndEr[k];
            int i = mb.iBndEr[k];
            
            outputfileEr << i << ","<< j 
                << ","<< gc.x[i] << ","<< gc.r[j] << ","<<mb.sBndEr[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileEr.close();
    }
    //********* Er境界条件 (下) *********
    //x0
    for (int i=gc.i_flc_bl[0][0]+1;i<=gc.i_flc_bl[0][1];i++){
        int j=gc.j_flc_bl[0][0]+1;
        mb.iBndSr.push_back(i);
        mb.jBndSr.push_back(j);
        mb.sBndSr.push_back(2); //ガウス
    }
    //x0
    {
        int i=gc.i_flc_bl[0][0];
        int j=gc.j_flc_bl[0][0]+1;
        mb.iBndSr.push_back(i);
        mb.jBndSr.push_back(j);
        mb.sBndSr.push_back(3); //凹角 左側壁
    }
    //x2
    for (int i=gc.i_flc_bl[5][0];i<=gc.i_flc_bl[2][1];i++){
        int j=gc.j_flc_bl[5][0]+1;
        mb.iBndSr.push_back(i);
        mb.jBndSr.push_back(j);
        mb.sBndSr.push_back(2); //ガウス
    }
    //出力
    if(pm.flag_chk == 1){
        std::ofstream outputfileSr("results/BCr2.csv");
        outputfileSr <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<mb.iBndSr.size();k++){
            int i = mb.iBndSr[k];
            int j = mb.jBndSr[k];
            
            outputfileSr << i << ","<< j 
                << ","<< gc.x[i] << ","<< gc.r[j] << ","<< mb.sBndSr[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileSr.close();
    }
    //********* Er境界条件 (上) *********
    //x1
    for (int i=gc.i_flc_bl[0][0]+1;i<=gc.i_flc_bl[1][1]-1;i++){
        int j=gc.j_flc_bl[0][1];
        mb.iBndNr.push_back(i);
        mb.jBndNr.push_back(j);
        mb.sBndNr.push_back(2); //ガウス
    }
    //x1
    {
        int i=gc.i_flc_bl[0][0];
        int j=gc.j_flc_bl[0][1];
        mb.iBndNr.push_back(i);
        mb.jBndNr.push_back(j);
        mb.sBndNr.push_back(3); //凹角 左側壁
    }
    //x1
    {
        int i=gc.i_flc_bl[1][1];
        int j=gc.j_flc_bl[0][1];
        mb.iBndNr.push_back(i);
        mb.jBndNr.push_back(j);
        mb.sBndNr.push_back(4); //凹角 右側壁
    }
    //x3
    for (int i=gc.i_flc_bl[5][0];i<=gc.i_flc_bl[5][1];i++){
        int j=gc.j_flc_bl[5][1];
        mb.iBndNr.push_back(i);
        mb.jBndNr.push_back(j);
        mb.sBndNr.push_back(2); //ガウス
    }
    //x4
    for (int i=gc.i_flc_bl[1][1]+1;i<=gc.i_flc_bl[3][1];i++){
        int j=gc.j_flc_bl[3][1];
        mb.iBndNr.push_back(i);
        mb.jBndNr.push_back(j);
        mb.sBndNr.push_back(2); //ガウス
    }
    //x5
    for (int i=gc.i_flc_bl[4][0]+1;i<=gc.i_flc_bl[4][1]-1;i++){
        int j=gc.j_flc_bl[4][1];
        mb.iBndNr.push_back(i);
        mb.jBndNr.push_back(j);
        mb.sBndNr.push_back(2); //ガウス
    }
    //x5
    {
        int i=gc.i_flc_bl[4][0];
        int j=gc.j_flc_bl[4][1];
        mb.iBndNr.push_back(i);
        mb.jBndNr.push_back(j);
        mb.sBndNr.push_back(3); //凹角 左側壁
    }
    //x5
    {
        int i=gc.i_flc_bl[4][1];
        int j=gc.j_flc_bl[4][1];
        mb.iBndNr.push_back(i);
        mb.jBndNr.push_back(j);
        mb.sBndNr.push_back(4); //凹角 右側壁
    }
    //出力
    if(pm.flag_chk == 1){
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

    //********* Ep境界条件 (左) *********
    //z0
    for (int j=gc.j_flc_bl[0][0];j<=gc.j_flc_bl[0][1];j++){
        int i=gc.i_flc_bl[0][0];
        mb.iBndWp.push_back(i);
        mb.jBndWp.push_back(j);
        mb.sBndWp.push_back(0); //ディリクレ
    }
    //z1
    for (int j=gc.j_flc_bl[5][1]+1;j<=gc.j_flc_bl[0][0]-1;j++){
        int i=gc.i_flc_bl[2][0];
        mb.iBndWp.push_back(i);
        mb.jBndWp.push_back(j);
        mb.sBndWp.push_back(0); //ディリクレ
    }
    //z2
    for (int j=gc.j_flc_bl[3][0];j<=gc.j_flc_bl[2][0]-1;j++){
        int i=gc.i_flc_bl[3][0];
        mb.iBndWp.push_back(i);
        mb.jBndWp.push_back(j);
        mb.sBndWp.push_back(0); //ディリクレ
    }
    //z4
    for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[4][1];j++){
        int i=gc.i_flc_bl[4][0];
        mb.iBndWp.push_back(i);
        mb.jBndWp.push_back(j);
        mb.sBndWp.push_back(0); //ディリクレ
    }
    //z6
    for (int j=gc.j_flc_bl[5][0];j<=gc.j_flc_bl[5][1];j++){
        int i=gc.i_flc_bl[5][0];
        mb.iBndWp.push_back(i);
        mb.jBndWp.push_back(j);
        mb.sBndWp.push_back(1); //開放
    }
    //出力
    if(pm.flag_chk == 1){
        std::ofstream outputfileWp("results/BCp0.csv");
        outputfileWp <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<mb.iBndWp.size();k++){
            int i = mb.iBndWp[k];
            int j = mb.jBndWp[k];
            
            outputfileWp << i << ","<< j 
                << ","<< gc.x[i] << ","<< gc.r[j] << ","<<mb.sBndWp[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileWp.close();
    }
    //********* Ep境界条件 (右) *********
    //z3
    for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[1][1];j++){
        int i=gc.i_flc_bl[1][1];
        mb.iBndEp.push_back(i);
        mb.jBndEp.push_back(j);
        mb.sBndEp.push_back(0); //ディリクレ
    }
    //z5
    for (int j=gc.j_flc_bl[4][0];j<=gc.j_flc_bl[4][1];j++){
        int i=gc.i_flc_bl[4][1];
        mb.iBndEp.push_back(i);
        mb.jBndEp.push_back(j);
        mb.sBndEp.push_back(0); //ディリクレ
    }
    //出力
    if(pm.flag_chk == 1){
        std::ofstream outputfileEp("results/BCp1.csv");
        outputfileEp <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<mb.iBndEp.size();k++){
            int i = mb.iBndEp[k];
            int j = mb.jBndEp[k];
            
            outputfileEp << i << ","<< j 
                << ","<< gc.x[i] << ","<< gc.r[j] << ","<<mb.sBndEp[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileEp.close();
    }
    //********* Ep境界条件 (下) *********
    //x0
    for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[0][1];i++){
        int j=gc.j_flc_bl[0][0];
        mb.iBndSp.push_back(i);
        mb.jBndSp.push_back(j);
        mb.sBndSp.push_back(0); //ディリクレ
    }
    //x2
    for (int i=gc.i_flc_bl[5][0];i<=gc.i_flc_bl[2][1];i++){
        int j=gc.j_flc_bl[5][0];
        mb.iBndSp.push_back(i);
        mb.jBndSp.push_back(j);
        mb.sBndSp.push_back(0); //ディリクレ
    }
    //出力
    if(pm.flag_chk == 1){
        std::ofstream outputfileSp("results/BCp2.csv");
        outputfileSp <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<mb.iBndSp.size();k++){
            int i = mb.iBndSp[k];
            int j = mb.jBndSp[k];
            
            outputfileSp << i << ","<< j 
                << ","<< gc.x[i] << ","<< gc.r[j] << ","<<mb.sBndSp[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileSp.close();
    }
    //********* Ep境界条件 (上) *********
    //x1
    for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[1][1];i++){
        int j=gc.j_flc_bl[0][1];
        mb.iBndNp.push_back(i);
        mb.jBndNp.push_back(j);
        mb.sBndNp.push_back(0); //ディリクレ
    }
    //x3
    for (int i=gc.i_flc_bl[5][0];i<=gc.i_flc_bl[5][1];i++){
        int j=gc.j_flc_bl[5][1];
        mb.iBndNp.push_back(i);
        mb.jBndNp.push_back(j);
        mb.sBndNp.push_back(0); //ディリクレ
    }
    //x4
    for (int i=gc.i_flc_bl[1][1]+1;i<=gc.i_flc_bl[3][1];i++){
        int j=gc.j_flc_bl[3][1];
        mb.iBndNp.push_back(i);
        mb.jBndNp.push_back(j);
        mb.sBndNp.push_back(0); //ディリクレ
    }
    //x5
    for (int i=gc.i_flc_bl[4][0];i<=gc.i_flc_bl[4][1];i++){
        int j=gc.j_flc_bl[4][1];
        mb.iBndNp.push_back(i);
        mb.jBndNp.push_back(j);
        mb.sBndNp.push_back(1); //開放
    }
    //x5 (仮)
    //for (int i=gc.i_flc_bl[4][0];i<=gc.i_flc_bl[4][1];i++){
    //    int j=gc.j_flc_bl[4][1];
    //    mb.iBndNp.push_back(i);
    //    mb.jBndNp.push_back(j);
    //    mb.sBndNp.push_back(0); //ディリクレ
    //}

    //出力
    if(pm.flag_chk == 1){
        std::ofstream outputfileNp("results/BCp3.csv");
        outputfileNp <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<mb.iBndNp.size();k++){
            int i = mb.iBndNp[k];
            int j = mb.jBndNp[k];
            
            outputfileNp << i << ","<< j 
                << ","<< gc.x[i] << ","<< gc.r[j] << ","<<mb.sBndNp[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileNp.close();
    }
}

//*****************************************************************
//**                                                             **
//**           void makeBoundary                                 **
//**                                                             **
//*****************************************************************
void InitialModule::makeBoundary_impedanceTest(Params &pm, GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr, MicrowaveBC &mb)
{
    double x_tmp = 0.0;
    double r_tmp = 0.0;

    pm.n_bl = 1;
    
    pm.x1 = pm.xL;
    pm.x2 = pm.xR;
   
    pm.r1 =  pm.rmin;
    pm.r2 =  pm.rmax;

    gc.i_flc_bl[0][0] = 1;
    gc.j_flc_bl[0][0] = 1;
    gc.i_flc_bl[0][1] = pm.ni;
    gc.j_flc_bl[0][1] = pm.nj;

    gx.i_flx_bl[0][0] = 2;
    gx.j_flx_bl[0][0] = 1;
    gx.i_flx_bl[0][1] = pm.ni;
    gx.j_flx_bl[0][1] = pm.nj;

    gr.i_flr_bl[0][0] = 1;
    gr.j_flr_bl[0][0] = 2;
    gr.i_flr_bl[0][1] = pm.ni;
    gr.j_flr_bl[0][1] = pm.nj;

    std::cout << "flc" << std::endl;
    std::cout << gc.i_flc_bl[0][0] << " , "<<gc.j_flc_bl[0][0] << " , "<< gc.i_flc_bl[0][1] << " , "<<gc.j_flc_bl[0][1] << std::endl;
    
    std::cout << "flx" << std::endl;
    std::cout << gx.i_flx_bl[0][0] << " , "<<gx.j_flx_bl[0][0] << " , "<< gx.i_flx_bl[0][1] << " , "<<gx.j_flx_bl[0][1] << std::endl;
    
    std::cout << "flr" << std::endl;
    std::cout << gr.i_flr_bl[0][0] << " , "<<gr.j_flr_bl[0][0] << " , "<< gr.i_flr_bl[0][1] << " , "<<gr.j_flr_bl[0][1] << std::endl;
    


    for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[0][1];i++){ 
        for (int j=gc.j_flc_bl[0][0];j<=gc.j_flc_bl[0][1];j++){
            gc.jdgBnd_Ep[i][j] = 1;
        }
    }

    for (int i=gx.i_flx_bl[0][0];i<=gx.i_flx_bl[0][1];i++){ 
        for (int j=gx.j_flx_bl[0][0];j<=gx.j_flx_bl[0][1];j++){
            gx.jdgBnd_Ex[i][j] = 1;
        }
    }

    for (int i=gr.i_flr_bl[0][0];i<=gr.i_flr_bl[0][1];i++){ 
        for (int j=gr.j_flr_bl[0][0];j<=gr.j_flr_bl[0][1];j++){
            gr.jdgBnd_Er[i][j] = 1;
        }
    }

    //Output boundary check file
    if(pm.flag_chk == 1){

        std::ofstream outputfile1("results/boundary0.csv");
        //outputfile1<<"i,j,x,r,jdgBnd_Ep,jdgBnd_Ex,jdgBnd_Er,jdgBnd_flc,jdgBnd_flx,jdgBnd_flr,kx,kr,kp,zero" << std::endl;
        outputfile1<<"i,j,x,r,jdgBnd_Ep,jdgBnd_Ex,jdgBnd_Er,zero" << std::endl;

        for(int i=0;i<=pm.ni+1;i++){
            for(int j=0;j<=pm.nj+1;j++){
                outputfile1<< i << ","<< j << "," << gc.x[i]<< ","<< gc.r[j]
                    << "," << gc.jdgBnd_Ep[i][j] << "," << gx.jdgBnd_Ex[i][j]<< "," << gr.jdgBnd_Er[i][j]
                    //<< "," << gx.kx[i][j]<< "," << gr.kr[i][j]<< "," << gc.kp[i][j]
                    << "," << 0.0<< std::endl;
            }
        }

        outputfile1.close();
    }
    

    //**************** 境界条件保持配列作成 ****************
    //********* Ex境界条件 (左) *********
    //z6
    for (int j=gx.j_flx_bl[0][0];j<=gx.j_flx_bl[0][1];j++){
        int i=gx.i_flx_bl[0][0];
        mb.iBndWx.push_back(i);
        mb.jBndWx.push_back(j);
        mb.sBndWx.push_back(0); //ディリクレ
    }
    //出力
    if(pm.flag_chk == 1){
        std::ofstream outputfileWx("results/BCx0.csv");
        outputfileWx <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<mb.iBndWx.size();k++){
            int i = mb.iBndWx[k];
            int j = mb.jBndWx[k];
            
            outputfileWx << i << ","<< j 
                << ","<< gc.x[i] << ","<<gc.r[j] << ","<<mb.sBndWx[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileWx.close();
    }
    //********* Ex境界条件 (右) *********
    //z5
    for (int j=gx.j_flx_bl[0][0];j<=gx.j_flx_bl[0][1];j++){
        int i=gx.i_flx_bl[0][1];
        mb.jBndEx.push_back(j);
        mb.iBndEx.push_back(i);
        mb.sBndEx.push_back(2); //ガウス
    }
    //出力
    if(pm.flag_chk == 1){
        std::ofstream outputfileEx("results/BCx1.csv");
        outputfileEx <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<mb.iBndEx.size();k++){
            int i = mb.iBndEx[k];
            int j = mb.jBndEx[k];
            
            outputfileEx << i << ","<< j 
                << ","<< gc.x[i] << ","<<gc.r[j] << ","<<mb.sBndEx[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileEx.close();
    }
    //********* Ex境界条件 (下) *********
    //x2
    for (int i=gx.i_flx_bl[0][0];i<=gx.i_flx_bl[0][1];i++){
        int j=gx.j_flx_bl[0][0];
        mb.iBndSx.push_back(i);
        mb.jBndSx.push_back(j);
        mb.sBndSx.push_back(0); //ディリクレ
    }
    //出力
    if(pm.flag_chk == 1){
        std::ofstream outputfileSx("results/BCx2.csv");
        outputfileSx <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<mb.iBndSx.size();k++){
            int i = mb.iBndSx[k];
            int j = mb.jBndSx[k];
            
            outputfileSx << i << ","<< j 
                << ","<< gc.x[i] << ","<<gc.r[j] << ","<<mb.sBndSx[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileSx.close();
    }
    //********* Ex境界条件 (上) *********
    //x3
    for (int i=gx.i_flx_bl[0][0];i<=gx.i_flx_bl[0][1];i++){
        int j=gx.j_flx_bl[0][1];
        mb.iBndNx.push_back(i);
        mb.jBndNx.push_back(j);
        mb.sBndNx.push_back(0); //ディリクレ
    }
    //出力
    if(pm.flag_chk == 1){
        std::ofstream outputfileNx("results/BCx3.csv");
        outputfileNx <<"i,j,x,r,BC,zero" << std::endl;
        for (int k=0;k<mb.iBndNx.size();k++){
            int i = mb.iBndNx[k];
            int j = mb.jBndNx[k];
            
            outputfileNx << i << ","<< j 
                << ","<< gc.x[i] << ","<<gc.r[j] << ","<<mb.sBndNx[k] 
                << ","<<0.0 << std::endl;
        }
        outputfileNx.close();
    }

    //Er
    //=============================================
    //********* Er境界条件 (左) *********
    for (int j=gr.j_flr_bl[0][0];j<=gr.j_flr_bl[0][1];j++){
        int i=gr.i_flr_bl[0][0];
        mb.iBndWr.push_back(i);
        mb.jBndWr.push_back(j);
        //mb.sBndWr.push_back(1); //開放
        mb.sBndWr.push_back(7); //励振
    }
    //出力
    if(pm.flag_chk == 1){
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
    for (int j=gr.j_flr_bl[0][0];j<=gr.j_flr_bl[0][1];j++){
        int i=gr.i_flr_bl[0][1];
        mb.iBndEr.push_back(i);
        mb.jBndEr.push_back(j);
        mb.sBndEr.push_back(1); //開放
    }
    //出力
    if(pm.flag_chk == 1){
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
    for (int i=gr.i_flr_bl[0][0];i<=gr.i_flr_bl[0][1];i++){
        int j=gr.j_flr_bl[0][0];
        mb.iBndSr.push_back(i);
        mb.jBndSr.push_back(j);
        mb.sBndSr.push_back(2); //ガウス
    }
    //出力
    if(pm.flag_chk == 1){
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
    for (int i=gr.i_flr_bl[0][0];i<=gr.i_flr_bl[0][1];i++){
        int j=gr.j_flr_bl[0][1];
        mb.iBndNr.push_back(i);
        mb.jBndNr.push_back(j);
        mb.sBndNr.push_back(2); //ガウス
    }
    //出力
    if(pm.flag_chk == 1){
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

    //本番用はmake profileだが
    for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[0][1];i++){ 
        for (int j=gc.j_flc_bl[0][0];j<=gc.j_flc_bl[0][1];j++){
            gc.epsr[i][j] = pm.epsr_diele;
        }
    }

}

//*****************************************************************
//**                                                             **
//**           void makeBoundary                                 **
//**                                                             **
//*****************************************************************
void makeBoundary_impedanceTest_old(Params &pm, GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr, MicrowaveBC &mb)
{
    double x_tmp = 0.0;
    double r_tmp = 0.0;
    
    pm.x1 = 0.0e-3;
    pm.x2 = pm.xR;
   
    pm.r1 =  0.8e-3;
    pm.r2 =  2.0e-3;

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

    //Output boundary check file
    if(pm.flag_chk == 1){

        std::ofstream outputfile1("results/boundary0.csv");
        //outputfile1<<"i,j,x,r,jdgBnd_Ep,jdgBnd_Ex,jdgBnd_Er,jdgBnd_flc,jdgBnd_flx,jdgBnd_flr,kx,kr,kp,zero" << std::endl;
        outputfile1<<"i,j,x,r,jdgBnd_Ep,jdgBnd_Ex,jdgBnd_Er,jdgBnd_flc,jdgBnd_flx,jdgBnd_flr,zero" << std::endl;

        for(int i=0;i<=pm.ni+1;i++){
            for(int j=0;j<=pm.nj+1;j++){
                outputfile1<< i << ","<< j << "," << gc.x[i]<< ","<< gc.r[j]
                    << "," << gc.jdgBnd_Ep[i][j] << "," << gx.jdgBnd_Ex[i][j]<< "," << gr.jdgBnd_Er[i][j]
                    << "," << gc.jdgBnd_flc[i][j]<< "," << gx.jdgBnd_flx[i][j]<< "," << gr.jdgBnd_flr[i][j]
                    //<< "," << gx.kx[i][j]<< "," << gr.kr[i][j]<< "," << gc.kp[i][j]
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
        //mb.sBndWr.push_back(1); //開放
        mb.sBndWr.push_back(7); //励振
    }
    //出力
    if(pm.flag_chk == 1){
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
    if(pm.flag_chk == 1){
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
    if(pm.flag_chk == 1){
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
    if(pm.flag_chk == 1){
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

//*****************************************************************
//**                                                             **
//**           void makeBoundary                                 **
//**                                                             **
//*****************************************************************
void InitialModule::makeProfile(Params &pm, GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr)
{

    ////プラズマ密度プロファイル作成
    double rho_max = 4.0e18*1.0; //最大値　(カットオフ密度は7.45E+16)
    double rho_min = 1.0e14*1.0; //最大値　(カットオフ密度は7.45E+16)
    double sigmax_rho = 0.04; //標準偏差 m
    double sigmar_rho = 0.01; //標準偏差 m 0.05
    double xCen_rho = 0.025; //xの中心 0.009 + 0.0078
    double rCen_rho = 0.0; //rの中心 0.007
    
    for (int i=0;i<=pm.ni+1;i++){
        for (int j=0;j<=pm.nj+1;j++){
            double x_tmp = gc.x[i];
            double r_tmp = gc.r[j];
            
            double z_x = (gc.x[i]-xCen_rho)/sigmax_rho;
            double z_r = (gc.r[j]-rCen_rho)/sigmar_rho;
            gc.rate_ionize[i][j] = 6.0e21*exp(-z_x*z_x/0.12-z_r*z_r/0.12)*gc.jdgBnd_flc[i][j];
            //gc.rhoi[i][j] = 1e17*exp(-z_x*z_x/0.12-z_r*z_r/0.12)*gc.jdgBnd_flc[i][j];
            gc.rhoi[i][j] = pm.rhoe_ini*gc.jdgBnd_flc[i][j];
            //Pabs[i][j] = 2e9*exp(-z_x*z_x/0.12-z_r*z_r/0.12); //nablaGのテスト用
            gc.rhom[i][j] = 1e-10*gc.jdgBnd_flc[i][j];
        }
    }

    //プラズマ密度プロファイル作成
    for (int i=0;i<=pm.ni+1;i++){
        for (int j=0;j<=pm.nj+1;j++){
            gc.rhoi_old[i][j] = gc.rhoi[i][j];
            gc.rhoe[i][j] = gc.rhoi[i][j];
            gc.rhoe_old[i][j] = gc.rhoi[i][j];
        }
    }

    //電子温度プロファイル作成
    double Te_max = 4.0*ph::e0/ph::Boltz; //最大値 4.0
    double Te_min = 1.0*ph::e0/ph::Boltz; //最大値 2.0
    double sigmax_Te = 0.04; //標準偏差 m
    double sigmar_Te = 0.01; //標準偏差 m 0.12
    double xCen_Te = 0.025; //xの中心 0.009 + 0.0078
    double rCen_Te = 0.0; //rの中心 0.025 0.06

    for (int i=0;i<=pm.ni+1;i++){
        for (int j=0;j<=pm.nj+1;j++){
            double x_tmp = gc.x[i];
            double r_tmp = gc.r[j];
            
            double z_x = (gc.x[i]-xCen_Te)/sigmax_Te;
            double z_r = (gc.r[j]-rCen_Te)/sigmar_Te;
            
            //gc.Te[i][j] = ((Te_max-Te_min)*exp(-z_x*z_x/0.12-z_r*z_r/0.12) + Te_min)*gc.jdgBnd_flc[i][j];
            gc.Te[i][j] = pm.Te_ini*gc.jdgBnd_flc[i][j];
            //Te[i][j] = (Te_max-Te_min)*pow(fmax(Ap[i][j],0.0)/0.0000466606,0.5) + Te_min;
            gc.rhoeps[i][j] = 3.0/2.0*gc.rhoe[i][j]*ph::Boltz*gc.Te[i][j]*gc.jdgBnd_flc[i][j];
        }
    }

    /*
    //左 壁2
    //std::cout << " Left2 "<< std::endl;
    //std::cout << " i = "  << i_flc_bl[1][0]   << std::endl;
    //std::cout << " j = "  << j_flc_bl[1][0] << " ~ "<< j_flc_bl[0][0]-1 << std::endl;
    for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[0][0]-1;j++){
        int i = gc.i_flc_bl[1][0];
        gc.Te[i-1][j] = gc.Te[i][j];
    }
    
    //下 壁
    //std::cout << " Bottom wall "<< std::endl;
    //std::cout << " j = "  << j_flc_bl[0][0] << std::endl;
    //std::cout << " i = "  << i_flc_bl[0][0] << " ~ "<< i_flc_bl[0][1] << std::endl;
    for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[0][1]-1;i++){ 
        int j = gc.j_flc_bl[0][0];
        gc.Te[i][j-1] = gc.Te[i][j];
    }
    */

    //比誘電率プロファイル作成
    for (int i=gc.i_flc_bl[pm.n_bl-1][0];i<=gc.i_flc_bl[pm.n_bl-1][1];i++){ 
        for (int j=gc.j_flc_bl[pm.n_bl-1][0];j<=gc.j_flc_bl[pm.n_bl-1][1];j++){
            gc.epsr[i][j] = pm.epsr_diele;
        }
    }
    
    //基底中性粒子密度プロファイル作成
    for (int i=0;i<=pm.ni+1;i++){
        for (int j=0;j<=pm.nj+1;j++){
            gc.rhon[i][j] = pm.rhon_ini;
            gc.rhon_old[i][j] = pm.rhon_ini;
        }
    }

    //異常衝突周波数の分布
    for (int i=0;i<=pm.ni+1;i++){
        for (int j=0;j<=pm.nj+1;j++){

            double Bmag = sqrt(gc.Bx[i][j]*gc.Bx[i][j] + gc.Br[i][j]*gc.Br[i][j]); // masse/e0;//
            
            //一定
            gc.nu_ano[i][j] = ph::e0*Bmag/pm.masse/pm.alpha_Bohm;

            /*
            //分布 2つ
            double dL = (x5-x4)/2.0; //遷移幅の片側
            double xCen = (x5+x4)/2.0; //遷移部の中心座標
            double alpha_Bohm_L = 16.0;
            double alpha_Bohm_R = 160;

            
            if(x[i]<= xCen - dL){ 
                nu_ano[i][j] = e0*Bmag/masse/alpha_Bohm_L;
            }else if(x[i]> xCen - dL && x[i]<= xCen + dL){ 
                double nu_ano_L = e0*Bmag/masse/alpha_Bohm_L;
                double nu_ano_R = e0*Bmag/masse/alpha_Bohm_R;

                double t1 = (x[i]-(xCen - dL))/(2.0*dL);
                nu_ano[i][j] = (1.0-t1)*nu_ano_L + t1*nu_ano_R;

            }else if(x[i]> xCen + dL){
                nu_ano[i][j] = e0*Bmag/masse/alpha_Bohm_R;
            }
        
            //分布 3つ
            double alpha_Bohm_L = 16.0;
            double alpha_Bohm_C = 16.0;
            double alpha_Bohm_R = 16.0;

            double dL1 = 1.0e-3; //遷移幅の片側
            double xCen1 = x4+dL1; //遷移部の中心座標
            
            double dL2 = 1.0e-3; //遷移幅の片側
            double xCen2 = x5-dL2; //遷移部の中心座標
            
            if(x[i]<= xCen1 - dL1){ 
                nu_ano[i][j] = e0*Bmag/masse/alpha_Bohm_L;
            }else if(x[i]> xCen1 - dL1 && x[i]<= xCen1 + dL1){ 
                double nu_ano_L = e0*Bmag/masse/alpha_Bohm_L;
                double nu_ano_C = e0*Bmag/masse/alpha_Bohm_C;

                double t1 = (x[i]-(xCen1 - dL1))/(2.0*dL1);
                nu_ano[i][j] = (1.0-t1)*nu_ano_L + t1*nu_ano_C;

            }else if(x[i]> xCen1 + dL1 && x[i]<= xCen2 - dL2){
                nu_ano[i][j] = e0*Bmag/masse/alpha_Bohm_C;
            }else if(x[i]> xCen2 - dL2 && x[i]<= xCen2 + dL2){ 
                double nu_ano_C = e0*Bmag/masse/alpha_Bohm_C;
                double nu_ano_R = e0*Bmag/masse/alpha_Bohm_R;

                double t2 = (x[i]-(xCen2 - dL2))/(2.0*dL2);
                nu_ano[i][j] = (1.0-t2)*nu_ano_C + t2*nu_ano_R;

            }else if(x[i]> xCen2 + dL2){
                nu_ano[i][j] = e0*Bmag/masse/alpha_Bohm_R;
            }     
            */

        }
    }
}


//===========================================================================
//                           private functions
//===========================================================================

//*****************************************************************
//**                                                             **
//**           CmpMat2 sParamToABCD()                            **
//**                                                             **
//*****************************************************************
void InitialModule::sParamToABCD(
    CmpMat2& T,
    const complex<double>& S11,
    const complex<double>& S21,
    const complex<double>& S12,
    const complex<double>& S22,
    const double Z0
){
    
    if (abs(S21) < 1.0e-14) {
        throw runtime_error("sParamToABCD: S21 is too small.");
    }

    const complex<double> one(1.0, 0.0);
    const complex<double> two(2.0, 0.0);
    const complex<double> denom = 2.0*S21;

    T.A = ((1.0 + S11)*(1.0 - S22) + S12*S21)/denom;
    T.B = Z0*((1.0 + S11)*(1.0 + S22) - S12*S21)/denom;
    T.C = 1.0/Z0*((1.0 - S11)*(1.0 - S22) - S12*S21)/denom;
    T.D = ((1.0 - S11)*(1.0 + S22) + S12*S21)/denom;
}

//*****************************************************************
//**                                                             **
//**           void ABCDToSParam()                               **
//**                                                             **
//*****************************************************************
void InitialModule::ABCDToSParam(
    complex<double>& S11,
    complex<double>& S21,
    complex<double>& S12,
    complex<double>& S22,
    const CmpMat2& T,
    const double Z0
){

    const complex<double> denom = T.A + T.B/Z0 + T.C*Z0 + T.D;

    if (abs(denom) < 1.0e-14) {
        throw runtime_error("ABCDToSParam: denominator is too small.");
    }

    S11 = (T.A + T.B/Z0 - T.C*Z0 - T.D)/denom;

    S21 = 2.0/denom;

    S12 = 2.0*(T.A*T.D - T.B*T.C)/denom;

    S22 = (-T.A + T.B/Z0 - T.C*Z0 + T.D)/denom;
}
