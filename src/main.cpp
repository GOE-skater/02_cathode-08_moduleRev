//*****************************************************************
//*****************************************************************
//**                                                             **
//**  Plasma fluid simulation of microwave neutralizer           **
//**          coded by Ryo Shirakawa.                            **
//*****************************************************************
//*****************************************************************

//general
#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <time.h>
#include <vector>
#include <sstream>
#include <mpi.h>

//variables
#include "constants.hpp"
#include "params.hpp"
#include "arrays.hpp"

//functions
#include "miniFuncs.hpp"

//modules
#include "inputModule.hpp"
#include "initialModule.hpp"
#include "outputModule.hpp"
#include "fluidModule.hpp"
#include "emfieldModule.hpp"

using namespace std;

//*****************************************************************
//**                                                             **
//**           int main                                          **
//**                                                             **
//*****************************************************************
int main(int argc, char *argv[])
{

    //input data
    //-------------------------------------
    Params pm;
    BolsigVec bo;
    GridK gk;
    MicrowaveBC mb;
    SeeVec se;
    //-------------------------------------

    //modules
    //-------------------------------------
    InputModule inpM;
    InitialModule iniM;
    OutputModule outM;
    FluidModule fluM;
    EmfieldModule emfM;
    //FieldModule fieldM;
    //ParticleModule pclM;
    //-------------------------------------

    int ndiv_out = 10;
    int icon_end = 0;
    double GnuFactor = 0.1;
    int nGnuOutTime = 1000*GnuFactor;
    int nGnuMaxTimeRange = 100000*GnuFactor; //Gnuplotで粒子の履歴を表示する最大幅
    int nGnuDivTime = fmax(10*GnuFactor,1); //Gnuplotで何ステップごとに出力するか
    
    int nErr = 13; //エラーを表示する数
    vector<vector<double> > error_history(nErr,vector<double>(0)); //エラーの履歴
    vector<vector<double> > current_history(3,vector<double>(0)); //エラーの履歴
    vector<int> itime_history(0,0); //時間ステップの履歴

    clock_t gstart = clock();
    clock_t time1  = clock();
    clock_t time2  = clock();

    FILE* error_plot;
    FILE* current_plot;
    FILE* rhon_plot;
    FILE* Te_plot;
    FILE* rhoi_plot;
    FILE* Uex_plot;

    if(pm.icon_gnuRes == 1){
        error_plot = popen("gnuplot", "w"); 
        fprintf(error_plot, "set title 'Error History'\n");

        current_plot = popen("gnuplot", "w"); 
        fprintf(current_plot, "set title 'Current History (mA)'\n");

        //rhon_plot = popen("gnuplot", "w"); 
        //fprintf(rhon_plot, "set title 'Neutral Density'\n");

        //rhoi_plot = popen("gnuplot", "w"); 
        //fprintf(rhoi_plot, "set title 'Ion Density'\n");

        //Uex_plot = popen("gnuplot", "w"); 
        //fprintf(Uex_plot, "set title 'Electron Velocity'\n");

        //Te_plot = popen("gnuplot", "w"); 
        //fprintf(Te_plot, "set title 'Electron Temperature'\n");

    }

    MPI_Init(&argc, &argv); //hypre用

    //parameter input
    //-------------------------------------
    inpM.inputParam(pm,"setup.yaml");
    //-------------------------------------

    //initialization of arrays
    //-------------------------------------
    GridCenter gc(pm.ni,pm.nj);
    GridInterfaceX gx(pm.ni,pm.nj);
    GridInterfaceR gr(pm.ni,pm.nj);
    //-------------------------------------

    //initialization of parameters
    //-------------------------------------
    iniM.iniParam(pm,gc,gx,gr); 
    //-------------------------------------
    
    
    if(pm.icon_impTest == 1){
        
        iniM.makeBoundary_impedanceTest(pm, gc, gx, gr, gk, mb);
        emfM.solve_Microwave_impedanceTest(pm, gc, gx, gr, gk, mb); //マイクロ波更新
        outM.output_phase(pm, gc, gx, gr);
        //outM.output(pm, gc, gx, gr, bo); 
        return 0;
    }

    iniM.makeBoundary(pm, gc, gx, gr, gk, mb);
    inpM.input_Bfield_data(pm,gc,"Bfield_data.csv"); //input B-field data
    inpM.input_SEE_data(pm, se, "coefEISEE.csv"); //input SEE data
    iniM.makeProfile(pm, gc, gx, gr); //set the initial profile
    inpM.input_restart_data(pm, gc, gx, gr, "restart.csv"); //input restart data
    inpM.input_BOLSIG_data(pm, bo, "rateCoef_e.csv"); //input Bolsig data
    
    fluM.update_transport_coef(pm, gc, gx, gr, bo); //update transport coefficients

    //マイクロ波計算
    {
        emfM.solve_Microwave(pm, gc, gx, gr, gk, mb); //マイクロ波更新
        emfM.update_energy_profile(pm, gc, gx, gr); //電力吸収プロファイル更新
    }

    //ファイルにアウトプット
    outM.output_phase(pm, gc, gx, gr);
    outM.output(pm, gc, gx, gr, bo); 

    ofstream outputfile1("results/residuals.csv");
    outputfile1 << "itime,time,rhoi,Uix,Uir,Uip,phi,rhoe,rhoUex,rhoUer,rhoeps,Gx,Gr,rhom,rhon" << endl;

    do {

       //if(itime%1 == 0) V_bias = fmin(V_bias + 1,100);
        
        fluM.check_CFL(pm, gc, gx, gr); //CFL条件を元にdtを決定する　あるいは　dtからCFLを調べる
        pm.gtime = pm.gtime + pm.dt;
        pm.itime = pm.itime + 1;

        if(pm.itime % ndiv_out == 0){
            time2 = clock();
            double timePerCycle = (double)(time2-time1)/CLOCKS_PER_SEC/ndiv_out;
            double restSec = double(pm.ntime - pm.itime)*timePerCycle;
            int restHour = int(restSec/3600);
            int restMin = int(restSec/60) - restHour*60;
            
            cout << "*********************************************" << endl;
            cout << "itime = " << pm.itime << "/" << pm.ntime << " CFL=" << pm.CFL << " dt=" << pm.dt << " t=" << pm.gtime << endl;
            cout << " The rest of time = "<<restHour<<" h "<<restMin<<" min | Progress = "<<(double)pm.itime/pm.ntime*100 <<" % "<<endl;
            cout <<  endl;

            time1 = clock();

        }
        
        //solve ion
        //------------------------------------
        if(pm.itime % pm.ndt_i == 0){
            //temporary velocity update
            //------------------------------------
            fluM.solve_Uix_tmp(pm, gc, gx, gr);
            fluM.solve_Uir_tmp(pm, gc, gx, gr);
            fluM.solve_Uip_tmp(pm, gc, gx, gr);
            //------------------------------------

            //solve density
            //------------------------------------
            fluM.solve_rhoi_constTe(pm, gc, gx, gr);
            //------------------------------------

            //correct velocity
            //------------------------------------
            fluM.correct_Ui_constTe(pm, gc, gx, gr);
            //------------------------------------
        }
        //------------------------------------


        if(pm.icon_PC == 1){
            
            fluM.solve_phi_couple_wdTe_wSEE_PC(pm, gc, gx, gr, se);

            if(pm.icon_inertia == 0){

                fluM.solve_rhoe_wdTe_wSEE_PC(pm, gc, gx, gr, se);

                fluM.solve_Te_wdTe_wSEE_PC(pm, gc, gx, gr, se);
            
            }else if(pm.icon_inertia == 1){
                //solve_rhoe_wdTe_wSEE_wInertia_PC();
                //solve_Te_wdTe_wSEE_wInertia_PC();
            }
            
        }else{
            //------------------------------------
            // no functions
            //------------------------------------
        }

        if(pm.itime % pm.ndt_n == 0){
            fluM.update_rhon(pm, gc, gx, gr); //Diffusion方程式
            //update_rhon_log(); //Diffusion方程式
        }

        if(pm.itime % pm.ndt_m == 0){
            //metastable更新
            fluM.update_rhom(pm, gc, gx, gr); //論文と同じ実装
            //update_rhom_kinetic(); //境界条件kinetic修正版
        }

        
        //輸送係数更新
        fluM.update_transport_coef(pm, gc, gx, gr, bo); //論文と同じ実装
        //update_transport_coef_mod(); //nu_effの範囲の修正適用

        if(pm.itime % pm.ndiv_MW == 0){
            emfM.solve_Microwave(pm, gc, gx, gr, gk, mb); //マイクロ波更新
            emfM.update_energy_profile(pm, gc, gx, gr); //電力吸収プロファイル更新
        }


        if(pm.itime % ndiv_out == 0){
            cout << "error_rhoi = "    << pm.error_rhoi 
                <<  " error_Uix = "    << pm.error_Uix
                <<  " error_Uir = "    << pm.error_Uir
                <<  " error_Uip = "    << pm.error_Uip
                <<  " error_phi = "    << pm.error_phi
                <<  " error_rhoe = "   << pm.error_rhoe
                <<  " error_rhoUex = " << pm.error_rhoUex
                <<  " error_rhoUer = " << pm.error_rhoUer
                <<  " error_rhoeps = " << pm.error_rhoeps
                <<  " error_Gx = "     << pm.error_Gx
                <<  " error_Gr = "     << pm.error_Gr
                <<  " error_rhom = "   << pm.error_rhom
                <<  " error_rhon = "   << pm.error_rhon
                <<  " error_max = "    << fmax(fmax(fmax(fmax(fmax(fmax(fmax(fmax(fmax(fmax(fmax(fmax(pm.error_rhoi,pm.error_Uix),pm.error_Uir),pm.error_Uip),pm.error_phi),pm.error_rhoe),pm.error_rhoUex),pm.error_rhoUer)
                                        ,pm.error_rhoeps),pm.error_Gx),pm.error_Gr),pm.error_rhom),pm.error_rhon)
                << endl;
        }

        outputfile1 << pm.itime
                <<","<< pm.gtime
                <<","<< pm.error_rhoi 
                <<","<< pm.error_Uix
                <<","<< pm.error_Uir
                <<","<< pm.error_Uip
                <<","<< pm.error_phi
                <<","<< pm.error_rhoe
                <<","<< pm.error_rhoUex
                <<","<< pm.error_rhoUer
                <<","<< pm.error_rhoeps
                <<","<< pm.error_Gx
                <<","<< pm.error_Gr
                <<"," << pm.error_rhom
                <<"," << pm.error_rhon
                << endl;

        //****************** Gnuplot結果出力 History ****************** 
        if(pm.icon_gnuRes == 1){
            
            if(pm.itime % nGnuDivTime ==0){
                vector<double> errorVec(nErr,0.0);
                errorVec[0] = pm.error_rhoi;
                errorVec[1] = pm.error_Uix;
                errorVec[2] = pm.error_Uir;
                errorVec[3] = pm.error_Uip;
                errorVec[4] = pm.error_phi;
                errorVec[5] = pm.error_rhoe;
                errorVec[6] = pm.error_rhoUex;
                errorVec[7] = pm.error_rhoUer;
                errorVec[8] = pm.error_rhoeps;
                errorVec[9] = pm.error_Gx;
                errorVec[10] = pm.error_Gr;
                errorVec[11] = pm.error_rhom;
                errorVec[12] = pm.error_rhon;

                vector<double> currentVec(3,0.0);
                currentVec[0] = pm.Ii_Anode*1000;
                currentVec[1] = pm.Ie_Anode*1000;
                currentVec[2] = pm.I_Anode*1000;
                
                //値の追加
                itime_history.push_back(pm.itime);
                for(int iSp = 0;iSp<nErr;iSp++){
                    error_history[iSp].push_back(errorVec[iSp]);
                }
                for(int iSp = 0;iSp<3;iSp++){
                    current_history[iSp].push_back(currentVec[iSp]);
                }

                //値の削除
                if(pm.itime > nGnuMaxTimeRange){
                    itime_history.erase(itime_history.begin());
                    for(int iSp = 0;iSp<nErr;iSp++){
                        error_history[iSp].erase(error_history[iSp].begin());
                    }
                    for(int iSp = 0;iSp<3;iSp++){
                        current_history[iSp].erase(current_history[iSp].begin());
                    }
                }

                if(fmod(pm.itime,nGnuOutTime) == 0){
                    plotHistory(error_plot,error_history,itime_history,nGnuMaxTimeRange);
                }

                if(fmod(pm.itime,nGnuOutTime) == 0){
                    plotCurrentHistory(current_plot,current_history,itime_history,nGnuMaxTimeRange);
                }
            }
        }

        if(pm.itime % pm.ndiv_fout == 0){
            outM.output(pm, gc, gx, gr, bo);
        }

        if(pm.itime >= pm.ntime){
            emfM.solve_Microwave(pm, gc, gx, gr, gk, mb);
            emfM.update_energy_profile(pm, gc, gx, gr);
            emfM.solve_Microwave(pm, gc, gx, gr, gk, mb);
            outM.output(pm, gc, gx, gr, bo);

            oneMoreTime:
            int tmp;

            if(pm.icon_autoFinish == 0){
                cout << "[Terminal Outpput] Continue Calculation? Yes = 1, No = 0"<<endl;
                cin >> tmp;
            }else{
                tmp = 0;
            }
            
            if(tmp == 0){
                int input1;
                if(pm.icon_autoFinish == 0){
                    cout << "[Terminal Outpput] Really? Yes = 1, No = 0"<<endl;
                    cin >> input1;
                }else{
                    icon_end = 1;
                }

                if(input1 == 1){
                    icon_end = 1;
                }else{
                    goto oneMoreTime;
                }

            }else{
                int tmp1;
                int input;
                cout << "[Terminal Outpput] Next Max Steps = "<<endl;
                cin >> pm.ntime;
                cout << "[Terminal Outpput] Same CFL? Yes = 1, No = 0"<<endl;
                cin >> tmp1;

                if(tmp1 == 0){
                    if(pm.icon_adp_dt == 1){
                        cout << "[Terminal Outpput] CFL = "<<endl;
                        cin >> pm.CFL;
                    }else{
                        cout << "[Terminal Outpput] dt = "<<endl;
                        cin >> pm.dt;
                    }
                }

                if(pm.icon_adp_dt == 1) cout << "[Terminal Outpput] Next Max Steps = "<< pm.ntime << " CFL = " << pm.CFL << endl;
                if(pm.icon_adp_dt == 0) cout << "[Terminal Outpput] Next Max Steps = "<< pm.ntime << " dt = " << pm.dt << endl;
                cout << "[Terminal Outpput] Is it OK? Yes = 1, No = 0"<< endl;
                cin >> input;
                if(input == 1){
                    icon_end = 0;
                }else{
                    goto oneMoreTime;
                }
            }
        }

        if(pm.itime==pm.ntime){
            icon_end = 1;
        }

    } while (icon_end == 0);

    outM.output_phase(pm, gc, gx, gr);
    outputfile1.close();

    clock_t gend = clock();
    double time_real = (double)(gstart-gend)/CLOCKS_PER_SEC;
    cout << "elapsed time = " << time_real << " sec" << endl;

    MPI_Finalize();


}