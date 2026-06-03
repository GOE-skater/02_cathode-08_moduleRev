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
#include <mpi.h>

//*****************************************************************
//**                                                             **
//**           int main                                          **
//**                                                             **
//*****************************************************************
int main(int argc, char *argv[])
{
    int ndiv_out = 10;
    int icon_end = 0;
    double GnuFactor = 0.1;
    int nGnuOutTime = 1000*GnuFactor;
    int nGnuMaxTimeRange = 100000*GnuFactor; //Gnuplotで粒子の履歴を表示する最大幅
    int nGnuDivTime = fmax(10*GnuFactor,1); //Gnuplotで何ステップごとに出力するか
    
    int nErr = 13; //エラーを表示する数
    std::vector<std::vector<double> > error_history(nErr,std::vector<double>(0)); //エラーの履歴
    std::vector<std::vector<double> > current_history(3,std::vector<double>(0)); //エラーの履歴
    std::vector<int> itime_history(0,0); //時間ステップの履歴

    clock_t gstart = clock();
    clock_t time1  = clock();
    clock_t time2  = clock();

    FILE* error_plot;
    FILE* current_plot;
    FILE* rhon_plot;
    FILE* Te_plot;
    FILE* rhoi_plot;
    FILE* Uex_plot;

    if(icon_gnuRes == 1){
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

    //初期化
    inputParam(); //setup.csvからインプットパラメータを読み込み
    iniparameter(); //パラメータ初期化
    
    if(icon_impTest == 1){
        
        makeBoundary_impedanceTest(); //形状の生成
    
        solve_Microwave_impedanceTest(); //マイクロ波更新
        output_phase();
        //output(); //ファイルにアウトプット
        return 0;
    }

    makeBoundary(); //形状の生成
    input_Bfield_data(); //磁場の読み込み
    input_SEE_data();
    makeProfile(); //事前定義プロファイルの作成
    input_restart_data(); //リスタートデータの読み込み
    input_BOLSIG_data(); //レート係数の読み込み
    update_transport_coef(); //輸送係数更新

    //マイクロ波計算
    {
        solve_Microwave(); //マイクロ波更新
        //update_energy_profile(); //電力吸収プロファイル更新
    }
    

    /*
    //初期の準安定種のプロファイルを生成
    if(icon_restart != 1){
        update_rhom_diff();
    }
    */
    output_phase();
    output(); //ファイルにアウトプット
    return 0;

    std::ofstream outputfile1("results/residuals.csv");
    outputfile1 << "itime,time,rhoi,Uix,Uir,Uip,phi,rhoe,rhoUex,rhoUer,rhoeps,Gx,Gr,rhom,rhon" << std::endl;


    do {

       //if(itime%1 == 0) V_bias = fmin(V_bias + 1,100);
        
        check_CFL(); //CFL条件を元にdtを決定する　あるいは　dtからCFLを調べる
        gtime = gtime+dt;
        itime=itime+1;

        if(itime % ndiv_out == 0){
            time2 = clock();
            double timePerCycle = (double)(time2-time1)/CLOCKS_PER_SEC/ndiv_out;
            double restSec = double(ntime - itime)*timePerCycle;
            int restHour = int(restSec/3600);
            int restMin = int(restSec/60) - restHour*60;
            
            std::cout << "*********************************************" << std::endl;
            std::cout << "itime = " << itime << "/" << ntime << " CFL=" << CFL << " dt=" << dt << " t=" << gtime << std::endl;
            std::cout << " The rest of time = "<<restHour<<" h "<<restMin<<" min | Progress = "<<(double)itime/ntime*100<<" % "<<std::endl;
            std::cout <<  std::endl;

            time1 = clock();

        }
        
        
        if(itime % ndt_i == 0){
            //イオン密度・速度更新
            solve_Uix_tmp(); //仮のUix計算
            solve_Uir_tmp(); //仮のUir計算
            solve_Uip_tmp(); //仮のUip計算
            
            //solve_rhoi(); //イオン密度計算
            solve_rhoi_constTe(); //イオン密度計算
            //solve_rhoi_my(); //イオン密度計算
            //update_Ui(); //速度場更新 Ui
            update_Ui_constTe(); //速度場更新 Ui
            //update_Ui_my(); //速度場更新 Ui
        }

        if(icon_PC == 0){
            //電場更新
            solve_phi_couple();
            //solve_phi_couple_PC();
            //solve_phi_couple_wdTe_PC();
            //solve_phi_couple_SG_PC();
            //solve_phi_couple_directional();
            //solve_phi_couple_dieleOpen();

            //電子密度更新
            //solve_rhoe_direct(); //論文と同じ実装
            //solve_rhoe_wdTe_direct();
            solve_rhoe_wdTe_wSEE_direct();
            //solve_rhoe_direct_dieleOpen();
            //solve_rhoe_PC(); //PC method
            //solve_rhoe_wdTe_PC();
            //solve_rhoe_TV(); //transverse flux method
            
            //電子温度更新
            //solve_Te_direct(); //論文と同じ実装
            //solve_Te_wdTe_direct(); //論文と同じ実装
            solve_Te_wdTe_wSEE_direct(); //論文と同じ実装
            //solve_Te_PC(); //PC method
            //solve_Te_wdTe_PC(); //PC method
        }else{
            solve_phi_couple_wdTe_wSEE_PC();

            if(icon_inertia == 0){
                //solve_rhoe_PC(); //PC method
                solve_rhoe_wdTe_wSEE_PC();

                //solve_Te_PC(); //PC method
                solve_Te_wdTe_wSEE_PC(); //PC method
            }else if(icon_inertia == 1){
                //solve_rhoe_PC(); //PC method
                solve_rhoe_wdTe_wSEE_wInertia_PC();

                //solve_Te_PC(); //PC method
                solve_Te_wdTe_wSEE_wInertia_PC(); //PC method
            }
        }

        //solve_rhoe_Te_PC(); //rho+Te PC method

        //基底の更新
        //reconstAll_for_neut(); //空間再構築
        //riemann_FVS_for_neut(); //リーマンソルバー
        //update_neut(); //更新

        if(itime % ndt_n == 0){
            update_rhon(); //Diffusion方程式
            //update_rhon_log(); //Diffusion方程式
        }

        if(itime % ndt_m == 0){
            //metastable更新
            update_rhom(); //論文と同じ実装
            //update_rhom_kinetic(); //境界条件kinetic修正版
        }
        
        //輸送係数更新
        update_transport_coef(); //論文と同じ実装
        //update_transport_coef_mod(); //nu_effの範囲の修正適用

        if(itime % ndiv_MW == 0){
            solve_Microwave(); //マイクロ波更新
            update_energy_profile(); //電力吸収プロファイル更新
        }


        if(itime % ndiv_out == 0){
            std::cout << "error_rhoi = " << error_rhoi 
                <<  " error_Uix = " << error_Uix
                <<  " error_Uir = " << error_Uir
                <<  " error_Uip = " << error_Uip
                <<  " error_phi = " << error_phi
                <<  " error_rhoe = " << error_rhoe
                <<  " error_rhoUex = " << error_rhoUex
                <<  " error_rhoUer = " << error_rhoUer
                <<  " error_rhoeps = " << error_rhoeps
                <<  " error_Gx = " << error_Gx
                <<  " error_Gr = " << error_Gr
                <<  " error_rhom = " << error_rhom
                <<  " error_rhon = " << error_rhon
                <<  " error_max = " << fmax(fmax(fmax(fmax(fmax(fmax(fmax(fmax(fmax(fmax(fmax(fmax(error_rhoi,error_Uix),error_Uir),error_Uip),error_phi),error_rhoe),error_rhoUex),error_rhoUer)
                                        ,error_rhoeps),error_Gx),error_Gr),error_rhom),error_rhon)
                << std::endl;
        }

        outputfile1 << itime
                <<","<< gtime
                <<","<< error_rhoi 
                <<","<< error_Uix
                <<","<< error_Uir
                <<","<< error_Uip
                <<","<< error_phi
                <<","<< error_rhoe
                <<","<< error_rhoUex
                <<","<< error_rhoUer
                <<","<< error_rhoeps
                <<","<< error_Gx
                <<","<< error_Gr
                <<"," << error_rhom
                <<"," << error_rhon
                << std::endl;

        //****************** Gnuplot結果出力 History ****************** 
        if(icon_gnuRes == 1){
            
            if(itime % nGnuDivTime ==0){
                std::vector<double> errorVec(nErr,0.0);
                errorVec[0] = error_rhoi;
                errorVec[1] = error_Uix;
                errorVec[2] = error_Uir;
                errorVec[3] = error_Uip;
                errorVec[4] = error_phi;
                errorVec[5] = error_rhoe;
                errorVec[6] = error_rhoUex;
                errorVec[7] = error_rhoUer;
                errorVec[8] = error_rhoeps;
                errorVec[9] = error_Gx;
                errorVec[10] = error_Gr;
                errorVec[11] = error_rhom;
                errorVec[12] = error_rhon;

                std::vector<double> currentVec(3,0.0);
                currentVec[0] = Ii_Anode*1000;
                currentVec[1] = Ie_Anode*1000;
                currentVec[2] = I_Anode*1000;
                
                //値の追加
                itime_history.push_back(itime);
                for(int iSp = 0;iSp<nErr;iSp++){
                    error_history[iSp].push_back(errorVec[iSp]);
                }
                for(int iSp = 0;iSp<3;iSp++){
                    current_history[iSp].push_back(currentVec[iSp]);
                }

                //値の削除
                if(itime > nGnuMaxTimeRange){
                    itime_history.erase(itime_history.begin());
                    for(int iSp = 0;iSp<nErr;iSp++){
                        error_history[iSp].erase(error_history[iSp].begin());
                    }
                    for(int iSp = 0;iSp<3;iSp++){
                        current_history[iSp].erase(current_history[iSp].begin());
                    }
                }

                if(fmod(itime,nGnuOutTime) == 0){
                    plotHistory(error_plot,error_history,itime_history,nGnuMaxTimeRange);
                }

                if(fmod(itime,nGnuOutTime) == 0){
                    plotCurrentHistory(current_plot,current_history,itime_history,nGnuMaxTimeRange);
                }
            }
        }

        if(itime % ndiv_fout == 0){
            output();
        }

        if(itime >= ntime){
            solve_Microwave();
            update_energy_profile();
            solve_Microwave();
            output();

            oneMoreTime:
            int tmp;

            if(icon_autoFinish == 0){
                std::cout << "[Terminal Outpput] Continue Calculation? Yes = 1, No = 0"<<std::endl;
                std::cin >> tmp;
            }else{
                tmp = 0;
            }
            
            if(tmp == 0){
                int input1;
                if(icon_autoFinish == 0){
                    std::cout << "[Terminal Outpput] Really? Yes = 1, No = 0"<<std::endl;
                    std::cin >> input1;
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
                std::cout << "[Terminal Outpput] Next Max Steps = "<<std::endl;
                std::cin >> ntime;
                std::cout << "[Terminal Outpput] Same CFL? Yes = 1, No = 0"<<std::endl;
                std::cin >> tmp1;

                if(tmp1 == 0){
                    if(icon_adp_dt == 1){
                        std::cout << "[Terminal Outpput] CFL = "<<std::endl;
                        std::cin >> CFL;
                    }else{
                        std::cout << "[Terminal Outpput] dt = "<<std::endl;
                        std::cin >> dt;
                    }
                }

                if(icon_adp_dt == 1) std::cout << "[Terminal Outpput] Next Max Steps = "<< ntime << " CFL = " << CFL << std::endl;
                if(icon_adp_dt == 0) std::cout << "[Terminal Outpput] Next Max Steps = "<< ntime << " dt = " << dt << std::endl;
                std::cout << "[Terminal Outpput] Is it OK? Yes = 1, No = 0"<< std::endl;
                std::cin >> input;
                if(input == 1){
                    icon_end = 0;
                }else{
                    goto oneMoreTime;
                }
            }
        }

        if(itime==ntime){
            icon_end = 1;
        }

    } while (icon_end == 0);

    output_phase();
    outputfile1.close();

    clock_t gend = clock();
    double time_real = (double)(gstart-gend)/CLOCKS_PER_SEC;
    std::cout << "elapsed time = " << time_real << " sec" << std::endl;

    MPI_Finalize();
}