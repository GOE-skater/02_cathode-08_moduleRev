#pragma once
#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <time.h>
#include <vector>
#include <sstream>

#include "constants.hpp"
#include "params.hpp"

using namespace std;

//*****************************************************************
//**                                                             **
//**                          split                              **
//**                                                             **
//*****************************************************************
inline vector<string> split(string& input, char delimiter){
    istringstream stream(input);
    string field;
    vector<string> result;
    while (getline(stream, field, delimiter)) {
        result.push_back(field);
    }
    return result;
}

//*****************************************************************
//**                                                             **
//**           double func_Te                                    **
//**                                                             **
//*****************************************************************
inline double func_Te(vector<double> value,double Te_tmp,double dTe_tmp)
{
    int iL = int(Te_tmp/dTe_tmp);
    int iR = int(Te_tmp/dTe_tmp)+1;
    double Te_L = iL*dTe_tmp;

    if(iR < value.size()-1){
        return (value[iR] - value[iL])/(dTe_tmp+1e-100)*(Te_tmp-Te_L) + value[iL];
    }else{
        return value[value.size()-1];
    }
}

//*****************************************************************
//**                                                             **
//**           double func_Te                                    **
//**                                                             **
//*****************************************************************
inline double func_dTe(vector<double> value,double Te_tmp,double dTe_tmp)
{
    int iL = int(Te_tmp/dTe_tmp);
    int iR = int(Te_tmp/dTe_tmp)+1;

    if(iR >= value.size()-1){
        iL = value.size()-2;
        iR = value.size()-1;
    }

    double Te_L = iL*dTe_tmp;
    double Te_R = iR*dTe_tmp;

    return (value[iR] - value[iL])/(Te_R-Te_L+1e-100);
}


//*****************************************************************
//**                                                             **
//**           double func_Ei                                    **
//**                                                             **
//*****************************************************************
inline double func_Ei(vector<double> value,double Ei_tmp,double dEi_tmp, double massi)
{
    int iL = int(Ei_tmp/dEi_tmp);
    int iR = int(Ei_tmp/dEi_tmp)+1;
    double Ei_L = iL*dEi_tmp;

    if(iR < value.size()-1){
        return (value[iR] - value[iL])/(dEi_tmp+1e-100)*(Ei_tmp - Ei_L) + value[iL];
    }else{
        double imax = value.size()-1;
        double v_max = sqrt(imax*dEi_tmp/massi);
        double v_tmp = sqrt(Ei_tmp/massi);
        double sigma_max = value[value.size()-1]/v_max;
        //cout << "Emax = " << imax*dEi/e0 << ", sigma_tmp = " << sigma_tmp <<endl;
        return sigma_max*v_tmp; //有効衝突断面積を一定に保ちながら速度の1乗に比例させる
    }
}


//*****************************************************************
//**                                                             **
//**           double truncate                                   **
//**                                                             **
//*****************************************************************
inline double truncate(double c){

    if(fabs(c) <1e-100){
        return 0.0;
    }else{
        return c;
    }

}

//*****************************************************************
//**                                                             **
//**           void plot                                         **
//**                                                             **
//*****************************************************************
inline void plotHistory(FILE* gnuplot_name,vector<vector<double> > value,vector<int> itime_history, int nGnuMaxTimeRange)
{
    //errorVec[0] = error_rhoi;
    //errorVec[1] = error_Uix;
    //errorVec[2] = error_Uir;
    //errorVec[3] = error_Uip;
    //errorVec[4] = error_phi;
    //errorVec[5] = error_rhoe;
    //errorVec[6] = error_rhoUex;
    //errorVec[7] = error_rhoUer;
    //errorVec[8] = error_rhoeps;
    //errorVec[9] = error_Gx;
    //errorVec[10] = error_Gr;
    //errorVec[11] = error_rhom;

    fprintf(gnuplot_name, "set logscale y\n");
    string graphLegend = "plot ";
    
    graphLegend = graphLegend+"'-' with lines title '";
    graphLegend = graphLegend+"rhoi";
    graphLegend = graphLegend+"'";
    
    graphLegend = graphLegend+", ";

    graphLegend = graphLegend+"'-' with lines title '";
    graphLegend = graphLegend+"Uix";
    graphLegend = graphLegend+"'";

    graphLegend = graphLegend+", ";

    graphLegend = graphLegend+"'-' with lines title '";
    graphLegend = graphLegend+"Uir";
    graphLegend = graphLegend+"'";

    graphLegend = graphLegend+", ";

    graphLegend = graphLegend+"'-' with lines title '";
    graphLegend = graphLegend+"Uip";
    graphLegend = graphLegend+"'";

    graphLegend = graphLegend+", ";

    graphLegend = graphLegend+"'-' with lines title '";
    graphLegend = graphLegend+"phi";
    graphLegend = graphLegend+"'";

    graphLegend = graphLegend+", ";

    graphLegend = graphLegend+"'-' with lines title '";
    graphLegend = graphLegend+"nUex";
    graphLegend = graphLegend+"'";

    graphLegend = graphLegend+", ";

    graphLegend = graphLegend+"'-' with lines title '";
    graphLegend = graphLegend+"nUer";
    graphLegend = graphLegend+"'";

    graphLegend = graphLegend+", ";

    graphLegend = graphLegend+"'-' with lines title '";
    graphLegend = graphLegend+"rhoe";
    graphLegend = graphLegend+"'";

    graphLegend = graphLegend+", ";

    graphLegend = graphLegend+"'-' with lines title '";
    graphLegend = graphLegend+"rhoeUex";
    graphLegend = graphLegend+"'";

    graphLegend = graphLegend+", ";

    graphLegend = graphLegend+"'-' with lines title '";
    graphLegend = graphLegend+"rhoeUer";
    graphLegend = graphLegend+"'";

    graphLegend = graphLegend+", ";

    graphLegend = graphLegend+"'-' with lines title '";
    graphLegend = graphLegend+"rhoeps";
    graphLegend = graphLegend+"'";

    graphLegend = graphLegend+", ";

    graphLegend = graphLegend+"'-' with lines title '";
    graphLegend = graphLegend+"Gx";
    graphLegend = graphLegend+"'";

    graphLegend = graphLegend+", ";

    graphLegend = graphLegend+"'-' with lines title '";
    graphLegend = graphLegend+"Gr";
    graphLegend = graphLegend+"'";

    graphLegend = graphLegend+", ";

    graphLegend = graphLegend+"'-' with lines title '";
    graphLegend = graphLegend+"rhom";
    graphLegend = graphLegend+"'";

    graphLegend = graphLegend+", ";

    graphLegend = graphLegend+"'-' with lines title '";
    graphLegend = graphLegend+"rhon";
    graphLegend = graphLegend+"'";
    

    fprintf(gnuplot_name,"%s\n",graphLegend.c_str());

    for(int iSp = 0;iSp<value.size();iSp++){
        //cout << value.size() << ","<<itime_history.size() << endl;
        for(int n=0;n<itime_history.size();n++){
            //cout << itime_history[n] << ","<<double(value[iSp][n]) << endl;
            fprintf(gnuplot_name, "%d, %lf\n",itime_history[n],value[iSp][n]);
        }
        fprintf(gnuplot_name, "e\n"); 
        fflush(gnuplot_name);
    }
}

//*****************************************************************
//**                                                             **
//**           void plot                                         **
//**                                                             **
//*****************************************************************
inline void plotCurrentHistory(FILE* gnuplot_name,vector<vector<double> > value,vector<int> itime_history, int nGnuMaxTimeRange)
{
    //errorVec[0] = error_rhoi;
    //errorVec[1] = error_Uix;
    //errorVec[2] = error_Uir;
    //errorVec[3] = error_Uip;
    //errorVec[4] = error_phi;
    //errorVec[5] = error_rhoe;
    //errorVec[6] = error_rhoUex;
    //errorVec[7] = error_rhoUer;
    //errorVec[8] = error_rhoeps;
    //errorVec[9] = error_Gx;
    //errorVec[10] = error_Gr;
    //errorVec[11] = error_rhom;

    //fprintf(gnuplot_name, "set logscale y\n");
    string graphLegend = "plot ";
    
    graphLegend = graphLegend+"'-' with lines title '";
    graphLegend = graphLegend+"I_{i,Anode}";
    graphLegend = graphLegend+"'";
    
    graphLegend = graphLegend+", ";

    graphLegend = graphLegend+"'-' with lines title '";
    graphLegend = graphLegend+"I_{e,Anode}";
    graphLegend = graphLegend+"'";

    graphLegend = graphLegend+", ";

    graphLegend = graphLegend+"'-' with lines title '";
    graphLegend = graphLegend+"I_{Anode}";
    graphLegend = graphLegend+"'";
    

    fprintf(gnuplot_name,"%s\n",graphLegend.c_str());

    for(int iSp = 0;iSp<value.size();iSp++){
        //cout << value.size() << ","<<itime_history.size() << endl;
        for(int n=0;n<itime_history.size();n++){
            //cout << itime_history[n] << ","<<double(value[iSp][n]) << endl;
            fprintf(gnuplot_name, "%d, %lf\n",itime_history[n],value[iSp][n]);
        }
        fprintf(gnuplot_name, "e\n"); 
        fflush(gnuplot_name);
    }
}

//*****************************************************************
//**                                                             **
//**           void calcRes()                                    **
//**                                                             **
//*****************************************************************
inline void calcRes(double &error, 
        std::vector<std::vector<double> > &val_new, 
        std::vector<std::vector<double> > &val_old,
        std::vector<std::vector<int> > i_fl_bl,
        std::vector<std::vector<int> > j_fl_bl,
        int n_bl,
        double EPS,
        int icon_error //0: cell-center, 1: cell-interface-x, 2: cell-interface-r
    ){

    error = 0.0;

    int ni = val_old.size()-2;
    int nj = val_old[0].size()-2;
    
    //0:Linf_個別正規化 (max)
    //------------------------------------
    if(icon_error == 0){
        for (int iblock=0;iblock<n_bl;iblock++){ 
            for (int i=i_fl_bl[iblock][0];i<=i_fl_bl[iblock][1];i++){ 
                for (int j=j_fl_bl[iblock][0];j<=j_fl_bl[iblock][1];j++){
                    double error_tmp = fabs(val_new[i][j]-val_old[i][j]);
                    if(error_tmp > error){
                        error = error_tmp;
                    }
                }
            }
        }
    //------------------------------------

    //1:L2_個別正規化
    //------------------------------------
    }else if(icon_error == 1){
        int ncount = 0;
        for (int iblock=0;iblock<n_bl;iblock++){ 
            for (int i=i_fl_bl[iblock][0];i<=i_fl_bl[iblock][1];i++){ 
                for (int j=j_fl_bl[iblock][0];j<=j_fl_bl[iblock][1];j++){
                    error += pow(val_new[i][j]-val_old[i][j],2)/(pow(val_old[i][j],2)+EPS);
                    ncount = ncount + 1;
                }
            }
        }
        error = sqrt(error/double(ncount));
    //------------------------------------

    //2:Linf (max)_一括正規化
    //------------------------------------
    }else if(icon_error == 2){
        double norm = 0.0;
        for (int iblock=0;iblock<n_bl;iblock++){ 
            for (int i=i_fl_bl[iblock][0];i<=i_fl_bl[iblock][1];i++){ 
                for (int j=j_fl_bl[iblock][0];j<=j_fl_bl[iblock][1];j++){
                    double error_tmp = fabs(val_new[i][j]-val_old[i][j]);
                    if(error_tmp > error){
                        error = error_tmp;
                    }
                    double norm_tmp = fabs(val_old[i][j]);
                    if(norm_tmp > norm){
                        norm = norm_tmp;
                    }
                }
            }
        }
        error = error/norm;
    //------------------------------------
    
    //3:L2_一括正規化
    //------------------------------------
    }else if(icon_error == 3){
        double norm = 0.0;
        for (int iblock=0;iblock<n_bl;iblock++){ 
            for (int i=i_fl_bl[iblock][0];i<=i_fl_bl[iblock][1];i++){ 
                for (int j=j_fl_bl[iblock][0];j<=j_fl_bl[iblock][1];j++){
                    error += pow(val_new[i][j]-val_old[i][j],2);
                    norm += pow(val_old[i][j],2);
                }
            }
        }
        error = sqrt(error/norm);
    }
    //------------------------------------

    //========================================
}

//*****************************************************************
//**                                                             **
//**           void output_residual_for_PC                       **
//**                                                             **
//*****************************************************************
void output_residual_for_PC(std::string name, double gtime, int itime, double error_val, double error_Fx, double error_Fr, int itime_PC, int ncount){
    
    std::ofstream outputfile1;
    if(itime_PC == 1){
        outputfile1.open("results/residuals_PC_" + name + ".csv");
        outputfile1 << "gtime,itime,ncount,itime_PC,err_val,err_Fx,err_Fr" << std::endl;
    }else{
        outputfile1.open("results/residuals_PC_" + name + ".csv",std::ios::app);
    }
    outputfile1 << gtime
                <<","<< itime
                <<","<< ncount
                <<","<< itime_PC
                <<","<< error_val
                <<","<< error_Fx
                <<","<< error_Fr
                << std::endl;
    
    outputfile1.close();      
}
