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

