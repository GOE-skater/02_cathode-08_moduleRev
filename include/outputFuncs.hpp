
#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <time.h>
#include <vector>
#include <sstream>

#include "constants.hpp"
#include "params.hpp"

//*****************************************************************
//**                                                             **
//**           void plot                                         **
//**                                                             **
//*****************************************************************
void plotCurrentHistory_old(FILE* gnuplot_name,std::vector<double> value,std::vector<int> itime_history, int nGnuMaxTimeRange)
{

    std::string graphLegend = "plot ";
    fprintf(gnuplot_name,"%s\n",graphLegend.c_str());
    
    graphLegend = graphLegend+"'-' with lines title '";
    graphLegend = graphLegend+"ion";
    graphLegend = graphLegend+"'";

    fprintf(gnuplot_name,"%s\n",graphLegend.c_str());

    //std::cout << value.size() << ","<<itime_history.size() << std::endl;
    for(int n=0;n<itime_history.size();n++){
        //std::cout << itime_history[n] << ","<<double(value[iSp][n]) << std::endl;
        fprintf(gnuplot_name, "%d, %lf\n",itime_history[n],value[n]);
    }
    fprintf(gnuplot_name, "e\n"); 
    fflush(gnuplot_name);

}

//*****************************************************************
//**                                                             **
//**           void plot                                         **
//**                                                             **
//*****************************************************************
void plotCurrentHistory(FILE* gnuplot_name,std::vector<std::vector<double> > value,std::vector<int> itime_history, int nGnuMaxTimeRange)
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
    std::string graphLegend = "plot ";
    
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
        //std::cout << value.size() << ","<<itime_history.size() << std::endl;
        for(int n=0;n<itime_history.size();n++){
            //std::cout << itime_history[n] << ","<<double(value[iSp][n]) << std::endl;
            fprintf(gnuplot_name, "%d, %lf\n",itime_history[n],value[iSp][n]);
        }
        fprintf(gnuplot_name, "e\n"); 
        fflush(gnuplot_name);
    }
}