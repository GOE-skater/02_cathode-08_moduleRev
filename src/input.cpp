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
#include <complex>
#include <yaml-cpp/yaml.h>
#include <stdexcept>
//*****************************************************************
//**                                                             **
//**           void inputParam                                   **
//**                                                             **
//*****************************************************************
void inputParam()
{
    const std::string filename = "input/setup.yaml";

    YAML::Node setup;

    try {
        setup = YAML::LoadFile(filename);
    }
    catch (const YAML::BadFile&) {
        std::cout << "[" << filename << "] Error! File not found." << std::endl;
        abort();
    }

    int count = 0;

    auto readDouble = [&](double& dst,
                          const std::string& category,
                          const std::string& key) {
        const YAML::Node categoryNode = setup[category];

        if (!categoryNode) {
            throw std::runtime_error("Missing YAML category: " + category);
        }
        if (!categoryNode[key]) {
            throw std::runtime_error("Missing YAML parameter: " + category + "." + key);
        }

        dst = categoryNode[key].as<double>();
        count += 1;
    };

    auto readInt = [&](int& dst,
                       const std::string& category,
                       const std::string& key) {
        const YAML::Node categoryNode = setup[category];

        if (!categoryNode) {
            throw std::runtime_error("Missing YAML category: " + category);
        }
        if (!categoryNode[key]) {
            throw std::runtime_error("Missing YAML parameter: " + category + "." + key);
        }

        dst = categoryNode[key].as<int>();
        count += 1;
    };

#define READ_DOUBLE(category, var) readDouble((var), (category), #var)
#define READ_INT(category, var)    readInt((var), (category), #var)

    try {
        std::string category;

        category = "microwave";
        READ_DOUBLE(category, Pmw);
        READ_DOUBLE(category, omegam);
        READ_DOUBLE(category, nu_eff);
        READ_DOUBLE(category, deltaECR);

        category = "plasma";
        READ_DOUBLE(category, Ti);
        READ_DOUBLE(category, Tn);
        READ_DOUBLE(category, rhon_ini);
        READ_DOUBLE(category, DmN);
        READ_DOUBLE(category, Te_rep_eV);

        category = "material";
        READ_DOUBLE(category, epsr_diele);

        category = "bias";
        READ_DOUBLE(category, V_bias);

        category = "neutral_inlet";
        READ_DOUBLE(category, Q_neutIn_mgs);
        READ_DOUBLE(category, width_neutIn);

        category = "transport_model";
        READ_INT(category, icon_Bohm);
        READ_INT(category, icon_Sagdeev);
        READ_DOUBLE(category, alpha_Bohm);
        READ_DOUBLE(category, scale_inertia);

        category = "see";
        READ_DOUBLE(category, coefIISEE_ts);
        READ_DOUBLE(category, coefMISEE_ts);
        READ_DOUBLE(category, Te_emitSEE_eV);
        READ_DOUBLE(category, ratioEngy_EISEE_rd);

        category = "scheme";
        READ_INT(category, icon_PC);
        READ_INT(category, icon_inertia);
        READ_INT(category, icon_adp_dt);
        READ_DOUBLE(category, dt_ini);
        READ_INT(category, ndt_i);
        READ_INT(category, ndt_m);
        READ_INT(category, ndt_n);
        READ_DOUBLE(category, CFL);

        category = "solver";
        READ_DOUBLE(category, error_cnv_SOR_Ui);
        READ_INT(category, maxITR_SOR_Ui);
        READ_INT(category, icon_iter_Ui);
        READ_DOUBLE(category, error_cnv_SOR_rhoi);
        READ_INT(category, maxITR_SOR_rhoi);
        READ_INT(category, icon_iter_rhoi);
        READ_DOUBLE(category, error_cnv_SOR_phi);
        READ_INT(category, maxITR_SOR_phi);
        READ_INT(category, icon_iter_phi);
        READ_DOUBLE(category, error_cnv_HES_phi);
        READ_INT(category, maxITR_HES_phi);
        READ_DOUBLE(category, error_cnv_SOR_rhoe);
        READ_INT(category, maxITR_SOR_rhoe);
        READ_INT(category, icon_iter_rhoe);
        READ_DOUBLE(category, error_cnv_HES_rhoe);
        READ_INT(category, maxITR_HES_rhoe);
        READ_DOUBLE(category, error_cnv_SOR_rhoeps);
        READ_INT(category, maxITR_SOR_rhoeps);
        READ_INT(category, icon_iter_rhoeps);
        READ_DOUBLE(category, error_cnv_HES_rhoeps);
        READ_INT(category, maxITR_HES_rhoeps);

        category = "microwave_coupling";
        READ_INT(category, ndiv_MW);

        category = "relaxation";
        READ_DOUBLE(category, beta_rhoe);
        READ_DOUBLE(category, beta_rhoUe);

        category = "simulation";
        READ_INT(category, ntime);
        READ_INT(category, icon_autoFinish);

        category = "output";
        READ_INT(category, icon_chk);
        READ_INT(category, icon_gnuRes);
        READ_INT(category, ndiv_fout);
    }
    catch (const std::exception& e) {
        std::cout << std::endl;
        std::cout << "[setup.yaml] Error while reading parameters." << std::endl;
        std::cout << e.what() << std::endl;
        std::cout << std::endl;
        abort();
    }

#undef READ_DOUBLE
#undef READ_INT

    std::cout << std::endl;
    std::cout << "[" << count << "/58 parameters are input to simulation]" << std::endl;

    if (count != 58) {
        std::cout << std::endl;
        std::cout << "[setup.yaml] Error! Number of input parameters is incorrect." << std::endl;
        abort();
    }

    std::cout << std::endl;
}


//*****************************************************************
//**                                                             **
//**           void iniparameter                                 **
//**                                                             **
//*****************************************************************
void inputParam_csv()
{

    //インプットファイル指定
    std::ifstream ifs("input/setup.csv");

    if (!ifs) {
        std::cout << "[setup.csv] Error! File not found." << std::endl;
        abort();
    }

    std::string line;
    //std::vector<double> input(10,0.0);
    int count=0;

    //読み込み
    while (getline(ifs, line)) {
        
        std::vector<std::string> strvec = split(line, ',');

        //パラメータ設定
        if(strvec[0] == "Pmw"){
            Pmw = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "omegam"){
            omegam = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "nu_eff"){
            nu_eff = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "deltaECR"){
            deltaECR = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "Ti"){
            Ti = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "Tn"){
            Tn = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "rhon_ini"){
            rhon_ini = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "DmN"){
            DmN = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "epsr_diele"){
            epsr_diele = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "V_bias"){
            V_bias = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "Q_neutIn_mgs"){
            Q_neutIn_mgs = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "width_neutIn"){
            width_neutIn = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "Te_rep_eV"){
            Te_rep_eV = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "icon_Bohm"){
            icon_Bohm = std::stoi(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "icon_Sagdeev"){
            icon_Sagdeev = std::stoi(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "alpha_Bohm"){
            alpha_Bohm = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "scale_inertia"){
            scale_inertia = std::stod(strvec[1]);
            count+=1;
        }
        //SEE
        if(strvec[0] == "coefIISEE_ts"){
            coefIISEE_ts = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "coefMISEE_ts"){
            coefMISEE_ts = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "Te_emitSEE_eV"){
            Te_emitSEE_eV = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "ratioEngy_EISEE_rd"){
            ratioEngy_EISEE_rd = std::stod(strvec[1]);
            count+=1;
        }
        //流体スキーム設定
        if(strvec[0] == "icon_PC"){
            icon_PC = std::stoi(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "icon_inertia"){
            icon_inertia = std::stoi(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "icon_adp_dt"){
            icon_adp_dt = std::stoi(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "dt_ini"){
            dt_ini = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "ndt_i"){
            ndt_i = std::stoi(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "ndt_m"){
            ndt_m = std::stoi(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "ndt_n"){
            ndt_n = std::stoi(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "CFL"){
            CFL = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "error_cnv_SOR_Ui"){
            error_cnv_SOR_Ui = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "maxITR_SOR_Ui"){
            maxITR_SOR_Ui = std::stoi(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "icon_iter_Ui"){
            icon_iter_Ui = std::stoi(strvec[1]);
            count+=1;
        }

        if(strvec[0] == "error_cnv_SOR_rhoi"){
            error_cnv_SOR_rhoi = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "maxITR_SOR_rhoi"){
            maxITR_SOR_rhoi = std::stoi(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "icon_iter_rhoi"){
            icon_iter_rhoi = std::stoi(strvec[1]);
            count+=1;
        }
        
        if(strvec[0] == "error_cnv_SOR_phi"){
            error_cnv_SOR_phi = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "maxITR_SOR_phi"){
            maxITR_SOR_phi = std::stoi(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "icon_iter_phi"){
            icon_iter_phi = std::stoi(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "error_cnv_HES_phi"){
            error_cnv_HES_phi = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "maxITR_HES_phi"){
            maxITR_HES_phi = std::stoi(strvec[1]);
            count+=1;
        }

        if(strvec[0] == "error_cnv_SOR_rhoe"){
            error_cnv_SOR_rhoe = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "maxITR_SOR_rhoe"){
            maxITR_SOR_rhoe = std::stoi(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "icon_iter_rhoe"){
            icon_iter_rhoe = std::stoi(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "error_cnv_HES_rhoe"){
            error_cnv_HES_rhoe = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "maxITR_HES_rhoe"){
            maxITR_HES_rhoe = std::stoi(strvec[1]);
            count+=1;
        }

        if(strvec[0] == "error_cnv_SOR_rhoeps"){
            error_cnv_SOR_rhoeps = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "maxITR_SOR_rhoeps"){
            maxITR_SOR_rhoeps = std::stoi(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "icon_iter_rhoeps"){
            icon_iter_rhoeps = std::stoi(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "error_cnv_HES_rhoeps"){
            error_cnv_HES_rhoeps = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "maxITR_HES_rhoeps"){
            maxITR_HES_rhoeps = std::stoi(strvec[1]);
            count+=1;
        }

        if(strvec[0] == "ndiv_MW"){
            ndiv_MW = std::stoi(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "beta_rhoe"){
            beta_rhoe = std::stod(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "beta_rhoUe"){
            beta_rhoUe = std::stod(strvec[1]);
            count+=1;
        }

        //全体の計算設定
        if(strvec[0] == "ntime"){
            ntime = std::stoi(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "icon_autoFinish"){
            icon_autoFinish = std::stoi(strvec[1]);
            count+=1;
        }

        //出力設定
        if(strvec[0] == "icon_chk"){
            icon_chk = std::stoi(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "icon_gnuRes"){
            icon_gnuRes = std::stoi(strvec[1]);
            count+=1;
        }
        if(strvec[0] == "ndiv_fout"){
            ndiv_fout = std::stoi(strvec[1]);
            count+=1;
        }
        
    }

    //std::cout<<"nPcli = "<<nPcli<<std::endl;

    std::cout<<std::endl;
    //std::cout << "icon_WF = " << icon_WF << std::endl;
    std::cout<<"["<<count<<"/58 parameters are input to simulation]"<<std::endl;
    
    if(count != 58){
        std::cout<<std::endl;
        abort();
    }

    std::cout<<std::endl;
}


//*****************************************************************
//**                                                             **
//**           void input_Bfield_data()                          **
//**                                                             **
//*****************************************************************
void input_Bfield_data()
{

    //****************** 磁場データ読み込み ******************
    //インプットファイル指定
    std::ifstream ifs("input/Bfield_data.csv");

    if (!ifs) {
        std::cout << "[Bfield_data.csv] Error! File not found." << std::endl;
        abort();
    }

    std::string line;

    //1行目読む
    getline(ifs, line);
    std::vector<std::string> strvec = split(line, ',');

    //2行目以降読む
    double delta = 0.0;
    while (getline(ifs, line)) {
        
        strvec = split(line, ',');
        int i = std::stoi(strvec[0]);
        int j = std::stoi(strvec[1]);
        Bx[i][j] = std::stod(strvec[4]);
        Br[i][j] = std::stod(strvec[5]);
        Ap[i][j] = std::stod(strvec[6]);
    }
    ifs.close();
    //補間
    for(int i=1;i<ni;i++){
        Bx[i][0] = Bx[i][1];
        Bx[i][nj+1] = Bx[i][nj];
        Br[i][0] = Br[i][1];
        Br[i][nj+1] = Br[i][nj];
        Ap[i][0] = Ap[i][1];
        Ap[i][nj+1] = Ap[i][nj];
    }
    for(int j=1;j<nj;j++){
        Bx[0][j] = Bx[1][j];
        Bx[ni+1][j] = Bx[ni][j];
        Br[0][j] = Br[1][j];
        Br[ni+1][j] = Br[ni][j];
        Ap[0][j] = Ap[1][j];
        Ap[ni+1][j] = Ap[ni][j];
    }

    /*
    //有効衝突周波数を与えるために必要な距離ファイル読み込み
    std::ifstream ifs1("distance.csv");
    if (!ifs1) {
        std::cout << "[distance.csv] Error! File not found." << std::endl;
        abort();
    }
    std::string line1;
    //1行目読む
    getline(ifs1, line1);
    std::vector<std::string> strvec1 = split(line1, ',');

    //2行目以降読む
    while (getline(ifs1, line1)) {
        strvec1 = split(line1, ',');
        int i = std::stoi(strvec1[0]);
        int j = std::stoi(strvec1[1]);
        distECR[i][j] = std::stod(strvec1[8]);
    }

    //補間
    for(int i=1;i<ni;i++){
        distECR[i][0] = distECR[i][1];
        distECR[i][nj+1] = distECR[i][nj];
    }
    for(int j=1;j<nj;j++){
        distECR[0][j] = distECR[1][j];
        distECR[ni+1][j] = distECR[ni][j];
    }
    ifs1.close();
    */

}

//*****************************************************************
//**                                                             **
//**           void input_Bfield_data()                          **
//**                                                             **
//*****************************************************************
void input_restart_data()
{

    //****************** 磁場データ読み込み ******************
    //インプットファイル指定
    std::ifstream ifs("restart/restart.csv");

    if (!ifs) {
        std::cout << "[restart.csv] File not found." << std::endl;
        return;
    }

    std::string line;
    //std::vector<double> input(10,0.0);

    //1行目読む
    getline(ifs, line);
    std::vector<std::string> strvec = split(line, ',');

    //2行目以降読む
    double delta = 0.0;
    while (getline(ifs, line)) {
        
        strvec = split(line, ',');
        int i        = std::stoi(strvec[0]);
        int j        = std::stoi(strvec[1]);

        rhoi[i][j]   = std::stod(strvec[4]);
        rhoi_old[i][j]   = std::stod(strvec[5]);
        Uix[i][j]    = std::stod(strvec[6]);
        Uix_old[i][j]    = std::stod(strvec[7]);
        Uir[i][j]    = std::stod(strvec[8]);
        Uir_old[i][j]    = std::stod(strvec[9]);
        Uip[i][j]    = std::stod(strvec[10]);
        Uip_old[i][j]    = std::stod(strvec[11]);
        rhoe[i][j]   = std::stod(strvec[12]);
        rhoe_old[i][j]   = std::stod(strvec[13]);
        rhoUex[i][j] = std::stod(strvec[14]);
        rhoUex_old[i][j] = std::stod(strvec[15]);
        rhoUer[i][j] = std::stod(strvec[16]);
        rhoUer_old[i][j] = std::stod(strvec[17]);
        Te[i][j] = std::stod(strvec[18])*e0/Boltz;
        Te_old[i][j] = std::stod(strvec[19])*e0/Boltz;
        Gx[i][j]     = std::stod(strvec[20]);
        Gx_old[i][j]     = std::stod(strvec[21]);
        Gr[i][j]     = std::stod(strvec[22]);
        Gr_old[i][j]     = std::stod(strvec[23]);
        phi[i][j]    = std::stod(strvec[24]);
        phi_old[i][j]    = std::stod(strvec[25]);
        rhom[i][j]   = std::stod(strvec[26]);
        
        rhom_old[i][j]   = std::stod(strvec[27]);
        
        rhoUmx[i][j]   = std::stod(strvec[28]);
        rhoUmr[i][j]   = std::stod(strvec[29]);
        rate_ionize[i][j]   = std::stod(strvec[30]);

        if(strvec.size() >= 33){
            //std::cout << "Bingo" << std::endl;
            scx[i][j]   = std::stod(strvec[31]);
            scr[i][j]   = std::stod(strvec[32]);
        }

        if(strvec.size() >= 37){
            //std::cout << "Bingo" << std::endl;
            nUex[i][j]     = std::stod(strvec[33]);
            nUex_old[i][j] = std::stod(strvec[34]);
            nUer[i][j]     = std::stod(strvec[35]);
            nUer_old[i][j] = std::stod(strvec[36]);
        }

        if(strvec.size() >= 41){
            //std::cout << "Bingo" << std::endl;
            rhon[i][j]     = std::stod(strvec[37]);
            rhon_old[i][j] = std::stod(strvec[38]);
            rhoUnx[i][j]   = std::stod(strvec[39]);
            rhoUnr[i][j]   = std::stod(strvec[40]);
        }

        if(strvec.size() >= 43){
            //std::cout << "Bingo" << std::endl;
            rhoeps[i][j] = std::stod(strvec[41]);
            rhoeps_old[i][j] = std::stod(strvec[42]);
        }else{
            rhoeps[i][j] = 3.0/2.0*rhoe[i][j]*Boltz*Te[i][j];
            rhoeps_old[i][j] = 3.0/2.0*rhoe_old[i][j]*Boltz*Te_old[i][j];
        }
    }

    ifs.close();

    //Update Ex
    //------------------------------------
    for (int iblock=0;iblock<n_bl-1;iblock++){
        for (int i=i_flx_bl[iblock][0];i<=i_flx_bl[iblock][1];i++){
            for (int j=j_flx_bl[iblock][0];j<=j_flx_bl[iblock][1];j++){
                Ex[i][j] = - (phi[i][j]-phi[i-1][j])/dx;
            }
        }
    }
    //------------------------------------

    //calculate electric field
    //=====================================================================
    //left-wall-BC (z0) potential = 0
    //------------------------------------
    for (int j=j_flc_bl[0][0];j<=j_flc_bl[0][1];j++){
        int i=i_flc_bl[0][0];
        phi[i-1][j] = - phi[i][j];
        Ex[i][j] =  - (phi[i][j]-phi[i-1][j])/dx;
    }
    //------------------------------------

    //left-wall-BC (z1) potential = 0
    //------------------------------------
    for (int j=j_flc_bl[2][0];j<=j_flc_bl[0][0]-1;j++){
        int i=i_flc_bl[2][0];
        phi[i-1][j] = - phi[i][j];
        Ex[i][j] =  - (phi[i][j]-phi[i-1][j])/dx;
    }
    //------------------------------------

    //left-wall-BC (z2) potential = 0
    //------------------------------------
    for (int j=j_flc_bl[3][0];j<=j_flc_bl[2][0]-1;j++){
        int i=i_flc_bl[3][0];
        phi[i-1][j] = - phi[i][j];
        Ex[i][j] =  - (phi[i][j]-phi[i-1][j])/dx;
    }
    //------------------------------------

    //left-wall-BC (z4) potential = 0
    //------------------------------------
    for (int j=j_flc_bl[1][0];j<=j_flc_bl[4][1];j++){
        int i=i_flc_bl[4][0];
        phi[i-1][j] = - phi[i][j];
        Ex[i][j] =  - (phi[i][j]-phi[i-1][j])/dx;
    }
    //------------------------------------

    //right-wall-BC (z3) potential = 0
    //------------------------------------
    for (int j=j_flc_bl[1][0];j<=j_flc_bl[1][1];j++){
        int i=i_flc_bl[1][1];
        phi[i+1][j] = -phi[i][j];
        Ex[i+1][j] =  - (phi[i+1][j]-phi[i][j])/dx;
    }
    //------------------------------------
    
    //right-wall-BC (z5) potential = V_bias
    //------------------------------------
    for (int j=j_flc_bl[4][0];j<=j_flc_bl[4][1];j++){
        int i=i_flc_bl[4][1];
        phi[i+1][j] = 2.0*V_bias - phi[i][j];
        Ex[i+1][j] =  - (phi[i+1][j]-phi[i][j])/dx;
    }
    //------------------------------------

    //Update Er
    //------------------------------------
    for (int iblock=0;iblock<n_bl-1;iblock++){
        for (int i=i_flr_bl[iblock][0];i<=i_flr_bl[iblock][1];i++){
            for (int j=j_flr_bl[iblock][0];j<=j_flr_bl[iblock][1];j++){
                Er[i][j] = - (phi[i][j]-phi[i][j-1])/dr;
            }
        }
    }
    //------------------------------------

    //lower-wall-BC (x0) potential = 0
    //------------------------------------
    for (int i=i_flc_bl[0][0];i<=i_flc_bl[0][1];i++){
        int j=j_flc_bl[0][0];
        phi[i][j-1] = - phi[i][j];
        Er[i][j] = - (phi[i][j]-phi[i][j-1])/dr;
    }
    //------------------------------------

    //lower-wall-BC (x2) potential = 0
    //------------------------------------
    for (int i=i_flc_bl[2][0];i<=i_flc_bl[2][1];i++){
        int j=j_flc_bl[2][0];
        phi[i][j-1] = - phi[i][j];
        Er[i][j] = - (phi[i][j]-phi[i][j-1])/dr;
    }
    //------------------------------------

    //centerline-normal-BC (x6) Er = 0
    //------------------------------------
    for (int i=i_flc_bl[3][0];i<=i_flc_bl[4][1];i++){ 
        int j = j_flc_bl[3][0];
        Er[i][j] = 0.0;
    }
    //------------------------------------

    //upper-wall-BC (x1) potential = 0
    //------------------------------------
    for (int i=i_flc_bl[0][0];i<=i_flc_bl[1][1];i++){
        int j=j_flc_bl[0][1];
        phi[i][j+1] = - phi[i][j];
        Er[i][j+1] = - (phi[i][j+1]-phi[i][j])/dr;
    }
    //------------------------------------

    //upper-wall-BC (x4) potential = 0
    //------------------------------------
    for (int i=i_flc_bl[1][1]+1;i<=i_flc_bl[3][1];i++){
        int j=j_flc_bl[3][1];
        phi[i][j+1] = - phi[i][j];
        Er[i][j+1] = - (phi[i][j+1]-phi[i][j])/dr;
    }
    //------------------------------------

    //upper-wall-BC (x5) potential = 0
    //------------------------------------
    for (int i=i_flc_bl[4][0];i<=i_flc_bl[4][1];i++){
        int j=j_flc_bl[4][1];
        phi[i][j+1] = phi[i][j];
        Er[i][j+1] = - (phi[i][j+1]-phi[i][j])/dr;
    }
    //------------------------------------

    std::cout << "[restart.csv] Imported" << std::endl;
    icon_restart = 1;



    /*
    ////プラズマ密度プロファイル作成
    double rho_max = 4.0e18*1.0; //最大値　(カットオフ密度は7.45E+16)
    double rho_min = 1.0e14*1.0; //最大値　(カットオフ密度は7.45E+16)
    double sigmax_rho = 0.04; //標準偏差 m
    double sigmar_rho = 0.01; //標準偏差 m 0.05
    double xCen_rho = 0.025; //xの中心 0.009 + 0.0078
    double rCen_rho = 0.0; //rの中心 0.007
    
    for (int i=0;i<=ni+1;i++){
        for (int j=0;j<=nj+1;j++){
            double x_tmp = x[i];
            double r_tmp = r[j];

            double z_x = (x[i]-xCen_rho)/sigmax_rho;
            double z_r = (r[j]-rCen_rho)/sigmar_rho;
        
            rhoi[i][j] = rhoi[i][j] + 1e17*exp(-z_x*z_x/0.12-z_r*z_r/0.12)*jdgBnd_flc[i][j];
            rhoi_old[i][j] = rhoi_old[i][j] + 1e17*exp(-z_x*z_x/0.12-z_r*z_r/0.12)*jdgBnd_flc[i][j];
            rhoe[i][j] = rhoe[i][j] + 1e17*exp(-z_x*z_x/0.12-z_r*z_r/0.12)*jdgBnd_flc[i][j];
            rhoe_old[i][j] = rhoe_old[i][j] + 1e17*exp(-z_x*z_x/0.12-z_r*z_r/0.12)*jdgBnd_flc[i][j];
        }
    }
        */
}


//*****************************************************************
//**                                                             **
//**           void input_BOLSIG_data()                          **
//**                                                             **
//*****************************************************************
void input_BOLSIG_data()
{

    //****************** データ読み込み ******************
    //電子衝突レートインプットファイル指定
    std::ifstream ifs("input/rateCoef_e.csv");

    if (!ifs) {
        std::cout << "[rateCoef_e.csv] Error! File not found." << std::endl;
        abort();
    }

    std::string line;
    //std::vector<double> input(10,0.0);

    //1行目読む
    getline(ifs, line);
    std::vector<std::string> strvec = split(line, ',');
    int nSp = strvec.size()-1; //衝突の種類数

    std::vector<int> col_type; 
    for (int iSp=0;iSp<nSp;iSp++){
        //std::cout << strvec[iSp+1] << std::endl;
        if(strvec[iSp+1] == std::string("Elastic")){
            col_type.push_back(0);
        }else if(strvec[iSp+1] == std::string("Excitation")){
            col_type.push_back(1);
        }else if(strvec[iSp+1] == std::string("Excitation_reso")){
            col_type.push_back(2);
        }else if(strvec[iSp+1] == std::string("Excitation_meta")){
            col_type.push_back(3);
        }else if(strvec[iSp+1] == std::string("Ionization")){
            col_type.push_back(4);
        }else if(strvec[iSp+1] == std::string("Ionization_step")){
            col_type.push_back(5);
        }else if(strvec[iSp+1] == std::string("Excitation_step")){
            col_type.push_back(6);
        }else{
            col_type.push_back(7); //super-elastic
        }
    }

    //2行目読む
    getline(ifs, line);
    strvec = split(line, ',');

    //閾値のエネルギーを書き込み
    std::vector<double> threshold(nSp,0.0); 
    for (int iSp=0;iSp<nSp;iSp++){
        threshold[iSp] = std::stod(strvec[iSp+1])*e0; //eV -> J
        //std::cout << "iSp = " << iSp << ","<<threshold[iSp]/e0 << ","<<threshold[iSp]  << std::endl;
    }

    std::vector<std::vector<double> > data(nSp,std::vector<double>());
    int ncount=0;

    //3行目以降読む
    double tmp0 = 0;
    double tmp1 = 0;

    while (getline(ifs, line)) {
        
        strvec = split(line, ',');

        //std::cout <<strvec[0]<< std::endl;
        for (int iSp=0;iSp<nSp;iSp++){
            //std::cout << ncount << ","<<iSp << std::endl;
            data[iSp].push_back(std::stod(strvec[iSp+1]));
        }

        if(ncount == 0){
            tmp0 = std::stod(strvec[0]);
        }else if(ncount == 1){
            tmp1 = std::stod(strvec[0]);
        }

        ncount = ncount + 1;
    }
    
    int ndata_BSG = ncount;
    double delta = tmp1-tmp0;

    //データ数カウント
    int nElas = 0;     //弾性衝突
    int nExcReso = 0;  //励起衝突 (1s4と1s2)
    int nExcMeta = 0;  //準安定衝突
    int nExc = 0;      //励起衝突
    int nIonz = 0;     //電離衝突
    int nIonzStep = 0; //段階電離衝突
    int nExcStep = 0;  //段階電離衝突
    int nSuper = 0;    //超弾性衝突
    
    for (int iSp=0;iSp<nSp;iSp++){
        //std::cout << col_name[iSp] << std::endl;
        if(col_type[iSp] == 0){
            nElas++;
        }else if(col_type[iSp] == 1){
            nExc++;
        }else if(col_type[iSp] == 2){
            nExcReso++;
        }else if(col_type[iSp] == 3){
            nExcMeta++;
        }else if(col_type[iSp] == 4){
            nIonz++;
        }else if(col_type[iSp] == 5){
            nIonzStep++;
        }else if(col_type[iSp] == 6){
            nExcStep++;
        }else if(col_type[iSp] == 7){
            nSuper++;
        }
    }

    std::cout<<"["<<nElas << " Elastic & " 
    << nExcReso << " Excitation_reso & "<<nExcMeta << " Excitation_meta & "<<nExc << " Excitation_upper & " 
    << nIonz << " Ionization & " 
    << nIonzStep << " Ionization_step " << nExcStep << " Excitation_step "  << nSuper << " Superelastic " 
    <<  "data are input to simulation]"<<std::endl;
    std::cout<< " Total data num = " << ndata_BSG << ", delta = " << delta << " eV, Maximum = "<< delta*(ndata_BSG-1) << " eV"<< std::endl;
    std::cout<<std::endl;
    ifs.close();

    //イオン衝突レートインプットファイル指定
    std::ifstream ifs1("input/rateCoef_i.csv");

    if (!ifs1) {
        std::cout << "[rateCoef_i.csv] Error! File not found." << std::endl;
        abort();
    }

    std::string line1;
    //std::vector<double> input(10,0.0);

    //1行目読む
    getline(ifs1, line1);
    std::vector<std::string> strvec1 = split(line1, ',');
    int nSp1 = strvec1.size()-1; //衝突の種類数

    //2行目読む
    getline(ifs1, line1);
    strvec1 = split(line1, ',');


    //閾値のエネルギーを書き込み
    std::vector<double> threshold1(nSp1,0.0); 
    for (int iSp=0;iSp<nSp1;iSp++){
        threshold1[iSp] = std::stod(strvec1[iSp+1])*e0; //eV -> J
        //std::cout << "iSp = " << iSp << ","<<threshold[iSp]/e0 << ","<<threshold[iSp]  << std::endl;
    }

    std::vector<std::vector<double> > data1(nSp1,std::vector<double>());
    ncount=0;

    //3行目以降読む
    
    while (getline(ifs1, line1)) {
        
        strvec1 = split(line1, ',');

        //std::cout <<strvec1[0]<< std::endl;
        for (int iSp=0;iSp<nSp1;iSp++){
            data1[iSp].push_back(std::stod(strvec1[iSp+1]));
        }

        if(ncount == 0){
            tmp0 = std::stod(strvec1[0]);
            //std::cout << "delta = " << delta << std::endl;
        }else if(ncount == 1){
            tmp1 = std::stod(strvec1[0]);
            //std::cout << "delta = " << delta << std::endl;
        }

        ncount = ncount + 1;
    }
    int ndata_File = ncount;
    double delta1 = tmp1-tmp0;

    std::cout<<"["<<nSp1 << " Ion collision data are input to simulation]"<<std::endl;
    std::cout<< " Total data num = " << ndata_File << ", delta = " << delta1 << " eV, Maximum = "<< delta1*(ndata_File-1) << " eV"<< std::endl;
    std::cout<<std::endl;


    //****************** プログラム内で扱う配列に変換する ******************
    //電子側
    nu_elas_N.resize(ndata_BSG,0.0); //弾性周波数/nn <σelas*vth>
    nu_exc_N.resize(ndata_BSG,0.0); //励起周波数/nn <σexc*vth>
    nu_excReso_N.resize(ndata_BSG,0.0); //励起周波数/nn <σexc*vth>
    nu_excMeta_N.resize(ndata_BSG,0.0); //励起周波数/nn <σexc*vth>
    nu_ionz_N.resize(ndata_BSG,0.0); //イオン化周波数/nn <σionz*vth>
    nu_ionzStep_N.resize(ndata_BSG,0.0); //イオン化周波数/nn <σionz*vth>
    nu_excStep_N.resize(ndata_BSG,0.0); //段階電離周波数/nn <σionz*vth>
    nu_super_N.resize(ndata_BSG,0.0); //超弾性周波数/nn <σsuper*vth>
    rate_eloss_n_N.resize(ndata_BSG,0.0); //エネルギーロスレート/nn
    rate_eloss_m_N.resize(ndata_BSG,0.0); //エネルギーロスレート/nm

    //std::cout << "OK1" << std::endl;

    for (int i=0;i<ndata_BSG;i++){
        for (int iSp=0;iSp<nSp;iSp++){
            //std::cout << i << ","<<iSp << std::endl;
            if(col_type[iSp] == 0){ //elastic
                nu_elas_N[i] = nu_elas_N[i] + data[iSp][i];
                rate_eloss_n_N[i] = rate_eloss_n_N[i] +  data[iSp][i]*threshold[iSp];
            }else if(col_type[iSp] == 1){ //excitation
                nu_exc_N[i] = nu_exc_N[i] + data[iSp][i];
                rate_eloss_n_N[i] = rate_eloss_n_N[i] +  data[iSp][i]*threshold[iSp];
            }else if(col_type[iSp] == 2){ //resonance-excitation
                nu_excReso_N[i] = nu_excReso_N[i] + data[iSp][i];
                rate_eloss_n_N[i] = rate_eloss_n_N[i] +  data[iSp][i]*threshold[iSp];
            }else if(col_type[iSp] == 3){ //metastable-excitation
                nu_excMeta_N[i] = nu_excMeta_N[i] + data[iSp][i];
                rate_eloss_n_N[i] = rate_eloss_n_N[i] +  data[iSp][i]*threshold[iSp];
            }else if(col_type[iSp] == 4){ //ionization
                nu_ionz_N[i] = nu_ionz_N[i] + data[iSp][i];
                rate_eloss_n_N[i] = rate_eloss_n_N[i] +  data[iSp][i]*threshold[iSp];
            }else if(col_type[iSp] == 5){ //stepwise ionization
                nu_ionzStep_N[i] = nu_ionzStep_N[i] + data[iSp][i];
                rate_eloss_m_N[i] = rate_eloss_m_N[i] +  data[iSp][i]*threshold[iSp];
            }else if(col_type[iSp] == 6){ //stepwise excitation
                nu_excStep_N[i] = nu_excStep_N[i] + data[iSp][i];
                rate_eloss_m_N[i] = rate_eloss_m_N[i] +  data[iSp][i]*threshold[iSp];
            }else if(col_type[iSp] == 7){ //superelastic
                nu_super_N[i] = nu_super_N[i] + data[iSp][i];
                rate_eloss_m_N[i] = rate_eloss_m_N[i] +  data[iSp][i]*threshold[iSp];
            }
        }
    }

    //イオン
    nui_m_N.resize(ndata_File,0.0); //イオン衝突周波数/nn <σexc*vth>
    for (int i=0;i<ndata_File;i++){
        for (int iSp=0;iSp<nSp1;iSp++){
            nui_m_N[i] = nui_m_N[i] + data1[iSp][i];
        }
    }

    //dTe,dEi計算
    dTe = delta*2.0/3.0*e0/Boltz;
    dEi = delta1*e0;

}


//*****************************************************************
//**                                                             **
//**           void input_BOLSIG_data()                          **
//**                                                             **
//*****************************************************************
void input_SEE_data()
{

    //****************** データ読み込み ******************
    //電子衝突レートインプットファイル指定
    std::ifstream ifs("input/coefEISEE.csv");

    if (!ifs) {
        std::cout << "[coefEISEE.csv] Error! File not found." << std::endl;
        abort();
    }

    std::string line;
    //std::vector<double> input(10,0.0);

    //1行目読む
    getline(ifs, line);
    std::vector<std::string> strvec = split(line, ',');
   
    std::vector<std::vector<double> > data(3,std::vector<double>());
    int ncount=0;

    //2行目以降読む
    double tmp0 = 0;
    double tmp1 = 0;

    while (getline(ifs, line)) {
        
        strvec = split(line, ',');

        coefEISEE_eb_table.push_back(std::stod(strvec[1])); //EB
        coefEISEE_rd_table.push_back(std::stod(strvec[2])); //RD
        coefEISEE_ts_table.push_back(std::stod(strvec[3])); //TS

        if(ncount == 0){
            tmp0 = std::stod(strvec[0]);
        }else if(ncount == 1){
            tmp1 = std::stod(strvec[0]);
        }

        ncount = ncount + 1;
    }

    //dTe
    double delta = tmp1-tmp0;
    dTe_SEE = delta*0.5*e0/Boltz;

    std::cout<<"[ 3 EI-SEE data are input to simulation]"<<std::endl;
    std::cout<< " Total data num = " << ncount << ", delta = " << delta << " eV (mean-energy-base) "
    << delta*0.5 << " eV (Te-base) "<< std::endl;
    std::cout<<std::endl;
    ifs.close();
}
