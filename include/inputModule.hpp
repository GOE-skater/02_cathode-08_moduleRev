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
#include "constants.hpp"
#include "params.hpp"
#include <sstream>
#include <mpi.h>
#include <yaml-cpp/yaml.h>
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <cstdlib>
#include <algorithm>

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
class InputModule
{
    private:

    public:
        void inputParam(Params &pm,string inputFileName);
        void input_Bfield_data(Params &pm, GridCenter &gc, string inputFileName);
        void input_SEE_data(Params &pm, SeeVec &se, string inputFileName);
        void input_restart_data(Params &pm, GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr, string inputFileName);
        void input_BOLSIG_data(Params &pm, BolsigVec &bo, string inputFileName);
};

//*****************************************************************
//**                                                             **
//**           void param                                        **
//**                                                             **
//*****************************************************************
void InputModule::inputParam(Params &pm, string inputFileName)
{

    const string filename = "input/"+inputFileName;

    YAML::Node setup;

    try {
        setup = YAML::LoadFile(filename);
    }
    catch (const YAML::BadFile&) {
        cout << "[" << filename << "] Error! File not found." << endl;
        abort();
    }

    vector<string> missingParams;

    // READ_DOUBLE(category, pm.Pmw) の #var は "pm.Pmw" になるため、
    // YAML のキーとして使う前に "pm." を外して "Pmw" に正規化する。
    // "pm. error_cnv_HES_phi" のようにドット後に空白がある場合も吸収する。
    auto normalizeYamlKey = [](const string& rawKey) {
        string key;
        key.reserve(rawKey.size());

        for (unsigned char c : rawKey) {
            if (!isspace(c)) {
                key.push_back(static_cast<char>(c));
            }
        }

        const size_t dotPos = key.find_last_of('.');
        if (dotPos != string::npos) {
            key = key.substr(dotPos + 1);
        }

        return key;
    };

    auto readDouble = [&](double& dst,
                          const string& category,
                          const string& key) {
        const YAML::Node categoryNode = setup[category];
        const string yamlKey = normalizeYamlKey(key);

        if (!categoryNode) {
            missingParams.push_back(category + "." + yamlKey + "  [missing category: " + category + "]");
            return;
        }

        if (!categoryNode[yamlKey]) {
            missingParams.push_back(category + "." + yamlKey + "  [C++ name: " + key + "]");
            return;
        }

        try {
            dst = categoryNode[yamlKey].as<double>();
        }
        catch (const exception& e) {
            throw runtime_error(
                "Type conversion error: " + category + "." + yamlKey + " should be double"
            );
        }
    };

    auto readInt = [&](int& dst,
                       const string& category,
                       const string& key) {
        const YAML::Node categoryNode = setup[category];
        const string yamlKey = normalizeYamlKey(key);

        if (!categoryNode) {
            missingParams.push_back(category + "." + yamlKey + "  [missing category: " + category + "]");
            return;
        }

        if (!categoryNode[yamlKey]) {
            missingParams.push_back(category + "." + yamlKey + "  [C++ name: " + key + "]");
            return;
        }

        try {
            dst = categoryNode[yamlKey].as<int>();
        }
        catch (const exception& e) {
            throw runtime_error(
                "Type conversion error: " + category + "." + yamlKey + " should be int"
            );
        }
    };

#define READ_DOUBLE(category, var) readDouble((var), (category), #var)
#define READ_INT(category, var)    readInt((var), (category), #var)

    try {
        string category;

        category = "microwave";
        READ_DOUBLE(category, pm.Pmw);
        READ_INT(category,    pm.icon_mwRef);
        READ_INT(category,    pm.icon_impTest);
        READ_DOUBLE(category, pm.omegam);
        READ_DOUBLE(category, pm.nu_eff);
        READ_DOUBLE(category, pm.deltaECR);
        READ_DOUBLE(category, pm.S11_mag);     // (double) Magnitude of S11
        READ_DOUBLE(category, pm.S11_arg_deg); // (double) Argument of S11 (deg)
        READ_DOUBLE(category, pm.S21_mag);     // (double) Magnitude of S21
        READ_DOUBLE(category, pm.S21_arg_deg); // (double) Argument of S21 (deg)
        READ_DOUBLE(category, pm.S12_mag);     // (double) Magnitude of S12
        READ_DOUBLE(category, pm.S12_arg_deg); // (double) Argument of S12 (deg)
        READ_DOUBLE(category, pm.S22_mag);     // (double) Magnitude of S22
        READ_DOUBLE(category, pm.S22_arg_deg); // (double) Argument of S22 (deg)
        READ_DOUBLE(category, pm.Z0_base); // (double) Argument of S22 (deg)

        category = "plasma";
        READ_DOUBLE(category, pm.Ti);
        READ_DOUBLE(category, pm.Tn);
        READ_DOUBLE(category, pm.rhon_ini);
        READ_DOUBLE(category, pm.DmN);
        READ_DOUBLE(category, pm.Te_rep_eV);

        category = "material";
        READ_DOUBLE(category, pm.epsr_diele);

        category = "bias";
        READ_DOUBLE(category, pm.V_bias);

        category = "neutral_inlet";
        READ_DOUBLE(category, pm.Q_neutIn_mgs);
        READ_DOUBLE(category, pm.width_neutIn);

        category = "transport_model";
        READ_INT(category,    pm.icon_Bohm);
        READ_INT(category,    pm.icon_Sagdeev);
        READ_DOUBLE(category, pm.alpha_Bohm);
        READ_DOUBLE(category, pm.scale_inertia);

        category = "see";
        READ_DOUBLE(category, pm.coefIISEE_ts);
        READ_DOUBLE(category, pm.coefMISEE_ts);
        READ_DOUBLE(category, pm.Te_emitSEE_eV);
        READ_DOUBLE(category, pm.ratioEngy_EISEE_rd);

        category = "scheme";
        READ_INT(category,    pm.icon_PC);
        READ_INT(category,    pm.icon_inertia);
        READ_INT(category,    pm.icon_adp_dt);
        READ_DOUBLE(category, pm.dt_ini);
        READ_INT(category,    pm.ndt_i);
        READ_INT(category,    pm.ndt_m);
        READ_INT(category,    pm.ndt_n);
        READ_DOUBLE(category, pm.CFL);

        category = "solver";
        READ_DOUBLE(category, pm.error_cnv_SOR_Ui);
        READ_INT(category,    pm.maxITR_SOR_Ui);
        READ_INT(category,    pm.icon_iter_Ui);

        READ_DOUBLE(category, pm.error_cnv_SOR_rhoi);
        READ_INT(category,    pm.maxITR_SOR_rhoi);
        READ_INT(category,    pm.icon_iter_rhoi);

        READ_DOUBLE(category, pm.error_cnv_SOR_phi);
        READ_INT(category,    pm.maxITR_SOR_phi);
        READ_INT(category,    pm.icon_iter_phi);
        READ_DOUBLE(category, pm. error_cnv_HES_phi);
        READ_INT(category,    pm.maxITR_HES_phi);

        READ_DOUBLE(category, pm.error_cnv_SOR_rhoe);
        READ_INT(category,    pm.maxITR_SOR_rhoe);
        READ_INT(category,    pm.icon_iter_rhoe);
        READ_DOUBLE(category, pm.error_cnv_HES_rhoe);
        READ_INT(category,    pm.maxITR_HES_rhoe);

        READ_DOUBLE(category, pm.error_cnv_SOR_rhoeps);
        READ_INT(category,    pm.maxITR_SOR_rhoeps);
        READ_INT(category,    pm.icon_iter_rhoeps);
        READ_DOUBLE(category, pm.error_cnv_HES_rhoeps);
        READ_INT(category,    pm.maxITR_HES_rhoeps);

        category = "microwave_coupling";
        READ_INT(category,    pm.ndiv_MW);

        category = "relaxation";
        READ_DOUBLE(category, pm.beta_rhoe);
        READ_DOUBLE(category, pm.beta_rhoUe);

        category = "simulation";
        READ_INT(category,    pm.ntime);
        READ_INT(category,    pm.icon_autoFinish);

        category = "output";
        READ_INT(category,    pm.icon_chk);
        READ_INT(category,    pm.icon_gnuRes);
        READ_INT(category,    pm.ndiv_fout);
    }
    catch (const exception& e) {
        cout << endl;
        cout << "[setup.yaml] Error while reading parameters." << endl;
        cout << e.what() << endl;
        cout << endl;
        abort();
    }

#undef READ_DOUBLE
#undef READ_INT

    if (!missingParams.empty()) {
        cout << endl;
        cout << "[setup.yaml] Error! Missing YAML parameters:" << endl;

        for (const auto& name : missingParams) {
            cout << "  - " << name << endl;
        }

        cout << endl;
        abort();
    }

    cout << endl;
    cout << "[All setup parameters are successfully loaded]" << endl;
    cout << endl;

}

//*****************************************************************
//**                                                             **
//**           void param                                        **
//**                                                             **
//*****************************************************************
void inputParam_old(Params &pm, string inputFileName)
{

    const string filename = "input/"+inputFileName;

    YAML::Node setup;

    try {
        setup = YAML::LoadFile(filename);
    }
    catch (const YAML::BadFile&) {
        cout << "[" << filename << "] Error! File not found." << endl;
        abort();
    }

    vector<string> missingParams;

    auto readDouble = [&](double& dst,
                          const string& category,
                          const string& key) {
        const YAML::Node categoryNode = setup[category];

        if (!categoryNode) {
            missingParams.push_back(category + "." + key + "  [missing category: " + category + "]");
            return;
        }

        if (!categoryNode[key]) {
            missingParams.push_back(category + "." + key);
            return;
        }

        try {
            dst = categoryNode[key].as<double>();
        }
        catch (const exception& e) {
            throw runtime_error(
                "Type conversion error: " + category + "." + key + " should be double"
            );
        }
    };

    auto readInt = [&](int& dst,
                       const string& category,
                       const string& key) {
        const YAML::Node categoryNode = setup[category];

        if (!categoryNode) {
            missingParams.push_back(category + "." + key + "  [missing category: " + category + "]");
            return;
        }

        if (!categoryNode[key]) {
            missingParams.push_back(category + "." + key);
            return;
        }

        try {
            dst = categoryNode[key].as<int>();
        }
        catch (const exception& e) {
            throw runtime_error(
                "Type conversion error: " + category + "." + key + " should be int"
            );
        }
    };

#define READ_DOUBLE(category, var) readDouble((var), (category), #var)
#define READ_INT(category, var)    readInt((var), (category), #var)

    try {
        string category;

        category = "microwave";
        READ_DOUBLE(category, pm.Pmw);
        READ_INT(category,    pm.icon_mwRef);
        READ_INT(category,    pm.icon_impTest);
        READ_DOUBLE(category, pm.omegam);
        READ_DOUBLE(category, pm.nu_eff);
        READ_DOUBLE(category, pm.deltaECR);
        READ_DOUBLE(category, pm.S11_mag);     // (double) Magnitude of S11
        READ_DOUBLE(category, pm.S11_arg_deg); // (double) Argument of S11 (deg)
        READ_DOUBLE(category, pm.S21_mag);     // (double) Magnitude of S21
        READ_DOUBLE(category, pm.S21_arg_deg); // (double) Argument of S21 (deg)
        READ_DOUBLE(category, pm.S12_mag);     // (double) Magnitude of S12
        READ_DOUBLE(category, pm.S12_arg_deg); // (double) Argument of S12 (deg)
        READ_DOUBLE(category, pm.S22_mag);     // (double) Magnitude of S22
        READ_DOUBLE(category, pm.S22_arg_deg); // (double) Argument of S22 (deg)
        READ_DOUBLE(category, pm.Z0_base); // (double) Argument of S22 (deg)

        category = "plasma";
        READ_DOUBLE(category, pm.Ti);
        READ_DOUBLE(category, pm.Tn);
        READ_DOUBLE(category, pm.rhon_ini);
        READ_DOUBLE(category, pm.DmN);
        READ_DOUBLE(category, pm.Te_rep_eV);

        category = "material";
        READ_DOUBLE(category, pm.epsr_diele);

        category = "bias";
        READ_DOUBLE(category, pm.V_bias);

        category = "neutral_inlet";
        READ_DOUBLE(category, pm.Q_neutIn_mgs);
        READ_DOUBLE(category, pm.width_neutIn);

        category = "transport_model";
        READ_INT(category,    pm.icon_Bohm);
        READ_INT(category,    pm.icon_Sagdeev);
        READ_DOUBLE(category, pm.alpha_Bohm);
        READ_DOUBLE(category, pm.scale_inertia);

        category = "see";
        READ_DOUBLE(category, pm.coefIISEE_ts);
        READ_DOUBLE(category, pm.coefMISEE_ts);
        READ_DOUBLE(category, pm.Te_emitSEE_eV);
        READ_DOUBLE(category, pm.ratioEngy_EISEE_rd);

        category = "scheme";
        READ_INT(category,    pm.icon_PC);
        READ_INT(category,    pm.icon_inertia);
        READ_INT(category,    pm.icon_adp_dt);
        READ_DOUBLE(category, pm.dt_ini);
        READ_INT(category,    pm.ndt_i);
        READ_INT(category,    pm.ndt_m);
        READ_INT(category,    pm.ndt_n);
        READ_DOUBLE(category, pm.CFL);

        category = "solver";
        READ_DOUBLE(category, pm.error_cnv_SOR_Ui);
        READ_INT(category,    pm.maxITR_SOR_Ui);
        READ_INT(category,    pm.icon_iter_Ui);

        READ_DOUBLE(category, pm.error_cnv_SOR_rhoi);
        READ_INT(category,    pm.maxITR_SOR_rhoi);
        READ_INT(category,    pm.icon_iter_rhoi);

        READ_DOUBLE(category, pm.error_cnv_SOR_phi);
        READ_INT(category,    pm.maxITR_SOR_phi);
        READ_INT(category,    pm.icon_iter_phi);
        READ_DOUBLE(category, pm. error_cnv_HES_phi);
        READ_INT(category,    pm.maxITR_HES_phi);

        READ_DOUBLE(category, pm.error_cnv_SOR_rhoe);
        READ_INT(category,    pm.maxITR_SOR_rhoe);
        READ_INT(category,    pm.icon_iter_rhoe);
        READ_DOUBLE(category, pm.error_cnv_HES_rhoe);
        READ_INT(category,    pm.maxITR_HES_rhoe);

        READ_DOUBLE(category, pm.error_cnv_SOR_rhoeps);
        READ_INT(category,    pm.maxITR_SOR_rhoeps);
        READ_INT(category,    pm.icon_iter_rhoeps);
        READ_DOUBLE(category, pm.error_cnv_HES_rhoeps);
        READ_INT(category,    pm.maxITR_HES_rhoeps);

        category = "microwave_coupling";
        READ_INT(category,    pm.ndiv_MW);

        category = "relaxation";
        READ_DOUBLE(category, pm.beta_rhoe);
        READ_DOUBLE(category, pm.beta_rhoUe);

        category = "simulation";
        READ_INT(category,    pm.ntime);
        READ_INT(category,    pm.icon_autoFinish);

        category = "output";
        READ_INT(category,    pm.icon_chk);
        READ_INT(category,    pm.icon_gnuRes);
        READ_INT(category,    pm.ndiv_fout);
    }
    catch (const exception& e) {
        cout << endl;
        cout << "[setup.yaml] Error while reading parameters." << endl;
        cout << e.what() << endl;
        cout << endl;
        abort();
    }

#undef READ_DOUBLE
#undef READ_INT

    if (!missingParams.empty()) {
        cout << endl;
        cout << "[setup.yaml] Error! Missing YAML parameters:" << endl;

        for (const auto& name : missingParams) {
            cout << "  - " << name << endl;
        }

        cout << endl;
        abort();
    }

    cout << endl;
    cout << "[All setup parameters are successfully loaded]" << endl;
    cout << endl;

}

//*****************************************************************
//**                                                             **
//**           void input_Bfield_data()                          **
//**                                                             **
//*****************************************************************
void InputModule::input_Bfield_data(Params &pm, GridCenter &gc, string inputFileName)
{

    //****************** 磁場データ読み込み ******************
    //インプットファイル指定
    //ifstream ifs("input/Bfield_data.csv");
    ifstream ifs("input/" + inputFileName);
    
    if (!ifs) {
        cout << "[Bfield_data.csv] Error! File not found." << endl;
        abort();
    }

    string line;

    //1行目読む
    getline(ifs, line);
    vector<string> strvec = split(line, ',');

    //2行目以降読む
    double delta = 0.0;
    while (getline(ifs, line)) {
        
        strvec = split(line, ',');
        int i = stoi(strvec[0]);
        int j = stoi(strvec[1]);
        gc.Bx[i][j] = stod(strvec[4]);
        gc.Br[i][j] = stod(strvec[5]);
        gc.Ap[i][j] = stod(strvec[6]);
    }
    ifs.close();
    //補間
    for(int i=1;i<pm.ni;i++){
        gc.Bx[i][0] = gc.Bx[i][1];
        gc.Bx[i][pm.nj+1] = gc.Bx[i][pm.nj];
        gc.Br[i][0] = gc.Br[i][1];
        gc.Br[i][pm.nj+1] = gc.Br[i][pm.nj];
        gc.Ap[i][0] = gc.Ap[i][1];
        gc.Ap[i][pm.nj+1] = gc.Ap[i][pm.nj];
    }
    for(int j=1;j<pm.nj;j++){
        gc.Bx[0][j] = gc.Bx[1][j];
        gc.Bx[pm.ni+1][j] = gc.Bx[pm.ni][j];
        gc.Br[0][j] = gc.Br[1][j];
        gc.Br[pm.ni+1][j] = gc.Br[pm.ni][j];
        gc.Ap[0][j] = gc.Ap[1][j];
        gc.Ap[pm.ni+1][j] = gc.Ap[pm.ni][j];
    }

    /*
    //有効衝突周波数を与えるために必要な距離ファイル読み込み
    ifstream ifs1("distance.csv");
    if (!ifs1) {
        cout << "[distance.csv] Error! File not found." << endl;
        abort();
    }
    string line1;
    //1行目読む
    getline(ifs1, line1);
    vector<string> strvec1 = split(line1, ',');

    //2行目以降読む
    while (getline(ifs1, line1)) {
        strvec1 = split(line1, ',');
        int i = stoi(strvec1[0]);
        int j = stoi(strvec1[1]);
        distECR[i][j] = stod(strvec1[8]);
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
//**           void input_SEE_data()                             **
//**                                                             **
//*****************************************************************
void InputModule::input_SEE_data(Params &pm, SeeVec &se, string inputFileName)
{

    //****************** データ読み込み ******************
    //電子衝突レートインプットファイル指定
    //ifstream ifs("input/coefEISEE.csv");
    ifstream ifs("input/" + inputFileName);

    if (!ifs) {
        cout << "[coefEISEE.csv] Error! File not found." << endl;
        abort();
    }

    string line;
    //vector<double> input(10,0.0);

    //1行目読む
    getline(ifs, line);
    vector<string> strvec = split(line, ',');
   
    vector<vector<double> > data(3,vector<double>());
    int ncount=0;

    //2行目以降読む
    double tmp0 = 0;
    double tmp1 = 0;

    while (getline(ifs, line)) {
        
        strvec = split(line, ',');

        se.coefEISEE_eb_table.push_back(stod(strvec[1])); //EB
        se.coefEISEE_rd_table.push_back(stod(strvec[2])); //RD
        se.coefEISEE_ts_table.push_back(stod(strvec[3])); //TS

        if(ncount == 0){
            tmp0 = stod(strvec[0]);
        }else if(ncount == 1){
            tmp1 = stod(strvec[0]);
        }

        ncount = ncount + 1;
    }

    //dTe
    double delta = tmp1-tmp0;
    pm.dTe_SEE = delta*0.5*ph::e0/ph::Boltz;

    cout<<"[ 3 EI-SEE data are input to simulation]"<<endl;
    cout<< " Total data num = " << ncount << ", delta = " << delta << " eV (mean-energy-base) "
    << delta*0.5 << " eV (Te-base) "<< endl;
    cout<<endl;
    ifs.close();
}


//*****************************************************************
//**                                                             **
//**           void input_restart_data()                         **
//**                                                             **
//*****************************************************************
void InputModule::input_restart_data(Params &pm, GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr, string inputFileName)
{

    //****************** Read restart data ******************

    ifstream ifs("restart/" + inputFileName);

    if (!ifs) {
        cout << "[restart.csv] File not found." << endl;
        return;
    }

    string line;

    // Read the header line.
    if (!getline(ifs, line)) {
        cerr << "[restart.csv] Empty file." << endl;
        abort();
    }

    vector<string> header = split(line, ',');

    // Remove spaces, tabs, and CR/LF from header names.
    // This is important because the last CSV header may contain '\r'.
    auto trim = [](string s) -> string {
        const string whitespace = " \t\r\n";

        size_t first = s.find_first_not_of(whitespace);
        if (first == string::npos) {
            return "";
        }

        size_t last = s.find_last_not_of(whitespace);
        return s.substr(first, last - first + 1);
    };

    for (auto &name : header) {
        name = trim(name);
    }

    // Create a map from column name to column index.
    unordered_map<string, int> col;
    for (int k = 0; k < static_cast<int>(header.size()); k++) {
        col[header[k]] = k;
    }
    
    // Helper function to check whether a column exists.
    auto has = [&](const string &name) -> bool {
        return col.find(name) != col.end();
    };

    // Helper function to abort if a required column is missing.
    auto require_column = [&](const string &name) {
        if (!has(name)) {
            cerr << "[restart.csv] Required column is missing: " << name << endl;
            abort();
        }
    };

    // i, j, x, and r define the cell arrangement and must exist.
    require_column("i");
    require_column("j");
    require_column("x");
    require_column("r");

    // A restart variable consists of the CSV header name and the assignment operation.
    struct RestartVariable {
        string name;
        function<void(int, int, double)> assign;
    };

    vector<RestartVariable> programVars;

    // Register a variable whose CSV header name is identical to the program variable name.
    #define READ_RESTART_VAR(obj, var) \
        programVars.push_back({#var, [&](int i, int j, double v) { obj.var[i][j] = v; }})

    // Register a variable that needs a custom conversion when reading from CSV.
    #define READ_RESTART_VAR_CONVERT(obj, var, expr) \
        programVars.push_back({#var, [&](int i, int j, double v) { obj.var[i][j] = (expr); }})

    // Register restart variables.
    // To add a new variable, add one line here.
    READ_RESTART_VAR(gc, rhoi);
    READ_RESTART_VAR(gc, rhoi_old);

    READ_RESTART_VAR(gx, Uix);
    READ_RESTART_VAR(gx, Uix_old);

    READ_RESTART_VAR(gr, Uir);
    READ_RESTART_VAR(gr, Uir_old);

    READ_RESTART_VAR(gc, Uip);
    READ_RESTART_VAR(gc, Uip_old);

    READ_RESTART_VAR(gc, rhoe);
    READ_RESTART_VAR(gc, rhoe_old);

    READ_RESTART_VAR(gx, rhoUex);
    READ_RESTART_VAR(gx, rhoUex_old);

    READ_RESTART_VAR(gr, rhoUer);
    READ_RESTART_VAR(gr, rhoUer_old);

    // Te in restart.csv is assumed to be written in eV.
    // The program stores Te in Kelvin.
    READ_RESTART_VAR_CONVERT(gc, Te,     v * ph::e0 / ph::Boltz);
    READ_RESTART_VAR_CONVERT(gc, Te_old, v * ph::e0 / ph::Boltz);

    READ_RESTART_VAR(gx, Gx);
    READ_RESTART_VAR(gx, Gx_old);

    READ_RESTART_VAR(gr, Gr);
    READ_RESTART_VAR(gr, Gr_old);

    READ_RESTART_VAR(gc, phi);
    READ_RESTART_VAR(gc, phi_old);

    READ_RESTART_VAR(gc, rhom);
    READ_RESTART_VAR(gc, rhom_old);

    READ_RESTART_VAR(gx, rhoUmx);
    READ_RESTART_VAR(gr, rhoUmr);

    READ_RESTART_VAR(gc, rate_ionize);

    READ_RESTART_VAR(gx, scx);
    READ_RESTART_VAR(gr, scr);

    READ_RESTART_VAR(gx, nUex);
    READ_RESTART_VAR(gx, nUex_old);

    READ_RESTART_VAR(gr, nUer);
    READ_RESTART_VAR(gr, nUer_old);

    READ_RESTART_VAR(gc, rhon);
    READ_RESTART_VAR(gc, rhon_old);

    READ_RESTART_VAR(gx, rhoUnx);
    READ_RESTART_VAR(gr, rhoUnr);

    READ_RESTART_VAR(gc, rhoeps);
    READ_RESTART_VAR(gc, rhoeps_old);

    #undef READ_RESTART_VAR
    #undef READ_RESTART_VAR_CONVERT

    // Create a set of variable names used by this program.
    unordered_set<string> programVarNames;
    for (const auto &var : programVars) {
        programVarNames.insert(var.name);
    }

    // Columns used only for checking the cell definition.
    unordered_set<string> coordinateVars = {
        "i", "j", "x", "r"
    };

    // Columns intentionally ignored even if they exist in restart.csv.
    unordered_set<string> ignoredCsvVars = {
        "zero"
    };

    for (const auto &name : header) {
        if (coordinateVars.count(name) == 0 &&
            programVarNames.count(name) == 0 &&
            ignoredCsvVars.count(name) == 0) {

            cerr << "[restart.csv] Warning: column exists in CSV but is not used by the program: "
                 << name << endl;
        }
    }

    // Warn if variables expected by this program do not exist in restart.csv.
    for (const auto &var : programVars) {
        if (!has(var.name)) {
            cerr << "[restart.csv] Warning: variable exists in the program but not in CSV: "
                 << var.name << endl;
        }
    }

    // Helper function to get a string value from the current row.
    auto get_string = [&](const vector<string> &row, const string &name, int lineNumber) -> string {
        auto it = col.find(name);

        if (it == col.end()) {
            cerr << "[restart.csv] Internal error: missing column: " << name << endl;
            abort();
        }

        int idx = it->second;

        if (idx < 0 || idx >= static_cast<int>(row.size())) {
            cerr << "[restart.csv] Malformed row at line " << lineNumber << "." << endl;
            cerr << "Column '" << name
                 << "' exists in the header, but this row does not have that field." << endl;
            abort();
        }

        return row[idx];
    };

    // Helper function to get a double value from the current row.
    auto get_double = [&](const vector<string> &row, const string &name, int lineNumber) -> double {
        try {
            return stod(get_string(row, name, lineNumber));
        }
        catch (const exception &e) {
            cerr << "[restart.csv] Failed to convert value to double." << endl;
            cerr << "Line: " << lineNumber << ", column: " << name << endl;
            cerr << "Error: " << e.what() << endl;
            abort();
        }
    };

    // Helper function to get an integer value from the current row.
    auto get_int = [&](const vector<string> &row, const string &name, int lineNumber) -> int {
        try {
            return stoi(get_string(row, name, lineNumber));
        }
        catch (const exception &e) {
            cerr << "[restart.csv] Failed to convert value to int." << endl;
            cerr << "Line: " << lineNumber << ", column: " << name << endl;
            cerr << "Error: " << e.what() << endl;
            abort();
        }
    };

    // Check whether two coordinates are inconsistent.
    // absTol controls the absolute tolerance.
    // relTol controls the tolerance relative to the coordinate scale.
    auto coordinate_mismatch = [](double a, double b) -> bool {
        const double absTol = 1.0e-12;
        const double relTol = 1.0e-10;

        double diff = fabs(a - b);
        double scale = max({1.0, fabs(a), fabs(b)});

        return diff > absTol + relTol * scale;
    };

    // found[i][j] is used to verify that all restart cells exist exactly once.
    vector<vector<bool>> found(pm.ni + 2, vector<bool>(pm.nj + 2, false));

    int lineNumber = 1;

    // Read data lines.
    while (getline(ifs, line)) {
        lineNumber++;

        if (line.empty()) {
            continue;
        }

        vector<string> strvec = split(line, ',');

        int i = get_int(strvec, "i", lineNumber);
        int j = get_int(strvec, "j", lineNumber);

        double x_restart = get_double(strvec, "x", lineNumber);
        double r_restart = get_double(strvec, "r", lineNumber);

        // Check the index range.
        if (i < 0 || i > pm.ni + 1 || j < 0 || j > pm.nj + 1) {
            cerr << "[restart.csv] Cell definition mismatch." << endl;
            cerr << "Index is out of range at line " << lineNumber << "." << endl;
            cerr << "i = " << i << ", valid range: 0 to " << pm.ni + 1 << endl;
            cerr << "j = " << j << ", valid range: 0 to " << pm.nj + 1 << endl;
            abort();
        }

        // Check duplicated cells.
        if (found[i][j]) {
            cerr << "[restart.csv] Cell definition mismatch." << endl;
            cerr << "Duplicated cell was found at line " << lineNumber << "." << endl;
            cerr << "i = " << i << ", j = " << j << endl;
            abort();
        }

        found[i][j] = true;

        // Check consistency between restart coordinates and current grid coordinates.
        if (coordinate_mismatch(gc.x[i], x_restart) || coordinate_mismatch(gc.r[j], r_restart)) {
            cerr << "[restart.csv] Cell definition mismatch." << endl;
            cerr << "Grid coordinate is inconsistent at line " << lineNumber << "." << endl;
            cerr << "i = " << i << ", j = " << j << endl;
            cerr << "restart x = " << x_restart << ", current x = " << gc.x[i] << endl;
            cerr << "restart r = " << r_restart << ", current r = " << gc.r[j] << endl;
            abort();
        }

        // Read only variables that exist in restart.csv.
        for (const auto &var : programVars) {
            if (has(var.name)) {
                double value = get_double(strvec, var.name, lineNumber);
                var.assign(i, j, value);
            }
        }
    }

    ifs.close();

    // Verify that all cells exist in restart.csv.
    for (int i = 0; i <= pm.ni + 1; i++) {
        for (int j = 0; j <= pm.nj + 1; j++) {
            if (!found[i][j]) {
                cerr << "[restart.csv] Cell definition mismatch." << endl;
                cerr << "Missing cell in restart.csv." << endl;
                cerr << "i = " << i << ", j = " << j << endl;
                abort();
            }
        }
    }

    cout << "[restart.csv] Restart data was successfully loaded." << endl;


    //Update Ex
    //------------------------------------
    for (int iblock=0;iblock<pm.n_bl-1;iblock++){
        for (int i=gx.i_flx_bl[iblock][0];i<=gx.i_flx_bl[iblock][1];i++){
            for (int j=gx.j_flx_bl[iblock][0];j<=gx.j_flx_bl[iblock][1];j++){
                gx.Ex[i][j] = - (gc.phi[i][j]-gc.phi[i-1][j])/pm.dx;
            }
        }
    }
    //------------------------------------

    //calculate electric field
    //=====================================================================
    //left-wall-BC (z0) potential = 0
    //------------------------------------
    for (int j=gc.j_flc_bl[0][0];j<=gc.j_flc_bl[0][1];j++){
        int i=gc.i_flc_bl[0][0];
        gc.phi[i-1][j] = - gc.phi[i][j];
        gx.Ex[i][j] =  - (gc.phi[i][j]-gc.phi[i-1][j])/pm.dx;
    }
    //------------------------------------

    //left-wall-BC (z1) potential = 0
    //------------------------------------
    for (int j=gc.j_flc_bl[2][0];j<=gc.j_flc_bl[0][0]-1;j++){
        int i=gc.i_flc_bl[2][0];
        gc.phi[i-1][j] = - gc.phi[i][j];
        gx.Ex[i][j] =  - (gc.phi[i][j]-gc.phi[i-1][j])/pm.dx;
    }
    //------------------------------------

    //left-wall-BC (z2) potential = 0
    //------------------------------------
    for (int j=gc.j_flc_bl[3][0];j<=gc.j_flc_bl[2][0]-1;j++){
        int i=gc.i_flc_bl[3][0];
        gc.phi[i-1][j] = - gc.phi[i][j];
        gx.Ex[i][j] =  - (gc.phi[i][j]-gc.phi[i-1][j])/pm.dx;
    }
    //------------------------------------

    //left-wall-BC (z4) potential = 0
    //------------------------------------
    for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[4][1];j++){
        int i=gc.i_flc_bl[4][0];
        gc.phi[i-1][j] = - gc.phi[i][j];
        gx.Ex[i][j] =  - (gc.phi[i][j]-gc.phi[i-1][j])/pm.dx;
    }
    //------------------------------------

    //right-wall-BC (z3) potential = 0
    //------------------------------------
    for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[1][1];j++){
        int i=gc.i_flc_bl[1][1];
        gc.phi[i+1][j] = -gc.phi[i][j];
        gx.Ex[i+1][j] =  - (gc.phi[i+1][j]-gc.phi[i][j])/pm.dx;
    }
    //------------------------------------
    
    //right-wall-BC (z5) potential = V_bias
    //------------------------------------
    for (int j=gc.j_flc_bl[4][0];j<=gc.j_flc_bl[4][1];j++){
        int i=gc.i_flc_bl[4][1];
        gc.phi[i+1][j] = 2.0*pm.V_bias - gc.phi[i][j];
        gx.Ex[i+1][j] =  - (gc.phi[i+1][j]-gc.phi[i][j])/pm.dx;
    }
    //------------------------------------

    //Update Er
    //------------------------------------
    for (int iblock=0;iblock<pm.n_bl-1;iblock++){
        for (int i=gr.i_flr_bl[iblock][0];i<=gr.i_flr_bl[iblock][1];i++){
            for (int j=gr.j_flr_bl[iblock][0];j<=gr.j_flr_bl[iblock][1];j++){
                gr.Er[i][j] = - (gc.phi[i][j]-gc.phi[i][j-1])/pm.dr;
            }
        }
    }
    //------------------------------------

    //lower-wall-BC (x0) potential = 0
    //------------------------------------
    for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[0][1];i++){
        int j=gc.j_flc_bl[0][0];
        gc.phi[i][j-1] = - gc.phi[i][j];
        gr.Er[i][j] = - (gc.phi[i][j]-gc.phi[i][j-1])/pm.dr;
    }
    //------------------------------------

    //lower-wall-BC (x2) potential = 0
    //------------------------------------
    for (int i=gc.i_flc_bl[2][0];i<=gc.i_flc_bl[2][1];i++){
        int j=gc.j_flc_bl[2][0];
        gc.phi[i][j-1] = - gc.phi[i][j];
        gr.Er[i][j] = - (gc.phi[i][j]-gc.phi[i][j-1])/pm.dr;
    }
    //------------------------------------

    //centerline-normal-BC (x6) Er = 0
    //------------------------------------
    for (int i=gc.i_flc_bl[3][0];i<=gc.i_flc_bl[4][1];i++){ 
        int j = gc.j_flc_bl[3][0];
        gr.Er[i][j] = 0.0;
    }
    //------------------------------------

    //upper-wall-BC (x1) potential = 0
    //------------------------------------
    for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[1][1];i++){
        int j=gc.j_flc_bl[0][1];
        gc.phi[i][j+1] = - gc.phi[i][j];
        gr.Er[i][j+1] = - (gc.phi[i][j+1]-gc.phi[i][j])/pm.dr;
    }
    //------------------------------------

    //upper-wall-BC (x4) potential = 0
    //------------------------------------
    for (int i=gc.i_flc_bl[1][1]+1;i<=gc.i_flc_bl[3][1];i++){
        int j=gc.j_flc_bl[3][1];
        gc.phi[i][j+1] = - gc.phi[i][j];
        gr.Er[i][j+1] = - (gc.phi[i][j+1]-gc.phi[i][j])/pm.dr;
    }
    //------------------------------------

    //upper-wall-BC (x5) potential = 0
    //------------------------------------
    for (int i=gc.i_flc_bl[4][0];i<=gc.i_flc_bl[4][1];i++){
        int j=gc.j_flc_bl[4][1];
        gc.phi[i][j+1] = gc.phi[i][j];
        gr.Er[i][j+1] = - (gc.phi[i][j+1]-gc.phi[i][j])/pm.dr;
    }
    //------------------------------------
    
    pm.icon_restart = 1;

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
//**           void input_restart_data()                         **
//**                                                             **
//*****************************************************************
void input_restart_data_old(Params &pm, GridCenter &gc, GridInterfaceX &gx, GridInterfaceR &gr, string inputFileName)
{

    //****************** 磁場データ読み込み ******************
    //インプットファイル指定
    //ifstream ifs("restart/restart.csv");
    ifstream ifs("restart/" + inputFileName);

    if (!ifs) {
        cout << "[restart.csv] File not found." << endl;
        return;
    }

    string line;
    //vector<double> input(10,0.0);

    //1行目読む
    getline(ifs, line);
    vector<string> strvec = split(line, ',');

    //2行目以降読む
    double delta = 0.0;
    while (getline(ifs, line)) {
        
        strvec = split(line, ',');
        int i        = stoi(strvec[0]);
        int j        = stoi(strvec[1]);

        gc.rhoi[i][j]   = stod(strvec[4]);
        gc.rhoi_old[i][j]   = stod(strvec[5]);
        gx.Uix[i][j]    = stod(strvec[6]);
        gx.Uix_old[i][j]    = stod(strvec[7]);
        gr.Uir[i][j]    = stod(strvec[8]);
        gr.Uir_old[i][j]    = stod(strvec[9]);
        gc.Uip[i][j]    = stod(strvec[10]);
        gc.Uip_old[i][j]    = stod(strvec[11]);
        gc.rhoe[i][j]   = stod(strvec[12]);
        gc.rhoe_old[i][j]   = stod(strvec[13]);
        gx.rhoUex[i][j] = stod(strvec[14]);
        gx.rhoUex_old[i][j] = stod(strvec[15]);
        gr.rhoUer[i][j] = stod(strvec[16]);
        gr.rhoUer_old[i][j] = stod(strvec[17]);
        gc.Te[i][j] = stod(strvec[18])*ph::e0/ph::Boltz;
        gc.Te_old[i][j] = stod(strvec[19])*ph::e0/ph::Boltz;
        gx.Gx[i][j]     = stod(strvec[20]);
        gx.Gx_old[i][j]     = stod(strvec[21]);
        gr.Gr[i][j]     = stod(strvec[22]);
        gr.Gr_old[i][j]     = stod(strvec[23]);
        gc.phi[i][j]    = stod(strvec[24]);
        gc.phi_old[i][j]    = stod(strvec[25]);
        gc.rhom[i][j]   = stod(strvec[26]);
        gc.rhom_old[i][j]   = stod(strvec[27]);
        gx.rhoUmx[i][j]   = stod(strvec[28]);
        gr.rhoUmr[i][j]   = stod(strvec[29]);
        gc.rate_ionize[i][j]   = stod(strvec[30]);

        if(strvec.size() >= 33){
            //cout << "Bingo" << endl;
            gx.scx[i][j]   = stod(strvec[31]);
            gr.scr[i][j]   = stod(strvec[32]);
        }

        if(strvec.size() >= 37){
            //cout << "Bingo" << endl;
            gx.nUex[i][j]     = stod(strvec[33]);
            gx.nUex_old[i][j] = stod(strvec[34]);
            gr.nUer[i][j]     = stod(strvec[35]);
            gr.nUer_old[i][j] = stod(strvec[36]);
        }

        if(strvec.size() >= 41){
            //cout << "Bingo" << endl;
            gc.rhon[i][j]     = stod(strvec[37]);
            gc.rhon_old[i][j] = stod(strvec[38]);
            gx.rhoUnx[i][j]   = stod(strvec[39]);
            gr.rhoUnr[i][j]   = stod(strvec[40]);
        }

        if(strvec.size() >= 43){
            //cout << "Bingo" << endl;
            gc.rhoeps[i][j] = stod(strvec[41]);
            gc.rhoeps_old[i][j] = stod(strvec[42]);
        }else{
            gc.rhoeps[i][j] = 3.0/2.0*gc.rhoe[i][j]*ph::Boltz*gc.Te[i][j];
            gc.rhoeps_old[i][j] = 3.0/2.0*gc.rhoe_old[i][j]*ph::Boltz*gc.Te_old[i][j];
        }
    }

    ifs.close();

    //Update Ex
    //------------------------------------
    for (int iblock=0;iblock<pm.n_bl-1;iblock++){
        for (int i=gx.i_flx_bl[iblock][0];i<=gx.i_flx_bl[iblock][1];i++){
            for (int j=gx.j_flx_bl[iblock][0];j<=gx.j_flx_bl[iblock][1];j++){
                gx.Ex[i][j] = - (gc.phi[i][j]-gc.phi[i-1][j])/pm.dx;
            }
        }
    }
    //------------------------------------

    //calculate electric field
    //=====================================================================
    //left-wall-BC (z0) potential = 0
    //------------------------------------
    for (int j=gc.j_flc_bl[0][0];j<=gc.j_flc_bl[0][1];j++){
        int i=gc.i_flc_bl[0][0];
        gc.phi[i-1][j] = - gc.phi[i][j];
        gx.Ex[i][j] =  - (gc.phi[i][j]-gc.phi[i-1][j])/pm.dx;
    }
    //------------------------------------

    //left-wall-BC (z1) potential = 0
    //------------------------------------
    for (int j=gc.j_flc_bl[2][0];j<=gc.j_flc_bl[0][0]-1;j++){
        int i=gc.i_flc_bl[2][0];
        gc.phi[i-1][j] = - gc.phi[i][j];
        gx.Ex[i][j] =  - (gc.phi[i][j]-gc.phi[i-1][j])/pm.dx;
    }
    //------------------------------------

    //left-wall-BC (z2) potential = 0
    //------------------------------------
    for (int j=gc.j_flc_bl[3][0];j<=gc.j_flc_bl[2][0]-1;j++){
        int i=gc.i_flc_bl[3][0];
        gc.phi[i-1][j] = - gc.phi[i][j];
        gx.Ex[i][j] =  - (gc.phi[i][j]-gc.phi[i-1][j])/pm.dx;
    }
    //------------------------------------

    //left-wall-BC (z4) potential = 0
    //------------------------------------
    for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[4][1];j++){
        int i=gc.i_flc_bl[4][0];
        gc.phi[i-1][j] = - gc.phi[i][j];
        gx.Ex[i][j] =  - (gc.phi[i][j]-gc.phi[i-1][j])/pm.dx;
    }
    //------------------------------------

    //right-wall-BC (z3) potential = 0
    //------------------------------------
    for (int j=gc.j_flc_bl[1][0];j<=gc.j_flc_bl[1][1];j++){
        int i=gc.i_flc_bl[1][1];
        gc.phi[i+1][j] = -gc.phi[i][j];
        gx.Ex[i+1][j] =  - (gc.phi[i+1][j]-gc.phi[i][j])/pm.dx;
    }
    //------------------------------------
    
    //right-wall-BC (z5) potential = V_bias
    //------------------------------------
    for (int j=gc.j_flc_bl[4][0];j<=gc.j_flc_bl[4][1];j++){
        int i=gc.i_flc_bl[4][1];
        gc.phi[i+1][j] = 2.0*pm.V_bias - gc.phi[i][j];
        gx.Ex[i+1][j] =  - (gc.phi[i+1][j]-gc.phi[i][j])/pm.dx;
    }
    //------------------------------------

    //Update Er
    //------------------------------------
    for (int iblock=0;iblock<pm.n_bl-1;iblock++){
        for (int i=gr.i_flr_bl[iblock][0];i<=gr.i_flr_bl[iblock][1];i++){
            for (int j=gr.j_flr_bl[iblock][0];j<=gr.j_flr_bl[iblock][1];j++){
                gr.Er[i][j] = - (gc.phi[i][j]-gc.phi[i][j-1])/pm.dr;
            }
        }
    }
    //------------------------------------

    //lower-wall-BC (x0) potential = 0
    //------------------------------------
    for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[0][1];i++){
        int j=gc.j_flc_bl[0][0];
        gc.phi[i][j-1] = - gc.phi[i][j];
        gr.Er[i][j] = - (gc.phi[i][j]-gc.phi[i][j-1])/pm.dr;
    }
    //------------------------------------

    //lower-wall-BC (x2) potential = 0
    //------------------------------------
    for (int i=gc.i_flc_bl[2][0];i<=gc.i_flc_bl[2][1];i++){
        int j=gc.j_flc_bl[2][0];
        gc.phi[i][j-1] = - gc.phi[i][j];
        gr.Er[i][j] = - (gc.phi[i][j]-gc.phi[i][j-1])/pm.dr;
    }
    //------------------------------------

    //centerline-normal-BC (x6) Er = 0
    //------------------------------------
    for (int i=gc.i_flc_bl[3][0];i<=gc.i_flc_bl[4][1];i++){ 
        int j = gc.j_flc_bl[3][0];
        gr.Er[i][j] = 0.0;
    }
    //------------------------------------

    //upper-wall-BC (x1) potential = 0
    //------------------------------------
    for (int i=gc.i_flc_bl[0][0];i<=gc.i_flc_bl[1][1];i++){
        int j=gc.j_flc_bl[0][1];
        gc.phi[i][j+1] = - gc.phi[i][j];
        gr.Er[i][j+1] = - (gc.phi[i][j+1]-gc.phi[i][j])/pm.dr;
    }
    //------------------------------------

    //upper-wall-BC (x4) potential = 0
    //------------------------------------
    for (int i=gc.i_flc_bl[1][1]+1;i<=gc.i_flc_bl[3][1];i++){
        int j=gc.j_flc_bl[3][1];
        gc.phi[i][j+1] = - gc.phi[i][j];
        gr.Er[i][j+1] = - (gc.phi[i][j+1]-gc.phi[i][j])/pm.dr;
    }
    //------------------------------------

    //upper-wall-BC (x5) potential = 0
    //------------------------------------
    for (int i=gc.i_flc_bl[4][0];i<=gc.i_flc_bl[4][1];i++){
        int j=gc.j_flc_bl[4][1];
        gc.phi[i][j+1] = gc.phi[i][j];
        gr.Er[i][j+1] = - (gc.phi[i][j+1]-gc.phi[i][j])/pm.dr;
    }
    //------------------------------------

    cout << "[restart.csv] Imported" << endl;
    pm.icon_restart = 1;

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
void InputModule::input_BOLSIG_data(Params &pm, BolsigVec &bo, string inputFileName)
{

    //****************** データ読み込み ******************
    //電子衝突レートインプットファイル指定
    //ifstream ifs("input/rateCoef_e.csv");
    ifstream ifs("input/" + inputFileName);

    if (!ifs) {
        cout << "[rateCoef_e.csv] Error! File not found." << endl;
        abort();
    }

    string line;
    //vector<double> input(10,0.0);

    //1行目読む
    getline(ifs, line);
    vector<string> strvec = split(line, ',');
    int nSp = strvec.size()-1; //衝突の種類数

    vector<int> col_type; 
    for (int iSp=0;iSp<nSp;iSp++){
        //cout << strvec[iSp+1] << endl;
        if(strvec[iSp+1] == string("Elastic")){
            col_type.push_back(0);
        }else if(strvec[iSp+1] == string("Excitation")){
            col_type.push_back(1);
        }else if(strvec[iSp+1] == string("Excitation_reso")){
            col_type.push_back(2);
        }else if(strvec[iSp+1] == string("Excitation_meta")){
            col_type.push_back(3);
        }else if(strvec[iSp+1] == string("Ionization")){
            col_type.push_back(4);
        }else if(strvec[iSp+1] == string("Ionization_step")){
            col_type.push_back(5);
        }else if(strvec[iSp+1] == string("Excitation_step")){
            col_type.push_back(6);
        }else{
            col_type.push_back(7); //super-elastic
        }
    }

    //2行目読む
    getline(ifs, line);
    strvec = split(line, ',');

    //閾値のエネルギーを書き込み
    vector<double> threshold(nSp,0.0); 
    for (int iSp=0;iSp<nSp;iSp++){
        threshold[iSp] = stod(strvec[iSp+1])*ph::e0; //eV -> J
        //cout << "iSp = " << iSp << ","<<threshold[iSp]/e0 << ","<<threshold[iSp]  << endl;
    }

    vector<vector<double> > data(nSp,vector<double>());
    int ncount=0;

    //3行目以降読む
    double tmp0 = 0;
    double tmp1 = 0;

    while (getline(ifs, line)) {
        
        strvec = split(line, ',');

        //cout <<strvec[0]<< endl;
        for (int iSp=0;iSp<nSp;iSp++){
            //cout << ncount << ","<<iSp << endl;
            data[iSp].push_back(stod(strvec[iSp+1]));
        }

        if(ncount == 0){
            tmp0 = stod(strvec[0]);
        }else if(ncount == 1){
            tmp1 = stod(strvec[0]);
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
        //cout << col_name[iSp] << endl;
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

    cout<<"["<<nElas << " Elastic & " 
    << nExcReso << " Excitation_reso & "<<nExcMeta << " Excitation_meta & "<<nExc << " Excitation_upper & " 
    << nIonz << " Ionization & " 
    << nIonzStep << " Ionization_step " << nExcStep << " Excitation_step "  << nSuper << " Superelastic " 
    <<  "data are input to simulation]"<<endl;
    cout<< " Total data num = " << ndata_BSG << ", delta = " << delta << " eV, Maximum = "<< delta*(ndata_BSG-1) << " eV"<< endl;
    cout<<endl;
    ifs.close();

    //イオン衝突レートインプットファイル指定
    ifstream ifs1("input/rateCoef_i.csv");

    if (!ifs1) {
        cout << "[rateCoef_i.csv] Error! File not found." << endl;
        abort();
    }

    string line1;
    //vector<double> input(10,0.0);

    //1行目読む
    getline(ifs1, line1);
    vector<string> strvec1 = split(line1, ',');
    int nSp1 = strvec1.size()-1; //衝突の種類数

    //2行目読む
    getline(ifs1, line1);
    strvec1 = split(line1, ',');


    //閾値のエネルギーを書き込み
    vector<double> threshold1(nSp1,0.0); 
    for (int iSp=0;iSp<nSp1;iSp++){
        threshold1[iSp] = stod(strvec1[iSp+1])*ph::e0; //eV -> J
        //cout << "iSp = " << iSp << ","<<threshold[iSp]/e0 << ","<<threshold[iSp]  << endl;
    }

    vector<vector<double> > data1(nSp1,vector<double>());
    ncount=0;

    //3行目以降読む
    
    while (getline(ifs1, line1)) {
        
        strvec1 = split(line1, ',');

        //cout <<strvec1[0]<< endl;
        for (int iSp=0;iSp<nSp1;iSp++){
            data1[iSp].push_back(stod(strvec1[iSp+1]));
        }

        if(ncount == 0){
            tmp0 = stod(strvec1[0]);
            //cout << "delta = " << delta << endl;
        }else if(ncount == 1){
            tmp1 = stod(strvec1[0]);
            //cout << "delta = " << delta << endl;
        }

        ncount = ncount + 1;
    }
    int ndata_File = ncount;
    double delta1 = tmp1-tmp0;

    cout<<"["<<nSp1 << " Ion collision data are input to simulation]"<<endl;
    cout<< " Total data num = " << ndata_File << ", delta = " << delta1 << " eV, Maximum = "<< delta1*(ndata_File-1) << " eV"<< endl;
    cout<<endl;


    //****************** プログラム内で扱う配列に変換する ******************
    //電子側
    bo.nu_elas_N.resize(ndata_BSG,0.0); //弾性周波数/nn <σelas*vth>
    bo.nu_exc_N.resize(ndata_BSG,0.0); //励起周波数/nn <σexc*vth>
    bo.nu_excReso_N.resize(ndata_BSG,0.0); //励起周波数/nn <σexc*vth>
    bo.nu_excMeta_N.resize(ndata_BSG,0.0); //励起周波数/nn <σexc*vth>
    bo.nu_ionz_N.resize(ndata_BSG,0.0); //イオン化周波数/nn <σionz*vth>
    bo.nu_ionzStep_N.resize(ndata_BSG,0.0); //イオン化周波数/nn <σionz*vth>
    bo.nu_excStep_N.resize(ndata_BSG,0.0); //段階電離周波数/nn <σionz*vth>
    bo.nu_super_N.resize(ndata_BSG,0.0); //超弾性周波数/nn <σsuper*vth>
    bo.rate_eloss_n_N.resize(ndata_BSG,0.0); //エネルギーロスレート/nn
    bo.rate_eloss_m_N.resize(ndata_BSG,0.0); //エネルギーロスレート/nm

    //cout << "OK1" << endl;

    for (int i=0;i<ndata_BSG;i++){
        for (int iSp=0;iSp<nSp;iSp++){
            //cout << i << ","<<iSp << endl;
            if(col_type[iSp] == 0){ //elastic
                bo.nu_elas_N[i] = bo.nu_elas_N[i] + data[iSp][i];
                bo.rate_eloss_n_N[i] = bo.rate_eloss_n_N[i] +  data[iSp][i]*threshold[iSp];
            }else if(col_type[iSp] == 1){ //excitation
                bo.nu_exc_N[i] = bo.nu_exc_N[i] + data[iSp][i];
                bo.rate_eloss_n_N[i] = bo.rate_eloss_n_N[i] +  data[iSp][i]*threshold[iSp];
            }else if(col_type[iSp] == 2){ //resonance-excitation
                bo.nu_excReso_N[i] = bo.nu_excReso_N[i] + data[iSp][i];
                bo.rate_eloss_n_N[i] = bo.rate_eloss_n_N[i] +  data[iSp][i]*threshold[iSp];
            }else if(col_type[iSp] == 3){ //metastable-excitation
                bo.nu_excMeta_N[i] = bo.nu_excMeta_N[i] + data[iSp][i];
                bo.rate_eloss_n_N[i] = bo.rate_eloss_n_N[i] +  data[iSp][i]*threshold[iSp];
            }else if(col_type[iSp] == 4){ //ionization
                bo.nu_ionz_N[i] = bo.nu_ionz_N[i] + data[iSp][i];
                bo.rate_eloss_n_N[i] = bo.rate_eloss_n_N[i] +  data[iSp][i]*threshold[iSp];
            }else if(col_type[iSp] == 5){ //stepwise ionization
                bo.nu_ionzStep_N[i] = bo.nu_ionzStep_N[i] + data[iSp][i];
                bo.rate_eloss_m_N[i] = bo.rate_eloss_m_N[i] +  data[iSp][i]*threshold[iSp];
            }else if(col_type[iSp] == 6){ //stepwise excitation
                bo.nu_excStep_N[i] = bo.nu_excStep_N[i] + data[iSp][i];
                bo.rate_eloss_m_N[i] = bo.rate_eloss_m_N[i] +  data[iSp][i]*threshold[iSp];
            }else if(col_type[iSp] == 7){ //superelastic
                bo.nu_super_N[i] = bo.nu_super_N[i] + data[iSp][i];
                bo.rate_eloss_m_N[i] = bo.rate_eloss_m_N[i] +  data[iSp][i]*threshold[iSp];
            }
        }
    }

    //イオン
    bo.nui_m_N.resize(ndata_File,0.0); //イオン衝突周波数/nn <σexc*vth>
    for (int i=0;i<ndata_File;i++){
        for (int iSp=0;iSp<nSp1;iSp++){
            bo.nui_m_N[i] = bo.nui_m_N[i] + data1[iSp][i];
        }
    }

    //dTe,dEi計算
    pm.dTe = delta*2.0/3.0*ph::e0/ph::Boltz;
    pm.dEi = delta1*ph::e0;

}

