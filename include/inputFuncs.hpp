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
class InputFuncs
{
    private:

    public:
        void inputParam(Params &prm,string inputFileName);


};



//*****************************************************************
//**                                                             **
//**           void param                                        **
//**                                                             **
//*****************************************************************
void InputFuncs::inputParam(Params &prm, string inputFileName)
{

    const std::string filename = "input/"+inputFileName;

    YAML::Node setup;

    try {
        setup = YAML::LoadFile(filename);
    }
    catch (const YAML::BadFile&) {
        std::cout << "[" << filename << "] Error! File not found." << std::endl;
        abort();
    }

    std::vector<std::string> missingParams;

    auto readDouble = [&](double& dst,
                          const std::string& category,
                          const std::string& key) {
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
        catch (const std::exception& e) {
            throw std::runtime_error(
                "Type conversion error: " + category + "." + key + " should be double"
            );
        }
    };

    auto readInt = [&](int& dst,
                       const std::string& category,
                       const std::string& key) {
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
        catch (const std::exception& e) {
            throw std::runtime_error(
                "Type conversion error: " + category + "." + key + " should be int"
            );
        }
    };

#define READ_DOUBLE(category, var) readDouble((var), (category), #var)
#define READ_INT(category, var)    readInt((var), (category), #var)

    try {
        std::string category;

        category = "microwave";
        READ_DOUBLE(category, prm.Pmw);
        READ_INT(category,    prm.icon_mwRef);
        READ_INT(category,    prm.icon_impTest);
        READ_DOUBLE(category, prm.omegam);
        READ_DOUBLE(category, prm.nu_eff);
        READ_DOUBLE(category, prm.deltaECR);
        READ_DOUBLE(category, prm.S11_mag);     // (double) Magnitude of S11
        READ_DOUBLE(category, prm.S11_arg_deg); // (double) Argument of S11 (deg)
        READ_DOUBLE(category, prm.S21_mag);     // (double) Magnitude of S21
        READ_DOUBLE(category, prm.S21_arg_deg); // (double) Argument of S21 (deg)
        READ_DOUBLE(category, prm.S12_mag);     // (double) Magnitude of S12
        READ_DOUBLE(category, prm.S12_arg_deg); // (double) Argument of S12 (deg)
        READ_DOUBLE(category, prm.S22_mag);     // (double) Magnitude of S22
        READ_DOUBLE(category, prm.S22_arg_deg); // (double) Argument of S22 (deg)
        READ_DOUBLE(category, prm.Z0_base); // (double) Argument of S22 (deg)

        category = "plasma";
        READ_DOUBLE(category, prm.Ti);
        READ_DOUBLE(category, prm.Tn);
        READ_DOUBLE(category, prm.rhon_ini);
        READ_DOUBLE(category, prm.DmN);
        READ_DOUBLE(category, prm.Te_rep_eV);

        category = "material";
        READ_DOUBLE(category, prm.epsr_diele);

        category = "bias";
        READ_DOUBLE(category, prm.V_bias);

        category = "neutral_inlet";
        READ_DOUBLE(category, prm.Q_neutIn_mgs);
        READ_DOUBLE(category, prm.width_neutIn);

        category = "transport_model";
        READ_INT(category,    prm.icon_Bohm);
        READ_INT(category,    prm.icon_Sagdeev);
        READ_DOUBLE(category, prm.alpha_Bohm);
        READ_DOUBLE(category, prm.scale_inertia);

        category = "see";
        READ_DOUBLE(category, prm.coefIISEE_ts);
        READ_DOUBLE(category, prm.coefMISEE_ts);
        READ_DOUBLE(category, prm.Te_emitSEE_eV);
        READ_DOUBLE(category, prm.ratioEngy_EISEE_rd);

        category = "scheme";
        READ_INT(category,    prm.icon_PC);
        READ_INT(category,    prm.icon_inertia);
        READ_INT(category,    prm.icon_adp_dt);
        READ_DOUBLE(category, prm.dt_ini);
        READ_INT(category,    prm.ndt_i);
        READ_INT(category,    prm.ndt_m);
        READ_INT(category,    prm.ndt_n);
        READ_DOUBLE(category, prm.CFL);

        category = "solver";
        READ_DOUBLE(category, prm.error_cnv_SOR_Ui);
        READ_INT(category,    prm.maxITR_SOR_Ui);
        READ_INT(category,    prm.icon_iter_Ui);

        READ_DOUBLE(category, prm.error_cnv_SOR_rhoi);
        READ_INT(category,    prm.maxITR_SOR_rhoi);
        READ_INT(category,    prm.icon_iter_rhoi);

        READ_DOUBLE(category, prm.error_cnv_SOR_phi);
        READ_INT(category,    prm.maxITR_SOR_phi);
        READ_INT(category,    prm.icon_iter_phi);
        READ_DOUBLE(category, prm. error_cnv_HES_phi);
        READ_INT(category,    prm.maxITR_HES_phi);

        READ_DOUBLE(category, prm.error_cnv_SOR_rhoe);
        READ_INT(category,    prm.maxITR_SOR_rhoe);
        READ_INT(category,    prm.icon_iter_rhoe);
        READ_DOUBLE(category, prm.error_cnv_HES_rhoe);
        READ_INT(category,    prm.maxITR_HES_rhoe);

        READ_DOUBLE(category, prm.error_cnv_SOR_rhoeps);
        READ_INT(category,    prm.maxITR_SOR_rhoeps);
        READ_INT(category,    prm.icon_iter_rhoeps);
        READ_DOUBLE(category, prm.error_cnv_HES_rhoeps);
        READ_INT(category,    prm.maxITR_HES_rhoeps);

        category = "microwave_coupling";
        READ_INT(category,    prm.ndiv_MW);

        category = "relaxation";
        READ_DOUBLE(category, prm.beta_rhoe);
        READ_DOUBLE(category, prm.beta_rhoUe);

        category = "simulation";
        READ_INT(category,    prm.ntime);
        READ_INT(category,    prm.icon_autoFinish);

        category = "output";
        READ_INT(category,    prm.icon_chk);
        READ_INT(category,    prm.icon_gnuRes);
        READ_INT(category,    prm.ndiv_fout);
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

    if (!missingParams.empty()) {
        std::cout << std::endl;
        std::cout << "[setup.yaml] Error! Missing YAML parameters:" << std::endl;

        for (const auto& name : missingParams) {
            std::cout << "  - " << name << std::endl;
        }

        std::cout << std::endl;
        abort();
    }

    std::cout << std::endl;
    std::cout << "[All setup parameters are successfully loaded]" << std::endl;
    std::cout << std::endl;

}
    