#pragma once
#include <vector>

/******************** support.cpp ********************/
extern std::vector<std::string> split(std::string& input, char delimiter);
extern double func_Te(std::vector<double> value,double Te_tmp,double dTe_tmp);
extern double func_dTe(std::vector<double> value,double Te_tmp,double dTe_tmp);
extern double func_Ei(std::vector<double> value,double Ei_tmp,double dEi_tmp);
extern void plotHistory(FILE* gnuplot_name,std::vector<std::vector<double> > value,std::vector<int> itime_history, int nGnuMaxTimeRange);
extern void plotCurrentHistory(FILE* gnuplot_name,std::vector<std::vector<double> > value,std::vector<int> itime_history, int nGnuMaxTimeRange);
extern void reconst0(std::vector<std::vector<double> > value
    ,std::vector<std::vector<double> > &value_Lx,std::vector<std::vector<double> > &value_Rx
    ,std::vector<std::vector<double> > &value_Lr,std::vector<std::vector<double> > &value_Rr);
extern double truncate(double c);
/******************** initial.cpp ********************/
extern void iniparameter();
extern void makeBoundary();
extern void makeProfile();

/******************** input.cpp ********************/
extern void inputParam();
extern void input_Bfield_data();
extern void input_restart_data();
extern void input_BOLSIG_data();
extern void input_SEE_data();
/******************** EMfield.cpp ********************/
extern void solve_Microwave();
extern void update_energy_profile();

/******************** fluid.cpp ********************/
extern void check_CFL();
extern void iniConditionMod();
extern void update_transport_coef();
extern void update_transport_coef_mod();
extern void solve_phi_couple();
extern void solve_phi_couple_PC();
extern void solve_phi_couple_wdTe_PC();
extern void solve_phi_couple_wdTe_wSEE_PC();
extern void solve_phi_couple_SG_PC();
extern void solve_phi_couple_directional();
extern void solve_phi_couple_dieleOpen();

extern void solve_Uix_tmp();
extern void solve_Uir_tmp();
extern void solve_Uip_tmp();
extern void solve_rhoi();
extern void solve_rhoi_constTe();
extern void solve_rhoi_my();
extern void solve_rhoi_kinetic();
extern void update_Ui();
extern void update_Ui_constTe();
extern void update_Ui_my();

extern void solve_rhoe_PC();
extern void solve_rhoe_wdTe_PC();
extern void solve_rhoe_wdTe_wSEE_PC();
extern void solve_rhoe_wdTe_wSEE_wInertia_PC();
extern void solve_rhoe_Te_PC();
extern void solve_rhoe_PC_pDif();
extern void solve_rhoe_TV();
extern void solve_rhoe_direct();
extern void solve_rhoe_wdTe_direct();
extern void solve_rhoe_wdTe_wSEE_direct();
extern void solve_rhoe_direct_dieleOpen();

extern void solve_Te_direct();
extern void solve_Te_wdTe_direct();
extern void solve_Te_wdTe_wSEE_direct();
extern void solve_Te_PC();
extern void solve_Te_wdTe_PC();
extern void solve_Te_wdTe_wSEE_PC();
extern void solve_Te_wdTe_wSEE_wInertia_PC();

extern void update_rhom_diff();
extern void update_rhom_explicit();
extern void update_rhom();
extern void update_rhom_kinetic();

extern void reconstAll_for_neut();
extern void riemann_FVS_for_neut();
extern void update_neut();
extern void update_rhon();
extern void update_rhon_log();


/******************** solver.cpp ********************/
extern void solver_SMG(std::vector<std::vector<double> > aP
    ,std::vector<std::vector<double> > aE
    ,std::vector<std::vector<double> > aW
    ,std::vector<std::vector<double> > aN
    ,std::vector<std::vector<double> > aS
    ,std::vector<std::vector<double> > b
    ,std::vector<std::vector<int> > i_bl
    ,std::vector<std::vector<int> > j_bl
    ,int maxITR,double error_cnv,int icon_msg
    ,std::vector<std::vector<double> > &solution
    );

extern void solver_SOR(std::vector<std::vector<double> > aP
    ,std::vector<std::vector<double> > aE
    ,std::vector<std::vector<double> > aW
    ,std::vector<std::vector<double> > aN
    ,std::vector<std::vector<double> > aS
    ,std::vector<std::vector<double> > b
    ,std::vector<std::vector<int> > i_bl
    ,std::vector<std::vector<int> > j_bl
    ,double alpha_SOR,int maxITR,double error_cnv,int icon_msg
    ,std::vector<std::vector<double> > &solution
    );

extern void solver_GMRES_9p(std::vector<std::vector<double> > aP
    ,std::vector<std::vector<double> > aE
    ,std::vector<std::vector<double> > aW
    ,std::vector<std::vector<double> > aN
    ,std::vector<std::vector<double> > aS
    ,std::vector<std::vector<double> > aEN
    ,std::vector<std::vector<double> > aES
    ,std::vector<std::vector<double> > aWN
    ,std::vector<std::vector<double> > aWS
    ,std::vector<std::vector<double> > b
    ,std::vector<std::vector<int> > i_bl
    ,std::vector<std::vector<int> > j_bl
    ,int maxITR,double error_cnv,int icon_msg
    ,std::vector<std::vector<double> > &solution
    );

extern void solver_LU_9p(std::vector<std::vector<double> > aP
    ,std::vector<std::vector<double> > aE
    ,std::vector<std::vector<double> > aW
    ,std::vector<std::vector<double> > aN
    ,std::vector<std::vector<double> > aS
    ,std::vector<std::vector<double> > aEN
    ,std::vector<std::vector<double> > aES
    ,std::vector<std::vector<double> > aWN
    ,std::vector<std::vector<double> > aWS
    ,std::vector<std::vector<double> > b
    ,std::vector<std::vector<int> > i_bl
    ,std::vector<std::vector<int> > j_bl
    ,int icon_msg
    ,std::vector<std::vector<double> > &solution
    );

/******************** output.cpp ********************/
extern void output_phase();
extern void output_simple();
extern void output();