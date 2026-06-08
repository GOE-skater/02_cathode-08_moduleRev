#pragma once
#include <complex>

using namespace std;

struct Params{ 

    // ============================================================
    // Read from input file
    // ============================================================
    
    //domain
    //--------------------------------
    int ni = 194; //x方向のセルの数 (実際はゴーストセルが+2) //155 + 39
    int nj = 70; //x方向のセルの数 (実際はゴーストセルが+2) //100
    double xL= 0.00; //計算領域x左端
    double xR= 0.0388; //計算領域x右端 0.031 + 0.0078
    double rmin= 0.0; //計算領域r左端
    double rmax= 0.014; //計算領域r右端 2.0e-3 0.014
    double width_neutIn = 0.0; //中性粒子のインレットの幅 m
    //--------------------------------

    //operating_condition
    //--------------------------------
    double V_bias = 0.0; //バイアス電圧
    double Q_neutIn_mgs = 0.0; //中性粒子の流量 mg/s
    double fmw = 2.45e9; //マイクロ波周波数 (Hz)
    double Pmw = 8.0; //マイクロ波電力 (W)
    int flag_mwRef = 1; //反射を解くか
    //--------------------------------

    //physical_parameter
    //--------------------------------
    double masse=9.10938356e-31; //電子質量
    double massi = 131.293e-3/6.02214076e23; //イオン質量 Xe 131.293
    double ri = 216e-12; //イオン半径 Xe = 216pm
    double Ti = 300; //イオン温度
    double Tn = 300; //中性粒子温度
    double DmN = 1.66e20; //準安定種の拡散係数 (Elis 1969) 1/ms
    double epsr_diele = 1.0; //誘電体の比電率
    //--------------------------------

    //transmission_line
    //--------------------------------
    double S11_mag = 0.0;     // (double) Magnitude of S11
    double S11_arg_deg = 0.0; // (double) Argument of S11 (deg)
    double S21_mag = 1.0;     // (double) Magnitude of S21
    double S21_arg_deg = 0.0; // (double) Argument of S21 (deg)
    double S12_mag = 1.0;     // (double) Magnitude of S12
    double S12_arg_deg = 0.0; // (double) Argument of S12 (deg)
    double S22_mag = 0.0;     // (double) Magnitude of S22
    double S22_arg_deg = 0.0; // (double) Argument of S22 (deg)
    double Z0_base = 0.0;     // (double) Base impedance (Ohm)
    //--------------------------------
    
    //secondary_electron_emission
    //--------------------------------
    double coefIISEE_ts = 0.08; //II true-secondary
    double coefMISEE_ts = 0.08; //MI true-secondary (Auger)
    double Te_emitSEE_eV = 2.5; //electron emission temperature (eV)
    double ratioEngy_EISEE_rd = 0.5; //energy ratio of EI rediffusion 
    //--------------------------------

    //governing_equation
    //--------------------------------
    int flag_inertia = 0; //慣性を考慮するか
    double scale_inertia = 1.0; //scaling factor of inertia term
    //--------------------------------

    //anomalous_scattering_model
    //--------------------------------
    int flag_Sagdeev = 0; //Sagdeevモデルを考慮するか
    int flag_Bohm = 0; //Bohm拡散を考慮するか
    double alpha_Bohm = 16.0; //ボーム拡散係数
    //--------------------------------

    //collisionless_heating_model
    //--------------------------------
    double nu_eff = 1.0e9; //有効衝突周波数 = sqrt(vth*omegam/deltaB)
    //--------------------------------

    //initial_condition
    //--------------------------------
    double rhon_ini = 2.41432e19; //中性粒子密度 (初期値)
    double rhoe_ini = 5e18; //電子密度 (初期値)
    double Te_ini = 3.0*ph::e0/ph::Boltz; //電子密度 (初期値)
    double J1r_exc = 0.212687; //励起電流
    //--------------------------------

    //numerical_scheme
    //--------------------------------
    int flag_PC = 0; //PC methodを使うか (0:直接離散化，1:PC method)
    double Te_rep_eV = 3.0; //代表電子温度 (イオンスキーム用)
    //--------------------------------

    //solver_setting
    //--------------------------------
    double error_cnv_SOR_Ui = 1e-10;  //Uiの収束基準
    int maxITR_SOR_Ui = 30;           //Uiの最大反復数
    int flag_iter_Ui = 1;             //Uiの反復法の選択 0:SOR 1:SMG

    double error_cnv_SOR_rhoi = 1e-8; //rhoiの収束基準
    int maxITR_SOR_rhoi = 2000;       //rhoiの最大反復数
    int flag_iter_rhoi = 1;           //rhoiの反復法の選択 0:SOR 1:SMG

    double error_cnv_SOR_phi = 1e-8;  //phiの収束基準
    int maxITR_SOR_phi = 2000;        //phiの最大反復数
    int flag_iter_phi = 1;            //phiの反復法の選択 0:SOR 1:SMG
    double error_cnv_HES_phi = 1e-6; //phiの擬似タイムステップの収束
    int maxITR_HES_phi = 1; //phiの擬似時間の最大反復数

    double error_cnv_SOR_rhoe = 1e-6;  //rhoeの収束基準
    int maxITR_SOR_rhoe = 30;        //rhoeの最大反復数
    int flag_iter_rhoe = 1;          //rhoeの反復法の選択 0:直接法 1:GMRES (0:SOR 1:SMG)
    double error_cnv_HES_rhoe = 1e-6; //rhoeの擬似タイムステップの収束
    int maxITR_HES_rhoe = 1; //rhoeの擬似時間の最大反復数
    double beta_rhoe = 1.0; //rhoeの不足緩和係数 (慣性考慮時)
    double beta_rhoUe = 1.0; //rhoUeの不足緩和係数 (慣性考慮時)

    double error_cnv_SOR_rhoeps = 1e-6;  //rhoepsの収束基準
    int maxITR_SOR_rhoeps = 30;        //rhoepsの最大反復数
    int flag_iter_rhoeps = 1;          //rhoepsの反復法の選択 0:GMRES 0:直接法
    double error_cnv_HES_rhoeps = 1e-6; //rhoepsの擬似タイムステップの収束
    int maxITR_HES_rhoeps = 1; //rhoepsの擬似時間の最大反復数
    //--------------------------------

    //simulation_control
    //--------------------------------
    int flag_adp_dt = 0;  //アダプティブなdtを用いるか 0:固定dt, 1:可変dt
    double dt_ini = 5.0e-12; //初期CFL 誘電緩和時間 3.38e-13 5.0e-13
    double CFL = 1.0e1; //CFL数 0.001
    int ntime = 2000; //3000000
    int ndt_i = 100; //イオンの電子に対するタイムステップ倍率
    int ndt_m = 1000; //中性 (準安定) の電子に対するタイムステップ倍率
    int ndt_n = 10000; //中性 (準安定) の電子に対するタイムステップ倍率
    int ndiv_MW = 200000; ///マイクロ波の計算間隔
    int flag_error = 1; //errorの種類 0:max 1:rms
    int flag_autoFinish = 0;
    //--------------------------------
    
    //output_control
    //--------------------------------
    int flag_chk = 0; //チェックファイルを出力するか
    int flag_gnuRes = 1; //gnuplot結果表示
    int ndiv_fout = 2000;
    //--------------------------------

    //test_mode
    //--------------------------------
    int flag_impTest= 1;      // (int) conduct impedance test
    //--------------------------------
    
    // ============================================================
    // Do not read
    // ============================================================
    int n_bl = 6; //計算ブロック数

    double deltaECR = 1e-3; //ECRの領域の幅

    double omegam = 0.0; //マイクロ波角周波数

    
    double Lx=0.0; //計算領域のx方向の大きさ
    double Lr=0.0; //計算領域のx方向の大きさ
    
    double gtime=0.0; //時間
    int itime=0; //時間(整数)
    int itime_PC_phi = 0; //時間(整数)
    int itime_PC_rhoe = 0; //時間(整数)
    int itime_PC_rhoeps = 0; //時間(整数)
    
    double dx=0.0; //x方向の刻み
    double dr=0.0; //r方向の刻み
    double dt=0.0; //時間刻み
    
    double I_Anode = 0.0; //アノード電流
    double Ii_Anode = 0.0; //アノード電流
    double Ie_Anode = 0.0; //アノード電流
    
    double I_Nozzle = 0.0; //ノズル電流
    double Ii_Nozzle = 0.0; //アノード電流
    double Ie_Nozzle = 0.0; //アノード電流
    
    int nOut = 0; //ファイルアウトプット回数
    
    int flag_restart = 0;
    
    double error_rhoi = 0.0;
    double error_Uix = 0.0;
    double error_Uir = 0.0;
    double error_Uip = 0.0;

    double error_rhoe = 0.0;
    double error_rhoUex = 0.0;
    double error_rhoUer = 0.0;

    double error_rhoeps = 0.0;
    double error_Gx = 0.0;
    double error_Gr = 0.0;
    
    double error_phi = 0.0;
    double error_nUex = 0.0;
    double error_nUer = 0.0;

    double error_rhom = 0.0;
    double error_rhon = 0.0;
    
    //反応レート係数データ
    double dTe = 0.0; //電子温度の刻み幅
    double dEi = 0.0; //イオンドリフトエネルギーの刻み幅
    double dTe_SEE = 0.0; //SEE用の電子温度の刻み幅
    
    //流入条件
    double fn_In = 0.0; //中性粒子の流入フラックス
    
    //ジオメトリ
    double x1 = 0.0;
    double x2 = 0.0;
    double x3 = 0.0;
    double x4 = 0.0;
    double x5 = 0.0;
    double x6 = 0.0;
    
    double r1 = 0.0;
    double r2 = 0.0;
    double r3 = 0.0;
    double r4 = 0.0;
    double r5 = 0.0;
    double r6 = 0.0;
    
};
