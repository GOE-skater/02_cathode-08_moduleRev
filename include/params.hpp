#pragma once
#include "constants.hpp"
#include <complex>

struct Params{ 


    
    //microwave
    //--------------------------------
    double Pmw = 8.0; //マイクロ波電力 (W)
    int icon_mwRef = 1; //反射を解くか
    int icon_impTest= 1;      // (int) conduct impedance test
    double omegam = 2.0*M_PI*2.45e9; //マイクロ波周波数
    double nu_eff = 1.0e9; //有効衝突周波数 = sqrt(vth*omegam/deltaB)
    double deltaECR = 1e-3; //ECRの領域の幅
    double S11_mag = 0.0;     // (double) Magnitude of S11
    double S11_arg_deg = 0.0; // (double) Argument of S11 (deg)
    double S21_mag = 1.0;     // (double) Magnitude of S21
    double S21_arg_deg = 0.0; // (double) Argument of S21 (deg)
    double S12_mag = 1.0;     // (double) Magnitude of S12
    double S12_arg_deg = 0.0; // (double) Argument of S12 (deg)
    double S22_mag = 0.0;     // (double) Magnitude of S22
    double S22_arg_deg = 0.0; // (double) Argument of S22 (deg)
    double Z0_base = 0.0;     // (double) Base impedance (Ohm)

    //plasma
    //--------------------------------
    double Ti = 300; //イオン温度
    double Tn = 300; //中性粒子温度
    double rhon_ini = 2.41432e19; //中性粒子密度 (一定値) m^-3 @1Pa
    double DmN = 1.66e20; //準安定種の拡散係数 (Elis 1969) 1/ms
    double Te_rep_eV = 3.0; //代表電子温度 (イオンスキーム用)

    //material
    //--------------------------------
    double epsr_diele = 1.0; //誘電体の比電率
    //--------------------------------

    //bias
    //--------------------------------
    double V_bias = 0.0; //バイアス電圧
    //--------------------------------

    double masse=9.10938356e-31; //電子質量
    double massi = 131.293e-3/6.02214076e23; //イオン質量 Xe 131.293
    double ri = 216e-12; //イオン半径 Xe = 216pm
    
    
   

    
    int ni = 194; //x方向のセルの数 (実際はゴーストセルが+2) //155 + 39
    int nj = 70; //x方向のセルの数 (実際はゴーストセルが+2) //100
    int n_bl = 6; //計算ブロック数
    double xL= 0.00; //計算領域x左端
    double xR= 0.0388; //計算領域x右端 0.031 + 0.0078
    double rmin= 0.0; //計算領域r左端
    double rmax= 0.014; //計算領域r右端 2.0e-3 0.014

    std::complex<double> J1r_exc = std::complex<double>(0.212687,0.0); //励起電流
    
    
    
    
    

    
    
    double Q_neutIn_mgs = 0.0; //中性粒子の流量 mg/s
    double width_neutIn = 0.0; //中性粒子のインレットの幅 m
    
    int icon_Bohm = 0; //Bohm拡散を考慮するか
    int icon_Sagdeev = 0; //Sagdeevモデルを考慮するか
    double alpha_Bohm = 16.0; //ボーム拡散係数
    double scale_inertia = 1.0; //scaling factor of inertia term

    double coefIISEE_ts = 0.08; //II true-secondary
    double coefMISEE_ts = 0.08; //MI true-secondary (Auger)
    double Te_emitSEE_eV = 2.5; //electron emission temperature (eV)
    double ratioEngy_EISEE_rd = 0.5; //energy ratio of EI rediffusion 

    int icon_PC = 0; //PC methodを使うか (0:直接離散化，1:PC method)
    int icon_inertia = 0; //慣性を考慮するか
    int icon_adp_dt = 0;  //アダプティブなdtを用いるか 0:固定dt, 1:可変dt
    double dt_ini = 5.0e-13*10.0; //初期CFL 誘電緩和時間 3.38e-13 5.0e-13
    int ndt_i = 100; //イオンの電子に対するタイムステップ倍率
    int ndt_m = 1000; //中性 (準安定) の電子に対するタイムステップ倍率
    int ndt_n = 10000; //中性 (準安定) の電子に対するタイムステップ倍率
    double CFL = 1.0e1; //CFL数 0.001

    double error_cnv_SOR_Ui = 1e-10;  //Uiの収束基準
    int maxITR_SOR_Ui = 30;           //Uiの最大反復数
    int icon_iter_Ui = 1;             //Uiの反復法の選択 0:SOR 1:SMG

    double error_cnv_SOR_rhoi = 1e-8; //rhoiの収束基準
    int maxITR_SOR_rhoi = 2000;       //rhoiの最大反復数
    int icon_iter_rhoi = 1;           //rhoiの反復法の選択 0:SOR 1:SMG

    double error_cnv_SOR_phi = 1e-8;  //phiの収束基準
    int maxITR_SOR_phi = 2000;        //phiの最大反復数
    int icon_iter_phi = 1;            //phiの反復法の選択 0:SOR 1:SMG
    double error_cnv_HES_phi = 1e-6; //phiの擬似タイムステップの収束
    int maxITR_HES_phi = 1; //phiの擬似時間の最大反復数

    double error_cnv_SOR_rhoe = 1e-6;  //rhoeの収束基準
    int maxITR_SOR_rhoe = 30;        //rhoeの最大反復数
    int icon_iter_rhoe = 1;          //rhoeの反復法の選択 0:直接法 1:GMRES (0:SOR 1:SMG)
    double error_cnv_HES_rhoe = 1e-6; //rhoeの擬似タイムステップの収束
    int maxITR_HES_rhoe = 1; //rhoeの擬似時間の最大反復数

    double error_cnv_SOR_rhoeps = 1e-6;  //rhoepsの収束基準
    int maxITR_SOR_rhoeps = 30;        //rhoepsの最大反復数
    int icon_iter_rhoeps = 1;          //rhoepsの反復法の選択 0:GMRES 0:直接法
    double error_cnv_HES_rhoeps = 1e-6; //rhoepsの擬似タイムステップの収束
    int maxITR_HES_rhoeps = 1; //rhoepsの擬似時間の最大反復数
    int icon_error = 1; //errorの種類 0:max 1:rms
    int ndiv_MW = 200000; ///マイクロ波の計算間隔

    double beta_rhoe = 1.0; //rhoeの不足緩和係数 (慣性考慮時)
    double beta_rhoUe = 1.0; //rhoUeの不足緩和係数 (慣性考慮時)

    //全体の計算設定
    int ntime = 2000; //3000000
    int icon_autoFinish = 0;

    //出力設定
    int icon_chk = 0; //チェックファイルを出力するか
    int icon_gnuRes = 1; //gnuplot結果表示
    int ndiv_fout = 2000;


    //domain parameter
    //--------------------------------
    struct Domain{ //domain parameter
        
        //read
        //--------------------------------
        int ni = 256; //cell number in x direction
        double xMin = 0.0; //minimum-x
        double xMax = 6.7e-2; //minimum-y
        //--------------------------------

        //do not read
        //--------------------------------
        double Lx = 0.0; //計算領域のx方向の大きさ
        double dx=0.0; //x方向の刻み
        //--------------------------------
        
    };
    Domain dm;
    //--------------------------------
    
    //plasma parameter
    //--------------------------------
    struct Plasma{ 
    
        //read
        //--------------------------------
        double masse = 9.109e-31; //electron mass
        double massi = 6.67e-27; //ion masss
        double rho_ini = 1.0e16;
        double Te_ini = 10.0*phys::e0/phys::Boltz;  //initial electron temperature
        double Ti_ini = 0.095*phys::e0/phys::Boltz; //initial ion temperature
        double rhon = 2.0e19; //neutral density
        //--------------------------------

        //do not read
        //--------------------------------
        double Debye_per_dx_min = 0.0;
        //--------------------------------
    
    };
    Plasma pl;
    //--------------------------------
    
    //operation parameter
    //--------------------------------
    struct Operation{ 
    
        //read
        //--------------------------------
        double V0 = 300.0; //bias voltage
        double f0 = 13.56e6; //Hz frequency
        //--------------------------------
    
    };
    Operation ope;
    //--------------------------------
    
    //fluid model parameter
    //--------------------------------
    struct FluidModel{ 
    
        //read
        //--------------------------------
        struct flag{ 
            int eModel = 1;  //electron model 0: DD, 1: FFM
            int iModel = 0; //electron model 0: DD, 1: FFM (const Ti) 2: FFM (variable Ti)
            
            int eKineBC = 1; //電子にkinetic BCを用いるか
            int iKineBC = 1; //イオンにkinetic BCを用いるか

            int rateCoef_i = 0; ///calculation method of rate coefficient (0: mobility, 1: shifted-Maxwellian)
            int rateCoef_e = 0; //calculation method of rate coefficient (0: Maxwellian LMEA, 1: 2-term LMEA, 2: 2-term LFA)

            int rateShift_e = 0; //rate coefでshifted-Maxwellianを考慮するか

            int dTe = 1; //温度勾配を考慮するか
            int heatF = 1; //heat flux closure (0: zero-flux, 1: Fourier-law)
        };
        flag flag;
        //--------------------------------
    
    };
    FluidModel flMd;
    //--------------------------------

    //fluid scheme parameter
    //--------------------------------
    struct FluidScheme{ 
    
        //read
        //--------------------------------
        struct flag{ 
            
            int Riemann_i = 0; //Riemann solver (0: HLL & HLLC, 1:LLF, 2:GLF)
            int Riemann_e = 0; //Riemann solver (0: HLL & HLLC, 1:LLF, 2:GLF)

            int rcnst_i = 5; //space reconstruction (0:1st upwind,1:lax-wendroff,2:Fromm,3:Beam-warming 4:TVD 5:MUSCL-TVD)
            int rcnst_e = 5; //space reconstruction (0:１次風上,1:lax-wendroff,2:Fromm,3:Beam-warming 4:TVD 5:MUSCL-TVD)

            int limit_i = 0; //slope limitter function (0:minod,1:superbee,2:van Leer,3:van Albada)
            int limit_e = 0; //slope limitter function (0:minod,1:superbee,2:van Leer,3:van Albada)

            int rhoGhost_i = 0; //how to give density ghost cell (0: 0th order extraporation, 1: zero-density)
            int rhoGhost_e = 0; //how to give density ghost cell (0: 0th order extraporation, 1: zero-density)

            int RK3 = 1; //3次精度SSPルンゲクッタで時間積分するか
            //int imp = 0;  //陰解法で計算するか
        };
        flag flag;

        double muscl_k_i = 1.0/2.0; //MUSCLスキームの定数κ (flag.rcnst=4の場合のみ)
        double muscl_k_e = 1.0/2.0; //MUSCLスキームの定数κ (flag.rcnst=4の場合のみ)
        double floor_rho = 1e-100;   //floor for zero density
        //--------------------------------
    
    };
    FluidScheme flSch;
    //--------------------------------
    
    //calculation parameter
    //--------------------------------
    struct Calculation{ 
    
        //read
        //--------------------------------
        struct flag{ 
            int autoFinish = 1;
        };
        flag flag;

        int ndtPerCycle = 400; //number of timesteps per RF cycle 400
        int nCycle = 1280; //number of RF cycle 1280
        int nCycleAve = 64;
        //--------------------------------

        //do not read
        //--------------------------------
        int ntime = 56000000; //40000000
        int ntime_ave = 16000000; //3000000
        double dt=0.0; //時間刻み
        int itime=0; //時間(整数)
        double gtime=0.0; //時間
         double CFL = 0.001; //CFL数 0.001
        //--------------------------------
    };
    Calculation cl;
    //--------------------------------
    
    //output parameter
    //--------------------------------
    struct Output{ 
    
        //read
        //--------------------------------
        struct flag{ 
            int gnuRes = 1; //gnuplot結果表示
            int colRate_e = 0; //output electron collision rate
            int wallFlux_e = 0; //output electron wall flux
            //int colRate_i = 0; //output ion collision rate
        };
        flag flag;

        int nGnuOutTime = 10000; //Gnuplotの出力間隔
        int nGnuMaxTimeRange = 1000000; //Gnuplotで粒子の履歴を表示する最大幅
        int nSample_gnu = 1000; //Gnuplotで何ステップごとにサンプルするか
        int ndiv_out_console = 10000; //コンソールへのアウトプット間隔
        int ndiv_out_file = 10000; //ファイルへのアウトプット間隔
        int ndtPerPhaseBin = 1; //phase bin size (number of dt for 1 bin)
        //--------------------------------

        //do not read
        //--------------------------------
        int nOutFile = 0; //ファイルアウトプット回数
        //int ifExist_wflux = 0; //wall fluxファイルの存在チェック
        //--------------------------------

    };
    Output out;
    //--------------------------------
    
    //******************* グローバル変数 **************************
    
    
    struct Error{ //error
        
        //do not read
        //--------------------------------
        double rhoi = 0.0;
        double rhoe = 0.0;
        double rhoUix = 0.0;
        double rhoUex = 0.0;
        double rhoepse = 0.0;
        double rhon = 0.0;
        double Te = 0.0;
        double phi = 0.0;
        //--------------------------------
    };
    Error err;

    struct Bolsig{ //BOLSIG+
        
        //do not read
        //--------------------------------
        int ndata_BSG = 0; //読み込むデータ数
        double dTe = 0.0; //電子温度の刻み幅
        double dTi = 0.0; //イオン温度の刻み幅
        //--------------------------------

    };
    Bolsig bol;

    //particle model parameter
    //--------------------------------
    struct Particle{ 
        
        //read
        //--------------------------------
        struct flag{ 
            int iPIC = 0; //0:fluid, 1:PIC
            int ePIC = 0; //0:fluid, 1:PIC
            int wgtFunc = 0; //weight function (0:NGP, 1:CIC)
        };
        flag flag;
        int nPclPerCell_ini = 512; //numper of partilces per cell
        //double ratioPclMax_ini = 3.0; //初期粒子数に対して何倍の粒子数分メモリ確保するか
        //--------------------------------

        //do not read
        //--------------------------------
        //int nPclMax = 0;
        //int nPcl_ini = 0;
        //double wgti = 0.0;
        //double wgte = 0.0;

        double sigmavMax_i = 0.0;
        double sigmavMax_e = 0.0;
        //--------------------------------
    };
    Particle pcl;
    //--------------------------------

    //make distribution function
    //--------------------------------
    struct DistributionFunction{ 
        
        //read
        //--------------------------------
        struct flag{ 
            int calcEEDF = 0;
        };
        flag flag;
        int nBinEngy_EEDF = 0; //numper of energy bin size for EDF
        double engyMax_EEDF = 0.0; //maximum energy for EDF
        //--------------------------------

    };
    DistributionFunction df;
    //--------------------------------

};
