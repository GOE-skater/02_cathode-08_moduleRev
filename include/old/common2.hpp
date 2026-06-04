#pragma once
#include <vector>
#include <complex>

//============================================================
//     Common file for MHD_1D_HLLD
//          coded by Ryo, Shirakawa.
//=============================================================
//******************* グローバル定数 **************************
//物理定数
extern const double Boltz; //ボルツマン定数
extern const double e0; //電気素量
extern double eps0; //真空の誘電率
extern const double mu0; //真空の透磁率
extern const double masse; //電子質量

//******************* 読み込みあり **************************
//物理定数
extern double massi; //イオン質量
extern double ri ; //イオン半径
extern double DmN; //準安定種の拡散係数 (Elis 1969)
extern double epsr_diele; //誘電体の比電率
extern double V_bias; //バイアス電圧
//extern double tanD_diele; //誘電体の誘電正接

//計算領域
extern int ni; //x方向のセルの数
extern int nj; //x方向のセルの数
extern int n_bl; //計算ブロック数
extern double xL; //計算領域x左端
extern double xR; //計算領域x右端
extern double rmin; //計算領域r左端
extern double rmax; //計算領域r右端

//パラメータ設定
extern std::complex<double> J1r_exc; //励起電流
extern double Pmw; //マイクロ波電力 (W)
extern int icon_mwRef; //反射を解くか
extern int icon_impTest;      // (int) conduct impedance test
extern double omegam ; //マイクロ波周波数
extern double nu_eff; //有効衝突周波数 = sqrt(vth*omegam/deltaB)
extern double deltaECR; //ECRの領域の幅
extern double S11_mag;     // (double) Magnitude of S11
extern double S11_arg_deg; // (double) Argument of S11 (deg)
extern double S21_mag;     // (double) Magnitude of S21
extern double S21_arg_deg; // (double) Argument of S21 (deg)
extern double S12_mag;     // (double) Magnitude of S12
extern double S12_arg_deg; // (double) Argument of S12 (deg)
extern double S22_mag;     // (double) Magnitude of S22
extern double S22_arg_deg; // (double) Argument of S22 (deg)
extern double Z0_base;     // (double) Base impedance (Ohm)

extern double Ti; //イオン温度
extern double Tn; //中性粒子温度
extern double rhon_ini; //中性粒子密度 (一定値)
extern double Q_neutIn_mgs; //中性粒子の流量 mg/s
extern double width_neutIn; //中性粒子のインレットの幅 m
extern double Te_rep_eV; //代表電子温度 (イオンスキーム用)
extern int icon_Bohm; //Bohm拡散を考慮するか
extern int icon_Sagdeev; //Sagdeevモデルを考慮するか
extern double alpha_Bohm; //ボーム拡散係数
extern double scale_inertia; //scaling factor of inertia term

//SEE
//extern double coefEISEE_eb;  //EI elastic-backscattering
//extern double coefEISEE_rd;  //EI rediffusion
//extern double coefEISEE_ts;  //EI true-secondary
extern double coefIISEE_ts; //II true-secondary
extern double coefMISEE_ts; //MI true-secondary (Auger)
extern double Te_emitSEE_eV; //electron emission temperature (eV)
extern double ratioEngy_EISEE_rd; //energy ratio of EI rediffusion 

//流体スキーム設定
extern int icon_PC; //PC methodを使うか (0:直接離散化，1:PC method)
extern int icon_inertia; //慣性を考慮するか
extern int icon_adp_dt; //アダプティブなdtを用いるか 0:固定dt, 1:可変dt
extern double dt_ini; //初期dt
extern int ndt_i; //イオンの電子に対するタイムステップ倍率
extern int ndt_m; //中性 (準安定) の電子に対するタイムステップ倍率
extern int ndt_n; //中性 (準安定) の電子に対するタイムステップ倍率
extern double CFL; //CFL数

extern double error_cnv_SOR_Ui;
extern int maxITR_SOR_Ui;
extern int icon_iter_Ui; //反復法の選択 0:SOR 1:SMG
extern double error_cnv_SOR_rhoi; //rhoiの収束基準
extern int maxITR_SOR_rhoi; //rhoiの最大反復数
extern int icon_iter_rhoi; //rhoiの反復法の選択 0:SOR 1:SMG
extern double error_cnv_SOR_phi;  //phiの収束基準
extern int maxITR_SOR_phi;        //phiの最大反復数
extern int icon_iter_phi;            //phiの反復法の選択 0:SOR 1:SMG
extern double error_cnv_HES_phi;
extern int maxITR_HES_phi;
extern double error_cnv_SOR_rhoe;  //rhoeの収束基準
extern int maxITR_SOR_rhoe;        //rhoeの最大反復数
extern int icon_iter_rhoe;            //rhoeの反復法の選択 0:SOR 1:SMG
extern double error_cnv_HES_rhoe;
extern int maxITR_HES_rhoe;
extern double error_cnv_SOR_rhoeps;  //rhoepsの収束基準
extern int maxITR_SOR_rhoeps;        //rhoepsの最大反復数
extern int icon_iter_rhoeps;            //rhoepsの反復法の選択 0:SOR 1:SMG
extern double error_cnv_HES_rhoeps;
extern int maxITR_HES_rhoeps;
extern int icon_error; //errorの種類 0:max 1:rms
extern int ndiv_MW; ///マイクロ波の計算間隔;

extern double beta_rhoe; //rhoeの不足緩和係数 (慣性考慮時)
extern double beta_rhoUe; //rhoUeの不足緩和係数 (慣性考慮時)

//全体の計算設定
extern int ntime;
extern int icon_autoFinish;

//出力設定
extern int icon_chk; //チェックファイルを出力するか
extern int icon_gnuRes;
extern int ndiv_fout;
//******************* グローバル変数 **************************
//変数
extern double Lx; //計算領域のx方向の大きさ
extern double Lr; //計算領域のx方向の大きさ

extern double gtime; //時間
extern int itime; //時間(整数)

extern double dx; //x方向の刻み
extern double dr; //r方向の刻み
extern double dt; //時間刻み

extern double I_Anode; //アノード電流
extern double Ii_Anode; //アノード電流
extern double Ie_Anode; //アノード電流

extern double I_Nozzle; //ノズル電流
extern double Ii_Nozzle; //ノズル電流
extern double Ie_Nozzle; //ノズル電流

extern int nOut; //ファイルアウトプット回数

extern int icon_restart;

extern double error_rhoi;
extern double error_rhoe;
extern double error_rhoUex;
extern double error_rhoUer;
extern double error_rhoeps;
extern double error_Gx;
extern double error_Gr;

extern double error_Uix;
extern double error_Uir;
extern double error_Uip;
extern double error_phi;
extern double error_rhom;
extern double error_rhon;

//反応レート係数データ
extern double dTe; //電子温度の刻み幅
extern double dEi; //イオンドリフトエネルギーの刻み幅
extern double dTe_SEE; //SEE用の電子温度の刻み幅

//BiCGSTAB用変数
extern int nk; //kの数 (ExとErの合計)
extern int nkx; //kの数 (Ex)
extern int nkr; //kの数 (Er)
extern int nkp; //kの数 (Er)
extern int nkfc; //kの数 (fluid)
extern int nkfc1; //kの数 (fluid)
extern int nkfx; //kの数 (fluid)
extern int nkfr; //kの数 (fluid)

//流入条件
extern double fn_In; //中性粒子の流入フラックス

//ジオメトリ
extern double x1;
extern double x2;
extern double x3;
extern double x4;
extern double x5;
extern double x6;

extern double r1;
extern double r2;
extern double r3;
extern double r4;
extern double r5;
extern double r6;

//計算領域分割変数
//extern std::vector<int> i_flc_bl1;
//extern std::vector<int> j_flc_bl1;
//extern std::vector<int> i_flc_bl2;
//extern std::vector<int> j_flc_bl2;
//extern std::vector<std::vector<int> > i_flc_bl;
//extern std::vector<std::vector<int> > j_flc_bl;
//extern std::vector<std::vector<int> > i_flx_bl;
//extern std::vector<std::vector<int> > j_flx_bl;
//extern std::vector<std::vector<int> > i_flr_bl;
//extern std::vector<std::vector<int> > j_flr_bl;

//extern std::vector<std::vector<int> > kx; //(i,j) → kへの変換 (Ex)
//extern std::vector<std::vector<int> > kr; //(i,j) → kへの変換 (Er)
//extern std::vector<std::vector<int> > kp; //(i,j) → kへの変換 (Er)
//extern std::vector<std::vector<int> > kfc; //(i,j) → kへの変換 (fluid)
//extern std::vector<std::vector<int> > kfc1; //(i,j) → kへの変換 (fluid)
//extern std::vector<std::vector<int> > kfx; //(i,j) → kへの変換 (fluid)
//extern std::vector<std::vector<int> > kfr; //(i,j) → kへの変換 (fluid)
//extern std::vector<int> ikx; //k→ iへの変換 (Ex)
//extern std::vector<int> ikr; //k→ iへの変換 (Er)
//extern std::vector<int> ikp; //k→ iへの変換 (Er)
extern std::vector<int> ikfc; //k→ iへの変換 (fluid)
extern std::vector<int> ikfc1; //k→ iへの変換 (fluid)
extern std::vector<int> ikfx; //k→ iへの変換 (fluid)
extern std::vector<int> ikfr; //k→ iへの変換 (fluid)
//extern std::vector<int> jkx; //k→ jへの変換 (Ex)
//extern std::vector<int> jkr; //k→ jへの変換 (Er)
//extern std::vector<int> jkp; //k→ jへの変換 (Er)
extern std::vector<int> jkfc; //k→ jへの変換 (fluid)
extern std::vector<int> jkfc1; //k→ jへの変換 (fluid)
extern std::vector<int> jkfx; //k→ jへの変換 (fluid)
extern std::vector<int> jkfr; //k→ jへの変換 (fluid)

//境界条件保持用変数
//extern std::vector<int> iBndWx;
//extern std::vector<int> jBndWx;
//extern std::vector<int> sBndWx;
//
//extern std::vector<int> iBndEx;
//extern std::vector<int> jBndEx;
//extern std::vector<int> sBndEx;
//
//extern std::vector<int> iBndSx;
//extern std::vector<int> jBndSx;
//extern std::vector<int> sBndSx;
//
//extern std::vector<int> iBndNx;
//extern std::vector<int> jBndNx;
//extern std::vector<int> sBndNx;
//
//extern std::vector<int> iBndWr;
//extern std::vector<int> jBndWr;
//extern std::vector<int> sBndWr;
//
//extern std::vector<int> iBndEr;
//extern std::vector<int> jBndEr;
//extern std::vector<int> sBndEr;
//
//extern std::vector<int> iBndSr;
//extern std::vector<int> jBndSr;
//extern std::vector<int> sBndSr;
//
//extern std::vector<int> iBndNr;
//extern std::vector<int> jBndNr;
//extern std::vector<int> sBndNr;
//
//extern std::vector<int> iBndWp;
//extern std::vector<int> jBndWp;
//extern std::vector<int> sBndWp;
//
//extern std::vector<int> iBndEp;
//extern std::vector<int> jBndEp;
//extern std::vector<int> sBndEp;
//
//extern std::vector<int> iBndSp;
//extern std::vector<int> jBndSp;
//extern std::vector<int> sBndSp;
//
//extern std::vector<int> iBndNp;
//extern std::vector<int> jBndNp;
//extern std::vector<int> sBndNp;

//配列 フィールドデータ
//extern std::vector<double> x; //x座標
//extern std::vector<double> r; //r座標

//境界判別 (解くべきところが1，ゴーストが2，解かないのが0)
//extern std::vector<std::vector<int> > jdgBnd_Ex;
//extern std::vector<std::vector<int> > jdgBnd_Er;
//extern std::vector<std::vector<int> > jdgBnd_Ep;
//extern std::vector<std::vector<int> > jdgBnd_flc;
//extern std::vector<std::vector<int> > jdgBnd_flx;
//extern std::vector<std::vector<int> > jdgBnd_flr;

//マイクロ波電界
//extern std::vector<std::vector<std::complex<double> > > E1r;
//extern std::vector<std::vector<std::complex<double> > > E1x;
//extern std::vector<std::vector<std::complex<double> > > E1p;

//励振電流密度
//extern std::vector<std::vector<std::complex<double> > > J1r;
//extern std::vector<std::vector<std::complex<double> > > J1x;
//extern std::vector<std::vector<std::complex<double> > > J1p;

//比誘電率
//extern std::vector<std::vector<double> > epsr;

//プラズマ変数
//extern std::vector<std::vector<double> > rhoi; //プラズマ密度
//extern std::vector<std::vector<double> > rhoi_old; //密度
//extern std::vector<std::vector<double> > rhoi_Lx; //密度
//extern std::vector<std::vector<double> > rhoi_Rx; //密度
//extern std::vector<std::vector<double> > rhoi_Lr; //密度
//extern std::vector<std::vector<double> > rhoi_Rr; //密度


//extern std::vector<std::vector<double> > rhon; //基底密度
//extern std::vector<std::vector<double> > rhon_old; //基底密度
extern std::vector<std::vector<double> > rhon_Lx; //基底密度
extern std::vector<std::vector<double> > rhon_Rx; //基底密度
extern std::vector<std::vector<double> > rhon_Lr; //基底密度
extern std::vector<std::vector<double> > rhon_Rr; //基底密度

//extern std::vector<std::vector<double> > rhoUnx; //フラックス
//extern std::vector<std::vector<double> > rhoUnx_old; //フラックス
extern std::vector<std::vector<double> > rhoUnx_Lx; //フラックス
extern std::vector<std::vector<double> > rhoUnx_Rx; //フラックス
extern std::vector<std::vector<double> > rhoUnx_Lr; //フラックス
extern std::vector<std::vector<double> > rhoUnx_Rr; //フラックス

//extern std::vector<std::vector<double> > rhoUnr; //フラックス
//extern std::vector<std::vector<double> > rhoUnr_old; //フラックス
extern std::vector<std::vector<double> > rhoUnr_Lx; //フラックス
extern std::vector<std::vector<double> > rhoUnr_Rx; //フラックス
extern std::vector<std::vector<double> > rhoUnr_Lr; //フラックス
extern std::vector<std::vector<double> > rhoUnr_Rr; //フラックス

//extern std::vector<std::vector<double> > rhom; //準安定密度
//extern std::vector<std::vector<double> > rhom_old; //準安定密度

//extern std::vector<std::vector<double> > rhoeps; //エネルギー密度
//extern std::vector<std::vector<double> > rhoeps_old; //エネルギー密度

extern std::vector<std::vector<double> > scx; //表面電荷
extern std::vector<std::vector<double> > scx_old; //表面電荷
extern std::vector<std::vector<double> > scr; //表面電荷
extern std::vector<std::vector<double> > scr_old; //表面電荷

//extern std::vector<std::vector<double> > rhoUix_wall; //イオンフラックス
//extern std::vector<std::vector<double> > rhoUir_wall; //イオンフラックス
//extern std::vector<std::vector<double> > rhoUex_wall; //電子フラックス
//extern std::vector<std::vector<double> > rhoUer_wall; //電子フラックス
//extern std::vector<std::vector<double> > rhoUmx_wall; //準安定フラック
//extern std::vector<std::vector<double> > rhoUmr_wall; //準安定フラック
//extern std::vector<std::vector<double> > rhoUnx_wall; //基底フラックス
//extern std::vector<std::vector<double> > rhoUnr_wall; //基底フラックス
//extern std::vector<std::vector<double> > Gx_wall; //エネルギーフラックス
//extern std::vector<std::vector<double> > Gr_wall; //エネルギーフラックス

//extern std::vector<std::vector<double> > Uix; //速度
//extern std::vector<std::vector<double> > Uix_old; //速度
//extern std::vector<std::vector<double> > rhoUix; //速度
//extern std::vector<std::vector<double> > Uir; //速度
//extern std::vector<std::vector<double> > Uir_old; //速度
//extern std::vector<std::vector<double> > rhoUir; //速度

//extern std::vector<std::vector<double> > Uip; //速度
//extern std::vector<std::vector<double> > Uip_old; //速度

extern std::vector<std::vector<double> > limit; //リミッター

//extern std::vector<std::vector<double> > fix1; //イオンフラックスベクトル
//extern std::vector<std::vector<double> > fix2; //イオンフラックスベクトル
//extern std::vector<std::vector<double> > fix3; //イオンフラックスベクトル
//extern std::vector<std::vector<double> > fix4; //イオンフラックスベクトル
//extern std::vector<std::vector<double> > fir1; //イオンフラックスベクトル
//extern std::vector<std::vector<double> > fir2; //イオンフラックスベクトル
//extern std::vector<std::vector<double> > fir3; //イオンフラックスベクトル
//extern std::vector<std::vector<double> > fir4; //イオンフラックスベクトル

extern std::vector<std::vector<double> > fnx1; //基底フラックスベクトル
extern std::vector<std::vector<double> > fnx2; //基底フラックスベクトル
extern std::vector<std::vector<double> > fnx3; //基底フラックスベクトル
extern std::vector<std::vector<double> > fnr1; //基底フラックスベクトル
extern std::vector<std::vector<double> > fnr2; //基底フラックスベクトル
extern std::vector<std::vector<double> > fnr3; //基底フラックスベクトル

//extern std::vector<std::vector<double> > Te; //電子温度
//extern std::vector<std::vector<double> > Te_old; //電子温度
extern std::vector<std::vector<double> > Te_Lx; //電子温度
extern std::vector<std::vector<double> > Te_Rx; //電子温度
extern std::vector<std::vector<double> > Te_Lr; //電子温度
extern std::vector<std::vector<double> > Te_Rr; //電子温度
//extern std::vector<std::vector<double> > Pabs; //マイクロ波電力吸収

extern std::vector<std::vector<double> > qex;    //熱流束
extern std::vector<std::vector<double> > qer;    //熱流束

//extern std::vector<std::vector<double> > Gx;  //エネルギーフラックス
//extern std::vector<std::vector<double> > Gr;  //エネルギーフラックス
//extern std::vector<std::vector<double> > Gx_old;  //エネルギーフラックス
//extern std::vector<std::vector<double> > Gr_old;  //エネルギーフラックス

//extern std::vector<std::vector<double> > phi;//両極性電位
//extern std::vector<std::vector<double> > phi_old;//両極性電位
//extern std::vector<std::vector<double> > Ex;//両極性電場
//extern std::vector<std::vector<double> > Er;//両極性電場

//extern std::vector<std::vector<double> > rhoUex; //電子フラックス
//extern std::vector<std::vector<double> > rhoUer; //電子フラックス
//extern std::vector<std::vector<double> > rhoUex_old; //電子フラックス
//extern std::vector<std::vector<double> > rhoUer_old; //電子フラックス

//extern std::vector<std::vector<double> > nUex; //電子フラックス (サブ)
//extern std::vector<std::vector<double> > nUer; //電子フラックス (サブ)
//extern std::vector<std::vector<double> > nUex_old; //電子フラックス (サブ)
//extern std::vector<std::vector<double> > nUer_old; //電子フラックス (サブ)

//extern std::vector<std::vector<double> > rhoUmx; //準安定フラックス
//extern std::vector<std::vector<double> > rhoUmr; //準安定フラックス
//extern std::vector<std::vector<double> > rhoUmx_old; //準安定フラックス
//extern std::vector<std::vector<double> > rhoUmr_old; //準安定フラックス

//extern std::vector<std::vector<double> > rate_ionize; //イオン化レート
//extern std::vector<std::vector<double> > Mache; //電子マッハ数

//衝突周波数
//extern std::vector<std::vector<double> > nu_m; //電子運動量移動周波数
//extern std::vector<std::vector<double> > nu_m1; //電子運動量移動周波数 (有効衝突追加)
//extern std::vector<std::vector<double> > nui_m; //イオン衝突周波数
//extern std::vector<std::vector<double> > nu_elas; //弾性周波数
//extern std::vector<std::vector<double> > nu_super; //超弾性周波数
//extern std::vector<std::vector<double> > nu_ionz; //電離衝突周波数
//extern std::vector<std::vector<double> > nu_ionzStep; //段階電離衝突周波数
//extern std::vector<std::vector<double> > nu_exc; //励起周波数
//extern std::vector<std::vector<double> > nu_excReso; //励起(共鳴)周波数
//extern std::vector<std::vector<double> > nu_excMeta; //励起(準安定)周波数
//extern std::vector<std::vector<double> > nu_excStep; //段階励起周波数
//extern std::vector<std::vector<double> > nu_ano; //異常衝突周波数
//extern std::vector<std::vector<double> > nu_ano_IAT; //異常衝突周波数

//エネルギー交換
//extern std::vector<std::vector<double> > rate_eloss; //電子エネルギー獲得量
//extern std::vector<std::vector<double> > delossdeps; //電子エネルギーロスの微係数

//移動度
//extern std::vector<std::vector<double> > mu_para;
//extern std::vector<std::vector<double> > mu_perp;
//extern std::vector<std::vector<double> > mu_min;
//extern std::vector<std::vector<double> > mu_xx;
//extern std::vector<std::vector<double> > mu_rr;
//extern std::vector<std::vector<double> > mu_xr;
//extern std::vector<std::vector<double> > mu_inv_xx;
//extern std::vector<std::vector<double> > mu_inv_rr;
//extern std::vector<std::vector<double> > mu_inv_xr;
//extern std::vector<std::vector<double> > kappa_para;
//extern std::vector<std::vector<double> > kappa_perp;
//extern std::vector<std::vector<double> > kappa_x;
//extern std::vector<std::vector<double> > kappa_r;
//extern std::vector<std::vector<double> > kappa_c;
//extern std::vector<std::vector<double> > kappa_x_inv;
//extern std::vector<std::vector<double> > kappa_r_inv;
//extern std::vector<std::vector<double> > kappa_c_inv;
//extern std::vector<std::vector<double> > cx; //係数保持
//extern std::vector<std::vector<double> > cr; //係数保持
//extern std::vector<std::vector<double> > sx; //係数保持
//extern std::vector<std::vector<double> > sr; //係数保持

//固定磁場
//extern std::vector<std::vector<double> > Bx;
//extern std::vector<std::vector<double> > Br;
//extern std::vector<std::vector<double> > Ap; //ベクトルポテンシャル
extern std::vector<std::vector<double> > distECR; //ECR領域からの距離

//テスト配列
extern std::vector<std::vector<double> > test1;
extern std::vector<std::vector<double> > test2;
extern std::vector<std::vector<double> > test3;

//配列 BOLSIG+
//extern std::vector<double> nu_elas_N; //弾性周波数/nn <σelas*vth>
//extern std::vector<double> nu_super_N; //超弾性周波数/nn <σsuper*vth>
//extern std::vector<double> nu_ionz_N; //イオン化周波数/nn <σionz*vth>
//extern std::vector<double> nu_ionzStep_N; //段階イオン化周波数/nn <σionz*vth>
//extern std::vector<double> nu_exc_N; //励起周波数/nn <σexc*vth>
//extern std::vector<double> nu_excReso_N; //共鳴励起周波数/nn <σexc*vth>
//extern std::vector<double> nu_excMeta_N; //準安定励起周波数/nn <σexc*vth>
//extern std::vector<double> nu_excStep_N; //段階励起周波数/nn <σexc*vth>
//extern std::vector<double> rate_eloss_n_N; //エネルギーロスレート/nn
//extern std::vector<double> rate_eloss_m_N; //エネルギーロスレート/nn
//extern std::vector<double> nui_m_N; //イオン運動量移動周波数/nn <σionz*vth>

//配列 EI-SEE table
//extern std::vector<double> coefEISEE_eb_table; //EI elastic-backscattering
//extern std::vector<double> coefEISEE_rd_table; //EI rediffusion
//extern std::vector<double> coefEISEE_ts_table; //EI true-secondary
