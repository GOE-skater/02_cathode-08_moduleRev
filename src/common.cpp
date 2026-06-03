//============================================================
//     Common file for MHD_1D_HLLD
//          coded by Ryo, Shirakawa.
//============================================================
#include <iostream>
#include <string>
#include <math.h>
#include "common.hpp"
#include <complex>

//******************* グローバル定数 **************************
//物理定数
const double Boltz=1.380649e-23; //ボルツマン定数
const double e0=1.60217662e-19; //電気素量
double eps0 = 8.85418782e-12*1e0; //真空の誘電率
const double mu0 = 1.2566370614e-6; //真空の透磁率
const double masse=9.10938356e-31; //電子質量

//******************* 読み込みあり **************************
//物理定数
double massi = 131.293e-3/6.02214076e23; //イオン質量 Xe 131.293
double ri = 216e-12; //イオン半径 Xe = 216pm
double DmN = 1.66e20; //準安定種の拡散係数 (Elis 1969) 1/ms
double epsr_diele = 1.0; //誘電体の比電率
double V_bias = 0.0; //バイアス電圧
//double tanD_diele = 0.0002; //誘電体の誘電正接

//計算領域
int ni = 194; //x方向のセルの数 (実際はゴーストセルが+2) //155 + 39
int nj = 70; //x方向のセルの数 (実際はゴーストセルが+2) //100
int n_bl = 6; //計算ブロック数
double xL= 0.00; //計算領域x左端
double xR= 0.0388; //計算領域x右端 0.031 + 0.0078
double rmin= 0.0; //計算領域r左端
double rmax= 0.014; //計算領域r右端 2.0e-3 0.014

//パラメータ設定
std::complex<double> J1r_exc(0.212687,0.0); //励起電流
double Pmw = 5.0; //マイクロ波電力 (W)
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
double Ti = 300; //イオン温度
double Tn = 300; //中性粒子温度
double rhon_ini = 1.0/(Boltz*Tn); //中性粒子密度 (一定値) m^-3 @1Pa
double Q_neutIn_mgs = 0.0; //中性粒子の流量 mg/s
double width_neutIn = 0.0; //中性粒子のインレットの幅 m
double Te_rep_eV = 3.0; //代表電子温度 (イオンスキーム用)
int icon_Bohm = 0; //Bohm拡散を考慮するか
int icon_Sagdeev = 0; //Sagdeevモデルを考慮するか
double alpha_Bohm = 16.0; //ボーム拡散係数
double scale_inertia = 1.0; //scaling factor of inertia term

//SEE
//double coefEISEE_eb = 0.40;  //EI elastic-backscattering
//double coefEISEE_rd = 0.40;  //EI rediffusion
//double coefEISEE_ts = 0.20;  //EI true-secondary
double coefIISEE_ts = 0.08; //II true-secondary
double coefMISEE_ts = 0.08; //MI true-secondary (Auger)
double Te_emitSEE_eV = 2.5; //electron emission temperature (eV)
double ratioEngy_EISEE_rd = 0.5; //energy ratio of EI rediffusion 

//流体スキーム設定
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
//******************* グローバル変数 **************************
//変数
double Lx=0.0; //計算領域のx方向の大きさ
double Lr=0.0; //計算領域のx方向の大きさ

double gtime=0.0; //時間
int itime=0; //時間(整数)

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

int icon_restart = 0;

double error_rhoi = 0.0;
double error_rhoe = 0.0;

double error_rhoUex = 0.0;
double error_rhoUer = 0.0;
double error_rhoeps = 0.0;
double error_Gx = 0.0;
double error_Gr = 0.0;


double error_Uix = 0.0;
double error_Uir = 0.0;
double error_Uip = 0.0;
double error_phi = 0.0;
double error_rhom = 0.0;
double error_rhon = 0.0;

//反応レート係数データ
double dTe = 0.0; //電子温度の刻み幅
double dEi = 0.0; //イオンドリフトエネルギーの刻み幅
double dTe_SEE = 0.0; //SEE用の電子温度の刻み幅

//BiCGSTAB用変数
int nk = 0; //kの数 (ExとErの合計)
int nkx = 0; //kの数 (Ex)
int nkr = 0; //kの数 (Er)
int nkp = 0; //kの数 (Ep)
int nkfc = 0; //kの数 (fluid)
int nkfc1 = 0; //kの数 (fluid)
int nkfx = 0; //kの数 (fluid)
int nkfr = 0; //kの数 (fluid)

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


//計算領域分割変数
//std::vector<int> i_flc_bl1(2,0);
//std::vector<int> j_flc_bl1(2,0);
//std::vector<int> i_flc_bl2(2,0);
//std::vector<int> j_flc_bl2(2,0);
std::vector<std::vector<int> > i_flc_bl(n_bl,std::vector<int>(2,0));
std::vector<std::vector<int> > j_flc_bl(n_bl,std::vector<int>(2,0));
std::vector<std::vector<int> > i_flx_bl(n_bl,std::vector<int>(2,0));
std::vector<std::vector<int> > j_flx_bl(n_bl,std::vector<int>(2,0));
std::vector<std::vector<int> > i_flr_bl(n_bl,std::vector<int>(2,0));
std::vector<std::vector<int> > j_flr_bl(n_bl,std::vector<int>(2,0));

std::vector<std::vector<int> > kx(ni+2,std::vector<int>(nj+2,-1)); //(i,j) → kへの変換 (Ex)
std::vector<std::vector<int> > kr(ni+2,std::vector<int>(nj+2,-1)); //(i,j) → kへの変換 (Er)
std::vector<std::vector<int> > kp(ni+2,std::vector<int>(nj+2,-1)); //(i,j) → kへの変換 (Er)
std::vector<std::vector<int> > kfc(ni+2,std::vector<int>(nj+2,-1)); //(i,j) → kへの変換 (fluid)
std::vector<std::vector<int> > kfc1(ni+2,std::vector<int>(nj+2,-1)); //(i,j) → kへの変換 (fluid)
std::vector<std::vector<int> > kfx(ni+2,std::vector<int>(nj+2,-1)); //(i,j) → kへの変換 (fluid)
std::vector<std::vector<int> > kfr(ni+2,std::vector<int>(nj+2,-1)); //(i,j) → kへの変換 (fluid)
std::vector<int> ikx; //k→ iへの変換 (Ex)
std::vector<int> ikr; //k→ iへの変換 (Er)
std::vector<int> ikp; //k→ iへの変換 (Ep)
std::vector<int> ikfc; //k→ iへの変換 (fluid)
std::vector<int> ikfc1; //k→ iへの変換 (fluid)
std::vector<int> ikfx; //k→ iへの変換 (fluid)
std::vector<int> ikfr; //k→ iへの変換 (fluid)
std::vector<int> jkx; //k→ jへの変換 (Ex)
std::vector<int> jkr; //k→ jへの変換 (Er)
std::vector<int> jkp; //k→ jへの変換 (Ep)
std::vector<int> jkfc; //k→ jへの変換 (fluid)
std::vector<int> jkfc1; //k→ jへの変換 (fluid)
std::vector<int> jkfx; //k→ jへの変換 (fluid)
std::vector<int> jkfr; //k→ jへの変換 (fluid)

//境界条件保持用変数
std::vector<int> iBndWx;
std::vector<int> jBndWx;
std::vector<int> sBndWx;

std::vector<int> iBndEx;
std::vector<int> jBndEx;
std::vector<int> sBndEx;

std::vector<int> iBndSx;
std::vector<int> jBndSx;
std::vector<int> sBndSx;

std::vector<int> iBndNx;
std::vector<int> jBndNx;
std::vector<int> sBndNx;

std::vector<int> iBndWr;
std::vector<int> jBndWr;
std::vector<int> sBndWr;

std::vector<int> iBndEr;
std::vector<int> jBndEr;
std::vector<int> sBndEr;

std::vector<int> iBndSr;
std::vector<int> jBndSr;
std::vector<int> sBndSr;

std::vector<int> iBndNr;
std::vector<int> jBndNr;
std::vector<int> sBndNr;

std::vector<int> iBndWp;
std::vector<int> jBndWp;
std::vector<int> sBndWp;

std::vector<int> iBndEp;
std::vector<int> jBndEp;
std::vector<int> sBndEp;

std::vector<int> iBndSp;
std::vector<int> jBndSp;
std::vector<int> sBndSp;

std::vector<int> iBndNp;
std::vector<int> jBndNp;
std::vector<int> sBndNp;

//  ディリクレ上側境界 (Ex)
std::vector<int> iBndDiNx;
std::vector<int> jBndDiNx;
//  ディリクレ下側境界 (Ex)
std::vector<int> iBndDiSx;
std::vector<int> jBndDiSx;
//  ディリクレ左側境界 (Ex)
std::vector<int> iBndDiWx;
std::vector<int> jBndDiWx;

//配列 フィールドデータ
std::vector<double> x(ni+2,0.0); //x座標
std::vector<double> r(nj+2,0.0); //r座標

//境界判別 (解くべきところが1，ゴーストが2，解かないのが0)
std::vector<std::vector<int> > jdgBnd_Ex(ni+2,std::vector<int>(nj+2,0));
std::vector<std::vector<int> > jdgBnd_Er(ni+2,std::vector<int>(nj+2,0));
std::vector<std::vector<int> > jdgBnd_Ep(ni+2,std::vector<int>(nj+2,0));
std::vector<std::vector<int> > jdgBnd_flc(ni+2,std::vector<int>(nj+2,0));
std::vector<std::vector<int> > jdgBnd_flx(ni+2,std::vector<int>(nj+2,0));
std::vector<std::vector<int> > jdgBnd_flr(ni+2,std::vector<int>(nj+2,0));

//マイクロ波電界
std::vector<std::vector<std::complex<double> > > E1r(ni+2,std::vector<std::complex<double> >(nj+2,0.0));
std::vector<std::vector<std::complex<double> > > E1x(ni+2,std::vector<std::complex<double> >(nj+2,0.0));
std::vector<std::vector<std::complex<double> > > E1p(ni+2,std::vector<std::complex<double> >(nj+2,0.0));

//励振電流密度
std::vector<std::vector<std::complex<double> > > J1r(ni+2,std::vector<std::complex<double> >(nj+2,0.0));
std::vector<std::vector<std::complex<double> > > J1x(ni+2,std::vector<std::complex<double> >(nj+2,0.0));
std::vector<std::vector<std::complex<double> > > J1p(ni+2,std::vector<std::complex<double> >(nj+2,0.0));

//比誘電率
std::vector<std::vector<double> > epsr(ni+2,std::vector<double>(nj+2,1.0));

//プラズマ変数
std::vector<std::vector<double> > rhoi(ni+2,std::vector<double>(nj+2,0.0)); //プラズマ密度
std::vector<std::vector<double> > rhoi_old(ni+2,std::vector<double>(nj+2,0.0)); //密度
//std::vector<std::vector<double> > rhoi_Lx(ni+2,std::vector<double>(nj+2,0.0)); //密度
//std::vector<std::vector<double> > rhoi_Rx(ni+2,std::vector<double>(nj+2,0.0)); //密度
//std::vector<std::vector<double> > rhoi_Lr(ni+2,std::vector<double>(nj+2,0.0)); //密度
//std::vector<std::vector<double> > rhoi_Rr(ni+2,std::vector<double>(nj+2,0.0)); //密度

std::vector<std::vector<double> > rhoe(ni+2,std::vector<double>(nj+2,0.0)); //電子密度
std::vector<std::vector<double> > rhoe_old(ni+2,std::vector<double>(nj+2,0.0)); //電子密度

std::vector<std::vector<double> > rhon(ni+2,std::vector<double>(nj+2,0.0)); //基底密度
std::vector<std::vector<double> > rhon_old(ni+2,std::vector<double>(nj+2,0.0)); //基底密度
std::vector<std::vector<double> > rhon_Lx(ni+2,std::vector<double>(nj+2,0.0)); //密度
std::vector<std::vector<double> > rhon_Rx(ni+2,std::vector<double>(nj+2,0.0)); //密度
std::vector<std::vector<double> > rhon_Lr(ni+2,std::vector<double>(nj+2,0.0)); //密度
std::vector<std::vector<double> > rhon_Rr(ni+2,std::vector<double>(nj+2,0.0)); //密度

std::vector<std::vector<double> > rhoUnx(ni+2,std::vector<double>(nj+2,0.0)); //フラックス
std::vector<std::vector<double> > rhoUnx_old(ni+2,std::vector<double>(nj+2,0.0)); //フラックス
std::vector<std::vector<double> > rhoUnx_Lx(ni+2,std::vector<double>(nj+2,0.0)); //フラックス
std::vector<std::vector<double> > rhoUnx_Rx(ni+2,std::vector<double>(nj+2,0.0)); //フラックス
std::vector<std::vector<double> > rhoUnx_Lr(ni+2,std::vector<double>(nj+2,0.0)); //フラックス
std::vector<std::vector<double> > rhoUnx_Rr(ni+2,std::vector<double>(nj+2,0.0)); //フラックス

std::vector<std::vector<double> > rhoUnr(ni+2,std::vector<double>(nj+2,0.0)); //フラックス
std::vector<std::vector<double> > rhoUnr_old(ni+2,std::vector<double>(nj+2,0.0)); //フラックス
std::vector<std::vector<double> > rhoUnr_Lx(ni+2,std::vector<double>(nj+2,0.0)); //フラックス
std::vector<std::vector<double> > rhoUnr_Rx(ni+2,std::vector<double>(nj+2,0.0)); //フラックス
std::vector<std::vector<double> > rhoUnr_Lr(ni+2,std::vector<double>(nj+2,0.0)); //フラックス
std::vector<std::vector<double> > rhoUnr_Rr(ni+2,std::vector<double>(nj+2,0.0)); //フラックス

std::vector<std::vector<double> > rhom(ni+2,std::vector<double>(nj+2,0.0)); //準安定密度
std::vector<std::vector<double> > rhom_old(ni+2,std::vector<double>(nj+2,0.0)); //準安定密度

std::vector<std::vector<double> > rhoeps(ni+2,std::vector<double>(nj+2,0.0)); //エネルギー密度
std::vector<std::vector<double> > rhoeps_old(ni+2,std::vector<double>(nj+2,0.0)); //エネルギー密度

std::vector<std::vector<double> > scx(ni+2,std::vector<double>(nj+2,0.0)); //表面電荷
std::vector<std::vector<double> > scx_old(ni+2,std::vector<double>(nj+2,0.0)); //表面電荷
std::vector<std::vector<double> > scr(ni+2,std::vector<double>(nj+2,0.0)); //表面電荷
std::vector<std::vector<double> > scr_old(ni+2,std::vector<double>(nj+2,0.0)); //表面電荷

std::vector<std::vector<double> > rhoUix_wall(ni+2,std::vector<double>(nj+2,0.0)); //イオンフラックス
std::vector<std::vector<double> > rhoUir_wall(ni+2,std::vector<double>(nj+2,0.0)); //イオンフラックス
std::vector<std::vector<double> > rhoUex_wall(ni+2,std::vector<double>(nj+2,0.0)); //電子フラックス
std::vector<std::vector<double> > rhoUer_wall(ni+2,std::vector<double>(nj+2,0.0)); //電子フラックス
std::vector<std::vector<double> > rhoUmx_wall(ni+2,std::vector<double>(nj+2,0.0)); //準安定フラックス
std::vector<std::vector<double> > rhoUmr_wall(ni+2,std::vector<double>(nj+2,0.0)); //準安定フラックス
std::vector<std::vector<double> > rhoUnx_wall(ni+2,std::vector<double>(nj+2,0.0)); //基底フラックス
std::vector<std::vector<double> > rhoUnr_wall(ni+2,std::vector<double>(nj+2,0.0)); //基底フラックス
std::vector<std::vector<double> > Gx_wall(ni+2,std::vector<double>(nj+2,0.0)); //エネルギーフラックス
std::vector<std::vector<double> > Gr_wall(ni+2,std::vector<double>(nj+2,0.0)); //エネルギーフラックス

std::vector<std::vector<double> > Uix(ni+2,std::vector<double>(nj+2,0.0)); //速度
std::vector<std::vector<double> > Uix_old(ni+2,std::vector<double>(nj+2,0.0)); //速度
std::vector<std::vector<double> > rhoUix(ni+2,std::vector<double>(nj+2,0.0)); //速度


std::vector<std::vector<double> > Uir(ni+2,std::vector<double>(nj+2,0.0)); //速度
std::vector<std::vector<double> > Uir_old(ni+2,std::vector<double>(nj+2,0.0)); //速度
std::vector<std::vector<double> > rhoUir(ni+2,std::vector<double>(nj+2,0.0)); //速度

std::vector<std::vector<double> > Uip(ni+2,std::vector<double>(nj+2,0.0)); //速度
std::vector<std::vector<double> > Uip_old(ni+2,std::vector<double>(nj+2,0.0)); //速度

std::vector<std::vector<double> > limit(ni+2,std::vector<double>(nj+2,0.0)); //リミッター

//std::vector<std::vector<double> > fix1(ni+2,std::vector<double>(nj+2,0.0)); //イオンフラックスベクトル
//std::vector<std::vector<double> > fix2(ni+2,std::vector<double>(nj+2,0.0)); //イオンフラックスベクトル
//std::vector<std::vector<double> > fix3(ni+2,std::vector<double>(nj+2,0.0)); //イオンフラックスベクトル
//std::vector<std::vector<double> > fix4(ni+2,std::vector<double>(nj+2,0.0)); //イオンフラックスベクトル
//std::vector<std::vector<double> > fir1(ni+2,std::vector<double>(nj+2,0.0)); //イオンフラックスベクトル
//std::vector<std::vector<double> > fir2(ni+2,std::vector<double>(nj+2,0.0)); //イオンフラックスベクトル
//std::vector<std::vector<double> > fir3(ni+2,std::vector<double>(nj+2,0.0)); //イオンフラックスベクトル
//std::vector<std::vector<double> > fir4(ni+2,std::vector<double>(nj+2,0.0)); //イオンフラックスベクトル

std::vector<std::vector<double> > fnx1(ni+2,std::vector<double>(nj+2,0.0)); //基底フラックスベクトル
std::vector<std::vector<double> > fnx2(ni+2,std::vector<double>(nj+2,0.0)); //基底フラックスベクトル
std::vector<std::vector<double> > fnx3(ni+2,std::vector<double>(nj+2,0.0)); //基底フラックスベクトル
std::vector<std::vector<double> > fnr1(ni+2,std::vector<double>(nj+2,0.0)); //基底フラックスベクトル
std::vector<std::vector<double> > fnr2(ni+2,std::vector<double>(nj+2,0.0)); //基底フラックスベクトル
std::vector<std::vector<double> > fnr3(ni+2,std::vector<double>(nj+2,0.0)); //基底フラックスベクトル

std::vector<std::vector<double> > Te(ni+2,std::vector<double>(nj+2,0.0));    //電子温度
std::vector<std::vector<double> > Te_old(ni+2,std::vector<double>(nj+2,0.0));    //電子温度
std::vector<std::vector<double> > Te_Lx(ni+2,std::vector<double>(nj+2,0.0));  //電子温度
std::vector<std::vector<double> > Te_Rx(ni+2,std::vector<double>(nj+2,0.0)); //電子温度
std::vector<std::vector<double> > Te_Lr(ni+2,std::vector<double>(nj+2,0.0)); //電子温度
std::vector<std::vector<double> > Te_Rr(ni+2,std::vector<double>(nj+2,0.0)); //電子温度

std::vector<std::vector<double> > qex(ni+2,std::vector<double>(nj+2,0.0));    //熱流束
std::vector<std::vector<double> > qer(ni+2,std::vector<double>(nj+2,0.0));    //熱流束

std::vector<std::vector<double> > Gx(ni+2,std::vector<double>(nj+2,0.0));  //エネルギーフラックス
std::vector<std::vector<double> > Gr(ni+2,std::vector<double>(nj+2,0.0));  //エネルギーフラックス
std::vector<std::vector<double> > Gx_old(ni+2,std::vector<double>(nj+2,0.0));  //エネルギーフラックス
std::vector<std::vector<double> > Gr_old(ni+2,std::vector<double>(nj+2,0.0));  //エネルギーフラックス

std::vector<std::vector<double> > Pabs(ni+2,std::vector<double>(nj+2,0.0)); //マイクロ波電力吸収

std::vector<std::vector<double> > phi(ni+2,std::vector<double>(nj+2,0.0));//両極性電位
std::vector<std::vector<double> > phi_old(ni+2,std::vector<double>(nj+2,0.0));//両極性電位
std::vector<std::vector<double> > Ex(ni+2,std::vector<double>(nj+2,0.0));//両極性電場
std::vector<std::vector<double> > Er(ni+2,std::vector<double>(nj+2,0.0));//両極性電場

std::vector<std::vector<double> > rhoUex(ni+2,std::vector<double>(nj+2,0.0)); //電子フラックス
std::vector<std::vector<double> > rhoUer(ni+2,std::vector<double>(nj+2,0.0)); //電子フラックス
std::vector<std::vector<double> > rhoUex_old(ni+2,std::vector<double>(nj+2,0.0)); //電子フラックス
std::vector<std::vector<double> > rhoUer_old(ni+2,std::vector<double>(nj+2,0.0)); //電子フラックス

std::vector<std::vector<double> > nUex(ni+2,std::vector<double>(nj+2,0.0)); //電子フラックス (サブ)
std::vector<std::vector<double> > nUer(ni+2,std::vector<double>(nj+2,0.0)); //電子フラックス (サブ)
std::vector<std::vector<double> > nUex_old(ni+2,std::vector<double>(nj+2,0.0)); //電子フラックス (サブ)
std::vector<std::vector<double> > nUer_old(ni+2,std::vector<double>(nj+2,0.0)); //電子フラックス (サブ)

std::vector<std::vector<double> > rhoUmx(ni+2,std::vector<double>(nj+2,0.0)); //準安定フラックス
std::vector<std::vector<double> > rhoUmr(ni+2,std::vector<double>(nj+2,0.0)); //準安定フラックス
std::vector<std::vector<double> > rhoUmx_old(ni+2,std::vector<double>(nj+2,0.0)); //準安定フラックス
std::vector<std::vector<double> > rhoUmr_old(ni+2,std::vector<double>(nj+2,0.0)); //準安定フラックス

std::vector<std::vector<double> > rate_ionize(ni+2,std::vector<double>(nj+2,0.0)); //イオン化レート
std::vector<std::vector<double> > Mache(ni+2,std::vector<double>(nj+2,0.0)); //電子マッハ数


//衝突周波数
std::vector<std::vector<double> > nu_m(ni+2,std::vector<double>(nj+2,0.0)); //電子運動量移動周波数
std::vector<std::vector<double> > nu_m1(ni+2,std::vector<double>(nj+2,0.0)); //電子運動量移動周波数 (有効衝突追加)
std::vector<std::vector<double> > nui_m(ni+2,std::vector<double>(nj+2,0.0)); //イオン衝突周波数
std::vector<std::vector<double> > nu_elas(ni+2,std::vector<double>(nj+2,0.0)); //弾性周波数
std::vector<std::vector<double> > nu_super(ni+2,std::vector<double>(nj+2,0.0)); //超弾性周波数
std::vector<std::vector<double> > nu_ionz(ni+2,std::vector<double>(nj+2,0.0)); //電離衝突周波数
std::vector<std::vector<double> > nu_ionzStep(ni+2,std::vector<double>(nj+2,0.0)); //段階電離衝突周波数
std::vector<std::vector<double> > nu_exc(ni+2,std::vector<double>(nj+2,0.0)); //励起周波数
std::vector<std::vector<double> > nu_excReso(ni+2,std::vector<double>(nj+2,0.0)); //励起(準安定)周波数
std::vector<std::vector<double> > nu_excMeta(ni+2,std::vector<double>(nj+2,0.0)); //励起(共鳴)周波数
std::vector<std::vector<double> > nu_excStep(ni+2,std::vector<double>(nj+2,0.0)); //段階励起周波数
std::vector<std::vector<double> > nu_ano(ni+2,std::vector<double>(nj+2,0.0)); //異常衝突周波数
std::vector<std::vector<double> > nu_ano_IAT(ni+2,std::vector<double>(nj+2,0.0)); //異常衝突周波数



//エネルギー交換
std::vector<std::vector<double> > rate_eloss(ni+2,std::vector<double>(nj+2,0.0)); //電子エネルギーロス
std::vector<std::vector<double> > delossdeps(ni+2,std::vector<double>(nj+2,0.0)); //電子エネルギーロスの微係数

//移動度
std::vector<std::vector<double> > mu_para(ni+2,std::vector<double>(nj+2,0.0));
std::vector<std::vector<double> > mu_perp(ni+2,std::vector<double>(nj+2,0.0));
std::vector<std::vector<double> > mu_min(ni+2,std::vector<double>(nj+2,0.0));
std::vector<std::vector<double> > mu_xx(ni+2,std::vector<double>(nj+2,0.0));
std::vector<std::vector<double> > mu_rr(ni+2,std::vector<double>(nj+2,0.0));
std::vector<std::vector<double> > mu_xr(ni+2,std::vector<double>(nj+2,0.0));
std::vector<std::vector<double> > mu_inv_xx(ni+2,std::vector<double>(nj+2,0.0));
std::vector<std::vector<double> > mu_inv_rr(ni+2,std::vector<double>(nj+2,0.0));
std::vector<std::vector<double> > mu_inv_xr(ni+2,std::vector<double>(nj+2,0.0));
std::vector<std::vector<double> > kappa_para(ni+2,std::vector<double>(nj+2,0.0));
std::vector<std::vector<double> > kappa_perp(ni+2,std::vector<double>(nj+2,0.0));
std::vector<std::vector<double> > kappa_x(ni+2,std::vector<double>(nj+2,0.0));
std::vector<std::vector<double> > kappa_r(ni+2,std::vector<double>(nj+2,0.0));
std::vector<std::vector<double> > kappa_c(ni+2,std::vector<double>(nj+2,0.0));
std::vector<std::vector<double> > kappa_x_inv(ni+2,std::vector<double>(nj+2,0.0));
std::vector<std::vector<double> > kappa_r_inv(ni+2,std::vector<double>(nj+2,0.0));
std::vector<std::vector<double> > kappa_c_inv(ni+2,std::vector<double>(nj+2,0.0));
std::vector<std::vector<double> > cx(ni+2,std::vector<double>(nj+2,0.0));
std::vector<std::vector<double> > cr(ni+2,std::vector<double>(nj+2,0.0));
std::vector<std::vector<double> > sx(ni+2,std::vector<double>(nj+2,0.0));
std::vector<std::vector<double> > sr(ni+2,std::vector<double>(nj+2,0.0));

//固定磁場
std::vector<std::vector<double> > Bx(ni+2,std::vector<double>(nj+2,0.0));
std::vector<std::vector<double> > Br(ni+2,std::vector<double>(nj+2,0.0));
std::vector<std::vector<double> > Ap(ni+2,std::vector<double>(nj+2,0.0)); //ベクトルポテンシャル
std::vector<std::vector<double> > distECR(ni+2,std::vector<double>(nj+2,0.0)); //ECR領域からの距離

//テスト配列
std::vector<std::vector<double> > test1(ni+2,std::vector<double>(nj+2,0.0));
std::vector<std::vector<double> > test2(ni+2,std::vector<double>(nj+2,0.0));
std::vector<std::vector<double> > test3(ni+2,std::vector<double>(nj+2,0.0));

//配列 BOLSIG+
std::vector<double> nu_elas_N; //弾性周波数/nn <σelas*vth>
std::vector<double> nu_super_N; //超弾性周波数/nn <σsuper*vth>
std::vector<double> nu_ionz_N; //イオン化周波数/nn <σionz*vth>
std::vector<double> nu_ionzStep_N; //段階イオン化周波数/nn <σionz*vth>
std::vector<double> nu_exc_N; //励起周波数/nn <σexc*vth>
std::vector<double> nu_excReso_N; //共鳴励起周波数/nn <σexc*vth>
std::vector<double> nu_excMeta_N; //準安定励起周波数/nn <σexc*vth>
std::vector<double> nu_excStep_N; //段階励起周波数/nn <σexc*vth>
std::vector<double> rate_eloss_n_N; //エネルギーロスレート/nn
std::vector<double> rate_eloss_m_N; //エネルギーロスレート/nn
std::vector<double> nui_m_N; //イオン運動量移動周波数/nn <σionz*vth>

//配列 EI-SEE table
std::vector<double> coefEISEE_eb_table; //EI elastic-backscattering
std::vector<double> coefEISEE_rd_table; //EI rediffusion
std::vector<double> coefEISEE_ts_table; //EI true-secondary
