/********************************************************************
Simulation of train running and feeder circuit operation

ver13_2
1 train, 4 substation

ver14
Simulate Takasaki line
(With Diagram function)

ver15
Simulate Takasaki line
(With Diagram function)
With Wayside ESD
並び替えのアルゴリズムを見直し（OFFのSSを片っ端から右端に移動させる方式）

ver16
Simulate Takasaki line
(With Diagram function)
With Wayside ESD & regenerative INV

ver16_1
Simulate Takasaki line
(Without Diagram function)
With Wayside ESD & regenerative INV

ver16_2
Simulate Takasaki line
(Without Diagram function)
With Wayside ESD & regenerative INV

ver16_3
Simulate Takasaki line
(Without Diagram function)
With Wayside ESD & regenerative INV
With Delay of departure

ver16_4
Simulate Takasaki line
(Without Diagram function)
With Wayside ESD(time limitation included) & regenerative INV
With Delay of departure
With INV-ESD interruction investigated


ver17
Simulate Takasaki line
(Without Diagram function)
With Wayside ESD(time limitation included) & regenerative INV
With Delay of departure
With INV-ESD interruction investigated
パラメータ自動変更機能追加？（）
ESD: V_charge_ESD, V_discharge_ESD
INV: V_charge_INV (V_fin)

ver18　小林先生のものを書き直す
小林先生のは元々，回生ブレーキ力b_reg_c[kN]は計算してないけど，回生ブレーキ力b_reg[N]は計算しててそれを回生エネルギーの計算に用いていた

ノッチ5力行時 v1_P=40,v2_P=63.82549701  i1d=85.7418508,iq=165.6086836,ws=11.28173426,Fmax=13.90011356

回生時v1_R=50.52,i1d=103.2188528,i1q=-137.5678441,ws=-7.784729348
Bmax=13.90011356


ノッチ4回生時v1_R=56.31,i1d=102.7443702,i1q=-111.2973885,ws=-6.32721445

ノッチ3
回生時i1d=102.2117125,i1q=-84.83142392,ws=-4.84776665

d軸電流とq軸電流計算を追加した
そこからモータパワーを求めることで，インバータ損失を無視し，インバータ電流を計算
き電電圧変動時，力行時の引張力も変わるようにした
回生時の回生時のブレーキ力を変更(軽負荷回生制御により性能が変化)

ver19　2022.03.27
小田急さんとの共同研究で小田急多摩線での走行シミュレーションをすることになった(多分)
とりあえず新百合ヶ丘～唐木田間のものに適用させてみる

********************************************************************/


//char dir[]= "C:\Users\宏泰\Desktop\共同研究\JRE\検討\2018_10\\";


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "mkl.h"
#include <string>
#include <sstream>
#include <iomanip>
#include <windows.h>
using namespace std;

/***  シミュレーション刻み時間  ***/
//#define dt 0.01				//刻み時間:10ms

#define dt 0.001			//刻み時間:1ms

//#define dt 0.0005			//刻み時間:500μs

//#define dt 0.005
//#define dt 0.0001			//刻み時間:100μs	//これより長井とエラーが起きる(CFCが小さいため)
//#define dt 0.00005		//刻み時間:50μs
//#define dt 0.00001		//刻み時間:10μs

double t, minute;

//#define NUM_tra_UP 1								//Number of UP trains(+1)
//#define NUM_tre_DOWN 1								//Number of DOWN trains(-1)
//#define NUM_tra (NUM_tra_UP+NUM_tre_DOWN)			//Number of trains

#define NUM_tra 3										//Number of trains
#define NUM_sub 2										//Number of substations
#define NUM_esd 0										//Number of ESDs
#define NUM_inv 0										//Number of regenerative inverters
#define NUM (NUM_tra + NUM_sub + NUM_esd + NUM_inv - 1)						//Number of elements
#define MAX_tra 21
#define N_node (NUM_tra+NUM_sub)					//Number of nodes
#define N_branch (NUM_sub +(NUM_sub-1)+NUM_tra)	//Number of branches 単線
#define Speed_change 0.0

#define NUM_station 2
#define NUM_final_station 2			//各駅停車　(駅数 - 1 が正解？)

#define NUM_station2 4				//急行(駅数 - 1 が正解？)
#define NUM_final_station2 4

#define NUM_station3 2				//急行(駅数 - 1 が正解？)
#define NUM_final_station3 2

#define Ratio_E	3600.0*1000.0							//kwh⇒Jへの変換比

#define T_charge 50.0
#define T_discharge 40.0
#define T_ope 180.0

#define Time_diagram 4
#define F_flag 2 //列車別ファイルを作成するかのフラグ (0なら作成, 1なら3次元グラフ作成(力行電力), 2なら回生電力、3なら力行ー回生)
#define ALL_Vss_Change 0 //SS電圧を独立に変更するかのフラグ (0ならSS電圧を独立に変更)
#define Consider_loop_flag 1 //1ならループしない





//#define Vss 1650.0						//無負荷時変電所出力電圧[V]
//#define Vss 1600.0						//無負荷時変電所出力電圧[V]
#define Vss 1550.0					//無負荷時変電所出力電圧[V]
#define Vss1 1580					//無負荷時変電所出力電圧[V]
#define Vss2 1570					//無負荷時変電所出力電圧[V]
#define Vss3 1580					//無負荷時変電所出力電圧[V]
//#define Vss 1500.0						//無負荷時変電所出力電圧[V]
#define Iss_rated 8000.0						//変電所定格電流[A]
#define Iss_rated2 5333.3

#define ess1 8.0						//変電所出力電圧変動率[%]
//#define ess1 5.89						//変電所出力電圧変動率[%]
#define ess2 6.0						//変電所出力電圧変動率[%]
//#define ess3 7.89						//変電所出力電圧変動率[%]

#define R	0.0345					//き電抵抗[Ω/km]
#define R1	0.0001				//き電抵抗[Ω/km]
//#define R_StoN_DOWN	0.03423					//き電抵抗[Ω/km]
//#define R_NtoK_UP	0.04018					//き電抵抗[Ω/km]
//#define R_NtoK_DOWN	0.04061					//き電抵抗[Ω/km]

//#define C_FC (10.05 / 1000.0)		//車両FC容量[F](5M分)
//#define L_FL 0.0080				//FL巻線[H](5M分)
//#define R_FL 0.0348				//FL巻線抵抗[Ω](4M分)
//#define C_FC 0.015		//車両FC容量[F](5M分)
//#define L_FL 0.019				//FL巻線[H](5M分)
#define R_FL 0.100				//FL巻線抵抗[Ω](5M分)


#define Vclim	1780.0
#define Vcmax	1830.0

#define g 9.80665						//重力加速度[m/s/s]

#define speed_limit 90.0				//制限速度[km/h]
#define speed_min   50.0                //再加速速度[km/h]

#define stoptime 20.0					//停車時間[s]

//#define P_SIV 142.0*3.0*1000.0/10.0				//補機電力[W]（1両当たり）
//#define P_SIV 300.0*1000.0/10.0				//補機電力[W]（1両当たり）
//#define P_SIV 20.0*1000.0/10.0
#define P_SIV 10*1000

#define P_teifuka 3000.0*1500.0
//#define P_SHINYURIGAOKA 3000.0*1500.0

#define count_lap 1

#define rownum (int)(speed_limit+15)*20

/****** 各種効率 ******/
#define e_gear 0.98		//ギア効率
#define e_motor 0.92	//モーター効率
#define e_inv 0.95		//インバータ効率
/**********************/


/////スイッチング素子の損失計算用パラメータ/////////////////
/*** CM1200HC-66H ***/
#define t_on 1.6/1000.0/1000.0
#define t_off 2.5/1000.0/1000.0
#define Vce 3.6
#define Vf 2.7
#define Irp 1200.0
#define t_rr 1.4/1000.0/1000.0
#define fc 1000.0


double t_out;		//データ間引き関係
double t_out_old;
double t_out2;		//シミュレーション進捗表示関係
double t_out2_old;
double t_ope;

#define vss_pattern 5
#define vel_pattern 3
#define vlb_pattern 1

double vss_dif[vss_pattern] = { -20.0, -10.0, 0.0, 10.0, 20.0 };
//double vss_dif[vss_pattern] = { 40.0, 20.0 };
//double vss_dif[vss_pattern] = { 0.0 };
//double vss_dif[vss_pattern] = { 100.0, 50.0, 0.0, -50.0, -100.0 };
//double vss_dif[vss_pattern] = { 90.0, 60.0, 30.0, 0.0, -30.0, -60.0, -90.0 };
//double vss_dif[vss_pattern] = { 120.0, 90.0, 60.0, 30.0, 0.0, -30.0, -60.0, -90.0, -120.0 };

//double vss_dif[vss_pattern] = { 120.0 };

double vel_dif[vel_pattern] = { 0.5, 0.0, -0.5 };
//double vel_dif[vel_pattern] = { 0.5, 0.0 };
//double vel_dif[vel_pattern] = { 0.0 };

//double vlb_dif[vlb_pattern] = { -60.0, -30.0, 0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0 };
double vlb_dif[vlb_pattern] = { 30.0 };

string vss_fol1[vss_pattern];
string vss_fol2[vss_pattern];
string vss_fol3[vss_pattern];
string vss_fol4[vss_pattern];
string vss_fol5[vss_pattern];
string vel_fol[vel_pattern];
string vlb_fol[vlb_pattern];

int vss_state1 = 0;
int vss_state2 = 0;
int vss_state3 = 0;
int vss_state4 = 0;
int vss_state5 = 0;
int vss_state6 = 0;
int vss_state7 = 0;

int vel_state = 0;
int vlb_state = 0;

int loop_endflag = 0;
int sim_endflag = 0;
int test_flag = 0;
int sim_loopcount = 0;
int end_flag;

int Flag_invest;		//INV-ESD間潮流検出

typedef struct {
	/**  定数  **/
	int name_train;			//列車の名前
	double mass;			//列車総重量[kg]
	double mass_c;			//列車総重量[t]
	double mass_M;			//M車重量[kg]
	double mass_T;			//T車重量[kg]
	double mass_P;			//荷重重量[kg]
	double num;				//両数
	double num_M;			//M車の両数
	double num_T;			//T車の両数
	double num_mot;			//モータ数
	double v0_P;           //(力行時)V/f終端速度
	double v1_P;			//(力行時)定出力領域に入る速度(5ノッチ)
	double v2_P;			//(力行時)特性領域に入る速度(5ノッチ)
	double v0_R;           //(回生時)V/f終端速度
	double v1_R;			//(回生時)定出力領域に入る速度(5ノッチ)
	double v2_R;			//(回生時)特性領域に入る速度(5ノッチ)
	double Fmax;			//定トルク領域の引張力(力行時)
	double Bmax;			//定トルク領域の引張力(回生時)
	int  direction;			//上り・下り判別(1:大宮⇒吹上方面，-1:吹上⇒大宮方面)
	//	double UP_DOWN;
	double B;				//常用最大減速度[m/s/s]
	double Bref;			//減速度指令[m/s/s]
	double BEref;			//ブレーキ力指令[N]
	double BEmax;			//最大ブレーキ力[N]
	double Cfc;//FC容量
	double Lfl;//FLインダクタンス
	double Rfl;//Fl抵抗

	//	double theta;
	//	double theta_old;

	double notch;				//ノッチ入力(Powering: 1～4)
	double brake;				//ノッチ入力(Regenerating: -1～-7)
	int Type;				//列車種別，0:各駅停車，1:急行
	int flag_direction;		//列車行先
	//	int t_count;			//試運転車運転方法
	//	int c_count;			//定速運転カウント
	//	double con_speed;		//定速運転速度

		/**  変数  **/
	int accelflag;			//accelflag=1:力行, accelflag=2:惰行, accelflag=3:減速, accelflag=4:停車, accelflag=5:速度調整
	int brakeflag;			//brakeflag=1:駅に止まるために減速中, brakeflag=0:それ以外
	int laststopflag;		//laststopflag=1:終点に到着, laststopflag=0:それ以外
	double X_recentstop;		//直前の停車位置[m]
	double X_nextstop;		//次の停車位置[m]
	double X_brake;			//ブレーキ開始位置[m]
	double x;				//位置[m]
	double x_c;				//位置[km]
	double v;				//速度[m/s]
	double v_c;				//速度[km/h]
	double v_new;
	double a;				//加速度[m/s/s]
	double a_c;				//加速度[km/h/s]
	double Fmot;			//引張力[N]
	double Fmot_c;			//引張力[kN]
	double Fmot_regM;		//最大回生ブレーキ力[N]
	double Fmot_regM_c;		//最大回生ブレーキ力[kN]
	double Ftot;			//トータル引張力[N]
	double Ftot_c;			//トータル引張力[kN]
	double b_reg;			//回生ブレーキ力[N]
	double b_reg_c;			//回生ブレーキ力[kN]
	double b_regF_c;
	double b_reg_judge;		//軽負荷回生制御パターンを満たしているかの判定用
	double judgement;		//繰り返すか否か
	double b_air;			//機械ブレーキ力[N]
	double b_air_c;			//機械ブレーキ力[kN]
	//double b_reg_loss;			//回生ブレーキ力[N]
	//double b_reg_loss_c;			//回生ブレーキ力[kN]

//	double i1d_P; //力行時，定トルク領域のd軸電流
//	double i1q_P; //力行時，定トルク領域のq軸電流
//	double i1d_R; //回生時，定トルク領域のd軸電流
//	double i1q_R; //回生時，定トルク領域のq軸電流
//	double i1q_R_1;//回生時，定出力領域がない特性で最大となる電流
//	double R1;//1次抵抗[Ω]
//	double R2;//2次抵抗[Ω]
//	double L1;//1次自己インダクタンス[H]
//	double L2;//2次自己インダクタンス[H]
//	double MORE1;//一次漏れインダクタンス[H]	
//	double MORE2;//二次漏れインダクタンス[H]
//	double g0;//励磁コンダクタンス[S]
//	double M_IM;//相互インダクタンス[H]
//	double SIGMA;//漏れ係数
//	double Gr;//ギア比
//	double rd;//車輪半径[m]
//	double POLE;//極対数
	double Tacr;//電流制御時定数

	//	double id_mot_P;       //力行時のd軸電流[A]
	//	double iq_mot_P;       //回生時のq軸電流[A]
	//	double id_mot_R;       //力行時のd軸電流[A]
	//	double iq_mot_R;       //回生時のq軸電流[A]
	//	double idtot;          //力行・惰行・回生時のd軸電流[A]
	//	double iqtot;          //力行・惰行・回生時のq軸電流[A]
	//	double idtotF;          //力行・惰行・回生時のd軸電流[A]
	//	double iqtotF;          //力行・惰行・回生時のq軸電流[A]
	//	double iq_reg;         //回生時のq軸電流(軽負荷回生制御込み)[A]
	//	double iq_regF;         //回生時のq軸電流(軽負荷回生制御込み)[A]
	double b_regF;
	double b_regF_old;
	double b_regF_old_L;
	double FmotF;
	//	double wr;//角周波数[rad/s]
	//	double ws;//すべり角周波数[rad/s]
	//	double we;//インバータ角周波数[rad/s]
	//	double v1d;//d軸電圧[V]
	//	double v1q;//q軸電圧[V]
	//	double Em;//インバータ出力電圧[V]

	double Pmot;			//モータパワー[W]
	double Pmot_air;		//機械ブレーキパワー[W]
	double Pmot_c;			//モータパワー[kW]
	double Pmot_air_c;		//機械ブレーキパワー[kW]

	//double Pmot_reg_loss;
	//double Pmot_reg_loss_c;

	double Ptot;			//車両パワー[W](モータパワー+機械ブレーキ分)
	double Ptot_c;			//車両パワー[kW](モータパワー+機械ブレーキ分)

	double Preg_nashi;
	double Preg_nashi_c;

	double Pres;			//走行抵抗での損失[W]
	double Pres_c;			//走行抵抗での損失[kW]

	double Pinv_on;			//スイッチング素子のオン損失[W]
	double Pinv_diode;		//ダイオードの導通損失[W]
	double Pinv_sw;			//スイッチング損失[W]
	double Pinv_rec;		//ダイオードのリカバリー損失[W]

	double Ploss_inv;		//インバータ損失[W]
	double Ploss_mot;		//モータ銅損[W]
	double Ploss_fe;		//モータ鉄損[W]
	double Ploss_fl;		//FL損失[W]
	double Ploss_all;		//総主回路損失[W]

	double Ploss_inv_c;
	double Ploss_mot_c;
	double Ploss_fl_c;
	double Ploss_all_c;

	double Pvh;				//車両入出力パワー[W] (Pvh_in + Pvh_out)
	double Pvh_in;			//車両入力パワー[W]
	double Pvh_out;			//車両出力パワー[W]
	double Pvh_cos;			//車両入力パワー(惰行時)[W]
	double Pvh_st;			//車両入力パワー(惰行、停車時)[W]
	double Pvh_c;			//車両入出力パワー[kW]
	double Pvh_in_c;		//車両入力パワー[kW]
	double Pvh_out_c;		//車両出力パワー[kW]
	double Pvh_st_c;			//車両入力パワー(惰行、停車時)[kW]

	double Pvh_in1;
	double Pvh_in2;
	double Pvh_in3;
	double Pvh_in4;
	double Pvh_in5;
	double Pvh_in6;
	double Pvh_in7;
	double Pvh_in8;

	double Pvh_out1;
	double Pvh_out2;
	double Pvh_out3;
	double Pvh_out4;
	double Pvh_out5;
	double Pvh_out6;
	double Pvh_out7;
	double Pvh_out8;

	double vfc;				//FC電圧[V]
	double iinv;				//インバータ電流[A]
	double ifc;				//FC電流[A]
	double ifl;				//FL電流[A]
	double ifl_sum;
	double vfc_old;				//FC電圧一周回前[A]
	double isiv;			//FL電流[A]
	double vp;				//車両パンタ点電圧[V]
	double vp_old;				//車両パンタ点電圧[V]
	double R_run;			//走行抵抗[N/t]
	double R_grade;			//勾配抵抗[N/t]
	double R_curve;			//曲線抵抗[N/t]
	double R_total;			//全走行抵抗[N]

	/** エネルギー計算関係 **/
	double Emot_pow;		//モータ力行エネルギー[J]
	double Emot_reg;		//モータ回生エネルギー[J]
	double Emot_air;		//機械ブレーキエネルギー[J]
	double Emot_reg_loss;		//モータ回生絞り込みエネルギー[J]

	double Emot_pow_c;		//モータ力行エネルギー[kWh]
	double Emot_reg_c;		//モータ回生エネルギー[kWh]
	double Emot_air_c;		//機械ブレーキエネルギー[kWh]
	double Emot_reg_loss_c;		//モータ回生絞り込みエネルギー[kWh]

	double Eloss_inv_pow;	//インバータ損失エネルギー(力行時)[J]
	double Eloss_mot_pow;	//モータ損失エネルギー(力行時)[J]
	double Eloss_fl_pow;	//FL損失エネルギー(力行時)[J]
	double Eloss_all_pow;	//総主回路損失エネルギー(力行時)[J]

	double Eloss_inv_reg;	//インバータ損失(回生時)[J]
	double Eloss_mot_reg;	//モータ損失(回生時)[J]
	double Eloss_fl_reg;	//FL損失(回生時)[J]
	double Eloss_all_reg;	//主回路損失(回生時)[J]

	double Eres_pow;		//走行抵抗での損失(力行時)[J]
	double Eres_reg;		//走行抵抗での損失(回生時)[J]
	double Eres_coa;		//走行抵抗での損失(惰行時)[J]

	double Esiv_pow;	//SIVの消費エネルギー(力行時)[J]
	double Esiv_reg;	//SIVの消費エネルギー(回生時)[J]
	double Esiv_coa;	//SIVの消費エネルギー(惰行時)[J]
	double Esiv_stp;	//SIVの消費エネルギー(停車時)[J]

	double Ereg_nashi;
	double Ereg_nashi_c;

	double Eloss_inv_pow_c;
	double Eloss_mot_pow_c;
	double Eloss_fl_pow_c;
	double Eloss_all_pow_c;

	double Eloss_inv_reg_c;
	double Eloss_mot_reg_c;
	double Eloss_fl_reg_c;
	double Eloss_all_reg_c;

	double Eres_pow_c;
	double Eres_reg_c;
	double Eres_coa_c;

	double Esiv_pow_c;		//SIVの消費エネルギー(力行時)[kWh]
	double Esiv_reg_c;		//SIVの消費エネルギー(回生時)[kWh]
	double Esiv_coa_c;		//SIVの消費エネルギー(惰行時)[kWh]
	double Esiv_stp_c;		//SIVの消費エネルギー(停車時)[kWh]

	double Evh;				//列車入出力エネルギー[J]
	double Evh_in;			//列車入力エネルギー[J]
	double Evh_out;			//列車出力エネルギー[J]
	double Evh_cos;
	double Evh_st;          //列車入力エネルギー(惰行、停車時)[J]

	double Evh_c;			//列車入出力エネルギー[kWh]
	double Evh_in_c;		//列車入力エネルギー[kWh]
	double Evh_out_c;		//列車出力エネルギー[kWh]
	double Evh_cos_c;
	double Evh_st_c;		//列車入力エネルギー(惰行、停車時)[kWh]


	double Evh_in1;
	double Evh_in2;
	double Evh_in3;
	double Evh_in4;
	double Evh_in5;
	double Evh_in6;
	double Evh_in7;
	double Evh_in8;
	double Evh_in1_c;
	double Evh_in2_c;
	double Evh_in3_c;
	double Evh_in4_c;
	double Evh_in5_c;
	double Evh_in6_c;
	double Evh_in7_c;
	double Evh_in8_c;

	double Evh_out1;
	double Evh_out2;
	double Evh_out3;
	double Evh_out4;
	double Evh_out5;
	double Evh_out6;
	double Evh_out7;
	double Evh_out8;
	double Evh_out1_c;
	double Evh_out2_c;
	double Evh_out3_c;
	double Evh_out4_c;
	double Evh_out5_c;
	double Evh_out6_c;
	double Evh_out7_c;
	double Evh_out8_c;


	/** 列車停止位置計算関係 **/
	int flag_station;			//停車駅数
	int startflag;				//発車指令
	/**  列車停車時間関係  **/
	double t_stop;				//
	double t_stop_old;			//駅到達時刻[s]
	//	double t_wait;
	double route[rownum][2];	//ブレーキルート{位置,速度}
	int diaflag;				//運行判定
	int lap;					//周回数
	int lapmax;					//目標周回数
	int wait_flag;				//待ち合わせフラグ
	int reg_flag;				//軽負荷回生制御フラグ

	double X_nextstop_old;		//次の停車位置[m]
	double Speedlimit;			//制限速度[km/h]
	double Reaccelspeed;
	//	double Reaccelspeed;			//制限速度[km/h]

	double T_delay;				//発車タイミング遅れ考慮[s]
	double Vst;					//速度引張力特性ベース電圧

	double run_time;
	//double ave_runtime;

} TRAIN;

typedef struct {
	/** 定数 **/
	int name_SS;		//変電所の名前
	double Vss_0;		//無負荷時変電所出力電圧
	double Iss_0;		//定格出力電流
	double e_ss;		//変電所出力電圧変動率
	double Xss;			//変電所位置
	double Jss_ini;		//変電所等価電流源[初期値]
	double Rss_ini;			//変電所等価出力抵抗[初期値]

	/** 変数 **/
	int diode;			//変電所ダイオードON/OFF判定フラグ(1:ON 0:OFF)
	double Rss;			//変電所等価出力抵抗[Ω]
	double Jss;			//変電所等価電流源[A](更新用)

	double vss;			//変電所出力電圧[V]
	double vss_e;		//変電所出力電圧[V]
	double iss;			//変電所出力電流[A]
	double Pss;			//変電所出力パワー[W]
	double Pss_c;		//変電所出力パワー[kW]
	double vout;		//変電所母線電圧[V]
	int flag;			//フラグ

	/** エネルギー計算関係 **/
	double Ess;			//変電所出力エネルギー[J]
	double Ess_c;		//変電所出力エネルギー[kWh]
	double Wss;			//変電所損失[J]
	double Wss_c;		//変電所損失[kWh]
} SUB;


typedef struct {
	/**  定数  **/
	char name_station[20];	//駅の名前
	double Xs;	//距離
} STATION;

typedef struct {
	double mass_M;
	double mass_T;
	double num_M;
	double num_T;
	double x;
	double v;
	int direction;
	double T_delay;
	int flag_direction;
	int Type;
	int wait_flag;
} INI_TRA;

typedef struct {
	double Xss;
	double e_ss;
	double Iss;
	double vss;
	double iss;
	int diode;
} INI_SUB;

typedef struct {
	int Number;			//ノード番号
	double X;			//位置[m]（基準点からの距離）
	double I;			//ノード電流[A]
	double V;			//ノード電圧[V]
	double r;			//ノード抵抗値[Ω](Only SS)
	int flag;		//ノードの分類（1：UP, -1：DOWN, 0：SS）
} NODE;

typedef struct {
	double X_start;		//branch始端[m]（基準点からの距離）
	double X_end;		//branch終端[m]（基準点からの距離）
	double r;			//Resistance [Ohm]
	int flag;		//ブランチの分類（1：UP, -1：DOWN, 0：SS）
	int Node_pos;	//正方向側に接続されているノード№
	int Node_neg;	//負方向側に接続されているノード№
} BRANCH;


void Make_file1(FILE** fp, int i)
{
	time_t start;
	errno_t error;
	struct tm today;
	char filename[128];

	time(&start);
	error = localtime_s(&today, &start);

	sprintf_s(filename, sizeof(filename), "Result\\錦糸町vss%.1f\\平井vss%.1f\\v%.1f\\%d月%d日%d時%d分%d秒_Train_%d.csv", vss_dif[vss_state1], vss_dif[vss_state2], vel_dif[vel_state], today.tm_mon + 1, today.tm_mday, today.tm_hour, today.tm_min, today.tm_sec, i);		//名前の作成
	//sprintf_s(filename, sizeof(filename), "Result\\Train_%d%d@%d.%d.%d_%d.csv", today.tm_mon + 1, today.tm_mday, today.tm_hour, today.tm_min, today.tm_sec, i);
	error = fopen_s(fp, filename, "w");
	//ファイルオープン
}

void Make_file2(FILE** fp)
{
	time_t start;
	errno_t error;
	struct tm today;
	char filename[128];

	time(&start);
	error = localtime_s(&today, &start);

	sprintf_s(filename, sizeof(filename), "Result\\錦糸町vss%.1f\\平井vss%.1f\\v%.1f\\%d月%d日%d時%d分%d秒_Train_consumption.csv", vss_dif[vss_state1], vss_dif[vss_state2], vel_dif[vel_state], today.tm_mon + 1, today.tm_mday, today.tm_hour, today.tm_min, today.tm_sec);		//名前の作成
	//sprintf_s(filename, sizeof(filename), "Result\\Train_consumption_%d%d@%d.%d.%d.csv", today.tm_mon + 1, today.tm_mday, today.tm_hour, today.tm_min, today.tm_sec);
	error = fopen_s(fp, filename, "w");						//ファイルオープン
}

void Make_file3(FILE** fp)
{
	time_t start;
	errno_t error;
	struct tm today;
	char filename[128];

	time(&start);
	error = localtime_s(&today, &start);

	sprintf_s(filename, sizeof(filename), "Result\\錦糸町vss%.1f\\平井vss%.1f\\v%.1f\\%d月%d日%d時%d分%d秒_substation.csv", vss_dif[vss_state1], vss_dif[vss_state2], vel_dif[vel_state], today.tm_mon + 1, today.tm_mday, today.tm_hour, today.tm_min, today.tm_sec);		//名前の作成
	//sprintf_s(filename, sizeof(filename), "Result\\substation_%d%d@%d.%d.%d.csv", today.tm_mon + 1, today.tm_mday, today.tm_hour, today.tm_min, today.tm_sec);
	error = fopen_s(fp, filename, "w");																														//ファイルオープン
}

void Make_file4(FILE** fp)
{
	time_t start;
	errno_t error;
	struct tm today;
	char filename[128];

	time(&start);
	error = localtime_s(&today, &start);

	sprintf_s(filename, sizeof(filename), "Result\\錦糸町vss%.1f\\平井vss%.1f\\v%.1f\\%d月%d日%d時%d分%d秒_substation_power.csv", vss_dif[vss_state1], vss_dif[vss_state2], vel_dif[vel_state], today.tm_mon + 1, today.tm_mday, today.tm_hour, today.tm_min, today.tm_sec);		//名前の作成
	//sprintf_s(filename, sizeof(filename), "Result\\substation_power_%d%d@%d.%d.%d.csv", today.tm_mon + 1, today.tm_mday, today.tm_hour, today.tm_min, today.tm_sec);
	error = fopen_s(fp, filename, "w");						//ファイルオープン
}

void Make_file5(FILE** fp)
{
	time_t start;
	errno_t error;
	struct tm today;
	char filename[128];

	time(&start);
	error = localtime_s(&today, &start);

	sprintf_s(filename, sizeof(filename), "Result\\錦糸町vss%.1f\\平井vss%.1f\\v%.1f\\%d月%d日%d時%d分%d秒_列車位置.csv", vss_dif[vss_state1], vss_dif[vss_state2], vel_dif[vel_state], today.tm_mon + 1, today.tm_mday, today.tm_hour, today.tm_min, today.tm_sec);		//名前の作成
	//sprintf_s(filename, sizeof(filename), "Result\\列車位置_%d月%d日%d時%d分%d秒.csv", today.tm_mon + 1, today.tm_mday, today.tm_hour, today.tm_min, today.tm_sec);
	error = fopen_s(fp, filename, "w");						//ファイルオープン
}

void Make_file6(FILE** fp)
{
	time_t start;
	errno_t error;
	struct tm today;
	char filename[128];

	time(&start);
	error = localtime_s(&today, &start);

	sprintf_s(filename, sizeof(filename), "Result\\%d月%d日%d時%d分%d秒_グラフ作成用.csv", today.tm_mon + 1, today.tm_mday, today.tm_hour, today.tm_min, today.tm_sec);		//名前の作成
	error = fopen_s(fp, filename, "w");						//ファイルオープン
}

void Make_file_x(FILE** fp)
{
	time_t start;
	errno_t error;
	struct tm today;
	char filename[128];

	time(&start);
	error = localtime_s(&today, &start);

	//sprintf_s(filename, sizeof(filename), "Result\\錦糸町vss%.1f\\v%.1f\\book_x.csv", vss_dif[vss_state3], vel_dif[vel_state]);		//名前の作成
	sprintf_s(filename, sizeof(filename), "Result\\book_x.csv"/*, vss_dif[vss_state1], vss_dif[vss_state2], vss_dif[vss_state3], vel_dif[vel_state]*/);
	error = fopen_s(fp, filename, "w");						//ファイルオープン
}

void Make_file_y(FILE** fp)
{
	time_t start;
	errno_t error;
	struct tm today;
	char filename[128];

	time(&start);
	error = localtime_s(&today, &start);

	//sprintf_s(filename, sizeof(filename), "Result\\錦糸町vss%.1f\\v%.1f\\book_y.csv", vss_dif[vss_state3], vel_dif[vel_state]);		//名前の作成
	sprintf_s(filename, sizeof(filename), "Result\\book_y.csv"/*, vss_dif[vss_state1], vss_dif[vss_state2], vss_dif[vss_state3], vel_dif[vel_state]*/);
	error = fopen_s(fp, filename, "w");						//ファイルオープン
}

void Make_file_z(FILE** fp)
{
	time_t start;
	errno_t error;
	struct tm today;
	char filename[128];

	time(&start);
	error = localtime_s(&today, &start);

	//sprintf_s(filename, sizeof(filename), "Result\\錦糸町vss%.1f\\v%.1f\\book_z.csv", vss_dif[vss_state3], vel_dif[vel_state]);		//名前の作成
	sprintf_s(filename, sizeof(filename), "Result\\book_z.csv"/*, vss_dif[vss_state1], vss_dif[vss_state2], vss_dif[vss_state3], vel_dif[vel_state]*/);
	error = fopen_s(fp, filename, "w");						//ファイルオープン
}


void Make_train(TRAIN** tra, INI_TRA* ini, int i)
{
	int m;
	int n;

	*tra = (TRAIN*)malloc(sizeof(TRAIN));
	if (*tra == NULL)
	{
		//エラー処理
	}

	(*tra)->name_train = i;
	(*tra)->Type = ini->Type;
	(*tra)->flag_direction = ini->flag_direction;
	//(*tra)->notch = 3.0;
	//(*tra)->brake = 5.0;

	(*tra)->mass_M = ini->mass_M * ini->num_M * 1.1 + 144.4666667 * 55.0 * 1.0 * ini->num_M;											//[kg](250%乗車)
	(*tra)->mass_T = ini->mass_T * ini->num_T * 1.05 + 144.4666667 * 55.0 * 1.0 * ini->num_T;											//[kg](250%乗車)
	//(*tra)->mass = ini->mass_M * ini->num_M + ini->mass_T * ini->num_T;				//[kg]
	//(*tra)->mass_c = (ini->mass_M * ini->num_M + ini->mass_T * ini->num_T) / 1000.0;										//[t](250%乗車)
	(*tra)->mass = ini->mass_M * ini->num_M * 1.1 + 144.4666667 * 55.0 * 1.0 * ini->num_M + ini->mass_T * ini->num_T * 1.05 + 144.4666667 * 55.0 * 1.0 * ini->num_T;				//[kg]
	(*tra)->mass_c = (ini->mass_M * ini->num_M * 1.1 + 144.4666667 * 55.0 * 1.0 * ini->num_M + ini->mass_T * ini->num_T * 1.05 + 144.4666667 * 55.0 * 1.0 * ini->num_T) / 1000.0;										//[t](250%乗車)
	(*tra)->num = ini->num_M + ini->num_T;
	(*tra)->num_mot = ini->num_M * 4.0;

	(*tra)->v1_P = 43.0 / 3.6;		//[m/s]
	(*tra)->v2_P = 87.72432417 / 3.6;		//[m/s]
	//(*tra)->v0_R = 52.0 / 3.6;        //[m/s]
	(*tra)->v1_R = 80.85520429 / 3.6;		//[m/s]
	(*tra)->v2_R = 80.85520429 / 3.6;		//[m/s]



	(*tra)->Fmax = 11.07585959 * 1000.0;         //[N/MM]  自分用(100%乗車時)(参考：起動加速度0.639[m/s/s])
	(*tra)->Bmax = 12.78144531 * 1000.0;		//[N/MM] (100%乗車時)
	(*tra)->Bref = 0.625;
	(*tra)->BEmax = 0.0;

	if ((*tra)->Type == 0) {
		(*tra)->mass_M = ini->mass_M * ini->num_M * 1.1 + 158.8 * 55.0 * 1.0 * ini->num_M;											//[kg](250%乗車)
		(*tra)->mass_T = ini->mass_T * ini->num_T * 1.05 + 158.8 * 55.0 * 1.0 * ini->num_T;											//[kg](250%乗車)
		(*tra)->mass = ini->mass_M * ini->num_M * 1.1 + 158.8 * 55.0 * 1.0 * ini->num_M + ini->mass_T * ini->num_T * 1.05 + 158.8 * 55.0 * 1.0 * ini->num_T;				//[kg]
		(*tra)->mass_c = (ini->mass_M * ini->num_M * 1.1 + 158.8 * 55.0 * 1.0 * ini->num_M + ini->mass_T * ini->num_T * 1.05 + 158.8 * 55.0 * 1.0 * ini->num_T) / 1000.0;										//[t](250%乗車)
		(*tra)->num = ini->num_M + ini->num_T;
		(*tra)->num_mot = ini->num_M * 4.0;

		(*tra)->v1_P = 38.0 / 3.6;		//[m/s]
		(*tra)->v2_P = 67.09631518 / 3.6;		//[m/s]
		//(*tra)->v0_R = 52.0 / 3.6;        //[m/s]
		(*tra)->v1_R = 86.75693301 / 3.6;		//[m/s]
		(*tra)->v2_R = 86.75693301 / 3.6;		//[m/s]

		(*tra)->Fmax = 13.67065539 * 1000.0;         //[N/MM]  自分用(100%乗車時)(参考：起動加速度0.639[m/s/s])
		(*tra)->Bmax = 9.261753472 * 1000.0;		//[N/MM] (100%乗車時)
		(*tra)->Bref = 0.583333333;
		(*tra)->BEmax = 0.0;

	}

	//	(*tra)->i1d_P = 89.61499;  //力行時，定トルク領域でのd軸電流
	//	(*tra)->i1q_P = 265.02112; //力行時，定トルク領域でのq軸電流
	//	(*tra)->i1d_R = 90.55028;  //回生時，定トルク領域でのd軸電流
	//	(*tra)->i1q_R = -216.2263; //回生時，定トルク領域でのq軸電流



	(*tra)->B = 0.533;

	(*tra)->direction = ini->direction;
	//	(*tra)->UP_DOWN = (double)(*tra)->direction;
	(*tra)->flag_station = 0;
	(*tra)->startflag = 0;
	(*tra)->lap = 0;
	(*tra)->wait_flag = ini->wait_flag;

	//	(*tra)->t_count = 0;
	//	(*tra)->c_count = 0;

	//	(*tra)->wr = 0.0;
	//	(*tra)->ws = 0.0;
	//	(*tra)->we = 0.0;

	//	(*tra)->con_speed = 0.0;

	//	(*tra)->Gr = 6.31;  //ギア比
	(*tra)->Tacr = 0.010;  //電流制御時定数
	//	(*tra)->POLE = 2.0;  //極対数
	//	(*tra)->rd = 0.41;   //車輪半径[m]
	//	(*tra)->R1 = 0.0970;
	//	(*tra)->R2 = 0.07327;  //2次抵抗[Ω]
	//	(*tra)->L1 = 0.030549;
	//	(*tra)->L2 = 0.030549;  //2次自己インダクタンス[H]
	//	(*tra)->M_IM = 0.029514;
	//	(*tra)->g0 = 0.001614;

		//(*tra)->Cfc = C_FC * ini->num_M;
		//(*tra)->Lfl = L_FL / ini->num_M;
	(*tra)->Rfl = R_FL / ini->num_M;

	//	(*tra)->theta = 1.0;
	//	(*tra)->theta_old = 1.0;


	(*tra)->Emot_pow = 0.0;			//モータ力行エネルギー[J]
	(*tra)->Emot_reg = 0.0;			//モータ回生エネルギー[J]
	(*tra)->Emot_air = 0.0;			//機械ブレーキエネルギー[J]

	(*tra)->Emot_pow_c = 0.0;		//モータ力行エネルギー[kWh]
	(*tra)->Emot_reg_c = 0.0;		//モータ回生エネルギー[kWh]
	(*tra)->Emot_air_c = 0.0;		//機械ブレーキエネルギー[kWh]

	(*tra)->Eloss_inv_pow = 0.0;
	(*tra)->Eloss_mot_pow = 0.0;
	(*tra)->Eloss_fl_pow = 0.0;
	(*tra)->Eloss_all_pow = 0.0;

	(*tra)->Eloss_inv_reg = 0.0;
	(*tra)->Eloss_mot_reg = 0.0;
	(*tra)->Eloss_fl_reg = 0.0;
	(*tra)->Eloss_all_reg = 0.0;

	(*tra)->Eres_pow = 0.0;
	(*tra)->Eres_reg = 0.0;
	(*tra)->Eres_coa = 0.0;

	(*tra)->Esiv_pow = 0.0;
	(*tra)->Esiv_reg = 0.0;
	(*tra)->Esiv_coa = 0.0;
	(*tra)->Esiv_stp = 0.0;

	(*tra)->Eloss_inv_pow_c = 0.0;
	(*tra)->Eloss_mot_pow_c = 0.0;
	(*tra)->Eloss_fl_pow_c = 0.0;
	(*tra)->Eloss_all_pow_c = 0.0;

	(*tra)->Eloss_inv_reg_c = 0.0;
	(*tra)->Eloss_mot_reg_c = 0.0;
	(*tra)->Eloss_fl_reg_c = 0.0;
	(*tra)->Eloss_all_reg_c = 0.0;

	(*tra)->Eres_pow_c = 0.0;
	(*tra)->Eres_reg_c = 0.0;
	(*tra)->Eres_coa_c = 0.0;

	(*tra)->Esiv_pow_c = 0.0;
	(*tra)->Esiv_reg_c = 0.0;
	(*tra)->Esiv_coa_c = 0.0;
	(*tra)->Esiv_stp_c = 0.0;

	(*tra)->Ereg_nashi = 0.0;
	(*tra)->Ereg_nashi_c = 0.0;



	(*tra)->Evh = 0.0;				//列車入出力エネルギー[J]
	(*tra)->Evh_in = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_out = 0.0;			//列車出力エネルギー[J]
	(*tra)->Evh_cos = 0.0;
	(*tra)->Evh_st = 0.0;


	(*tra)->Evh_in1 = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_in2 = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_in3 = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_in4 = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_in5 = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_in6 = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_in7 = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_in8 = 0.0;			//列車入力エネルギー[J]

	(*tra)->Evh_c = 0.0;			//列車入出力エネルギー[kWh]
	(*tra)->Evh_in_c = 0.0;			//列車入力エネルギー[kWh]
	(*tra)->Evh_out_c = 0.0;		//列車出力エネルギー[kWh]
	(*tra)->Evh_cos_c = 0.0;
	(*tra)->Evh_st_c = 0.0;

	(*tra)->Evh_in1_c = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_in2_c = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_in3_c = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_in4_c = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_in5_c = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_in6_c = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_in7_c = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_in8_c = 0.0;			//列車入力エネルギー[J]

	(*tra)->Evh_out1 = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_out2 = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_out3 = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_out4 = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_out5 = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_out6 = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_out7 = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_out8 = 0.0;			//列車入力エネルギー[J]

	(*tra)->Evh_out1_c = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_out2_c = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_out3_c = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_out4_c = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_out5_c = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_out6_c = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_out7_c = 0.0;			//列車入力エネルギー[J]
	(*tra)->Evh_out8_c = 0.0;			//列車入力エネルギー[J]

	(*tra)->accelflag = 4;
	(*tra)->brakeflag = 0;
	(*tra)->laststopflag = 0;
	(*tra)->X_recentstop = 0;
	(*tra)->X_nextstop = 0.0;
	(*tra)->X_brake = 0.0;
	(*tra)->x = ini->x;
	(*tra)->v = ini->v;
	(*tra)->v_new = 0.0;
	(*tra)->a = 0.0;
	(*tra)->Fmot = 0.0;
	(*tra)->Ftot = 0.0;

	//	(*tra)->id_mot_P = 0.0;
	//	(*tra)->id_mot_R = 0.0;
	//	(*tra)->iq_mot_P = 0.0;
	//	(*tra)->iq_mot_R = 0.0;

	//	(*tra)->v1d = 0.0;
	//	(*tra)->v1q = 0.0;
	//	(*tra)->Em = 0.0;
	//	(*tra)->idtot = 0.0;
	//	(*tra)->iqtot = 0.0;
	//	(*tra)->idtotF = 0.0;
	//	(*tra)->iqtotF = 0.0;
	//	(*tra)->iq_reg = 0.0;
	//	(*tra)->iq_regF = 0.0;
	(*tra)->b_regF = 0.0;
	//(*tra)->b_regF_old = 0.0;
	//(*tra)->b_regF_old_L = 0.0;
	(*tra)->FmotF = 0.0;
	(*tra)->Preg_nashi = 0.0;
	(*tra)->Preg_nashi_c = 0.0;

	(*tra)->x_c = ini->x / 1000.0;
	(*tra)->v_c = ini->v / 1000.0;

	(*tra)->Pmot = 0.0;
	(*tra)->Ploss_mot = 0.0;
	(*tra)->Ploss_fe = 0.0;
	(*tra)->Ploss_inv = 0.0;
	(*tra)->Ploss_fl = 0.0;
	(*tra)->iinv = 0.0;
	(*tra)->ifc = 0.0;
	(*tra)->ifl = 0.0;
	(*tra)->ifl_sum = 0.0;
	(*tra)->vfc_old = Vss1;
	(*tra)->isiv = 0.0;
	(*tra)->vfc = Vss1;
	(*tra)->vp = Vss1;
	(*tra)->vp_old = Vss1;
	(*tra)->R_grade = 0.0;
	(*tra)->R_curve = 0.0;
	(*tra)->R_run = 0.0;
	(*tra)->R_total = 0.0;
	(*tra)->t_stop = 0.0;
	(*tra)->t_stop_old = 0.0;
	//	(*tra)->t_wait = 0.0;
	(*tra)->diaflag = 0;

	(*tra)->Speedlimit = 0.0;
	(*tra)->Reaccelspeed = 0.0;
	//	(*tra)->Reaccelspeed = 0.0;
	(*tra)->T_delay = ini->T_delay;
	(*tra)->Vst = Vss1;

	(*tra)->run_time = 0.0;
	//(*tra)->ave_runtime = 0.0;

	//for (m = 0; m < rownum; m++)
	//{
	//	for (n = 0; n < 2; n++)
	//	{
	//		(*tra)->route[m][n] = 0;
	//	}
	//}
}

void Make_substation(SUB** sub, INI_SUB* ini, int i)
{

	*sub = (SUB*)malloc(sizeof(SUB));
	if (*sub == NULL)
	{
		//エラー処理
	}

	(*sub)->name_SS = i;
	(*sub)->Vss_0 = Vss;
	(*sub)->Iss_0 = Iss_rated;
	(*sub)->e_ss = ini->e_ss;
	(*sub)->Xss = ini->Xss;

	(*sub)->diode = ini->diode;
	(*sub)->vss = ini->vss;
	(*sub)->vss_e = ini->vss;
	(*sub)->iss = ini->iss;
	(*sub)->vout = ini->vss;
	(*sub)->flag = 1;

	(*sub)->Rss = (ini->e_ss * ini->vss) / (100.0 * ini->Iss);
	(*sub)->Rss_ini = (ini->e_ss * ini->vss) / (100.0 * ini->Iss);
	(*sub)->Jss_ini = -ini->vss / ((ini->e_ss * ini->vss) / (100.0 * ini->Iss));	//Jssの初期値
	(*sub)->Jss = -ini->vss / ((ini->e_ss * ini->vss) / (100.0 * ini->Iss));		//更新用のJss

	(*sub)->Ess = 0.0;			//変電所出力エネルギー[J]
	(*sub)->Ess_c = 0.0;		//変電所出力エネルギー[kWh]
	(*sub)->Wss = 0.0;
	(*sub)->Wss_c = 0.0;
}

void Make_NODE_TRAIN(NODE** node, TRAIN** tra)
{
	*node = (NODE*)malloc(sizeof(NODE));
	if (*node == NULL)
	{
		printf("MEMORY?? What is that????\n");
		exit(EXIT_FAILURE);//エラー処理
	}

	(*node)->Number = (*tra)->name_train;
	(*node)->X = (*tra)->x;
	(*node)->V = (*tra)->vp;
	(*node)->I = (*tra)->ifl;
	(*node)->r = 0.0;
	(*node)->flag = (*tra)->direction;
}

void Make_NODE_SS(NODE** node, SUB** sub)
{
	*node = (NODE*)malloc(sizeof(NODE));
	if (*node == NULL)
	{
		printf("MEMORY?? What is that????\n");
		exit(EXIT_FAILURE);//エラー処理
	}

	(*node)->Number = (*sub)->name_SS;
	(*node)->X = (*sub)->Xss;
	(*node)->V = (*sub)->vss;
	(*node)->I = (*sub)->Jss;
	(*node)->r = (*sub)->Rss;
	(*node)->flag = 0.0;
}

void Error_Detection_SS(SUB* sub)
{
	if (sub->vout > sub->vss_e + 0.1) {
		sub->Jss = sub->iss - sub->vout / sub->Rss;
	}
	else if (sub->vout < sub->vss_e - 0.1) sub->Jss = sub->Jss_ini;

	if (sub->iss < -0.1)		//変電所に電流が流入している⇒flag=0（エラー）
	{
		sub->flag = 0;
		sub->Jss = -sub->vout / sub->Rss;
	}
	else
	{
		sub->flag = 1;			//変電所に電流が流入していない⇒flag=1（OK）
	}

	/*	if (sub->vout >(Vcmax))	//電圧の制限
		{
			sub->flag = 0;
			sub->Jss = -(Vcmax) / sub->Rss;
		}*/
}

//void Error_Detection_REG(TRAIN* tra) {

	/********** 軽負荷回生制御パターンを満たすかどうかの判定 **********/
/*	if (tra->vfc > Vclim && tra->vfc < Vcmax)						//回生絞り込み中
	{
		tra->b_reg_judge = tra->Fmot * (Vcmax - tra->vfc) / (Vcmax - Vclim);
	}
	else if (tra->vfc >= Vcmax)										//回生絞り込み終了
	{
		tra->b_reg_judge = 0.0;
	}
	else															//回生絞り込みなし
	{
		tra->b_reg_judge = tra->Fmot;
	}

	tra->judgement = fabs(tra->b_reg - tra->b_reg_judge);

	if (tra->judgement > 1.0) {
		tra->reg_flag = 0;
	}
	else {
		tra->reg_flag = 1;
	}
}*/

void Calculate_R_run(TRAIN* tra)		//走行抵抗の算出
{
	tra->R_run = (tra->mass_M / 1000.0 * (1.65 + 0.0247 * tra->v_c) + tra->mass_T / 1000.0 * (0.78 + 0.0028 * tra->v_c) + (0.028 + 0.0078 * (tra->num - 1)) * tra->v_c * tra->v_c) * g;		//rr=g×W×(a+bv+cv^2/W)
	//	tra->R_run = g*(1.32+0.0614 * tra->v_c)*tra->mass_c+(0.028 + 0.0078 * (tra->num - 1)) * tra->v_c * tra->v_c;
}

void Calculate_R_total(TRAIN* tra)		//合計抵抗の算出
{
	Calculate_R_run(tra);

	tra->R_grade = 0.0;
	tra->R_curve = 0.0;

	tra->R_total = (tra->R_grade + tra->R_curve) * (tra->mass_c) + tra->R_run;
}

/*** 一次遅れフィルタ ***/
double funcdelay(double aTs, double aTf, double ayold, double au)
{
	double fudely;
	fudely = (aTs * au + aTf * ayold) / (aTf + aTs);
	return(fudely);
}
/************************/

void Traction_force(TRAIN* tra)		//力行時引張力の算出
{
	tra->v1_P = 43.0 / 3.6 * tra->vfc / 1500.0;
	tra->v2_P = 87.72432417 / 3.6 * tra->vfc / 1500.0;

	if (tra->Type == 0) {
		tra->v1_P = 38.0 / 3.6 * tra->vfc / 1500.0;
		tra->v2_P = 67.09631518 / 3.6 * tra->vfc / 1500.0;
	}
	if (tra->Type == 2) {
		tra->v1_P = 38.0 / 3.6 * tra->vfc / 1500.0;
		tra->v2_P = 123.3043194 / 3.6 * tra->vfc / 1500.0;
	}

	if (tra->v < tra->v1_P)		//定トルク領域
	{
		tra->Fmot = tra->Fmax * tra->num_mot;
	}
	else if (tra->v < tra->v2_P)		//定電力領域
	{
		tra->Fmot = tra->Fmax * tra->num_mot * tra->v1_P / tra->v;
	}
	else
	{
		tra->Fmot = tra->Fmax * tra->num_mot * tra->v1_P * tra->v2_P / (tra->v * tra->v);	//特性領域
	}
	tra->Ftot = tra->Fmot;
	tra->FmotF = funcdelay(dt, tra->Tacr, tra->FmotF, tra->Fmot);
}

void Regenerative_force(TRAIN* tra)		//回生時引張力の算出
{
	tra->v1_R = 80.85520429 / 3.6 * tra->vfc / 1650.0;
	tra->v2_R = 80.85520429 / 3.6 * tra->vfc / 1650.0;

	if (tra->Type == 0) {
		tra->v1_R = 86.75693301 / 3.6 * tra->vfc / 1650.0;
		tra->v2_R = 86.75693301 / 3.6 * tra->vfc / 1650.0;
	}
	if (tra->Type == 2) {
		tra->v1_R = 90.5139581 / 3.6 * tra->vfc / 1650.0;
		tra->v2_R = 90.5139581 / 3.6 * tra->vfc / 1650.0;
	}

	if (tra->v < tra->v1_R)		//定トルク領域
	{
		tra->Fmot_regM = -tra->Bmax * tra->num_mot;
	}
	/*else if (tra->v < tra->v2_R)		//定電力領域
	{
		tra->Fmot_regM = -tra->Bmax * tra->num_mot * tra->v1_R / tra->v;
	}*/
	else
	{
		tra->Fmot_regM = -tra->Bmax * tra->num_mot * tra->v1_R * tra->v2_R / (tra->v * tra->v);		//特性領域
	}

	tra->Fmot = tra->Fmot_regM;
	//if (tra->Fmot < -tra->BEmax) tra->Fmot = -tra->BEmax;
	tra->Ftot = -tra->BEmax;
}

void Solve_next_state(TRAIN* tra)
{
	tra->a = (tra->Ftot - tra->R_total) / (tra->mass);		//運動方程式ma=F
	tra->a_c = tra->a * 3600.0 / 1000.0;
	/*if (tra->accelflag == 3)
	{
		tra->a = -0.533;		//運動方程式ma=F
		tra->a_c = tra->a * 3600.0 / 1000.0;
	}
	else
	{
		tra->a = (tra->Ftot - tra->R_total) / (tra->mass);		//運動方程式ma=F
		tra->a_c = tra->a * 3600.0 / 1000.0;
	}*/
	tra->v_new = tra->v + tra->a * dt;

	tra->v = tra->v_new;
	tra->v_c = tra->v * 3.6;


	if (tra->direction == 1)
	{
		tra->x += tra->v * dt;
	}
	else
	{
		tra->x -= tra->v * dt;
	}

	tra->x_c = tra->x / 1000.0;
}

void Calculate_BEmax(TRAIN* tra)
{
	tra->BEmax = 409006.2;
	if (tra->Type == 0) { tra->BEmax = 222282.0833; }
}

void Run_pattern(TRAIN* tra, STATION* sta, STATION* sta2, STATION* sta3, double t)
{
	int i;

	tra->t_stop = t;

	/*** 次の駅を探索 ***/
	tra->X_nextstop_old = tra->X_nextstop;
	if (tra->Type == 0)
	{
		for (i = 0; i < NUM_station + NUM_final_station; i++)
		{
			if (tra->direction == 1.0)
			{
				if (sta[i].Xs <= tra->x && sta[i + 1].Xs + 40.0 >= tra->x)
				{
					tra->X_nextstop = sta[i + 1].Xs;
					break;
				}
			}
			else
			{
				if (sta[i].Xs <= tra->x && sta[i + 1].Xs - 40.0 >= tra->x)
				{
					tra->X_nextstop = sta[i].Xs;
					break;
				}
			}

		}

		/*終着駅到着判定*/
		if (tra->direction == 1.0)			//新横浜⇒橋本方面
		{
			if (tra->flag_direction == 1)		//野洲止まりの場合
			{
				if (i == 7)
				{
					tra->laststopflag = 1;
				}
			}
			else
			{
				if (i == NUM_station - 1)
				{
					tra->laststopflag = 1;
				}
			}

		}
		else
		{
			if (i == 1)						//橋本⇒新横浜方面
			{
				tra->laststopflag = 1;
			}
		}
	}
	else if (tra->Type == 1)                                           //Type快速
	{
		for (i = 0; i < NUM_station2 + NUM_final_station2; i++)
		{
			/*if (tra->direction == 2.0)
			{
				if (sta2[i].Xs <= tra->x && sta2[i + 1].Xs + 40.0 >= tra->x)
				{
					tra->X_nextstop = sta[i + 1].Xs;
					break;
				}
			}
			else
			{
				if (sta[i].Xs <= tra->x && sta[i + 1].Xs - 40.0 >= tra->x)
				{
					tra->X_nextstop = sta[i].Xs;
					break;
				}
			}*/
			if (sta2[i].Xs <= tra->x && sta2[i + 1].Xs >= tra->x)
			{
				if (tra->direction == 1.0)			//新横浜⇒橋本方面
				{
					tra->X_nextstop = sta2[i + 1].Xs;
					break;
				}
				else								//橋本⇒新横浜方面
				{
					tra->X_nextstop = sta2[i].Xs;
					break;
				}
			}
		}

		/*終着駅到着判定*/
		if (tra->direction == 1.0)			//新横浜⇒橋本方面
		{
			if (tra->flag_direction == 1)		//町田止まりの場合
			{
				if (i == 4)
				{
					tra->laststopflag = 1;
				}
			}
			else
			{
				if (i == NUM_station2 - 1)
				{
					tra->laststopflag = 1;
				}
			}

		}
		else
		{
			if (i == 1)						//橋本⇒新横浜方面
			{
				tra->laststopflag = 1;
			}
		}

	}
	else                                         //Type特急
	{
		for (i = 0; i < NUM_station3 + NUM_final_station3; i++)
		{
			if (sta3[i].Xs <= tra->x && sta3[i + 1].Xs >= tra->x)
			{
				if (tra->direction == 1.0)			//新横浜⇒橋本方面
				{
					tra->X_nextstop = sta3[i + 1].Xs;
					break;
				}
				else								//橋本⇒新横浜方面
				{
					tra->X_nextstop = sta3[i].Xs;
					break;
				}
			}
		}

		/*終着駅到着判定*/
		if (tra->direction == 1.0)			//新横浜⇒橋本方面
		{
			if (tra->flag_direction == 1)		//町田止まりの場合
			{
				if (i == 2)
				{
					tra->laststopflag = 1;
				}
			}
			else
			{
				if (i == NUM_station3 - 1)
				{
					tra->laststopflag = 1;
				}
			}

		}
		else
		{
			if (i == 0)						//橋本⇒新横浜方面
			{
				tra->laststopflag = 1;
			}
		}

	}

	/*速度制限設定*/
	if (tra->Type == 0)//|| tra->Type == 1)
	{
		if (tra->direction == 1.0)
		{
			/*if (sta[1].Xs < tra->x && tra->x < sta[2].Xs)
			{
				tra->Speedlimit = 88.0 - Speed_change;
				tra->Reaccelspeed = 0.0;
			}
			else if (sta[2].Xs < tra->x && tra->x < sta[3].Xs)
			{
				tra->Speedlimit = 93.0 - Speed_change;
				tra->Reaccelspeed = 0.0;
			}
			else if (sta[4].Xs < tra->x && tra->x < sta[5].Xs)
			{
				tra->Speedlimit = 93.0 - Speed_change;
				tra->Reaccelspeed = 0.0;
			}
			else if (sta[5].Xs < tra->x && tra->x < sta[6].Xs)
			{
				tra->Speedlimit = 85.0 - Speed_change;
				tra->Reaccelspeed = 0.0;
			}
			else if (sta[6].Xs < tra->x && tra->x < sta[7].Xs)
			{
				tra->Speedlimit = 90.0 - Speed_change;
				tra->Reaccelspeed = 0.0;
			}
			else if (sta[7].Xs < tra->x && tra->x < sta[8].Xs)
			{
				tra->Speedlimit = 85.0 - Speed_change;
				tra->Reaccelspeed = 0.0;
			}
			else if (sta[10].Xs < tra->x && tra->x < sta[11].Xs)
			{
				tra->Speedlimit = 72.0 - Speed_change;
				tra->Reaccelspeed = 0.0;
			}
			else if (sta[11].Xs < tra->x && tra->x < sta[12].Xs)
			{
				tra->Speedlimit = 88.0 - Speed_change;
				tra->Reaccelspeed = 0.0;
			}
			else
			{*/
			tra->Speedlimit = 80.0 + vel_dif[vel_state];
			tra->Reaccelspeed = 0.0;
			//}
		}
		else
		{
			/*if (sta[1].Xs < tra->x && tra->x < sta[2].Xs)
			{
				tra->Speedlimit = 90.0 - Speed_change;
				tra->Reaccelspeed = 0.0;
			}
			else if (sta[10].Xs < tra->x && tra->x < sta[11].Xs)
			{
				tra->Speedlimit = 75.0 - Speed_change;
				tra->Reaccelspeed = 0.0;
			}
			else if (sta[11].Xs < tra->x && tra->x < sta[12].Xs)
			{
				tra->Speedlimit = 93.0 - Speed_change;
				tra->Reaccelspeed = 0.0;
			}
			else
			{*/
			tra->Speedlimit = 80.0 + vel_dif[vel_state];
			tra->Reaccelspeed = 0.0;
			//}
		}


	}
	//快速用
	if (tra->Type == 1)
	{
		if (tra->direction == 1.0)
		{/*
			if (sta2[1].Xs < tra->x && tra->x < sta2[2].Xs)
			{
				tra->Speedlimit = 93.0 - Speed_change;
			}
			else if (sta2[3].Xs < tra->x && tra->x < sta2[4].Xs)
			{
				tra->Speedlimit = 93.0 - Speed_change;
			}
			else if (sta2[4].Xs < tra->x && tra->x < sta2[5].Xs)
			{
				tra->Speedlimit = 90.0 - Speed_change;
			}
			else
			{*/
			tra->Speedlimit = 105.0 + vel_dif[vel_state];
			//	}
		}
		else
		{
			tra->Speedlimit = 105.0 + vel_dif[vel_state];
		}

		/*
		if (tra->direction == 1.0)
		{
			if (sta[1].Xs < tra->x && tra->x < sta[2].Xs)
			{
				tra->Speedlimit = 88.0 - Speed_change;
				tra->Reaccelspeed = 0.0;
			}
			else if (sta[2].Xs < tra->x && tra->x < sta[3].Xs)
			{
				tra->Speedlimit = 93.0 - Speed_change;
				tra->Reaccelspeed = 0.0;
			}
			else if (sta[4].Xs < tra->x && tra->x < sta[5].Xs)
			{
				tra->Speedlimit = 93.0 - Speed_change;
				tra->Reaccelspeed = 0.0;
			}
			else if (sta[5].Xs < tra->x && tra->x < sta[6].Xs)
			{
				tra->Speedlimit = 85.0 - Speed_change;
				tra->Reaccelspeed = 0.0;
			}
			else if (sta[6].Xs < tra->x && tra->x < sta[7].Xs)
			{
				tra->Speedlimit = 90.0 - Speed_change;
				tra->Reaccelspeed = 0.0;
			}
			else if (sta[7].Xs < tra->x && tra->x < sta[8].Xs)
			{
				tra->Speedlimit = 85.0 - Speed_change;
				tra->Reaccelspeed = 0.0;
			}
			else if (sta[10].Xs < tra->x && tra->x < sta[11].Xs)
			{
				tra->Speedlimit = 72.0 - Speed_change;
				tra->Reaccelspeed = 0.0;
			}
			else if (sta[11].Xs < tra->x && tra->x < sta[12].Xs)
			{
				tra->Speedlimit = 88.0 - Speed_change;
				tra->Reaccelspeed = 0.0;
			}
			else
			{
				tra->Speedlimit = 95.0 - Speed_change;
				tra->Reaccelspeed = 0.0;
			}
		}
		else
		{
			if (sta[1].Xs < tra->x && tra->x < sta[2].Xs)
			{
				tra->Speedlimit = 95.0 - Speed_change;
				tra->Reaccelspeed = 0.0;
			}
			else if (sta[10].Xs < tra->x && tra->x < sta[11].Xs)
			{
				tra->Speedlimit = 95.0 - Speed_change;
				tra->Reaccelspeed = 0.0;
			}
			else if (sta[11].Xs < tra->x && tra->x < sta[12].Xs)
			{
				tra->Speedlimit = 95.0 - Speed_change;
				tra->Reaccelspeed = 0.0;
			}
			else
			{
				tra->Speedlimit = 95.0 - Speed_change;
				tra->Reaccelspeed = 0.0;
			}
		}*/
	}



	/**** 初期化をするかどうか判定&初期化 ****/
	if (tra->accelflag == 4)
	{
		(tra->flag_station)++;

		/*停車位置修正*/
//		if(tra->direction == 1.0)			//大宮⇒吹上方面
//		{
//			tra->x = sta[i].Xs;
//		}
//		else								//吹上⇒大宮方面
//		{
//			tra->x = sta[i+1].Xs;
//		}
	}

	/**** 走行抵抗計算 ****/
	Calculate_R_total(tra);



	/*** 走行モードを決定 ***/
	if (tra->accelflag == 4)		/*停車中の動作を決定*/
	{
		if (tra->laststopflag == 1)
		{
			tra->accelflag = 4;
		}
		else if (t < tra->T_delay)
		{
			tra->accelflag = 4;
		}
		else if (t == 0)
		{
			tra->accelflag = 1;
		}

		else if (tra->wait_flag == 1)
		{
			if (tra->direction == -1.0)
			{
				if (tra->x < sta2[4].Xs && tra->x > sta2[3].Xs)
				{
					if ((tra->t_stop - tra->t_stop_old) >= 180.0)
					{
						tra->accelflag = 1;
					}
				}
				else
				{
					if ((tra->t_stop - tra->t_stop_old) >= 30.0)
					{
						tra->accelflag = 1;
					}
				}
			}
			else if (tra->direction == 1.0)
			{
				if (tra->x < sta2[3].Xs && tra->x > sta2[2].Xs)
				{
					if ((tra->t_stop - tra->t_stop_old) >= 180.0)
					{
						tra->accelflag = 1;
					}
				}
				else
				{
					if ((tra->t_stop - tra->t_stop_old) >= 30.0)
					{
						tra->accelflag = 1;
					}
				}
			}
		}
		else if ((tra->t_stop - tra->t_stop_old) >= 30.0)
		{

			tra->accelflag = 1;

			//			tra->x_dummy = tra->X_nextstop;
		}
		else
		{
			tra->accelflag = 4;
		}

	}
	else if (tra->accelflag == 2)		/*惰行中の動作を決定*/
	{
		if (tra->direction == 1.0)		//大宮⇒吹上方面
		{
			if (tra->Type == 1)
			{
				tra->Reaccelspeed = 95.0 + vel_dif[vel_state];/*
				if (tra->x > sta[2].Xs && tra->x < sta[3].Xs) {
					tra->Speedlimit = 93.0 - Speed_change;
					tra->Reaccelspeed = 90.0;
				}
				if (tra->x > sta[5].Xs && tra->x < sta[6].Xs) {
					tra->Speedlimit = 85.0 - Speed_change;
					tra->Reaccelspeed = 87.0;
				}
				if (tra->x > sta[7].Xs && tra->x < sta[8].Xs) {
					tra->Speedlimit = 85.0 - Speed_change;
					tra->Reaccelspeed = 85.0;
				}
				if (tra->x > sta[9].Xs && tra->x < sta[10].Xs) {
					tra->Speedlimit = 95.0 - Speed_change;
					tra->Reaccelspeed = 85.0;
				}
				if (tra->x > sta[11].Xs && tra->x < sta[12].Xs) {
					tra->Speedlimit = 88.0 - Speed_change;
					tra->Reaccelspeed = 75.0;
				}*/

				//10kmマイナス
				/*if (tra->x > sta[1].Xs && tra->x < sta[3].Xs) {
					tra->Reaccelspeed = 85.0 - Speed_change;
				}
				if (tra->x > sta[4].Xs && tra->x < sta[6].Xs) {
					tra->Reaccelspeed = 85.0 - Speed_change;
				}
				if (tra->x > sta[6].Xs && tra->x < sta[8].Xs) {
					tra->Reaccelspeed = 85.0 - Speed_change;
				}
				if (tra->x > sta[8].Xs && tra->x < sta[12].Xs) {
					tra->Reaccelspeed = 85.0 - Speed_change;
				}*/
			}
			if (tra->Type == 2)
			{
				tra->Reaccelspeed = 110.0 + vel_dif[vel_state];
				if (tra->x > tra->X_nextstop) {
					tra->accelflag = 4;
					tra->v_c = 0.0;
				}
			}


			if (tra->x + tra->v * tra->v / 2.0 / tra->Bref < tra->X_nextstop)
			{
				if (tra->v_c > tra->Reaccelspeed || tra->v_c > tra->Speedlimit)
				{
					tra->accelflag = 2;
				}
				else
				{
					tra->accelflag = 1;
				}
			}
			else
			{
				if (tra->Type == 0 || tra->Type == 1)
				{
					tra->accelflag = 3;
					tra->brakeflag = 1;
				}
			}


		}
		else							//大宮⇒吹上方面
		{
			if (tra->Type == 1)
			{
				tra->Reaccelspeed = 95.0 + vel_dif[vel_state];/*
				if (tra->x > sta[9].Xs && tra->x < sta[10].Xs) {
					tra->Speedlimit = 95.0 - Speed_change;
					tra->Reaccelspeed = 75.0;
				}
				if (tra->x > sta[8].Xs && tra->x < sta[9].Xs) {
					tra->Speedlimit = 95.0 - Speed_change;
					tra->Reaccelspeed = 75.0;
				}
				if (tra->x > sta[6].Xs && tra->x < sta[7].Xs) {
					tra->Speedlimit = 95.0 - Speed_change;
					tra->Reaccelspeed = 87.0;
				}
				if (tra->x > sta[4].Xs && tra->x < sta[5].Xs) {
					tra->Speedlimit = 95.0 - Speed_change;
					tra->Reaccelspeed = 90.0;
				}
				if (tra->x > sta[1].Xs && tra->x < sta[2].Xs) {
					tra->Speedlimit = 90.0 - Speed_change;
					tra->Reaccelspeed = 85.0;
				}*/
				//10kmマイナス
				/*if (tra->x > sta[1].Xs && tra->x < sta[3].Xs) {
					tra->Reaccelspeed = 85.0 - Speed_change;
				}
				if (tra->x > sta[4].Xs && tra->x < sta[6].Xs) {
					tra->Reaccelspeed = 85.0 - Speed_change;
				}
				if (tra->x > sta[6].Xs && tra->x < sta[8].Xs) {
					tra->Reaccelspeed = 85.0 - Speed_change;
				}
				if (tra->x > sta[8].Xs && tra->x < sta[12].Xs) {
					tra->Reaccelspeed = 85.0 - Speed_change;
				}*/
			}
			if (tra->Type == 2)
			{
				tra->Reaccelspeed = 110.0 + vel_dif[vel_state];
				if (tra->x < tra->X_nextstop) {
					tra->accelflag = 4;
					tra->v_c = 0.0;
				}

			}
			if (tra->x - tra->v * tra->v / 2.0 / tra->Bref > tra->X_nextstop)
			{
				if (tra->v_c > tra->Reaccelspeed || tra->v_c > tra->Speedlimit)
				{
					tra->accelflag = 2;
				}
				else
				{
					tra->accelflag = 1;
				}
			}
			else
			{
				if (tra->Type == 0 || tra->Type == 1)
				{
					tra->accelflag = 3;
					tra->brakeflag = 1;
				}
			}

		}
	}
	else		/*加速中or減速中の動作を決定*/
	{
		if (tra->direction == 1.0)								//大宮⇒吹上方面
		{
			if (tra->brakeflag == 1)		/*減速中の動作*/
			{
				if (tra->v_c >= 0.0)
				{
					tra->accelflag = 3;		//減速区間
				}
				else
				{
					tra->accelflag = 4;		//停車
					tra->brakeflag = 0;
					tra->t_stop_old = tra->t_stop;
				}
			}
			else		/*加速中の動作*/
			{
				if (tra->x + tra->v * tra->v / 2.0 / tra->Bref > tra->X_nextstop && tra->v_c > 50)
				{
					tra->accelflag = 3;
					tra->brakeflag = 1;
				}
				else if (tra->v_c <= tra->Speedlimit) //&& (tra->a >= 0))
				{
					tra->accelflag = 1;		//加速区間
				}
				else
				{
					tra->accelflag = 2;		//惰行区間
				}
			}
		}
		else													//吹上⇒大宮方面
		{
			if (tra->brakeflag == 1)		/*減速中の動作*/
			{
				if (tra->v_c >= 0.0)
				{
					tra->accelflag = 3;		//減速区間
				}
				else
				{
					tra->accelflag = 4;		//停車
					tra->brakeflag = 0;
					tra->t_stop_old = tra->t_stop;
				}
			}
			else		/*加速中の動作*/
			{
				if (tra->x - tra->v * tra->v / 2.0 / tra->Bref < tra->X_nextstop && tra->v_c > 50)
				{
					tra->accelflag = 3;
					tra->brakeflag = 1;
				}
				else if ((tra->v_c <= tra->Speedlimit)) //&& (tra->a >= 0))
				{
					tra->accelflag = 1;		//加速区間
				}
				else
				{
					tra->accelflag = 2;		//惰行区間
				}
			}
		}
	}

}


void Calculation_traction_force_brake_force(TRAIN* tra) {

	/*** 各走行モードに応じた動作を指定 ***/
	if (tra->accelflag == 1)				//加速中
	{
		/** ブレーキ力計算(加速中なので0) **/
		tra->b_reg = 0.0;
		tra->b_regF = 0.0;
		//tra->b_regF_old = 0.0;
		//tra->b_regF_old_L = 0.0;
		tra->b_air = 0.0;
		//tra->b_reg_loss = 0.0;

		/** 車両引張力計算 **/
		Traction_force(tra);

		/** 車両加速度・速度・位置計算 **/
	//	Solve_next_state(tra);
	}
	else if (tra->accelflag == 2)			//惰行中
	{
		/** ブレーキ力計算(惰行中なので0) **/
		tra->b_reg = 0.0;
		tra->b_regF = 0.0;
		//tra->b_regF_old = 0.0;
		//tra->b_regF_old_L = 0.0;
		tra->b_air = 0.0;
		//tra->b_reg_loss = 0.0;

		/** 車両引張力計算 **/
		tra->Fmot = 0.0;
		tra->Ftot = 0.0;


	}
	else if (tra->accelflag == 3)			//減速中
	{
		/** 車両ブレーキ力計算 **/
		Regenerative_force(tra);


		/********** 軽負荷回生中のトルク絞り込み模擬部（軽負荷回生制御下の電気ブレーキ力を計算） **********/
		if (tra->vfc > Vclim && tra->vfc < Vcmax)						//回生絞り込み中
		{
			tra->b_reg = tra->Fmot_regM * (Vcmax - tra->vfc) / (Vcmax - Vclim);
			//	tra->b_regF_old_L = tra->b_regF_old * (Vcmax - tra->vfc) / (Vcmax - Vclim);

		}
		else if (tra->vfc >= Vcmax)										//回生絞り込み終了
		{
			tra->b_reg = 0.0;
			//	tra->b_regF_old_L = 0.0;

		}
		else															//回生絞り込みなし
		{
			tra->b_reg = tra->Fmot_regM;
			//	tra->b_regF_old_L = tra->b_regF_old;

		}

		tra->b_regF = funcdelay(dt, tra->Tacr, tra->b_regF, tra->b_reg);

		//tra->b_reg_loss = tra->Fmot_regM - tra->b_regF;		//回生絞り込み量[N]

		tra->b_air = -tra->BEmax - tra->b_regF;			//足りないブレーキ力は空気ブレーキで補う

		/** 車両加速度・速度・位置計算 **/
		//Solve_next_state(tra);


	}
	else if (tra->accelflag == 4)						//停車中
	{
		tra->Fmot = 0.0;
		tra->Ftot = 0.0;
		tra->v = 0.0;
		tra->v_c = 0.0;
		tra->a = 0.0;
		tra->b_reg = 0.0;
		tra->b_regF = 0.0;
		//tra->b_regF_old = 0.0;
		//tra->b_regF_old_L = 0.0;
		tra->b_air = 0.0;
		//tra->b_reg_loss = 0.0;


	}


	tra->Fmot_c = tra->Fmot / 1000.0;
	tra->b_reg_c = tra->b_reg / 1000.0;
	tra->b_reg_c = tra->b_regF / 1000.0;
	tra->Ftot_c = tra->Ftot / 1000.0;
	tra->b_air_c = tra->b_air / 1000.0;
	//tra->b_reg_loss_c = tra->b_reg_loss / 1000.0;
}

void Calculation_traction_circuit(TRAIN* tra)
{
	tra->isiv = P_SIV * tra->num / tra->vfc;
	if (tra->Type == 4) tra->isiv = P_teifuka / tra->vp;
	/*
	if (tra->laststopflag == 0)
	{
		if (t < tra->T_delay)
		{
			tra->isiv = 0.0;
		}
		else
		{
			tra->isiv = P_SIV * tra->num / tra->vfc;
		}
	}
	else if(tra->accelflag == 4)
	{
		tra->isiv = 0.0;
	}
	*/
	tra->ifl = tra->isiv + (tra->Pmot + tra->Ploss_mot + tra->Ploss_inv) / tra->vfc;
	//	tra->ifl = (tra->Pmot) / tra->vfc;



	tra->Ploss_mot = fabs((1 - e_motor) * tra->Pmot);
	tra->Ploss_fl = R_FL * tra->ifl * tra->ifl;
	/***インバータ各部損失計算***/
	tra->Pinv_on = 6.0 * (Vce * fabs(tra->ifl - tra->isiv) * 0.5);
	tra->Pinv_diode = 6.0 * (Vf * fabs(tra->ifl - tra->isiv) * 0.5);
	tra->Pinv_sw = 6.0 * 1000.0 * tra->vfc * fabs(tra->ifl - tra->isiv) / 6.0 * (t_on + t_off);		//スイッチング周波数1[kHz]

	if (fabs(tra->ifl - tra->isiv) > 5.0)
	{
		tra->Pinv_rec = 6.0 * 1000.0 * tra->vfc * Irp / 6.0 * t_rr;						//スイッチング周波数1[kHz]
	}
	else
	{
		tra->Pinv_rec = 0.0;
	}

	tra->Ploss_inv = (tra->Pinv_on + tra->Pinv_diode + tra->Pinv_sw + tra->Pinv_rec);			//インバータの損失[W]
	tra->Ploss_inv_c = tra->Ploss_inv / 1000.0;													//インバータの損失[kW]
	//tra->Ploss_inv = 0.0;

	tra->Ploss_mot_c = tra->Ploss_mot / 1000.0;				//モータの銅損[kW]

	tra->Ploss_all = tra->Ploss_mot + tra->Ploss_fl + tra->Ploss_inv;
	tra->Ploss_all_c = tra->Ploss_all / 1000.0;
}

void Calculation_power_motor(TRAIN* tra)
{
	tra->Ptot = tra->v * tra->Ftot;

	if (tra->accelflag == 1)				/*力行中*/
	{
		tra->Pmot = tra->v * tra->FmotF;
		tra->Pmot_air = 0.0;
		tra->Preg_nashi = 0.0;

	}
	else if (tra->accelflag == 3)		/*回生中*/
	{
		tra->Pmot = tra->v * tra->b_regF;
		tra->Pmot_air = tra->v * tra->b_air;
		tra->Preg_nashi = tra->v * tra->Fmot_regM;
	}
	else
	{
		tra->Pmot = 0.0;
		tra->Pmot_air = 0.0;
		//tra->Pmot_reg_loss = 0.0;
		tra->Preg_nashi = 0.0;
	}

	tra->Pres = tra->v * tra->R_total;		//走行抵抗でのロス

	tra->Ptot_c = tra->Ptot / 1000.0;
	tra->Pmot_c = tra->Pmot / 1000.0;
	tra->Pmot_air_c = tra->Pmot_air / 1000.0;
	tra->Pres_c = tra->Pres / 1000.0;
	//tra->Pmot_reg_loss_c = tra->Pmot_reg_loss / 1000.0;
	tra->Preg_nashi_c = tra->Preg_nashi / 1000.0;

}

void Calculation_power_train(TRAIN* tra)
{
	tra->Pvh = tra->vp * tra->ifl;

	if (tra->accelflag == 1)				//力行中,惰行中
	{
		tra->Pvh_in = tra->vp * tra->ifl;
		tra->Pvh_out = 0.0;
		tra->Pvh_st = 0.0;

		/*if (tra->x >= -1100.0 && tra->x <= 900.0) { tra->Pvh_in1 = tra->vp * tra->ifl; }
		else if (tra->x > 900.0 && tra->x <= 2850.0) { tra->Pvh_in2 = tra->vp * tra->ifl; }
		else if (tra->x > 2850.0 && tra->x <= 6050.0) { tra->Pvh_in3 = tra->vp * tra->ifl; }
		else if (tra->x > 6050.0 && tra->x <= 8960.0) { tra->Pvh_in4 = tra->vp * tra->ifl; }
		else if (tra->x > 8960.0 && tra->x <= 12608.0) { tra->Pvh_in5 = tra->vp * tra->ifl; }
		else if (tra->x > 12608.0 && tra->x <= 16918.0) { tra->Pvh_in6 = tra->vp * tra->ifl; }
		else if (tra->x > 16918.0 && tra->x <= 20008.0) { tra->Pvh_in7 = tra->vp * tra->ifl; }
		else if (tra->x > 20008.0 && tra->x <= 21700.0) { tra->Pvh_in8 = tra->vp * tra->ifl; }
		*/
	}
	else if (tra->accelflag == 3)		//回生中
	{
		tra->Pvh_in = 0.0;
		tra->Pvh_out = tra->vp * tra->ifl;
		tra->Pvh_st = 0.0;
		tra->Pvh_cos = 0.0;

		/*tra->Pvh_in1 = 0.0;
		tra->Pvh_in2 = 0.0;
		tra->Pvh_in3 = 0.0;
		tra->Pvh_in4 = 0.0;
		tra->Pvh_in5 = 0.0;
		tra->Pvh_in6 = 0.0;
		tra->Pvh_in7 = 0.0;
		tra->Pvh_in8 = 0.0;
		*/
	}
	else if (tra->accelflag == 2)
	{
		tra->Pvh_cos = tra->vp * tra->ifl;
		tra->Pvh_st = tra->vp * tra->ifl;
		tra->Pvh_in = 0.0;
		tra->Pvh_out = 0.0;
	}
	else
	{
		tra->Pvh_in = 0.0;
		tra->Pvh_out = 0.0;
		tra->Pvh_st = tra->vp * tra->ifl;
		tra->Pvh_cos = 0.0;
		/*
		tra->Pvh_in1 = 0.0;
		tra->Pvh_in2 = 0.0;
		tra->Pvh_in3 = 0.0;
		tra->Pvh_in4 = 0.0;
		tra->Pvh_in5 = 0.0;
		tra->Pvh_in6 = 0.0;
		tra->Pvh_in7 = 0.0;
		tra->Pvh_in8 = 0.0;
		*/
	}
	tra->Pvh_c = tra->Pvh / 1000.0;
	tra->Pvh_in_c = tra->Pvh_in / 1000.0;
	tra->Pvh_out_c = tra->Pvh_out / 1000.0;
	tra->Pvh_st_c = tra->Pvh_st / 1000;
}

void Calculation_energy_train(TRAIN* tra)
{
	if (tra->accelflag == 1)				//力行中
	{
		//tra->Emot_pow += tra->Pmot * dt;

		tra->Eloss_inv_pow += tra->Ploss_inv * dt;
		tra->Eloss_mot_pow += tra->Ploss_mot * dt;
		tra->Eloss_fl_pow += tra->Ploss_fl * dt;
		tra->Eloss_all_pow += tra->Ploss_all * dt;
		tra->Esiv_pow += P_SIV * tra->num * dt;

		tra->Eres_pow += tra->Pres * dt;
		tra->Emot_pow += (tra->Pvh_in - tra->Ploss_all - P_SIV * tra->num - tra->Pres) * dt;
	}
	else if (tra->accelflag == 3)		//回生中
	{
		//tra->Emot_reg += tra->Pmot * dt;
		//tra->Emot_reg_loss += tra->Pmot_reg_loss * dt;

		tra->Eloss_inv_reg += tra->Ploss_inv * dt;
		tra->Eloss_mot_reg += tra->Ploss_mot * dt;
		tra->Eloss_fl_reg += tra->Ploss_fl * dt;
		tra->Eloss_all_reg += tra->Ploss_all * dt;
		tra->Esiv_reg += P_SIV * tra->num * dt;

		tra->Eres_reg += tra->Pres * dt;
		tra->Emot_reg += (tra->Pvh_out - tra->Ploss_all - P_SIV * tra->num - tra->Pres) * dt;
		tra->Ereg_nashi += tra->Preg_nashi * dt;
	}
	else if (tra->accelflag == 2)
	{
		tra->Eres_coa += tra->Pres * dt;
		tra->Esiv_coa += P_SIV * tra->num * dt;
	}
	else if (tra->accelflag == 4)
	{
		tra->Esiv_stp += P_SIV * tra->num * dt;

	}
	tra->Emot_air += tra->Pmot_air * dt;


	tra->Emot_pow_c = tra->Emot_pow / 1000.0 / 3600.0;
	tra->Emot_reg_c = tra->Emot_reg / 1000.0 / 3600.0;
	tra->Emot_air_c = tra->Emot_air / 1000.0 / 3600.0;
	tra->Emot_reg_loss_c = tra->Emot_reg_loss / 1000.0 / 3600.0;
	tra->Ereg_nashi_c = tra->Ereg_nashi / 1000.0 / 3600.0;

	tra->Eloss_inv_pow_c = tra->Eloss_inv_pow / 1000.0 / 3600.0;
	tra->Eloss_mot_pow_c = tra->Eloss_mot_pow / 1000.0 / 3600.0;
	tra->Eloss_fl_pow_c = tra->Eloss_fl_pow / 1000.0 / 3600.0;
	tra->Eloss_all_pow_c = tra->Eloss_all_pow / 1000.0 / 3600.0;

	tra->Eloss_inv_reg_c = tra->Eloss_inv_reg / 1000.0 / 3600.0;
	tra->Eloss_mot_reg_c = tra->Eloss_mot_reg / 1000.0 / 3600.0;
	tra->Eloss_fl_reg_c = tra->Eloss_fl_reg / 1000.0 / 3600.0;
	tra->Eloss_all_reg_c = tra->Eloss_all_reg / 1000.0 / 3600.0;

	tra->Eres_pow_c = tra->Eres_pow / 1000.0 / 3600.0;
	tra->Eres_reg_c = tra->Eres_reg / 1000.0 / 3600.0;
	tra->Eres_coa_c = tra->Eres_coa / 1000.0 / 3600.0;

	tra->Esiv_pow_c = tra->Esiv_pow / 1000.0 / 3600.0;
	tra->Esiv_reg_c = tra->Esiv_reg / 1000.0 / 3600.0;
	tra->Esiv_coa_c = tra->Esiv_coa / 1000.0 / 3600.0;
	tra->Esiv_stp_c = tra->Esiv_stp / 1000.0 / 3600.0;

	tra->Evh += tra->Pvh * dt;
	tra->Evh_in += tra->Pvh_in * dt;
	tra->Evh_out += tra->Pvh_out * dt;
	tra->Evh_st += tra->Pvh_st * dt;
	tra->Evh_cos += tra->Pvh_cos * dt;

	/*tra->Evh_in1 += tra->Pvh_in1 * dt;
	tra->Evh_in2 += tra->Pvh_in2 * dt;
	tra->Evh_in3 += tra->Pvh_in3 * dt;
	tra->Evh_in4 += tra->Pvh_in4 * dt;
	tra->Evh_in5 += tra->Pvh_in5 * dt;
	tra->Evh_in6 += tra->Pvh_in6 * dt;
	tra->Evh_in7 += tra->Pvh_in7 * dt;
	tra->Evh_in8 += tra->Pvh_in8 * dt;
	*/
	tra->Evh_c = tra->Evh / 1000.0 / 3600.0;
	tra->Evh_in_c = tra->Evh_in / 1000.0 / 3600.0;
	tra->Evh_out_c = tra->Evh_out / 1000.0 / 3600.0;
	tra->Evh_cos_c = tra->Evh_cos / 1000.0 / 3600.0;
	tra->Evh_st_c = tra->Evh_st / 1000.0 / 3600.0;

	/*tra->Evh_in1_c = tra->Evh_in1 / 1000.0 / 3600.0;
	tra->Evh_in2_c = tra->Evh_in2 / 1000.0 / 3600.0;
	tra->Evh_in3_c = tra->Evh_in3 / 1000.0 / 3600.0;
	tra->Evh_in4_c = tra->Evh_in4 / 1000.0 / 3600.0;
	tra->Evh_in5_c = tra->Evh_in5 / 1000.0 / 3600.0;
	tra->Evh_in6_c = tra->Evh_in6 / 1000.0 / 3600.0;
	tra->Evh_in7_c = tra->Evh_in7 / 1000.0 / 3600.0;
	tra->Evh_in8_c = tra->Evh_in8 / 1000.0 / 3600.0;
	*/
}

void Calculation_power_sub(SUB* sub)
{
	sub->Pss = sub->iss * sub->vout;
	sub->Pss_c = sub->Pss / 1000.0;
}

void Calculation_energy_sub(SUB* sub)
{
	sub->Ess += sub->Pss * dt;						//[J]
	sub->Ess_c = sub->Ess / 1000.0 / 3600.0;		//[kWh]
	if (sub->diode == 1)
	{
		sub->Wss += sub->Rss * sub->iss * sub->iss * dt;		//[J]
		sub->Wss_c = sub->Wss / 1000 / 3600;				//[kWh]
	}
}

/*** xおよびyが指す要素を交換 ***/
void swap(NODE* x, NODE* y)
{
	NODE temp = { 0.0 };
	temp = *x;
	*x = *y;
	*y = temp;
}


/*** 配列data[]の先頭n個の要素を距離の昇順にソート ***/
void sort(NODE* data[], int n)
{
	int k = n - 1;
	while (k >= 0)
	{
		int i, j;
		for (i = 1, j = -1; i <= k; i++)
			if (data[i - 1]->X > data[i]->X) {
				j = i - 1;
				swap(data[i], data[j]);
			}
		k = j;
	}
}

void Initialize_BRANCH(BRANCH** branch)
{
	*branch = (BRANCH*)malloc(sizeof(BRANCH));
	if (*branch == NULL)
	{
		//エラー処理
	}
}

int Count_Trains_Direction(TRAIN* tra[], int direction) {
	int count_d = 0;
	for (int i = 0; i < NUM_tra; i++) {
		if (tra[i]->direction == direction) count_d++;
	}
	return count_d;
}

void Make_branch(BRANCH* branch[], NODE* temp[], TRAIN* tra[])
{
	int i, j;						//j：正方向に接続されているノード探索用，k：負方向に接続されているノード探索用
	int count_SS, count_UP, count_DOWN, NUM_tra_UP;	//ブランチ作成用カウント変数

	count_SS = 0;
	count_UP = NUM_sub;
	NUM_tra_UP = Count_Trains_Direction(tra, 1);

	//printf("%d\n", NUM_tra_UP);

	count_DOWN = 2 * NUM_sub + NUM_tra_UP - 1;

	/*①SSのブランチを作成*/
	//printf("-----------SS--------------\n");
	for (i = 0; i < N_node; i++)
	{
		//printf("[%d]flag = %d\n", temp[i]->Number, temp[i]->flag);
		if (temp[i]->flag == 0)				//In case of SS( i Starts with 0 to NUM_sub-1 )
		{
			//puts("SS");
			branch[count_SS]->X_start = temp[i]->X;		//位置[m]（基準点からの距離）
			branch[count_SS]->X_end = temp[i]->X;		//位置[m]（基準点からの距離）
			branch[count_SS]->r = temp[i]->r;			//Resistance [Ohm]
			branch[count_SS]->flag = temp[i]->flag;		//ブランチの分類（1：UP, -1：DOWN, 0：SS）
			branch[count_SS]->Node_pos = -1;			//正方向側に接続されているノード№
			branch[count_SS]->Node_neg = -1;			//負方向側に接続されているノード№
			count_SS = count_SS + 1;
		}
	}

	/*②上りのブランチを作成*/
	//printf("-----------UP--------------\n");
	for (i = 0; i < N_node; i++)
	{
		if (temp[i]->flag >= 0)			//In case of +1 direction
		{
			//printf("[%d]flag = %d\n", temp[i]->Number, temp[i]->flag);
			if (count_UP < NUM_tra_UP + NUM_sub + 1)		//終端以外
			{

				for (j = 0; temp[i + (j + 1)]->flag == -1; j++) {}				//temp[i]に対して正方向側に接続されているノードを探索
				//printf("\nj = %d\n", j);

				/*temp[i]に対して正方向に接続されているブランチBを作成*/
				branch[count_UP]->Node_pos = temp[i + (j + 1)]->Number;			//ブランチBについて，正方向側に接続されているノード№
				branch[count_UP]->Node_neg = temp[i]->Number;					//ブランチBについて，負方向側に接続されているノード№

				branch[count_UP]->X_start = temp[i]->X;							//[m]
				branch[count_UP]->X_end = temp[i + (j + 1)]->X;					//[m]


				//
				if (temp[i]->X >= 500.0 && temp[i + (j + 1)]->X <= 3500.0) {
					branch[count_UP]->r = R * fabs(temp[i + (j + 1)]->X - temp[i]->X) / 1000.0;		//Resistance [Ohm]
				}
				else {
					branch[count_UP]->r = R1 * fabs(temp[i + (j + 1)]->X - temp[i]->X) / 1000.0;		//Resistance [Ohm]
				}


				branch[count_UP]->flag = 1;										//ブランチの分類（1：UP, -1：DOWN, 0：SS）

				/*次のループのための準備*/
				j = 0;
				count_UP = count_UP + 1;
			}
			else {}		//終端のノードでは何もしない
		}
	}


	/*③下りのブランチを作成*/
	//printf("-----------DOWN--------------\n");
	/*for (i = 0; i < N_node; i++)
	{
		if (temp[i]->flag <= 0)			//In case of +1 direction
		{
			//printf("[%d]flag = %d\n", temp[i]->Number, temp[i]->flag);
			if (count_DOWN < N_branch)		//終端以外
			{
				for (j = 0; temp[i + (j + 1)]->flag == 1; j++) {}				//temp[i]に対して正方向側に接続されているノードを探索
				//printf("\nj = %d\n", j);

				/*正方向に接続されているブランチを作成*/
				/*			branch[count_DOWN]->Node_pos = temp[i + (j + 1)]->Number;		//正方向側に接続されているノード№
							branch[count_DOWN]->Node_neg = temp[i]->Number;					//負方向側に接続されているノード№

							branch[count_DOWN]->X_start = temp[i]->X;						//[m]
							branch[count_DOWN]->X_end = temp[i + (j + 1)]->X;				//[m]



							if (temp[i + (j + 1)]->X <= 2850.0) {
								branch[count_DOWN]->r = R * fabs(temp[i + (j + 1)]->X - temp[i]->X) / 1000.0;		//Resistance [Ohm]
							}
							else {
								branch[count_DOWN]->r = R * fabs(temp[i + (j + 1)]->X - temp[i]->X) / 1000.0;		//Resistance [Ohm]
							}


							branch[count_DOWN]->flag = -1;									//ブランチの分類（1：UP, -1：DOWN, 0：SS）

							/*次のループのための準備*/
							/*			j = 0;
										count_DOWN = count_DOWN + 1;
									}
									else {}		//終端のノードでは何もしない
								}
							}*/
}

void Initialize_matrix1(double Hm[], int M, int N)		/*M行N列の２次元配列を初期化する関数*/
{
	int i, j;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			Hm[i * N + j] = 0.0;
		}
	}
}

void Make_H_matrix(BRANCH* branch[], double* Hm)		/*H matrix作成用関数*/
{
	int i;	/*行*/
	int j;	/*列*/

	for (i = 0; i < N_node; i++)
	{
		for (j = 0; j < NUM_sub; j++)
		{
			//printf("\ni = %d\nj = %d\n",i,j);
			if (i == j)
			{
				//printf("\ni = %d\nj = %d\n",i,j);
				Hm[i + j * N_node] = -1.0;
			}
			else
			{
				Hm[i + j * N_node] = 0.0;
			}
		}
	}
	for (i = 0; i < N_node; i++)
	{
		for (j = NUM_sub; j < N_branch; j++)
		{
			if (i == branch[j]->Node_pos)
			{
				Hm[i + j * N_node] = 1.0;
			}
			else if (i == branch[j]->Node_neg)
			{
				Hm[i + j * N_node] = -1.0;
			}
			else
			{
				Hm[i + j * N_node] = 0.0;
			}
		}
	}
}

void Make_Y_matrix(BRANCH* branch[], double* Y)		/*Y matrix作成用関数*/
{
	int i;	/*行*/
	int j;	/*列*/

	for (i = 0; i < N_branch; i++)
	{
		for (j = 0; j < N_branch; j++)
		{
			//printf("\ni = %d\nj = %d\n",i,j);
			if (i == j)
			{
				//printf("\ni = %d\nj = %d\n",i,j);
				if (branch[j]->r < 0.0001)
				{
					Y[i + j * N_branch] = 10000.0;
				}
				else
				{
					Y[i + j * N_branch] = 1.0 / branch[j]->r;		/*ここでコンダクタンス[S]に変換*/
				}
			}
			else
			{
				Y[i + j * N_branch] = 0.0;
			}
		}
	}
}

void Make_In_vector(NODE* data[], double In[])		/*In vector作成用関数*/
{
	int i;

	for (i = 0; i < N_node; i++)
	{
		In[i] = data[i]->I;
	}
}

void Transpose(double* trans, double* X, int row, int column)		/*転置行列計算用関数*/
{
	int i, j;
	double temp;

	for (i = 0; i < row; i++)
	{
		for (j = 0; j < column; j++)
		{
			temp = X[i + j * row];
			trans[j + i * column] = temp;
		}
	}
}

double Calculate_Loss1(SUB* sub[], BRANCH* branch[], double X[], int n)
{
	double loss;
	int i;

	loss = 0.0;
	for (i = NUM_sub; i < n; i++)
		if (branch[i]->X_start >= 490.0 && branch[i]->X_end <= sub[0]->Xss)
		{
			{
				loss += branch[i]->r * X[i] * X[i] * dt;
			}
		}
	return loss;
}

double Calculate_Loss2(TRAIN* tra[], SUB* sub[], BRANCH* branch[], double X[], int n)
{
	double loss;
	int i;

	loss = 0.0;
	for (i = NUM_sub; i < n; i++)
		if (branch[i]->X_start >= sub[0]->Xss && branch[i]->X_end <= tra[0]->x)
		{
			{
				loss += branch[i]->r * X[i] * X[i] * dt;
			}
		}
	return loss;
}

double Calculate_Loss3(TRAIN* tra[], SUB* sub[], BRANCH* branch[], double X[], int n)
{
	double loss;
	int i;

	loss = 0.0;
	for (i = NUM_sub; i < n; i++)
		if (branch[i]->X_start >= tra[0]->x && branch[i]->X_end <= sub[1]->Xss)
		{
			{
				loss += branch[i]->r * X[i] * X[i] * dt;
			}
		}
	return loss;
}

double Calculate_Loss4(SUB* sub[], BRANCH* branch[], double X[], int n)
{
	double loss;
	int i;

	loss = 0.0;
	for (i = NUM_sub; i < n; i++)
		if (branch[i]->X_start >= sub[1]->Xss && branch[i]->X_end <= 3600.0)
		{
			{
				loss += branch[i]->r * X[i] * X[i] * dt;
			}
		}
	return loss;
}

/*
double Calculate_Loss5(SUB* sub[], BRANCH* branch[], double X[], int n)
{
	double loss;
	int i;

	loss = 0.0;
	for (i = NUM_sub; i < n; i++)
		if (branch[i]->X_start >= sub[3]->Xss && branch[i]->X_end <= sub[4]->Xss)
		{
			{
				loss += branch[i]->r * X[i] * X[i] * dt;
			}
		}
	return loss;
}

double Calculate_Loss6(SUB* sub[], BRANCH* branch[], double X[], int n)
{
	double loss;
	int i;

	loss = 0.0;
	for (i = NUM_sub; i < n; i++)
		if (branch[i]->X_start >= sub[4]->Xss && branch[i]->X_end <= sub[5]->Xss)
		{
			{
				loss += branch[i]->r * X[i] * X[i] * dt;
			}
		}
	return loss;
}

double Calculate_Loss7(SUB* sub[], BRANCH* branch[], double X[], int n)
{
	double loss;
	int i;

	loss = 0.0;
	for (i = NUM_sub; i < n; i++)
		if (branch[i]->X_start >= sub[5]->Xss && branch[i]->X_end <= sub[6]->Xss)
		{
			{
				loss += branch[i]->r * X[i] * X[i] * dt;
			}
		}
	return loss;
}

double Calculate_Loss8(SUB* sub[], BRANCH* branch[], double X[], int n)
{
	double loss;
	int i;

	loss = 0.0;
	for (i = NUM_sub; i < n; i++)
		if (branch[i]->X_start >= sub[6]->Xss && branch[i]->X_end <= 21700.0)
		{
			{
				loss += branch[i]->r * X[i] * X[i] * dt;
			}
		}
	return loss;
}
*/

double _CalcLoss(BRANCH* _bra, double _Ib, size_t i) {
	auto ls = (_bra->r * _Ib * _Ib * dt);
	return ls;

}

double Calculate_Loss(BRANCH* branch[], double X[], int n)
{
	double loss;
	int i;

	loss = 0.0;
	for (i = NUM_sub; i < n; i++)
	{
		loss += branch[i]->r * X[i] * X[i] * dt;
	}
	return loss;
}

double _CalculateLoss(BRANCH* _bra[], double _Ib[])
{
	auto loss = 0.0;

	for (size_t i = 0; i < N_branch; i++)
	{
		if (_bra[i]->flag != 0) // SSブランチ出なければき電損失をカウントする
		{
			loss += _bra[i]->r * _Ib[i] * _Ib[i] * dt;
		}
	}
	return loss;
}

void lsm(double x[], double y[], int N, double* a0, double* a1)
{
	int i;
	double A00 = 0.0, A01 = 0.0, A02 = 0.0, A11 = 0.0, A12 = 0.0;

	for (i = 0; i < N; i++) {
		A00 += 1.0;
		A01 += x[i];
		A02 += y[i];
		A11 += x[i] * x[i];
		A12 += x[i] * y[i];
	}

	*a0 = (A02 * A11 - A01 * A12) / (A00 * A11 - A01 * A01);
	*a1 = (A00 * A12 - A01 * A02) / (A00 * A11 - A01 * A01);
}

bool createFolder(const char* path) {
	if (CreateDirectoryA(path, NULL)) {
		return true;
	}
	else {
		return false;
	}
}


void Calculation_power(TRAIN* tra) {
	if (tra->accelflag == 1) {
		if (tra->x >= 0.0 && tra->x <= 500.0) { tra->Pvh_in1 = tra->vp * tra->ifl; tra->Evh_in1 += tra->Pvh_in1 * dt; tra->Evh_in1_c = tra->Evh_in1 / 1000.0 / 3600.0; }
		else if (tra->x > 500.0 && tra->x <= 1000.0) { tra->Pvh_in2 = tra->vp * tra->ifl; tra->Evh_in2 += tra->Pvh_in2 * dt; tra->Evh_in2_c = tra->Evh_in2 / 1000.0 / 3600.0; }
		else if (tra->x > 1000.0 && tra->x <= 1500.0) { tra->Pvh_in3 = tra->vp * tra->ifl; tra->Evh_in3 += tra->Pvh_in3 * dt; tra->Evh_in3_c = tra->Evh_in3 / 1000.0 / 3600.0; }
		else if (tra->x > 1500.0 && tra->x <= 2000.0) { tra->Pvh_in4 = tra->vp * tra->ifl; tra->Evh_in4 += tra->Pvh_in4 * dt; tra->Evh_in4_c = tra->Evh_in4 / 1000.0 / 3600.0; }
		else if (tra->x > 2000.0 && tra->x <= 2500.0) { tra->Pvh_in5 = tra->vp * tra->ifl; tra->Evh_in5 += tra->Pvh_in5 * dt; tra->Evh_in5_c = tra->Evh_in5 / 1000.0 / 3600.0; }
		else if (tra->x > 2500.0 && tra->x <= 3000.0) { tra->Pvh_in6 = tra->vp * tra->ifl; tra->Evh_in6 += tra->Pvh_in6 * dt; tra->Evh_in6_c = tra->Evh_in6 / 1000.0 / 3600.0; }
		else if (tra->x > 3000.0 && tra->x <= 3500.0) { tra->Pvh_in7 = tra->vp * tra->ifl; tra->Evh_in7 += tra->Pvh_in7 * dt; tra->Evh_in7_c = tra->Evh_in7 / 1000.0 / 3600.0; }
		else if (tra->x > 3500.0 && tra->x <= 4000.0) { tra->Pvh_in8 = tra->vp * tra->ifl; tra->Evh_in8 += tra->Pvh_in8 * dt; tra->Evh_in8_c = tra->Evh_in8 / 1000.0 / 3600.0; }
	}
	if (tra->accelflag == 3) {
		if (tra->x >= 0.0 && tra->x <= 500.0) { tra->Pvh_out1 = tra->vp * tra->ifl; tra->Evh_out1 += tra->Pvh_out1 * dt; tra->Evh_out1_c = tra->Evh_out1 / 1000.0 / 3600.0; }
		else if (tra->x > 500.0 && tra->x <= 1000.0) { tra->Pvh_out2 = tra->vp * tra->ifl; tra->Evh_out2 += tra->Pvh_out2 * dt; tra->Evh_out2_c = tra->Evh_out2 / 1000.0 / 3600.0; }
		else if (tra->x > 1000.0 && tra->x <= 1500.0) { tra->Pvh_out3 = tra->vp * tra->ifl; tra->Evh_out3 += tra->Pvh_out3 * dt; tra->Evh_out3_c = tra->Evh_out3 / 1000.0 / 3600.0; }
		else if (tra->x > 1500.0 && tra->x <= 2000.0) { tra->Pvh_out4 = tra->vp * tra->ifl; tra->Evh_out4 += tra->Pvh_out4 * dt; tra->Evh_out4_c = tra->Evh_out4 / 1000.0 / 3600.0; }
		else if (tra->x > 2000.0 && tra->x <= 2500.0) { tra->Pvh_out5 = tra->vp * tra->ifl; tra->Evh_out5 += tra->Pvh_out5 * dt; tra->Evh_out5_c = tra->Evh_out5 / 1000.0 / 3600.0; }
		else if (tra->x > 2500.0 && tra->x <= 3000.0) { tra->Pvh_out6 = tra->vp * tra->ifl; tra->Evh_out6 += tra->Pvh_out6 * dt; tra->Evh_out6_c = tra->Evh_out6 / 1000.0 / 3600.0; }
		else if (tra->x > 3000.0 && tra->x <= 3500.0) { tra->Pvh_out7 = tra->vp * tra->ifl; tra->Evh_out7 += tra->Pvh_out7 * dt; tra->Evh_out7_c = tra->Evh_out7 / 1000.0 / 3600.0; }
		else if (tra->x > 3500.0 && tra->x <= 4000.0) { tra->Pvh_out8 = tra->vp * tra->ifl; tra->Evh_out8 += tra->Pvh_out8 * dt; tra->Evh_out8_c = tra->Evh_out8 / 1000.0 / 3600.0; }
	}
	if (tra->accelflag == 3) {
		tra->ifl_sum += tra->ifl * dt;
	}

}

int main(void)
{
	int i, l;
	l = 0;

	FILE* fp6;				//ダイヤ作成用
	Make_file6(&fp6);

	double pow_energy[vel_pattern];
	double reg_energy[vel_pattern];
	double cos_energy[vel_pattern];
	double ave_runtime[vel_pattern];
	double sub_output[vel_pattern];
	double fed_loss[vel_pattern];
	double fed_loss1[vel_pattern];
	double fed_loss2[vel_pattern];
	double fed_loss3[vel_pattern];
	double fed_loss4[vel_pattern];
	/*
	double a_pow[vss_pattern * vss_pattern * vss_pattern];
	double a_reg[vss_pattern * vss_pattern * vss_pattern];
	double a_sub[vss_pattern * vss_pattern * vss_pattern];
	double a_fed[vss_pattern * vss_pattern * vss_pattern];

	double b_pow[vss_pattern * vss_pattern * vss_pattern];
	double b_reg[vss_pattern * vss_pattern * vss_pattern];
	double b_sub[vss_pattern * vss_pattern * vss_pattern];
	double b_fed[vss_pattern * vss_pattern * vss_pattern];

	double cal_pow[vss_pattern * vss_pattern * vss_pattern];
	double cal_reg[vss_pattern * vss_pattern * vss_pattern];
	double cal_sub[vss_pattern * vss_pattern * vss_pattern];
	double cal_fed[vss_pattern * vss_pattern * vss_pattern];

	Initialize_matrix1(a_pow, vss_pattern * vss_pattern * vss_pattern, 1);
	Initialize_matrix1(b_pow, vss_pattern * vss_pattern * vss_pattern, 1);
	Initialize_matrix1(a_reg, vss_pattern * vss_pattern * vss_pattern, 1);
	Initialize_matrix1(b_reg, vss_pattern * vss_pattern * vss_pattern, 1);
	Initialize_matrix1(a_sub, vss_pattern * vss_pattern * vss_pattern, 1);
	Initialize_matrix1(b_sub, vss_pattern * vss_pattern * vss_pattern, 1);
	Initialize_matrix1(a_fed, vss_pattern * vss_pattern * vss_pattern, 1);
	Initialize_matrix1(b_fed, vss_pattern * vss_pattern * vss_pattern, 1);
	Initialize_matrix1(cal_pow, vss_pattern * vss_pattern * vss_pattern, 1);
	Initialize_matrix1(cal_reg, vss_pattern * vss_pattern * vss_pattern, 1);
	Initialize_matrix1(cal_sub, vss_pattern * vss_pattern * vss_pattern, 1);
	Initialize_matrix1(cal_fed, vss_pattern * vss_pattern * vss_pattern, 1);
	*/
	double a_pow[vss_pattern * vss_pattern];
	double a_reg[vss_pattern * vss_pattern];
	double a_cos[vss_pattern * vss_pattern];
	double a_sub[vss_pattern * vss_pattern];
	double a_fed[vss_pattern * vss_pattern];
	double a_fed1[vss_pattern * vss_pattern];
	double a_fed2[vss_pattern * vss_pattern];
	double a_fed3[vss_pattern * vss_pattern];
	double a_fed4[vss_pattern * vss_pattern];

	double b_pow[vss_pattern * vss_pattern];
	double b_reg[vss_pattern * vss_pattern];
	double b_cos[vss_pattern * vss_pattern];
	double b_sub[vss_pattern * vss_pattern];
	double b_fed[vss_pattern * vss_pattern];
	double b_fed1[vss_pattern * vss_pattern];
	double b_fed2[vss_pattern * vss_pattern];
	double b_fed3[vss_pattern * vss_pattern];
	double b_fed4[vss_pattern * vss_pattern];

	double cal_pow[vss_pattern * vss_pattern];
	double cal_reg[vss_pattern * vss_pattern];
	double cal_st[vss_pattern * vss_pattern];
	double cal_sub[vss_pattern * vss_pattern];
	double cal_fed[vss_pattern * vss_pattern];
	double cal_fed1[vss_pattern * vss_pattern];
	double cal_fed2[vss_pattern * vss_pattern];
	double cal_fed3[vss_pattern * vss_pattern];
	double cal_fed4[vss_pattern * vss_pattern];

	Initialize_matrix1(a_pow, vss_pattern * vss_pattern, 1);
	Initialize_matrix1(b_pow, vss_pattern * vss_pattern, 1);
	Initialize_matrix1(a_reg, vss_pattern * vss_pattern, 1);
	Initialize_matrix1(b_reg, vss_pattern * vss_pattern, 1);
	Initialize_matrix1(a_cos, vss_pattern * vss_pattern, 1);
	Initialize_matrix1(b_cos, vss_pattern * vss_pattern, 1);
	Initialize_matrix1(a_sub, vss_pattern * vss_pattern, 1);
	Initialize_matrix1(b_sub, vss_pattern * vss_pattern, 1);
	Initialize_matrix1(a_fed, vss_pattern * vss_pattern, 1);
	Initialize_matrix1(b_fed, vss_pattern * vss_pattern, 1);
	Initialize_matrix1(a_fed1, vss_pattern * vss_pattern, 1);
	Initialize_matrix1(b_fed1, vss_pattern * vss_pattern, 1);
	Initialize_matrix1(a_fed2, vss_pattern * vss_pattern, 1);
	Initialize_matrix1(b_fed2, vss_pattern * vss_pattern, 1);
	Initialize_matrix1(a_fed3, vss_pattern * vss_pattern, 1);
	Initialize_matrix1(b_fed3, vss_pattern * vss_pattern, 1);
	Initialize_matrix1(a_fed4, vss_pattern * vss_pattern, 1);
	Initialize_matrix1(b_fed4, vss_pattern * vss_pattern, 1);
	Initialize_matrix1(cal_pow, vss_pattern * vss_pattern, 1);
	Initialize_matrix1(cal_reg, vss_pattern * vss_pattern, 1);
	Initialize_matrix1(cal_st, vss_pattern * vss_pattern, 1);
	Initialize_matrix1(cal_sub, vss_pattern * vss_pattern, 1);
	Initialize_matrix1(cal_fed, vss_pattern * vss_pattern, 1);
	Initialize_matrix1(cal_fed1, vss_pattern * vss_pattern, 1);
	Initialize_matrix1(cal_fed2, vss_pattern * vss_pattern, 1);
	Initialize_matrix1(cal_fed3, vss_pattern * vss_pattern, 1);
	Initialize_matrix1(cal_fed4, vss_pattern * vss_pattern, 1);

	/*** 列車構造体配列初期化用の構造体配列宣言部 ***/
	/*オフピークダイヤ*/
	INI_TRA ini_tra[NUM_tra] = {	//M車重量[kg/両]，T車重量[kg/両]，M車両数，T車両数，初期位置[m]，初期速度[km/h]，上下判定，発車時刻, 行先,列車種別,待ち合わせフラグ
		/*		{28.8 * 1000.0, 24.675 * 1000.0, 6.0, 4.0, 4000.0, 0.0, -1.0, 0.0 * 60.0, 0, 0, 0},     //0
				{28.8 * 1000.0, 24.675 * 1000.0, 6.0, 4.0, 4000.0, 0.0, -1.0, 5.0 * 60.0, 0, 0, 0},    //1
				{28.8 * 1000.0, 24.675 * 1000.0, 6.0, 4.0, 4000.0, 0.0, -1.0, 10.0 * 60.0, 0, 0, 0},    //2
				{28.8 * 1000.0, 24.675 * 1000.0, 6.0, 4.0, 4000.0, 0.0, -1.0, 15.0 * 60.0, 0, 0, 0},    //3
				{28.8 * 1000.0, 24.675 * 1000.0, 6.0, 4.0, 4000.0, 0.0, -1.0, 20.0 * 60.0, 0, 0, 0},    //4
				{28.8 * 1000.0, 24.675 * 1000.0, 6.0, 4.0, 4000.0, 0.0, -1.0, 25.0 * 60.0, 0, 0, 0},    //5
				{28.8 * 1000.0, 24.675 * 1000.0, 6.0, 4.0, 4000.0, 0.0, -1.0, 30.0 * 60.0, 0, 0, 0},    //6
				{28.8 * 1000.0, 24.675 * 1000.0, 6.0, 4.0, 4000.0, 0.0, -1.0, 35.0 * 60.0, 0, 0, 0},    //7
				{28.8 * 1000.0, 24.675 * 1000.0, 6.0, 4.0, 4000.0, 0.0, -1.0, 40.0 * 60.0, 0, 0, 0},    //8
				{28.8 * 1000.0, 24.675 * 1000.0, 6.0, 4.0, 4000.0, 0.0, -1.0, 45.0 * 60.0, 0, 0, 0},    //9
				{28.8 * 1000.0, 24.675 * 1000.0, 6.0, 4.0, 4000.0, 0.0, -1.0, 50.0 * 60.0, 0, 0, 0},     //10
				{28.8 * 1000.0, 24.675 * 1000.0, 6.0, 4.0, 4000.0, 0.0, -1.0, 55.0 * 60.0, 0, 0, 0},     //11
		*/



		/*下り各停　56分～に注意！　実際+4分*/
		{28.8 * 1000.0, 24.675 * 1000.0, 6.0, 4.0, 1000.0, 0.0, 1.0, 0.0 * 60.0, 0, 0, 0},			//12
		/*		{28.8 * 1000.0, 24.675 * 1000.0, 6.0, 4.0, 0.0, 0.0, 1.0, 5.0 * 60.0, 0, 0, 0},			//13
				{28.8 * 1000.0, 24.675 * 1000.0, 6.0, 4.0, 0.0, 0.0, 1.0, 10.0 * 60.0, 0, 0, 0},		//14
				{28.8 * 1000.0, 24.675 * 1000.0, 6.0, 4.0, 0.0, 0.0, 1.0, 15.0 * 60.0, 0, 0, 0},		//15
				{28.8 * 1000.0, 24.675 * 1000.0, 6.0, 4.0, 0.0, 0.0, 1.0, 20.0 * 60.0, 0, 0, 0},		//16
				{28.8 * 1000.0, 24.675 * 1000.0, 6.0, 4.0, 0.0, 0.0, 1.0, 25.0 * 60.0, 0, 0, 0},		//17
				{28.8 * 1000.0, 24.675 * 1000.0, 6.0, 4.0, 0.0, 0.0, 1.0, 30.0 * 60.0, 0, 0, 0},		//18
				{28.8 * 1000.0, 24.675 * 1000.0, 6.0, 4.0, 0.0, 0.0, 1.0, 35.0 * 60.0, 0, 0, 0},		//19
				{28.8 * 1000.0, 24.675 * 1000.0, 6.0, 4.0, 0.0, 0.0, 1.0, 40.0 * 60.0, 0, 0, 0},		//20
				{28.8 * 1000.0, 24.675 * 1000.0, 6.0, 4.0, 0.0, 0.0, 1.0, 45.0 * 60.0, 0, 0, 0},		//21
				{28.8 * 1000.0, 24.675 * 1000.0, 6.0, 4.0, 0.0, 0.0, 1.0, 50.0 * 60.0, 0, 0, 0},		//22
				{28.8 * 1000.0, 24.675 * 1000.0, 6.0, 4.0, 0.0, 0.0, 1.0, 55.0 * 60.0, 0, 0, 0},		//23
		*/

		/*上り快速　56分～に注意！　実際+4分*/
/*	{31.8625 * 1000.0, 34.67142857 * 1000.0, 8.0, 7.0, 20508.0, 0.0, -1.0, 10.0 * 60.0, 0, 1, 0},     //31
	{31.8625 * 1000.0, 34.67142857 * 1000.0, 8.0, 7.0, 20508.0, 0.0, -1.0, 21.0 * 60.0, 0, 1, 0},     //32
	{31.8625 * 1000.0, 34.67142857 * 1000.0, 8.0, 7.0, 20508.0, 0.0, -1.0, 30.0 * 60.0, 0, 1, 1},     //33
	{31.8625 * 1000.0, 34.67142857 * 1000.0, 8.0, 7.0, 20508.0, 0.0, -1.0, 44.0 * 60.0, 0, 1, 0},     //34
	{31.8625 * 1000.0, 34.67142857 * 1000.0, 8.0, 7.0, 20508.0, 0.0, -1.0, 59.0 * 60.0, 0, 1, 0},     //35
	{31.8625 * 1000.0, 34.67142857 * 1000.0, 8.0, 7.0, 12698.0, 0.0, -1.0, 8.0 * 60.0, 0, 1, 0},      //36 通過待ち便
	*/

	/*下り快速　56分～に注意！　実際+4分*/
/*	{31.8625 * 1000.0, 34.67142857 * 1000.0, 8.0, 7.0, 2158.0, 0.0, 1.0, 6.0 * 60.0, 0, 1, 1},      //37
	{31.8625 * 1000.0, 34.67142857 * 1000.0, 8.0, 7.0, 2158.0, 0.0, 1.0, 18.0 * 60.0, 0, 1, 0},     //38
	{31.8625 * 1000.0, 34.67142857 * 1000.0, 8.0, 7.0, 2158.0, 0.0, 1.0, 28.0 * 60.0, 0, 1, 0},     //39
	{31.8625 * 1000.0, 34.67142857 * 1000.0, 8.0, 7.0, 2158.0, 0.0, 1.0, 37.0 * 60.0, 0, 1, 1},     //40
	{31.8625 * 1000.0, 34.67142857 * 1000.0, 8.0, 7.0, 2158.0, 0.0, 1.0, 53.0 * 60.0, 0, 1, 0},     //41
	{31.8625 * 1000.0, 34.67142857 * 1000.0, 8.0, 7.0, 12698.0, 0.0, 1.0, 7.0 * 60.0, 0, 1, 0},     //42
*/

	{28.8 * 1000.0, 24.675 * 1000.0, 6.0, 4.0, 490.0, 0.0, 1.0, 9999.0 * 60.0, 0, 4, 0},			//定負荷
	{28.8 * 1000.0, 24.675 * 1000.0, 6.0, 4.0, 3510.0, 0.0, 1.0, 9999.0 * 60.0, 0, 4, 0},			//定負荷

	};

	/***  駅の構造体配列宣言部 (新百合ヶ丘駅基準) ***/
	STATION sta[] = {	//各停
		{"basement1", -100.0},			//駅名, 駅位置	0
		{"ryogoku", 1000.0},					//駅名, 駅位置	1
		{"Kinshicho", 3000.0},			//駅名, 駅位置	2
		//	{"Kameido",	4000.0},				//駅名, 駅位置	3
			{"basement2", 3100.0},			//駅名, 駅位置	12
	};

	STATION sta2[] = {	//快速
		{"basement1",-100.0},			//駅名, 駅位置	0
		{"Kinshicho",2158.0},					//駅名, 駅位置	1
		{"Shin-koiwa",	7313.0},				//駅名, 駅位置	2
		{"Ichikawa",12698.0},			//駅名, 駅位置	3
		{"Funabashi",20508.0},			//駅名, 駅位置	4
		{"basement2",	20608.0},			//駅名, 駅位置	5
	};

	STATION sta3[] = {
		{"basement1",-1100.0},			//駅名, 駅位置	0
		{"ryogoku", -1000.0},					//駅名, 駅位置	1
		{"Funabashi",21508.0},			//駅名, 駅位置	4
		{"basement2",	21608.0},			//駅名, 駅位置	5
	};

	for (int i = 0; i < vss_pattern; i++) {
		stringstream ss1;
		ss1 << "錦糸町vss" << fixed << setprecision(1) << vss_dif[i];
		vss_fol1[i] = ss1.str();
	}

	for (int i = 0; i < vss_pattern; i++) {
		stringstream ss2;
		ss2 << "平井vss" << fixed << setprecision(1) << vss_dif[i];
		vss_fol2[i] = ss2.str();
	}
	/*
	for (int i = 0; i < vss_pattern; i++) {
		stringstream ss3;
		ss3 << "小岩vss" << fixed << setprecision(1) << vss_dif[i];
		vss_fol3[i] = ss3.str();
	}
	for (int i = 0; i < vss_pattern; i++) {
		stringstream ss4;
		ss4 << "市川vss" << fixed << setprecision(1) << vss_dif[i];
		vss_fol4[i] = ss4.str();
	}

	for (int i = 0; i < vss_pattern; i++) {
		stringstream ss5;
		ss5 << "西船橋vss" << fixed << setprecision(1) << vss_dif[i];
		vss_fol5[i] = ss5.str();
	}
	*/
	for (int i = 0; i < vel_pattern; i++) {
		stringstream ssu;
		ssu << "v" << fixed << setprecision(1) << vel_dif[i];
		vel_fol[i] = ssu.str();
	}


	if (ALL_Vss_Change == 0) {

		// Create folder for each element in A
		for (int h = 0; h < vss_pattern; h++) {
			string folderNameA = "Result\\" + vss_fol1[h];

			// Create folder
			if (createFolder(folderNameA.c_str())) {
				l++;
			}
			else {
				l--;
			}
			// Create folder for each element in B inside folder A
			for (int i = 0; i < vss_pattern; i++) {
				string folderNameB = folderNameA + "\\" + vss_fol2[i];

				// Create folder
				if (createFolder(folderNameB.c_str())) {
					l++;
				}
				else {
					l--;
				}
				// Create folder for each element in B inside folder A
				for (int j = 0; j < vel_pattern; j++) {
					string folderNameC = folderNameB + "\\" + vel_fol[j];

					// Create folder
					if (createFolder(folderNameC.c_str())) {
						l++;
					}
					else {
						l--;
					}
				}
			}



		}

	}
	else {

		// Create folder for each element in A
		for (int h = 0; h < vss_pattern; h++) {
			string folderNameA = "Result\\" + vss_fol1[h];

			// Create folder
			if (createFolder(folderNameA.c_str())) {
				l++;
			}
			else {
				l--;
			}

			string folderNameB = folderNameA + "\\" + vss_fol2[h];

			// Create folder
			if (createFolder(folderNameB.c_str())) {
				l++;
			}
			else {
				l--;
			}

			string folderNameC = folderNameB + "\\" + vss_fol3[h];

			// Create folder
			if (createFolder(folderNameC.c_str())) {
				l++;
			}
			else {
				l--;
			}

			for (int j = 0; j < vlb_pattern; j++) {
				string folderNameD = folderNameC + "\\" + vlb_fol[j];

				// Create folder
				if (createFolder(folderNameD.c_str())) {
					l++;
				}
				else {
					l--;
				}

				// Create folder for each element in B inside folder A
				for (int j = 0; j < vel_pattern; j++) {
					string folderNameE = folderNameD + "\\" + vel_fol[j];

					// Create folder
					if (createFolder(folderNameE.c_str())) {
						l++;
					}
					else {
						l--;
					}
				}
			}

		}
	}
	while (sim_endflag == 0) {

		/***  各行列宣言部  ***/
		double H[N_node * N_branch];		//Connection Matrix
		double H_tra[N_branch * N_node];	//Transpose matrix of H
		double Y[N_branch * N_branch];		//Conductance Matrix
		double A[N_node * N_node];			//A matrix(H*Y*H_tra)
		double A_inv[N_node * N_node];		//Inverse of A matrix

		double TEMP[N_node * N_branch];	//Tast

		double In[N_node];					//Node Current Matrix
		double In_cpy[N_node];				//Copy of Node Current Matrix
		double Vn[N_node];					//Node Voltage Matrix
		double Ib[N_branch];				//Branch Current Matrix
		double Vb[N_branch];				//Branch Voltage Matrix
		//double ave_runtime = 0.0;

		int i, j, k, l;
		k = 0;
		l = 0;

		int ERROR_SS = 0;						//回路計算の再計算判定用フラグ（flagを全て加算してNUM_sub以上なら回路計算ループを抜ける。NUM_sub以下なら再計算。）
		int ERROR_REG = 0;
		int count_loop = 0;
		//	int test_flag = 0;

			/***  マトリックス計算関係  ***/
		const char trans = 'N'; //Normalな場合。これは規定値'N','T','C'から選択
		int K_node = N_node;
		int K_branch = N_branch;
		int ld_node = N_node;
		int ld_branch = N_branch;
		double alpha = 1.0;
		double beta = 0.0;

		int* ipiv_Vn;
		ipiv_Vn = (int*)calloc(N_node, sizeof(int));
		int info_Vn;
		int nrhs = 1;

		int incx = 1;
		int incy = 1;

		end_flag = 0;

		/***  構造体宣言部  ***/
		TRAIN* tra[NUM_tra];				//ポインタの配列として定義
		SUB* sub[NUM_sub];					//ポインタの配列として定義
		NODE* node[N_node];
		NODE* node_order[N_node];
		BRANCH* branch[N_branch];
		//	ESD *esd[NUM_esd];					//ポインタの配列として定義
		//	INV *inv[NUM_inv];					//ポインタの配列として定義

			/*** き電回路送電損失格納変数　***/
		double Lcir1 = 0.0;							//[J]
		double Lcir1_c = 0.0;						//[kWh]
		double Lcir2 = 0.0;							//[J]
		double Lcir2_c = 0.0;						//[kWh]
		double Lcir3 = 0.0;							//[J]
		double Lcir3_c = 0.0;						//[kWh]
		double Lcir4 = 0.0;							//[J]
		double Lcir4_c = 0.0;						//[kWh]
		double Lcir5 = 0.0;							//[J]
		double Lcir5_c = 0.0;						//[kWh]
		double Lcir6 = 0.0;							//[J]
		double Lcir6_c = 0.0;						//[kWh]
		double Lcir7 = 0.0;							//[J]
		double Lcir7_c = 0.0;						//[kWh]
		double Lcir8 = 0.0;							//[J]
		double Lcir8_c = 0.0;						//[kWh]
		double Lcir = 0.0;							//[J]
		double Lcir_c = 0.0;						//[kWh]

		/*** 変電所構造体配列初期化用の構造体配列宣言 ***/
		INI_SUB ini_sub[NUM_sub] = {
			{500.0, ess1, Iss_rated, Vss1 + vss_dif[vss_state1], 0.0, 1},						//[両国]	変電所位置，電圧変動率，等価出力抵抗，出力電流，出力電圧, 架線電圧，ダイオード
			{3500.0, ess1, Iss_rated, Vss1 + vss_dif[vss_state2], 0.0, 1},						//[錦糸町]	変電所位置，電圧変動率，等価出力抵抗，出力電流，出力電圧, 架線電圧，ダイオード
			//	{4500.0, ess1, Iss_rated, Vss1 + vss_dif[vss_state3], 0.0, 1},
		};

		/** 出力データファイル **/
		FILE* fp1[NUM_tra];		//車両走行計算結果出力用ファイル
		FILE* fp2;				//車両関係エネルギー計算結果出力用ファイル
		FILE* fp3;				//変電所計算結果出力用ファイル
		FILE* fp4;				//変電所出力エネルギー計算結果出力用ファイル
		FILE* fp5;				//ダイヤ作成用
		FILE* fpx;				//x
		FILE* fpy;				//y
		FILE* fpz;				//z


		Make_file2(&fp2);
		fprintf(fp2, "<力行時>,力行エネルギー[kWh],走行ロスエネルギー(力行時)[kWh],インバータ損失エネルギー[kWh],モータ損失エネルギー[kWh],FL損失エネルギー[kWh],補機消費エネルギー[kWh],<回生時>,モータ回生エネルギー[kWh],機械ブレーキ損失エネルギー[kWh],走行ロスエネルギー(回生時)[kWh],インバータ損失エネルギー[kWh],モータ損失エネルギー[kWh],FL損失エネルギー[kWh],補機消費エネルギー[kWh],<惰行時>,走行ロスエネルギー(惰行時)[kWh],補機消費エネルギー[kWh],<停車時>,補器消費エネルギー[kWh],,流入(力行),流出(回生),流入(惰行停止),全体,回生(軽負荷なし),力行～1[kWh],力行1～2[kWh],力行2～3[kWh],力行3～4[kWh],力行4～5[kWh],力行5～6[kWh],力行6～7[kWh],力行7～[kWh],回生～1[kWh],回生1～2[kWh],回生2～3[kWh],回生3～4[kWh],回生4～5[kWh],回生5～6[kWh],回生6～7[kWh],回生7～[kWh]\n");
		Make_file4(&fp4);
		fprintf(fp4, "<総入力エネルギー>,変電所1出力エネルギー[kWh],変電所1損失[kWh],変電所2出力エネルギー[kWh],変電所2損失[kWh],変電所3出力エネルギー[kWh],変電所3損失[kWh],変電所4出力エネルギー[kWh],変電所4損失[kWh],変電所5出力エネルギー[kWh],変電所5損失[kWh],変電所6出力エネルギー[kWh],変電所6損失[kWh],変電所7出力エネルギー[kWh],変電所7損失[kWh], 合計変電所出力エネルギー[kWh], き電損失～1[kWh],き電損失1～2[kWh],き電損失2～3[kWh],き電損失3～4[kWh],き電損失4～5[kWh],き電損失5～6[kWh],き電損失6～7[kWh],き電損失7～[kWh],合計き電損失[kWh]\n");
		printf("Definition completed\n");

		/***  CSVファイル1オープン&一行目書き込み  ***/

		for (i = 0; i < NUM_tra; i++)
		{
			Make_file1(&fp1[i], i);
			fprintf(fp1[i], "時間[sec],時間[min],ノッチ,列車位置1[m],速度1[km/h],加速度1[m/s/s],停止距離1[m],列車引張力1[kN],回生ブレーキ力1[kN],回生ブレーキ力(フィルタ)1[kN],機械ブレーキ力1[kN],走行抵抗1[kN],架線電圧[V],FC電圧[V],FL電流[A],補機電流[A],モータ電力1[kW],回生電力(軽負荷なし)1[kW]\n");
		}

		//Make_file_x(&fpx);
		//Make_file_y(&fpy);
		//Make_file_z(&fpz);
		Make_file3(&fp3);
		fprintf(fp3, "時間[sec],時間[min],変電所1出力電圧[V],変電所1出力電流[A],変電所1定電流源[A],変電所1出力パワー[kW],変電所2出力電圧[V],変電所2出力電流[A],変電所2定電流源[A],変電所2出力パワー[kW],　列車～SS1への電流, 列車～SS2への電流\n");

		Make_file5(&fp5);

		if (Time_diagram == 1) {
			//fprintf(fp5, "時間[min], Tra_0, Tra_1, Tra_2, Tra_3, Tra_4, Tra_5, Tra_6, Tra_7, Tra_8, Tra_9, Tra_10, Tra_11, Tra_12, Tra_13, Tra_14, Tra_15, Tra_16, Tra_17, Tra_18, Tra_19, Tra_20, Tra_21, Tra_22, Tra_23, Tra_24, Tra_25, Tra_26, Tra_27, Tra_28, Tra_29, Tra_30, Tra_31, Tra_32, Tra_33, Tra_34, Tra_35, Tra_36, Tra_37, Tra_38, Tra_39, Tra_40, Tra_41, Tra_42, Tra_43, Tra_44, Tra_45, Tra_46\n");  //各停、快速、特急
		}
		else if (Time_diagram == 2) {
			//fprintf(fp5, "時間[min], Tra_0, Tra_1, Tra_2, Tra_3, Tra_4, Tra_5, Tra_6, Tra_7, Tra_8, Tra_9, Tra_10, Tra_11, Tra_12, Tra_13, Tra_14, Tra_15, Tra_16, Tra_17, Tra_18, Tra_19, Tra_20, Tra_21, Tra_22, Tra_23, Tra_24, Tra_25, Tra_26, Tra_27, Tra_28, Tra_29, Tra_30, Tra_31, Tra_32, Tra_33, Tra_34, Tra_35, Tra_36, Tra_37, Tra_38, Tra_39, Tra_40, Tra_41, Tra_42\n");       //各停、快速
		}
		else if (Time_diagram == 3) {
			//fprintf(fp5, "時間[min], Tra_0, Tra_1, Tra_2, Tra_3, Tra_4, Tra_5, Tra_6, Tra_7, Tra_8, Tra_9, Tra_10, Tra_11, Tra_12, Tra_13, Tra_14, Tra_15, Tra_16, Tra_17, Tra_18, Tra_19, Tra_20, Tra_21, Tra_22, Tra_23, Tra_24, Tra_25, Tra_26, Tra_27, Tra_28, Tra_29, Tra_30\n");     //各停
		}

		if (vel_state == 0) {
			fprintf(fp6, "\nSS.A%.1f, SS.B%.1f, 走行時分[s], 力行電力量[kWh], 惰行電力量, 回生電力量[kWh], %s, SS.A[kWh],SS.B[kWh], 変電所出力電力量[kWh] , き電損失[kWh]\n", vss_dif[vss_state1], vss_dif[vss_state2], "");

			Initialize_matrix1(pow_energy, vel_pattern, 1);
			Initialize_matrix1(reg_energy, vel_pattern, 1);
			Initialize_matrix1(cos_energy, vel_pattern, 1);
			Initialize_matrix1(ave_runtime, vel_pattern, 1);
			Initialize_matrix1(sub_output, vel_pattern, 1);
			Initialize_matrix1(fed_loss, vel_pattern, 1);
			Initialize_matrix1(fed_loss1, vel_pattern, 1);
			Initialize_matrix1(fed_loss2, vel_pattern, 1);
			Initialize_matrix1(fed_loss3, vel_pattern, 1);

		}

		/***  列車構造体・変電所構造体・ESD構造体初期化  ***/
		for (i = 0; i < NUM_sub; i++)
		{
			Make_substation(&(sub[i]), &(ini_sub[i]), i);				//ナンバリングは通し番号
		}
		for (i = 0; i < NUM_tra; i++)
		{
			Make_train(&(tra[i]), &(ini_tra[i]), i + NUM_sub);							//ナンバリングは通し番号
		}

		for (i = 0; i < NUM_tra; i++)
		{
			Calculate_BEmax(tra[i]);
		}

		for (i = 0; i < NUM_tra; i++)
		{
			printf("%d : x=%f\n", tra[i]->name_train, tra[i]->x);
		}
		for (i = 0; i < NUM_sub; i++)
		{
			printf("%d : x=%f\n", sub[i]->name_SS, sub[i]->Xss);
		}

		////////////////////////// 数値計算部 /////////////////////////
		for (t = 0.0; t <= 60.0 * 2.5; t += dt)	/***** 1時間シミュレーション *****/
		{
			minute = t / 60.0;


			////////////////////////// 車両運動方程式計算部 /////////////////////
				/*** 車両走行パターン模擬 & 運動方程式計算 ***/
				/*** 車両パワー計算 ***/

			for (i = 0; i < NUM_tra; i++)
			{
				Calculate_BEmax(tra[i]);

				Run_pattern(tra[i], sta, sta2, sta3, t);

			}

			/////////////////////////////////////////////////////////////////////

		//	for (i = 0; i < NUM_sub; i++)
		//	{
		//		sub[i]->Jss = sub[i]->Jss_ini;
		//	}

			/*	if (t > 1020.05 && t < 1020.09) {
					printf("-----------t = %f--------------\n", t);
					printf("-----------accelflag vehicle--------------\n");
					for (i = 0; i < NUM_tra; i++)
					{
						printf("%d\n", tra[i]->accelflag);
					}
					printf("-----------Fmot vehicle--------------\n");
					for (i = 0; i < NUM_tra; i++)
					{
						printf("%f\n", tra[i]->Fmot);
					}
					printf("\n");
					printf("\n");
					printf("\n");
				}*/

			while (1)
			{
			start_while:				//ダイオード判定が間違っていた場合，ここから再計算する

				//if (count_loop < 10 && t > 1900 && tra[0]->accelflag == 2)
				if (count_loop > 300000)
				{
					printf("\n");
					printf("-----------t = %f--------------\n", t);
					printf("-----------accelflag & position vehicle--------------\n");
					for (i = 0; i < NUM_tra; i++)
					{
						printf("%d, %f\n", tra[i]->accelflag, tra[i]->x);
					}
					/*				printf("-----------Fmot vehicle--------------\n");
									for (i = 0; i < NUM_tra; i++)
									{
										printf("%f\n", tra[i]->Fmot);
									}*/
					printf("-----------Pmot vehicle--------------\n");
					for (i = 0; i < NUM_tra; i++)
					{
						printf("%f\n", tra[i]->Pmot);
					}
					printf("-----------Vp vehicle--------------\n");
					for (i = 0; i < NUM_tra; i++)
					{
						printf("%f\n", tra[i]->vp);
					}
					printf("-----------Vfc vehicle--------------\n");
					for (i = 0; i < NUM_tra; i++)
					{
						printf("%f\n", tra[i]->vfc);
					}
					printf("-----------ifl vehicle--------------\n");
					for (i = 0; i < NUM_tra; i++)
					{
						printf("%f\n", tra[i]->ifl);
					}
					/*		printf("-----------ifc vehicle--------------\n");
					for (i = 0; i < NUM_tra; i++)
					{
						printf("%f\n", tra[i]->ifc);
					}
					printf("-----------iinv vehicle--------------\n");
					for (i = 0; i < NUM_tra; i++)
					{
						printf("%f\n", tra[i]->iinv);
					}
									printf("-----------vss substation--------------\n");
									for (i = 0; i < NUM_sub; i++)
									{
										printf("%f\n", sub[i]->vout);
									}
									printf("-----------Jss substation--------------\n");
									for (i = 0; i < NUM_sub; i++)
									{
										printf("%f\n", sub[i]->Jss);
									}
									printf("-----------Jss_ini substation--------------\n");
									for (i = 0; i < NUM_sub; i++)
									{
										printf("%f\n", sub[i]->Jss_ini);
									}
									printf("-----------iss substation--------------\n");
									for (i = 0; i < NUM_sub; i++)
									{
										printf("%f\n", sub[i]->iss);
									}*/
					printf("\n");
					printf("\n");
					printf("\n");
				}

				if (count_loop > 300030) {
					printf("Holy shit!!\n");
					exit(EXIT_FAILURE);
				}

				/***  node[]の初期化  ***/
				for (i = 0; i < NUM_sub; i++)
				{
					Make_NODE_SS(&(node[i]), &(sub[i]));
					Make_NODE_SS(&(node_order[i]), &(sub[i]));
				}

				for (i = 0; i < NUM_tra; i++)
				{
					Make_NODE_TRAIN(&(node[i + NUM_sub]), &(tra[i]));
					Make_NODE_TRAIN(&(node_order[i + NUM_sub]), &(tra[i]));
				}

				for (i = 0; i < N_branch; i++)
				{
					Initialize_BRANCH(&(branch[i]));
				}

				/*** 各Matrixの初期化 ***/
				Initialize_matrix1(H, N_node, N_branch);
				Initialize_matrix1(H_tra, N_branch, N_node);
				Initialize_matrix1(Y, N_branch, N_branch);
				Initialize_matrix1(A, N_node, N_node);
				Initialize_matrix1(A_inv, N_node, N_node);
				Initialize_matrix1(TEMP, N_node, N_branch);
				Initialize_matrix1(In, N_node, 1);
				Initialize_matrix1(In_cpy, N_node, 1);
				Initialize_matrix1(Vn, N_node, 1);
				Initialize_matrix1(Ib, N_branch, 1);
				Initialize_matrix1(Vb, N_branch, 1);

				/***  node_order[]を距離順にソート  ***/
				sort(node_order, N_node);

				Make_branch(branch, node_order, tra);

				//H行列の生成(「回生車を含むき電システムの現状とあり方」p.19参照)
				Make_H_matrix(branch, H);

				//アドミタンス行列の生成
				Make_Y_matrix(branch, Y);

				//ノード電流行列の生成
				Make_In_vector(node, In);
				Make_In_vector(node, In_cpy);		//"dgetrs"では右辺ベクトルが方程式の解で上書きされてしまうので，計算用にInをコピーしたIn_cpyを作成

				//H行列を転置する
				Transpose(H_tra, H, N_node, N_branch);

				//H行列×Y行列 (dgemmはintel MKLを用いて行列の掛け算を計算)
				dgemm(&trans, &trans, &K_node, &K_branch, &K_branch, &alpha, H, &ld_node, Y, &ld_branch, &beta, TEMP, &ld_node);

				//上記行列×転置H行列 (dgemmはintel MKLを用いて行列の掛け算を計算)
				dgemm(&trans, &trans, &K_node, &K_node, &K_branch, &alpha, TEMP, &ld_node, H_tra, &ld_branch, &beta, A, &ld_node);

				//連立1次方程式を解くためにH行列×Y行列×転置H行列の解A(N行N列)をLU分解
				dgetrf(&ld_node, &K_branch, A, &ld_node, ipiv_Vn, &info_Vn);

				//dgetrsで連立1次方程式「A*Vn=In」を解き，Inに解(ノード電圧行列)を上書き
				dgetrs(&trans, &K_node, &nrhs, A, &ld_node, ipiv_Vn, In_cpy, &ld_node, &info_Vn);	//In_cpyに方程式の解Vnが格納されていることに注意

				Make_Y_matrix(branch, Y);		//なぜかVnを求める際の"dgetrf"でY-matrixが書き換えられているため，ここで再計算させてる（実際は無駄な計算なのでしたくないが…）

				for (i = 0; i < N_node; i++)	/*Vnベクトルはマイナスをつける*/
				{
					Vn[i] = -In_cpy[i];
				}

				//転置H行列×Vnベクトルでブランチ電圧ベクトルを計算
				dgemv(&trans, &K_branch, &K_node, &alpha, H_tra, &ld_branch, Vn, &incx, &beta, Vb, &incy);

				//Y行列×ブランチ電圧ベクトルでブランチ電流ベクトルを計算
				dgemv(&trans, &K_branch, &K_branch, &alpha, Y, &ld_branch, Vb, &incx, &beta, Ib, &incy);


				for (i = 0; i < N_node; i++)			//ノード電圧の計算結果をノード構造体へ返す（結果的には，Vn配列はノード構造体の順番通りに並んでいる？）
				{
					node[i]->V = Vn[i];
				}

				for (i = 0; i < NUM_sub; i++)			//ノード電圧・電流とブランチ電流から，変電所出力電圧・電流を計算
				{
					sub[i]->vout = node[i]->V;
					sub[i]->iss = -In[i] + Ib[i];
					Error_Detection_SS(sub[i]);
					ERROR_SS += sub[i]->flag;
				}

				for (i = 0; i < NUM_tra; i++)
				{
					for (j = 0; j < N_node; j++)
					{
						if (node[j]->Number == tra[i]->name_train)	//ノード番号と列車識別番号が一致していたらノード電圧を列車パンタ点電圧に返す
						{
							tra[i]->vp = node[j]->V;
							//	if (tra[i]->vp > Vcmax) tra[i]->vp = Vcmax;
							tra[i]->vfc = tra[i]->vp - tra[i]->Rfl * tra[i]->ifl;

							//引張力・ブレーキ力・モータ電流・モータ電圧を計算
							Calculation_traction_force_brake_force(tra[i]);

							//dq軸電流・電圧からモータパワーを計算
							Calculation_power_motor(tra[i]);

							//主回路計算
							Calculation_traction_circuit(tra[i]);
						}
					}
				}

				if (t >= 17.6 && t < 17.61) {
					printf("-----------t = %f--------------\n", t);
					printf("-----------H matrix--------------\n");
					for (i = 0; i < N_node; i++)
					{
						for (j = 0; j < N_branch; j++)
						{
							if (j < N_branch - 1)
							{
								printf("%f,", H[i + j * N_node]);
							}
							else
							{
								printf("%f\n", H[i + j * N_node]);
							}
						}
					}
					printf("-----------H_tra matrix--------------\n");
					for (i = 0; i < N_branch; i++)
					{
						for (j = 0; j < N_node; j++)
						{
							if (j < N_node - 1)
							{
								printf("%f,", H_tra[i + j * N_branch]);
							}
							else
							{
								printf("%f\n", H_tra[i + j * N_branch]);
							}
						}
					}
					printf("-----------Y matrix--------------\n");
					for (i = 0; i < N_branch; i++)
					{
						for (j = 0; j < N_branch; j++)
						{
							if (j < N_branch - 1)
							{
								printf("%f,", Y[i + j * N_branch]);
							}
							else
							{
								printf("%f\n", Y[i + j * N_branch]);
							}
						}
					}
					printf("-----------A matrix--------------\n");
					for (i = 0; i < N_node; i++)
					{
						for (j = 0; j < N_node; j++)
						{
							if (j < N_node - 1)
							{
								printf("%f,", A[i + j * N_node]);
							}
							else
							{
								printf("%f\n", A[i + j * N_node]);
							}
						}
					}
					printf("\n");
					printf("\n");
					printf("\n");
					test_flag++;
				}

				for (i = 0; i < NUM_tra; i++)
				{
					tra[i]->vfc_old = tra[i]->vfc;

					tra[i]->judgement = fabs(tra[i]->vfc - tra[i]->vfc_old) / tra[i]->vfc_old;

					if (tra[i]->judgement < 0.001) {
						tra[i]->reg_flag = 1;
					}
					else {
						tra[i]->reg_flag = 0;
					}

					if (tra[i]->vfc > Vcmax) tra[i]->reg_flag = 0;

					ERROR_REG += tra[i]->reg_flag;
				}

				/*		for (i = 0; i < NUM_tra; i++)
						{
							if (tra[i]->accelflag == 3 || (tra[i]->accelflag == 5 && tra[i]->v_c > tra[i]->con_speed)) {
								Error_Detection_REG(tra[i]);
							}
							else {
								tra[i]->reg_flag = 1;
							}

							ERROR_REG += tra[i]->reg_flag;
						}*/

						/*		if (t > 1020.05 && t < 1020.09) {
									printf("-----------t = %f--------------\n", t);
									printf("-----------accelflag vehicle--------------\n");
									for (i = 0; i < NUM_tra; i++)
									{
										printf("%d\n", tra[i]->accelflag);
									}
									printf("-----------Fmot vehicle--------------\n");
									for (i = 0; i < NUM_tra; i++)
									{
										printf("%f\n", tra[i]->Fmot);
									}
									printf("-----------Pmot vehicle--------------\n");
									for (i = 0; i < NUM_tra; i++)
									{
										printf("%f\n", tra[i]->Pmot);
									}
									printf("-----------Vp vehicle--------------\n");
									for (i = 0; i < NUM_tra; i++)
									{
										printf("%f\n", tra[i]->vp);
									}
									printf("-----------ifl vehicle--------------\n");
									for (i = 0; i < NUM_tra; i++)
									{
										printf("%f\n", tra[i]->ifl);
									}
									printf("-----------vss substation--------------\n");
									for (i = 0; i < NUM_sub; i++)
									{
										printf("%f\n", sub[i]->vout);
									}
									printf("-----------iss substation--------------\n");
									for (i = 0; i < NUM_sub; i++)
									{
										printf("%f\n", sub[i]->iss);
									}
									printf("\n");
									printf("\n");
									printf("\n");
								}*/

				if (ERROR_SS < NUM_sub || ERROR_REG < NUM_tra)
				{
					for (i = 0; i < N_node; i++)
					{
						free(node[i]);
						free(node_order[i]);
					}
					for (i = 0; i < N_branch; i++)
					{
						free(branch[i]);
					}

					ERROR_SS = 0;
					ERROR_REG = 0;
					count_loop++;
					goto start_while;					//while文の先頭に戻る
				}
				else
				{
					ERROR_SS = 0;
					ERROR_REG = 0;
					count_loop = 0;
					break;
				}
			}

			if (t >= 17.6 && t < 17.61) {

				printf("-----------Vn vector--------------\n");
				for (i = 0; i < N_node; i++)
				{
					printf("%f\n", Vn[i]);
				}
				printf("-----------In vector--------------\n");
				for (i = 0; i < N_node; i++)
				{
					printf("%f\n", In[i]);
				}
				printf("-----------Vp vehicle--------------\n");
				for (i = 0; i < NUM_tra; i++)
				{
					printf("%f\n", tra[i]->vp);
				}
				printf("-----------vfc vehicle--------------\n");
				for (i = 0; i < NUM_tra; i++)
				{
					printf("%f\n", tra[i]->vfc);
				}
				printf("-----------ifl vehicle--------------\n");
				for (i = 0; i < NUM_tra; i++)
				{
					printf("%f\n", tra[i]->ifl);
				}
				printf("-----------x vehicle--------------\n");
				for (i = 0; i < NUM_tra; i++)
				{
					printf("%f\n", tra[i]->x);
				}
				printf("-----------vss substation--------------\n");
				for (i = 0; i < NUM_sub; i++)
				{
					printf("%f\n", sub[i]->vout);
				}
				printf("-----------iss substation--------------\n");
				for (i = 0; i < NUM_sub; i++)
				{
					printf("%f\n", sub[i]->iss);
				}
			}
			/***  次の状態を計算  ***/
			for (i = 0; i < NUM_tra; i++)
			{
				//tra[i]->b_regF_old = tra[i]->b_regF;
				if (tra[i]->accelflag != 4) {
					Solve_next_state(tra[i]);
					tra[i]->run_time += dt;
				}
			}

			/***  車両エネルギー計算  ***/
			for (i = 0; i < NUM_tra; i++)
			{
				Calculation_power_train(tra[i]);
				Calculation_energy_train(tra[i]);
				Calculation_power(tra[i]);
			}

			for (i = 0; i < NUM_sub; i++)
			{
				Calculation_power_sub(sub[i]);
				Calculation_energy_sub(sub[i]);
			}

			//Lcir1 += Calculate_Loss1(sub, branch, Ib, N_branch);
			//Lcir1_c = Lcir1 / 1000.0 / 3600.0;
			//fed_loss1[vel_state] = Lcir1_c;
			//Lcir2 += Calculate_Loss2(tra, sub, branch, Ib, N_branch);
			//Lcir2_c = Lcir2 / 1000.0 / 3600.0;
			//fed_loss2[vel_state] = Lcir2_c;
			//Lcir3 += Calculate_Loss3(tra, sub, branch, Ib, N_branch);
			//Lcir3_c = Lcir3 / 1000.0 / 3600.0;
			//fed_loss3[vel_state] = Lcir3_c;
			//Lcir4 += Calculate_Loss4(sub, branch, Ib, N_branch);
			//Lcir4_c = Lcir4 / 1000.0 / 3600.0;
			//fed_loss4[vel_state] = Lcir4_c;
			///*	Lcir5 += Calculate_Loss5(sub, branch, Ib, N_branch);
			//	Lcir5_c = Lcir5 / 1000.0 / 3600.0;
			//	Lcir6 += Calculate_Loss6(sub, branch, Ib, N_branch);
			//	Lcir6_c = Lcir6 / 1000.0 / 3600.0;
			//	Lcir7 += Calculate_Loss7(sub, branch, Ib, N_branch);
			//	Lcir7_c = Lcir7 / 1000.0 / 3600.0;
			//	Lcir8 += Calculate_Loss8(sub, branch, Ib, N_branch);
			//	Lcir8_c = Lcir8 / 1000.0 / 3600.0;*/
			//Lcir += Calculate_Loss(branch, Ib, N_branch);
			//Lcir_c = Lcir / 1000.0 / 3600.0;
			//fed_loss[vel_state] = Lcir_c;

			// 川﨑追記
			// branch 0, 1 はグランドー架線を繋ぐブランチなので、2,3,4,5番を見る
			Lcir1 += _CalcLoss(branch[2], Ib[2], 2);
			Lcir1_c = Lcir1 / 1000.0 / 3600.0;
			fed_loss1[vel_state] = Lcir1_c;
			Lcir2 += _CalcLoss(branch[3], Ib[3], 3);
			Lcir2_c = Lcir2 / 1000.0 / 3600.0;
			fed_loss2[vel_state] = Lcir2_c;
			Lcir3 += _CalcLoss(branch[4], Ib[4], 4);
			Lcir3_c = Lcir3 / 1000.0 / 3600.0;
			fed_loss3[vel_state] = Lcir3_c;
			Lcir4 += _CalcLoss(branch[5], Ib[5], 5);
			Lcir4_c = Lcir4 / 1000.0 / 3600.0;
			fed_loss4[vel_state] = Lcir4_c;
			/*	Lcir5 += Calculate_Loss5(sub, branch, Ib, N_branch);
				Lcir5_c = Lcir5 / 1000.0 / 3600.0;
				Lcir6 += Calculate_Loss6(sub, branch, Ib, N_branch);
				Lcir6_c = Lcir6 / 1000.0 / 3600.0;
				Lcir7 += Calculate_Loss7(sub, branch, Ib, N_branch);
				Lcir7_c = Lcir7 / 1000.0 / 3600.0;
				Lcir8 += Calculate_Loss8(sub, branch, Ib, N_branch);
				Lcir8_c = Lcir8 / 1000.0 / 3600.0;*/
			Lcir += _CalculateLoss(branch, Ib);
			Lcir_c = Lcir / 1000.0 / 3600.0;
			fed_loss[vel_state] = Lcir_c;

			//

			for (i = 0; i < N_node; i++)
			{
				free(node[i]);
				free(node_order[i]);
			}

			for (i = 0; i < N_branch; i++)
			{
				free(branch[i]);
			}

			////////////////////////// CSVファイル出力部 /////////////////////////
			t_out = t;//[s]
			if (t_out - t_out_old >= (1.00 - dt) * 0.1)
			{
				t_out_old = t_out;
				for (i = 0; i < NUM_tra; i++)
				{
					//fprintf(fp1[i], "%f,%f,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", t, minute, tra[i]->accelflag, tra[i]->x, tra[i]->v_c, tra[i]->a, tra[i]->X_brake, tra[i]->Ftot_c, tra[i]->b_reg_c, tra[i]->b_air_c, tra[i]->R_total / 1000.0, tra[i]->idtot,
					fprintf(fp1[i], "%f,%f,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", t, minute, tra[i]->accelflag, tra[i]->x, tra[i]->v_c, tra[i]->a, tra[i]->X_brake, tra[i]->Ftot_c, tra[i]->b_reg_c, tra[i]->b_regF_c, tra[i]->b_air_c, tra[i]->R_total / 1000.0, tra[i]->vp, tra[i]->vfc, tra[i]->ifl, tra[i]->isiv, tra[i]->Pmot_c, tra[i]->Preg_nashi_c);
				}
				fprintf(fp3, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", t, minute, sub[0]->vout, sub[0]->iss, sub[0]->Jss, sub[0]->Pss_c, sub[1]->vout, sub[1]->iss, sub[1]->Jss, sub[1]->Pss_c, sub[0]->iss - tra[1]->ifl, sub[1]->iss - tra[2]->ifl);

				if (F_flag == 1) {
					if (Time_diagram == 1) {
						//	fprintf(fpx, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute);
						//	fprintf(fpy, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", tra[0]->x, tra[1]->x, tra[2]->x, tra[3]->x, tra[4]->x, tra[5]->x, tra[6]->x, tra[7]->x, tra[8]->x, tra[9]->x, tra[10]->x, tra[11]->x, tra[12]->x, tra[13]->x, tra[14]->x, tra[15]->x, tra[16]->x, tra[17]->x, tra[18]->x, tra[19]->x, tra[20]->x, tra[21]->x, tra[22]->x, tra[23]->x, tra[24]->x, tra[25]->x, tra[26]->x, tra[27]->x, tra[28]->x, tra[29]->x, tra[30]->x, tra[31]->x, tra[32]->x, tra[33]->x, tra[34]->x, tra[35]->x, tra[36]->x, tra[37]->x, tra[38]->x, tra[39]->x, tra[40]->x, tra[41]->x, tra[42]->x, tra[43]->x, tra[44]->x, tra[45]->x, tra[46]->x);
						//	fprintf(fpz, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", tra[0]->Pvh_in, tra[1]->Pvh_in, tra[2]->Pvh_in, tra[3]->Pvh_in, tra[4]->Pvh_in, tra[5]->Pvh_in, tra[6]->Pvh_in, tra[7]->Pvh_in, tra[8]->Pvh_in, tra[9]->Pvh_in, tra[10]->Pvh_in, tra[11]->Pvh_in, tra[12]->Pvh_in, tra[13]->Pvh_in, tra[14]->Pvh_in, tra[15]->Pvh_in, tra[16]->Pvh_in, tra[17]->Pvh_in, tra[18]->Pvh_in, tra[19]->Pvh_in, tra[20]->Pvh_in, tra[21]->Pvh_in, tra[22]->Pvh_in, tra[23]->Pvh_in, tra[24]->Pvh_in, tra[25]->Pvh_in, tra[26]->Pvh_in, tra[27]->Pvh_in, tra[28]->Pvh_in, tra[29]->Pvh_in, tra[30]->Pvh_in, tra[31]->Pvh_in, tra[32]->Pvh_in, tra[33]->Pvh_in, tra[34]->Pvh_in, tra[35]->Pvh_in, tra[36]->Pvh_in, tra[37]->Pvh_in, tra[38]->Pvh_in, tra[39]->Pvh_in, tra[40]->Pvh_in, tra[41]->Pvh_in, tra[42]->Pvh_in, tra[43]->Pvh_in, tra[44]->Pvh_in, tra[45]->Pvh_in, tra[46]->Pvh_in);
					}
					else if (Time_diagram == 2) {
						//	fprintf(fpx, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute);
						//	fprintf(fpy, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", tra[0]->x, tra[1]->x, tra[2]->x, tra[3]->x, tra[4]->x, tra[5]->x, tra[6]->x, tra[7]->x, tra[8]->x, tra[9]->x, tra[10]->x, tra[11]->x, tra[12]->x, tra[13]->x, tra[14]->x, tra[15]->x, tra[16]->x, tra[17]->x, tra[18]->x, tra[19]->x, tra[20]->x, tra[21]->x, tra[22]->x, tra[23]->x, tra[24]->x, tra[25]->x, tra[26]->x, tra[27]->x, tra[28]->x, tra[29]->x, tra[30]->x, tra[31]->x, tra[32]->x, tra[33]->x, tra[34]->x, tra[35]->x, tra[36]->x, tra[37]->x, tra[38]->x, tra[39]->x, tra[40]->x, tra[41]->x, tra[42]->x);
						//	fprintf(fpz, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", tra[0]->Pvh_in, tra[1]->Pvh_in, tra[2]->Pvh_in, tra[3]->Pvh_in, tra[4]->Pvh_in, tra[5]->Pvh_in, tra[6]->Pvh_in, tra[7]->Pvh_in, tra[8]->Pvh_in, tra[9]->Pvh_in, tra[10]->Pvh_in, tra[11]->Pvh_in, tra[12]->Pvh_in, tra[13]->Pvh_in, tra[14]->Pvh_in, tra[15]->Pvh_in, tra[16]->Pvh_in, tra[17]->Pvh_in, tra[18]->Pvh_in, tra[19]->Pvh_in, tra[20]->Pvh_in, tra[21]->Pvh_in, tra[22]->Pvh_in, tra[23]->Pvh_in, tra[24]->Pvh_in, tra[25]->Pvh_in, tra[26]->Pvh_in, tra[27]->Pvh_in, tra[28]->Pvh_in, tra[29]->Pvh_in, tra[30]->Pvh_in, tra[31]->Pvh_in, tra[32]->Pvh_in, tra[33]->Pvh_in, tra[34]->Pvh_in, tra[35]->Pvh_in, tra[36]->Pvh_in, tra[37]->Pvh_in, tra[38]->Pvh_in, tra[39]->Pvh_in, tra[40]->Pvh_in, tra[41]->Pvh_in, tra[42]->Pvh_in);
					}
					else if (Time_diagram == 3) {
						//	fprintf(fpx, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute);
						//	fprintf(fpy, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", tra[0]->x, tra[1]->x, tra[2]->x, tra[3]->x, tra[4]->x, tra[5]->x, tra[6]->x, tra[7]->x, tra[8]->x, tra[9]->x, tra[10]->x, tra[11]->x, tra[12]->x, tra[13]->x, tra[14]->x, tra[15]->x, tra[16]->x, tra[17]->x, tra[18]->x, tra[19]->x, tra[20]->x, tra[21]->x, tra[22]->x, tra[23]->x, tra[24]->x, tra[25]->x, tra[26]->x, tra[27]->x, tra[28]->x, tra[29]->x, tra[30]->x);
						//	fprintf(fpz, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", tra[0]->Pvh_in, tra[1]->Pvh_in, tra[2]->Pvh_in, tra[3]->Pvh_in, tra[4]->Pvh_in, tra[5]->Pvh_in, tra[6]->Pvh_in, tra[7]->Pvh_in, tra[8]->Pvh_in, tra[9]->Pvh_in, tra[10]->Pvh_in, tra[11]->Pvh_in, tra[12]->Pvh_in, tra[13]->Pvh_in, tra[14]->Pvh_in, tra[15]->Pvh_in, tra[16]->Pvh_in, tra[17]->Pvh_in, tra[18]->Pvh_in, tra[19]->Pvh_in, tra[20]->Pvh_in, tra[21]->Pvh_in, tra[22]->Pvh_in, tra[23]->Pvh_in, tra[24]->Pvh_in, tra[25]->Pvh_in, tra[26]->Pvh_in, tra[27]->Pvh_in, tra[28]->Pvh_in, tra[29]->Pvh_in, tra[30]->Pvh_in);
					}
				}
				else if (F_flag == 2) {
					if (Time_diagram == 1) {
						//	fprintf(fpx, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute);
						//	fprintf(fpy, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", tra[0]->x, tra[1]->x, tra[2]->x, tra[3]->x, tra[4]->x, tra[5]->x, tra[6]->x, tra[7]->x, tra[8]->x, tra[9]->x, tra[10]->x, tra[11]->x, tra[12]->x, tra[13]->x, tra[14]->x, tra[15]->x, tra[16]->x, tra[17]->x, tra[18]->x, tra[19]->x, tra[20]->x, tra[21]->x, tra[22]->x, tra[23]->x, tra[24]->x, tra[25]->x, tra[26]->x, tra[27]->x, tra[28]->x, tra[29]->x, tra[30]->x, tra[31]->x, tra[32]->x, tra[33]->x, tra[34]->x, tra[35]->x, tra[36]->x, tra[37]->x, tra[38]->x, tra[39]->x, tra[40]->x, tra[41]->x, tra[42]->x, tra[43]->x, tra[44]->x, tra[45]->x, tra[46]->x);
						//  fprintf(fpz, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", tra[0]->Pvh_out, tra[1]->Pvh_out, tra[2]->Pvh_out, tra[3]->Pvh_out, tra[4]->Pvh_out, tra[5]->Pvh_out, tra[6]->Pvh_out, tra[7]->Pvh_out, tra[8]->Pvh_out, tra[9]->Pvh_out, tra[10]->Pvh_out, tra[11]->Pvh_out, tra[12]->Pvh_out, tra[13]->Pvh_out, tra[14]->Pvh_out, tra[15]->Pvh_out, tra[16]->Pvh_out, tra[17]->Pvh_out, tra[18]->Pvh_out, tra[19]->Pvh_out, tra[20]->Pvh_out, tra[21]->Pvh_out, tra[22]->Pvh_out, tra[23]->Pvh_out, tra[24]->Pvh_out, tra[25]->Pvh_out, tra[26]->Pvh_out, tra[27]->Pvh_out, tra[28]->Pvh_out, tra[29]->Pvh_out, tra[30]->Pvh_out, tra[31]->Pvh_out, tra[32]->Pvh_out, tra[33]->Pvh_out, tra[34]->Pvh_out, tra[35]->Pvh_out, tra[36]->Pvh_out, tra[37]->Pvh_out, tra[38]->Pvh_out, tra[39]->Pvh_out, tra[40]->Pvh_out, tra[41]->Pvh_out, tra[42]->Pvh_out, tra[43]->Pvh_out, tra[44]->Pvh_out, tra[45]->Pvh_out, tra[46]->Pvh_out);
					}
					else if (Time_diagram == 2) {
						//	fprintf(fpx, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute);
						//	fprintf(fpy, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", tra[0]->x, tra[1]->x, tra[2]->x, tra[3]->x, tra[4]->x, tra[5]->x, tra[6]->x, tra[7]->x, tra[8]->x, tra[9]->x, tra[10]->x, tra[11]->x, tra[12]->x, tra[13]->x, tra[14]->x, tra[15]->x, tra[16]->x, tra[17]->x, tra[18]->x, tra[19]->x, tra[20]->x, tra[21]->x, tra[22]->x, tra[23]->x, tra[24]->x, tra[25]->x, tra[26]->x, tra[27]->x, tra[28]->x, tra[29]->x, tra[30]->x, tra[31]->x, tra[32]->x, tra[33]->x, tra[34]->x, tra[35]->x, tra[36]->x, tra[37]->x, tra[38]->x, tra[39]->x, tra[40]->x, tra[41]->x, tra[42]->x);
						//  fprintf(fpz, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", tra[0]->Pvh_out, tra[1]->Pvh_out, tra[2]->Pvh_out, tra[3]->Pvh_out, tra[4]->Pvh_out, tra[5]->Pvh_out, tra[6]->Pvh_out, tra[7]->Pvh_out, tra[8]->Pvh_out, tra[9]->Pvh_out, tra[10]->Pvh_out, tra[11]->Pvh_out, tra[12]->Pvh_out, tra[13]->Pvh_out, tra[14]->Pvh_out, tra[15]->Pvh_out, tra[16]->Pvh_out, tra[17]->Pvh_out, tra[18]->Pvh_out, tra[19]->Pvh_out, tra[20]->Pvh_out, tra[21]->Pvh_out, tra[22]->Pvh_out, tra[23]->Pvh_out, tra[24]->Pvh_out, tra[25]->Pvh_out, tra[26]->Pvh_out, tra[27]->Pvh_out, tra[28]->Pvh_out, tra[29]->Pvh_out, tra[30]->Pvh_out, tra[31]->Pvh_out, tra[32]->Pvh_out, tra[33]->Pvh_out, tra[34]->Pvh_out, tra[35]->Pvh_out, tra[36]->Pvh_out, tra[37]->Pvh_out, tra[38]->Pvh_out, tra[39]->Pvh_out, tra[40]->Pvh_out, tra[41]->Pvh_out, tra[42]->Pvh_out);
					}
					else if (Time_diagram == 3) {
						//	fprintf(fpx, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute, minute);
						//	fprintf(fpy, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", tra[0]->x, tra[1]->x, tra[2]->x, tra[3]->x, tra[4]->x, tra[5]->x, tra[6]->x, tra[7]->x, tra[8]->x, tra[9]->x, tra[10]->x, tra[11]->x, tra[12]->x, tra[13]->x, tra[14]->x, tra[15]->x, tra[16]->x, tra[17]->x, tra[18]->x, tra[19]->x, tra[20]->x, tra[21]->x, tra[22]->x, tra[23]->x, tra[24]->x, tra[25]->x, tra[26]->x, tra[27]->x, tra[28]->x, tra[29]->x, tra[30]->x);
						//  fprintf(fpz, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", tra[0]->Pvh_out, tra[1]->Pvh_out, tra[2]->Pvh_out, tra[3]->Pvh_out, tra[4]->Pvh_out, tra[5]->Pvh_out, tra[6]->Pvh_out, tra[7]->Pvh_out, tra[8]->Pvh_out, tra[9]->Pvh_out, tra[10]->Pvh_out, tra[11]->Pvh_out, tra[12]->Pvh_out, tra[13]->Pvh_out, tra[14]->Pvh_out, tra[15]->Pvh_out, tra[16]->Pvh_out, tra[17]->Pvh_out, tra[18]->Pvh_out, tra[19]->Pvh_out, tra[20]->Pvh_out, tra[21]->Pvh_out, tra[22]->Pvh_out, tra[23]->Pvh_out, tra[24]->Pvh_out, tra[25]->Pvh_out, tra[26]->Pvh_out, tra[27]->Pvh_out, tra[28]->Pvh_out, tra[29]->Pvh_out, tra[30]->Pvh_out);
					}
				}

				if (Time_diagram == 1) {
					//fprintf(fp5, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", minute, tra[0]->x, tra[1]->x, tra[2]->x, tra[3]->x, tra[4]->x, tra[5]->x, tra[6]->x, tra[7]->x, tra[8]->x, tra[9]->x, tra[10]->x, tra[11]->x, tra[12]->x, tra[13]->x, tra[14]->x, tra[15]->x, tra[16]->x, tra[17]->x, tra[18]->x, tra[19]->x, tra[20]->x, tra[21]->x, tra[22]->x, tra[23]->x, tra[24]->x, tra[25]->x, tra[26]->x, tra[27]->x, tra[28]->x, tra[29]->x, tra[30]->x, tra[31]->x, tra[32]->x, tra[33]->x, tra[34]->x, tra[35]->x, tra[36]->x, tra[37]->x, tra[38]->x, tra[39]->x, tra[40]->x, tra[41]->x, tra[42]->x, tra[43]->x, tra[44]->x, tra[45]->x, tra[46]->x);
				}
				else if (Time_diagram == 2) {
					//fprintf(fp5, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", minute, tra[0]->x, tra[1]->x, tra[2]->x, tra[3]->x, tra[4]->x, tra[5]->x, tra[6]->x, tra[7]->x, tra[8]->x, tra[9]->x, tra[10]->x, tra[11]->x, tra[12]->x, tra[13]->x, tra[14]->x, tra[15]->x, tra[16]->x, tra[17]->x, tra[18]->x, tra[19]->x, tra[20]->x, tra[21]->x, tra[22]->x, tra[23]->x, tra[24]->x, tra[25]->x, tra[26]->x, tra[27]->x, tra[28]->x, tra[29]->x, tra[30]->x, tra[31]->x, tra[32]->x, tra[33]->x, tra[34]->x, tra[35]->x, tra[36]->x, tra[37]->x, tra[38]->x, tra[39]->x, tra[40]->x, tra[41]->x, tra[42]->x);
				}
				else if (Time_diagram == 3) {
					//fprintf(fp5, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", minute, tra[0]->x, tra[1]->x, tra[2]->x, tra[3]->x, tra[4]->x, tra[5]->x, tra[6]->x, tra[7]->x, tra[8]->x, tra[9]->x, tra[10]->x, tra[11]->x, tra[12]->x, tra[13]->x, tra[14]->x, tra[15]->x, tra[16]->x, tra[17]->x, tra[18]->x, tra[19]->x, tra[20]->x, tra[21]->x, tra[22]->x, tra[23]->x, tra[24]->x, tra[25]->x, tra[26]->x, tra[27]->x, tra[28]->x, tra[29]->x, tra[30]->x);
				}
				else if (Time_diagram == 4) {
					//fprintf(fp5, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", minute, tra[0]->x, tra[1]->x, tra[2]->x, tra[3]->x, tra[4]->x, tra[5]->x, tra[6]->x, tra[7]->x, tra[8]->x, tra[9]->x, tra[10]->x, tra[11]->x, tra[12]->x, tra[13]->x, tra[14]->x, tra[15]->x, tra[16]->x, tra[17]->x, tra[18]->x, tra[19]->x, tra[20]->x, tra[21]->x, tra[22]->x, tra[23]->x);
				}
			}
			//////////////////////////////////////////////////////////////////////


			/////////////////////////////  cmd表示部  ////////////////////////////

			t_out2 = t;
			if (t_out2 - t_out2_old >= (1.00 - dt) * 1.0)
			{
				t_out2_old = t_out2;
				printf("State=%d || Position=%f || t=%f || Flag=%d ||   \r", tra[0]->accelflag, tra[0]->x, t, Flag_invest);
			}

			//////////////////////////////////////////////////////////////////////
		}
		/********** 数値計算ループ終了 **********/
		for (i = 0; i < NUM_tra; i++)
		{
			pow_energy[vel_state] += tra[i]->Evh_in_c;
			reg_energy[vel_state] += tra[i]->Evh_out_c;
			cos_energy[vel_state] += tra[i]->Esiv_coa_c;
		}
		for (i = 0; i < NUM_tra; i++) {
			ave_runtime[vel_state] += (tra[i]->run_time / NUM_tra);
		}
		for (i = 0; i < NUM_sub; i++)
		{
			sub_output[vel_state] += sub[i]->Ess_c;
		}

		fprintf(fp6, "%s, %.1f, %f, %f, %f, %f, %s, %f, %f, %f, %f, %f, %f, %f, %f\n", "", vel_dif[vel_state], ave_runtime[vel_state], pow_energy[vel_state], cos_energy[vel_state], reg_energy[vel_state], "", sub[0]->Ess_c, sub[1]->Ess_c, sub_output[vel_state], fed_loss[vel_state], fed_loss1[vel_state], fed_loss2[vel_state], fed_loss3[vel_state], fed_loss4[vel_state]);

		for (i = 0; i < NUM_tra; i++)
		{
			fprintf(fp2, "%s,%f,%f,%f,%f,%f,%f,%s,%f,%f,%f,%f,%f,%f,%f,%s,%f,%f,%s,%f,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", "", tra[i]->Emot_pow_c, tra[i]->Eres_pow_c, tra[i]->Eloss_inv_pow_c, tra[i]->Eloss_mot_pow_c, tra[i]->Eloss_fl_pow_c, tra[i]->Esiv_pow_c, "", -tra[i]->Emot_reg_c, -tra[i]->Emot_air_c, tra[i]->Eres_reg_c, tra[i]->Eloss_inv_reg_c, tra[i]->Eloss_mot_reg_c, tra[i]->Eloss_fl_reg_c, tra[i]->Esiv_reg_c, "", tra[i]->Eres_coa_c, tra[i]->Esiv_coa_c, "", tra[i]->Esiv_stp_c, "", tra[i]->Evh_in_c, tra[i]->Evh_out_c, tra[i]->Evh_st_c, tra[i]->Evh_c, tra[i]->Ereg_nashi_c, tra[i]->Evh_in1_c, tra[i]->Evh_in2_c, tra[i]->Evh_in3_c, tra[i]->Evh_in4_c, tra[i]->Evh_in5_c, tra[i]->Evh_in6_c, tra[i]->Evh_in7_c, tra[i]->Evh_in8_c, tra[i]->Evh_out1_c, tra[i]->Evh_out2_c, tra[i]->Evh_out3_c, tra[i]->Evh_out4_c, tra[i]->Evh_out5_c, tra[i]->Evh_out6_c, tra[i]->Evh_out7_c, tra[i]->Evh_out8_c, tra[i]->ifl_sum, tra[i]->run_time);
		}

		fprintf(fp4, "%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", "", sub[0]->Ess_c, sub[0]->Wss_c, sub[1]->Ess_c, sub[1]->Wss_c, sub[0]->Ess_c + sub[1]->Ess_c, Lcir1_c, Lcir2_c, Lcir3_c, Lcir4_c, Lcir_c);

		for (i = 0; i < NUM_tra; i++)
		{
			fclose(fp1[i]);
		}


		fclose(fp2);
		fclose(fp3);
		fclose(fp4);
		fclose(fp5);
		//fclose(fpx);
		//fclose(fpy);
		//fclose(fpz);

		/*****  メモリ開放  *****/
		for (i = 0; i < NUM_tra; i++)
		{
			free(tra[i]);
		}
		for (i = 0; i < NUM_sub; i++)
		{
			free(sub[i]);
		}

		if (vel_state == (vel_pattern - 1)) {
			lsm(ave_runtime, pow_energy, vel_pattern, &a_pow[vss_state2 + vss_state1 * vss_pattern], &b_pow[vss_state2 + vss_state1 * vss_pattern]);
			lsm(ave_runtime, reg_energy, vel_pattern, &a_reg[vss_state2 + vss_state1 * vss_pattern], &b_reg[vss_state2 + vss_state1 * vss_pattern]);
			lsm(ave_runtime, cos_energy, vel_pattern, &a_cos[vss_state2 + vss_state1 * vss_pattern], &b_cos[vss_state2 + vss_state1 * vss_pattern]);
			lsm(ave_runtime, sub_output, vel_pattern, &a_sub[vss_state2 + vss_state1 * vss_pattern], &b_sub[vss_state2 + vss_state1 * vss_pattern]);
			lsm(ave_runtime, fed_loss, vel_pattern, &a_fed[vss_state2 + vss_state1 * vss_pattern], &b_fed[vss_state2 + vss_state1 * vss_pattern]);
			lsm(ave_runtime, fed_loss1, vel_pattern, &a_fed1[vss_state2 + vss_state1 * vss_pattern], &b_fed1[vss_state2 + vss_state1 * vss_pattern]);
			lsm(ave_runtime, fed_loss2, vel_pattern, &a_fed2[vss_state2 + vss_state1 * vss_pattern], &b_fed2[vss_state2 + vss_state1 * vss_pattern]);
			lsm(ave_runtime, fed_loss3, vel_pattern, &a_fed3[vss_state2 + vss_state1 * vss_pattern], &b_fed3[vss_state2 + vss_state1 * vss_pattern]);
			lsm(ave_runtime, fed_loss4, vel_pattern, &a_fed4[vss_state2 + vss_state1 * vss_pattern], &b_fed4[vss_state2 + vss_state1 * vss_pattern]);

		}
		if (vss_state1 == (vss_pattern - 1) && vss_state2 == (vss_pattern - 1) /* && vss_state3 == (vss_pattern - 1) && vss_state4 == (vss_pattern - 1) && vss_state5 == (vss_pattern - 1) && vss_state6 == (vss_pattern - 1) && vss_state7 == (vss_pattern - 1) */ && vel_state == (vel_pattern - 1)) { sim_endflag++; break; }

		if (Consider_loop_flag == 1) break;

		vel_state++;

		if (vel_state == vel_pattern) {


			if (ALL_Vss_Change == 0) {




				if (vss_state2 == (vss_pattern - 1)) {
					vss_state2 = 0;
					vss_state1++;
				}
				else {
					vss_state2++;
				}

			}
			else {
				vss_state1++;
				vss_state2++;

			}

			vel_state = 0;

		}
		loop_endflag++;

		t_out_old = 0.0;

		sim_loopcount++;
	}

	fprintf(fp6, "\n");

	for (i = 0; i < vss_pattern * vss_pattern; i++) {
		fprintf(fp6, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", a_pow[i], b_pow[i], a_reg[i], b_reg[i], a_cos[i], b_cos[i], a_sub[i], b_sub[i], a_fed[i], b_fed[i], a_fed1[i], b_fed1[i], a_fed2[i], b_fed2[i], a_fed3[i], b_fed3[i], a_fed4[i], b_fed4[i]);
	}

	fprintf(fp6, "\n力行電力量[kWh], き電損失[kWh], 惰行停止電力量[kWh], 消費電力量(力行+き電)[kWh], 全消費電力量[kWh], 回生電力量[kWh], 変電所出力電力量[kWh], 変電所出力電力量(1列車分),き電損失1, き電損失2, き電損失3,き電損失4\n");

	for (i = 0; i < vss_pattern * vss_pattern; i++) {

		cal_pow[i] = b_pow[i] * 41.7 + a_pow[i];
		cal_reg[i] = b_reg[i] * 41.7 + a_reg[i];
		cal_st[i] = b_cos[i] * 41.7 + a_cos[i];
		cal_sub[i] = b_sub[i] * 41.7 + a_sub[i];
		cal_fed[i] = b_fed[i] * 41.7 + a_fed[i];
		cal_fed1[i] = b_fed1[i] * 41.7 + a_fed1[i];
		cal_fed2[i] = b_fed2[i] * 41.7 + a_fed2[i];
		cal_fed3[i] = b_fed3[i] * 41.7 + a_fed3[i];
		cal_fed4[i] = b_fed4[i] * 41.7 + a_fed4[i];

		fprintf(fp6, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", cal_pow[i], cal_fed[i], cal_st[i], cal_pow[i] + cal_fed[i], cal_pow[i] + cal_fed[i] + cal_st[i], cal_reg[i], cal_sub[i], cal_pow[i] + cal_fed[i] + cal_st[i] + cal_reg[i], cal_fed1[i], cal_fed2[i], cal_fed3[i], cal_fed4[i]);
	}

	fclose(fp6);

	end_flag = 1;
	exit(0);
}
