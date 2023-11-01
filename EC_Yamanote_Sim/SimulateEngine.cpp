#include "SimulateEngine.hpp"

SimulateEngine::SimulateEngine()
{
	Init();
}

bool SimulateEngine::Init()
{
	/***  各行列宣言部  ***/
	/***  各行列宣言部  ***/
	H = std::make_unique<double[]>(N_node * N_branch);		//Connection Matrix
	H_tra = std::make_unique<double[]>(N_branch * N_node);	//Transpose matrix of H
	Y = std::make_unique<double[]>(N_branch * N_branch);		//Conductance Matrix
	A = std::make_unique<double[]>(N_node * N_node);			//A matrix(H*Y*H_tra)
	TEMP = std::make_unique<double[]>(N_node * N_branch);	//Test
	In = std::make_unique<double[]>(N_node);					//Node Current Matrix
	In_cpy = std::make_unique<double[]>(N_node);				//Copy of Node Current Matrix
	Vn = std::make_unique<double[]>(N_node);					//Node Voltage Matrix
	Ib = std::make_unique<double[]>(N_branch);				//Branch Current Matrix
	Vb = std::make_unique<double[]>(N_branch);
	ipiv_Vn = std::make_unique<int[]>(N_node);
	//i = 0;
	//j = 0;

	ERROR = 0;						//回路計算の再計算判定用フラグ（flagを全て加算してNUM_sub以上なら回路計算ループを抜ける。NUM_sub以下なら再計算。）
	ERROR_REG = 0;
	count_loop = 0;
	count_loop_rec = 0;
	//	test_flag = 0;

		/***  マトリックス計算関係  ***/
	K_node = N_node;
	K_branch = N_branch;
	ld_node = N_node;
	ld_branch = N_branch;
	alpha = 1.0;
	beta = 0.0;
	//int* ipiv_Vn;

	info_Vn = 0;
	nrhs = 1;
	incx = 1;
	incy = 1;
	/***  構造体宣言部  ***/
	tra = std::make_unique<TRAIN[]>(NUM_tra);	//ポインタの配列として定義
	sub = std::make_unique<SUB[]>(NUM_sub);
	node = std::make_unique<NODE[]>(N_node);
	node_order = std::make_unique<NODE[]>(N_node);
	branch = std::make_unique<BRANCH[]>(N_branch);

	t = 0.0;
	minute = 0.0;

	sub_total = 0.0;

	return true;
}

void SimulateEngine::Make_train(TRAIN& tra, const INI_TRA& ini, int i)
{
	int m;
	int n;

	//*tra = (TRAIN*)malloc(sizeof(TRAIN));
	//if (*tra == NULL)
	//{
	//	//エラー処理
	//}

	tra.name_train = i;
	tra.Type = ini.Type;
	tra.notch = 3.0;
	tra.brake = 5.0;

	tra.mass_M = ini.mass_M * ini.num_M;			//[kg]自分用(100%乗車)
	tra.mass_T = ini.mass_T * ini.num_T;			//[kg]自分用(100%乗車)
	tra.mass_P = (ini.num_M * 152.0 + (ini.num_T - 2.0) * 152.0 + 2.0 * 144.0) * 55.0 / 2.0;

	if (tra.Type == 2) tra.mass_P = 0.0;

	tra.mass = tra.mass_M + tra.mass_T + tra.mass_P + 15.8 * 1000.0 * (ini.num_M + ini.num_T) / 10.0;				//[kg](250%乗車)
	tra.mass_c = tra.mass / 1000.0;										//[t](250%乗車)

	tra.num = ini.num_M + ini.num_T;
	tra.num_mot = ini.num_M * 4.0;

	tra.v0_P = 38.0 / 3.6;        //[m/s]
	tra.v1_P = 44.0 / 3.6;		//[m/s]
	tra.v2_P = 46.0 / 3.6;		//[m/s]
	tra.v0_R = 52.0 / 3.6;        //[m/s]
	tra.v1_R = 86.0 / 3.6;		//[m/s]
	tra.v2_R = 85.0 / 3.6;		//[m/s]

	if (tra.Type == 2) {
		tra.v0_P = 40.0 / 3.6;        //[m/s]
		tra.v1_P = 60.0 / 3.6;		//[m/s]
		tra.v2_P = 60.0 / 3.6;		//[m/s]
		tra.v0_R = 52.0 / 3.6;        //[m/s]
		tra.v1_R = 96.0 / 3.6;		//[m/s]
		tra.v2_R = 96.0 / 3.6;		//[m/s]
	}

	tra.Fmax = 18.818 * 1000.0;         //[N/MM]  自分用(100%乗車時)(参考：起動加速度0.639[m/s/s])
	tra.Bmax = 15.556 * 1000.0;		//[N/MM] (100%乗車時)
	tra.Bref = 1.111 * 5.0 / 7.0;
	tra.BEmax = 0.0;

	tra.i1d_P = 90.42780281;  //力行時，定トルク領域でのd軸電流
	tra.i1q_P = 237.1149987; //力行時，定トルク領域でのq軸電流
	tra.i1d_R = 91.09231795;  //回生時，定トルク領域でのd軸電流
	tra.i1q_R = -211.1226486; //回生時，定トルク領域でのq軸電流

	if (tra.Type == 2) {
		tra.notch = 4.0;
		tra.Fmax = 14.776 * 1000.0;         //[N/MM]  自分用(100%乗車時)(参考：起動加速度0.639[m/s/s])
		tra.Bmax = 13.910 * 1000.0;		//[N/MM] (100%乗車時)
		tra.Bref = 1.111 * 5.0 / 7.0;
		tra.BEmax = 0.0;

		tra.i1d_P = 91.877;  //力行時，定トルク領域でのd軸電流
		tra.i1q_P = 183.246; //力行時，定トルク領域でのq軸電流
		tra.i1d_R = 91.60374;  //回生時，定トルク領域でのd軸電流
		tra.i1q_R = -173.0142; //回生時，定トルク領域でのq軸電流
	}

	tra.B = 1.111 * 5.0 / 7.0;

	tra.direction = ini.direction;
	tra.UP_DOWN = (double)tra.direction;
	tra.flag_station = 0;
	tra.startflag = 0;
	tra.lap = 0;
	tra.lapmax = ini.lapmax;

	tra.t_count = 0;
	tra.c_count = 0;

	tra.wr = 0.0;
	tra.ws = 0.0;
	tra.we = 0.0;

	tra.con_speed = 0.0;

	tra.Gr = 6.31;  //ギア比
	tra.Ts = 0.00001;  //収束刻み
	tra.Tf = 0.0001;  //収束定数
	tra.POLE = 2.0;  //極対数
	tra.rd = 0.41;   //車輪半径[m]
	tra.R1 = 0.0970;
	tra.R2 = 0.07327;  //2次抵抗[Ω]
	tra.L1 = 0.030549;
	tra.L2 = 0.030549;  //2次自己インダクタンス[H]
	tra.M_IM = 0.029514;
	tra.g0 = 0.001614;

	tra.Cfc = C_FC * ini.num_M;
	tra.Lfl = L_FL / ini.num_M;
	tra.Rfl = R_FL / ini.num_M;

	tra.theta = 1.0;
	tra.theta_old = 1.0;


	tra.Emot_pow = 0.0;			//モータ力行エネルギー[J]
	tra.Emot_reg = 0.0;			//モータ回生エネルギー[J]
	tra.Emot_air = 0.0;			//機械ブレーキエネルギー[J]

	tra.Emot_pow_c = 0.0;		//モータ力行エネルギー[kWh]
	tra.Emot_reg_c = 0.0;		//モータ回生エネルギー[kWh]
	tra.Emot_air_c = 0.0;		//機械ブレーキエネルギー[kWh]

	tra.Eloss_inv_pow = 0.0;
	tra.Eloss_mot_pow = 0.0;
	tra.Eloss_fl_pow = 0.0;
	tra.Eloss_all_pow = 0.0;

	tra.Eloss_inv_reg = 0.0;
	tra.Eloss_mot_reg = 0.0;
	tra.Eloss_fl_reg = 0.0;
	tra.Eloss_all_reg = 0.0;

	tra.Eres_pow = 0.0;
	tra.Eres_reg = 0.0;
	tra.Eres_coa = 0.0;

	tra.Esiv_pow = 0.0;
	tra.Esiv_reg = 0.0;
	tra.Esiv_coa = 0.0;
	tra.Esiv_stp = 0.0;

	tra.Eloss_inv_pow_c = 0.0;
	tra.Eloss_mot_pow_c = 0.0;
	tra.Eloss_fl_pow_c = 0.0;
	tra.Eloss_all_pow_c = 0.0;

	tra.Eloss_inv_reg_c = 0.0;
	tra.Eloss_mot_reg_c = 0.0;
	tra.Eloss_fl_reg_c = 0.0;
	tra.Eloss_all_reg_c = 0.0;

	tra.Eres_pow_c = 0.0;
	tra.Eres_reg_c = 0.0;
	tra.Eres_coa_c = 0.0;

	tra.Esiv_pow_c = 0.0;
	tra.Esiv_reg_c = 0.0;
	tra.Esiv_coa_c = 0.0;
	tra.Esiv_stp_c = 0.0;



	tra.Evh = 0.0;				//列車入出力エネルギー[J]
	tra.Evh_in = 0.0;			//列車入力エネルギー[J]
	tra.Evh_out = 0.0;			//列車出力エネルギー[J]
	tra.Evh_st = 0.0;

	tra.Evh_c = 0.0;			//列車入出力エネルギー[kWh]
	tra.Evh_in_c = 0.0;			//列車入力エネルギー[kWh]
	tra.Evh_out_c = 0.0;		//列車出力エネルギー[kWh]
	tra.Evh_st_c = 0.0;

	tra.accelflag = 4;
	tra.brakeflag = 0;
	tra.laststopflag = 0;
	tra.X_recentstop = 0;
	tra.X_nextstop = 0.0;
	tra.X_brake = 0.0;
	tra.x = ini.x;
	tra.v = ini.v;
	tra.v_new = 0.0;
	tra.a = 0.0;
	tra.Fmot = 0.0;
	tra.Ftot = 0.0;

	tra.id_mot_P = 0.0;
	tra.id_mot_R = 0.0;
	tra.iq_mot_P = 0.0;
	tra.iq_mot_R = 0.0;

	tra.v1d = 0.0;
	tra.v1q = 0.0;
	tra.Em = 0.0;
	tra.idtot = 0.0;
	tra.iqtot = 0.0;
	tra.idtotF = 0.0;
	tra.iqtotF = 0.0;
	tra.iq_reg = 0.0;
	tra.iq_regF = 0.0;

	tra.x_c = ini.x / 1000.0;
	tra.v_c = ini.v / 1000.0;

	tra.Pmot = 0.0;
	tra.Ploss_mot = 0.0;
	tra.Ploss_fe = 0.0;
	tra.Ploss_inv = 0.0;
	tra.Ploss_fl = 0.0;
	tra.iinv = 0.0;
	tra.ifc = 0.0;
	tra.ifl = 0.0;
	tra.vfc_old = Vss1;
	tra.isiv = 0.0;
	tra.vfc = Vss1;
	tra.vp = Vss1;
	tra.vp_old = Vss1;
	tra.R_grade = 0.0;
	tra.R_curve = 0.0;
	tra.R_run = 0.0;
	tra.R_total = 0.0;
	tra.t_stop = 0.0;
	tra.t_stop_old = 0.0;
	tra.t_wait = 0.0;
	tra.diaflag = 0;

	tra.Speedlimit = 0.0;
	tra.Reaccelspeed = 0.0;
	tra.T_delay = ini.T_delay;
	tra.Vst = Vss1;

	for (m = 0; m < rownum; m++)
	{
		for (n = 0; n < 2; n++)
		{
			tra.route[m][n] = 0;
		}
	}
}

void SimulateEngine::Make_substation(SUB& sub, const INI_SUB& ini, int i)
{

	//*sub = (SUB*)malloc(sizeof(SUB));
	//if (*sub == NULL)
	//{
	//	//エラー処理
	//}

	sub.name_SS = i;
	sub.Vss_0 = Vss1;
	sub.Iss_0 = Iss_rated1;
	sub.e_ss = ini.e_ss;
	sub.Xss = ini.Xss;

	sub.diode = ini.diode;
	sub.vss = ini.vss;
	sub.vss_e = ini.vss;
	sub.iss = ini.iss;
	sub.vout = ini.vss;
	sub.flag = 1;

	sub.Rss = (ini.e_ss * ini.vss) / (100.0 * ini.Iss);
	sub.Rss_ini = (ini.e_ss * ini.vss) / (100.0 * ini.Iss);
	sub.Jss_ini = -ini.vss / ((ini.e_ss * ini.vss) / (100.0 * ini.Iss));	//Jssの初期値
	sub.Jss = -ini.vss / ((ini.e_ss * ini.vss) / (100.0 * ini.Iss));		//更新用のJss

	sub.Ess = 0.0;			//変電所出力エネルギー[J]
	sub.Ess_c = 0.0;		//変電所出力エネルギー[kWh]
	sub.Wss = 0.0;
	sub.Wss_c = 0.0;
}

void SimulateEngine::Make_NODE_TRAIN(NODE& node, TRAIN& tra)
{
	//*node = (NODE*)malloc(sizeof(NODE));
	//if (*node == NULL)
	//{
	//	printf("MEMORY?? What is that????\n");
	//	exit(EXIT_FAILURE);//エラー処理
	//}

	node.Number = tra.name_train;
	node.X = tra.x;
	node.V = tra.vp;
	node.I = tra.ifl;
	node.r = 0.0;
	node.flag = tra.direction;
}

void SimulateEngine::Make_NODE_SS(NODE& node, SUB& sub)
{
	//*node = (NODE*)malloc(sizeof(NODE));
	//if (*node == NULL)
	//{
	//	printf("MEMORY?? What is that????\n");
	//	exit(EXIT_FAILURE);//エラー処理
	//}

	node.Number = sub.name_SS;
	node.X = sub.Xss;
	node.V = sub.vss;
	node.I = sub.Jss;
	node.r = sub.Rss;
	node.flag = 0;
}

void SimulateEngine::Error_Detection_SS(SUB& sub)
{
	if (sub.vout > sub.vss_e + 0.1) {
		sub.Jss = sub.iss - sub.vout / sub.Rss;
	}
	else if (sub.vout < sub.vss_e - 0.1) sub.Jss = sub.Jss_ini;

	if (sub.iss < -0.1)		//変電所に電流が流入している⇒flag=0（エラー）
	{
		sub.flag = 0;
		sub.Jss = -sub.vout / sub.Rss;
	}
	else
	{
		sub.flag = 1;			//変電所に電流が流入していない⇒flag=1（OK）
	}

	/*	if (sub.vout >(Vcmax))	//電圧の制限
		{
			sub.flag = 0;
			sub.Jss = -(Vcmax) / sub.Rss;
		}*/
}

void SimulateEngine::Error_Detection_REG(TRAIN& tra) {

	/********** 軽負荷回生制御パターンを満たすかどうかの判定 **********/
	if (tra.vfc > Vclim && tra.vfc < Vcmax)						//回生絞り込み中
	{
		tra.b_reg_judge = tra.Fmot * (Vcmax - tra.vfc) / (Vcmax - Vclim);
	}
	else if (tra.vfc >= Vcmax)										//回生絞り込み終了
	{
		tra.b_reg_judge = 0.0;
	}
	else															//回生絞り込みなし
	{
		tra.b_reg_judge = tra.Fmot;
	}

	tra.judgement = std::abs(tra.b_reg - tra.b_reg_judge);

	if (tra.judgement > 1.0) {
		tra.reg_flag = 0;
	}
	else {
		tra.reg_flag = 1;
	}
}

void SimulateEngine::Calculate_R_run(TRAIN& tra)		//走行抵抗の算出
{
	tra.R_run = (tra.mass_M / 1000.0 * (1.65 + 0.0247 * tra.v_c) + tra.mass_T / 1000.0 * (0.78 + 0.0028 * tra.v_c) + (0.028 + 0.0078 * (tra.num - 1)) * tra.v_c * tra.v_c) * g;		//rr=g×W×(a+bv+cv^2/W)
	//	tra.R_run = g*(1.32+0.0614 * tra.v_c)*tra.mass_c+(0.028 + 0.0078 * (tra.num - 1)) * tra.v_c * tra.v_c;
}

void SimulateEngine::Calculate_R_total(TRAIN& tra)		//合計抵抗の算出
{
	Calculate_R_run(tra);

	tra.R_total = tra.R_run;
}

void SimulateEngine::Traction_force(TRAIN& tra)		//力行時引張力の算出
{

	double v1_pp = tra.v1_P * tra.vfc / 1350.0;
	double v2_pp = tra.v2_P * tra.vfc / 1350.0;

	if (tra.Type == 2) tra.Fmax = 14.776 * 1000.0 * (tra.notch / 4.0);

	if (tra.v < v1_pp)		//定トルク領域
	{
		tra.Fmot = tra.Fmax * tra.num_mot;
	}
	else if (tra.v < v2_pp)		//定電力領域
	{
		tra.Fmot = tra.Fmax * tra.num_mot * v1_pp / tra.v;
	}
	else     //特性領域
	{
		tra.Fmot = tra.Fmax * tra.num_mot * v1_pp * v2_pp / (tra.v * tra.v);	//特性領域
	}

	if (tra.Fmot > tra.Fmax * tra.num_mot) tra.Fmot = tra.Fmax * tra.num_mot;

	tra.Ftot = tra.Fmot;
}

void SimulateEngine::Regenerative_force(TRAIN& tra)		//回生時引張力の算出
{
	double v1_rr = tra.v1_R * tra.vfc / 1650.0;
	double v2_rr = tra.v2_R * tra.vfc / 1650.0;

	if (tra.Type == 2) {
		tra.Bmax = 13.910 * 1000.0 * tra.brake / 5.0;
		if (tra.brake > 5.0) {
			tra.Bmax = 15.092 * 1000.0;
		}
	}

	if (tra.v < v1_rr)		//定トルク領域
	{
		tra.Fmot_regM = -tra.Bmax * tra.num_mot;
	}
	else if (tra.v < v2_rr)		//定電力領域
	{
		tra.Fmot_regM = -tra.Bmax * tra.num_mot * v1_rr / tra.v;
	}
	else	//特性領域
	{
		tra.Fmot_regM = -tra.Bmax * tra.num_mot * v1_rr * v2_rr / (tra.v * tra.v);
	}

	tra.Fmot = tra.Fmot_regM;
	if (tra.Fmot < -tra.BEmax) tra.Fmot = -tra.BEmax;
	tra.Ftot = -tra.BEmax;
}

/*** 一次遅れフィルタ ***/
double SimulateEngine::funcdelay(double aTs, double aTf, double ayold, double au)
{
	double fudely;
	fudely = (aTs * au + aTf * ayold) / (aTf + aTs);
	return(fudely);
}

void SimulateEngine::d_current_power(TRAIN& tra)//力行時のd軸電流算出
{
	double v0_PP = tra.v0_P * tra.vfc / 1350.0;
	tra.id_mot_P = tra.i1d_P;

	if (tra.v > v0_PP)		//弱め磁束領域
	{
		tra.id_mot_P = tra.i1d_P * v0_PP / tra.v;
	}

	tra.idtot = tra.id_mot_P;
}

void SimulateEngine::d_current_regenerative(TRAIN& tra)//回生時のd軸電流算出
{
	double v0_RR = tra.v0_R * tra.vfc / 1650.0;
	tra.id_mot_R = tra.i1d_R;

	if (tra.v > v0_RR)		//弱め磁束領域
	{
		tra.id_mot_R = tra.i1d_R * v0_RR / tra.v;
	}
	tra.idtot = tra.id_mot_R;
}

void SimulateEngine::q_current_power(TRAIN& tra)//力行時のq軸電流算出
{
	tra.iq_mot_P = (tra.Fmot * tra.rd / tra.Gr) * tra.L2 / tra.POLE / (tra.M_IM * tra.M_IM) / tra.idtot / tra.num_mot;
	tra.iqtot = tra.iq_mot_P;
	tra.iqtotF = funcdelay(tra.Ts, tra.Tf, tra.iqtotF, tra.iqtot);
}

void SimulateEngine::q_current_regenerative(TRAIN& tra)//回生時のq軸電流算出
{
	tra.iq_mot_R = (tra.Fmot * tra.rd / tra.Gr) * tra.L2 / tra.POLE / (tra.M_IM * tra.M_IM) / tra.idtot / tra.num_mot;
	tra.iqtot = tra.iq_mot_R;
	tra.iqtotF = funcdelay(tra.Ts, tra.Tf, tra.iqtotF, tra.iqtot);
}

void SimulateEngine::d_q_voltage_calculation(TRAIN& tra) {
	if (tra.accelflag == 2) {
		//角周波数計算
		tra.ws = 0.0;
		tra.we = tra.ws + tra.wr;

		/*d軸電圧*/
		tra.v1d = 0.0;

		/*q軸電圧*/
		tra.v1q = 0.0;
		tra.Em = 0.0;
	}
	else if (tra.accelflag == 4) {

		tra.ws = 0.0;
		tra.we = 0.0;
		/*d軸電圧*/
		tra.v1d = 0.0;

		/*q軸電圧*/
		tra.v1q = 0.0;
		tra.Em = 0.0;
	}
	else {
		tra.ws = tra.R2 / tra.L2 * tra.iqtot / tra.idtot;
		tra.we = tra.ws + tra.wr;

		/*d軸電圧*/
		tra.v1d = tra.R1 * tra.idtot - (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.iqtotF;
		//		tra.v1d = - (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.iqtot;

				/*q軸電圧*/
		tra.v1q = tra.R1 * tra.iqtotF + (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.idtot + (tra.M_IM / tra.L2) * tra.M_IM * tra.we * tra.idtot;
		//		tra.v1q = (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.idtot + (tra.M_IM / tra.L2) * tra.M_IM * tra.we * tra.idtot;

				/*インバータ入力電圧*/
		tra.Em = sqrt(tra.v1d * tra.v1d + tra.v1q * tra.v1q);
	}
}

void SimulateEngine::Solve_next_state(TRAIN& tra)
{
	tra.a = (tra.Ftot - tra.R_total) / (tra.mass);		//運動方程式ma=F
	tra.a_c = tra.a * 3600.0 / 1000.0;

	tra.v_new = tra.v + tra.a * dt;

	tra.v = tra.v_new;
	tra.v_c = tra.v * 3.6;

	tra.wr = (tra.v / tra.rd) * tra.Gr * tra.POLE;

	if (tra.direction == 1)
	{
		tra.x += tra.v * dt;
	}
	else
	{
		tra.x -= tra.v * dt;
	}

	tra.x_c = tra.x / 1000.0;
}

void SimulateEngine::Calculate_BEmax(TRAIN& tra)
{
	if (tra.Type == 2) tra.B = 1.111 * tra.brake / 7.0;
	tra.BEmax = tra.B * tra.mass;
}

void SimulateEngine::change_direction(TRAIN& tra)
{
	tra.t_count = 0;
	tra.c_count = 0;

	switch (tra.direction) {
	case 1:
		tra.direction = -1;
		tra.UP_DOWN = -1;
		break;
	case -1:
		tra.direction = 1;
		tra.UP_DOWN = 1;
		break;
	}
}

void SimulateEngine::decide_final_station(TRAIN& tra, int i)
{
	if (tra.Type == 0 || tra.Type == 9 || tra.Type == 10 || tra.Type == 11 || tra.Type == 12 || tra.Type == 13 || tra.Type == 14) {
		if (tra.direction == 1.0)			//池袋⇒渋谷
		{
			if (i == NUM_station)
			{
				if (tra.lap != tra.lapmax) {
					change_direction(tra);
					if (tra.name_train == 7) tra.T_delay = 27.0 * 60.0;
					if (tra.name_train == 8) tra.T_delay = 37.0 * 60.0;
					tra.lap++;
				}
				else {
					tra.laststopflag = 1;
					tra.T_delay = 9999.0;
				}
			}
		}
		else
		{
			if (i == 0)						//渋谷⇒池袋
			{
				if (tra.lap != tra.lapmax) {
					change_direction(tra);
					if (tra.name_train == 3) tra.T_delay = 27.0 * 60.0 + 50.0;
					if (tra.name_train == 4) tra.T_delay = 37.0 * 60.0 + 30.0;
					tra.lap++;
				}
				else {
					tra.laststopflag = 1;
					tra.T_delay = 9999.0;
				}
			}
		}
	}
	else if (tra.Type == 1) {
		if (tra.direction == 1.0)			//大宮⇒吹上方面
		{
			if (i == NUM_station2)
			{
				if (tra.lap != tra.lapmax) {
					change_direction(tra);
					tra.lap++;
				}
				else {
					tra.laststopflag = 1;
					tra.T_delay = 9999.0;
				}
			}
		}
		else
		{
			if (i == 0)						//吹上⇒大宮方面
			{
				if (tra.lap != tra.lapmax) {
					change_direction(tra);
					if (tra.name_train == 5) tra.T_delay = 41.0 * 60.0 + 50.0;
					tra.lap++;
				}
				else {
					tra.laststopflag = 1;
					tra.T_delay = 9999.0;
				}
			}
		}
	}
}

//ここを変えたいby Kasai1012
void SimulateEngine::Run_pattern(TRAIN& tra, const STATION* sta, double t, std::vector<double> wait_time)
{
	int i = 0;

	tra.t_stop = t;

	/*** 次の駅を探索 ***/
	tra.X_nextstop_old = tra.X_nextstop;
	if (tra.Type == 0 || tra.Type == 9 || tra.Type == 10 || tra.Type == 11 || tra.Type == 12 || tra.Type == 13 || tra.Type == 14) {
		for (size_t i = 0; i < NUM_station + NUM_final_station; i++)
		{
			if (tra.direction == 1.0)			//池袋⇒渋谷
			{
				if (sta[i].Xs <= tra.x && sta[i + 1].Xs >= tra.x)
				{
					tra.X_nextstop = sta[i + 1].Xs;
					break;
				}
			}
			else								//渋谷⇒池袋
			{
				if (sta[i].Xs <= tra.x && sta[i + 1].Xs >= tra.x)
				{
					tra.X_nextstop = sta[i].Xs;
					break;
				}
			}
		}

		decide_final_station(tra, i);

	}

	/*速度制限設定*/
	if (sta[1].Xs < tra.x && tra.x < sta[2].Xs)
	{
		tra.Speedlimit = 65.0;
	}
	else if (sta[2].Xs < tra.x && tra.x < sta[3].Xs)
	{
		tra.Speedlimit = 65.0;
	}
	else if (sta[3].Xs < tra.x && tra.x < sta[4].Xs)
	{
		tra.Speedlimit = 65.0;
	}
	else if (sta[4].Xs < tra.x && tra.x < sta[5].Xs)
	{
		tra.Speedlimit = 65.0;
	}
	else if (sta[5].Xs < tra.x && tra.x < sta[6].Xs)
	{
		tra.Speedlimit = 65.0;
	}
	else if (sta[6].Xs < tra.x && tra.x < sta[7].Xs)
	{
		tra.Speedlimit = 65.0;
	}
	else if (sta[7].Xs < tra.x && tra.x < sta[8].Xs)
	{
		tra.Speedlimit = 65.0;
	}
	else
	{
		tra.Speedlimit = 65.0;
	}


	/**** 走行抵抗計算 ****/
	Calculate_R_total(tra);

	tra.Bref = 1.111 * tra.brake / 7.0;



	/*** 走行モードを決定 ***/
	if (tra.accelflag == 4)		/*停車中の動作を決定*/
	{
		if (tra.laststopflag == 1)
		{
			tra.accelflag = 4;
		}
		else if (t < tra.T_delay)
		{
			tra.accelflag = 4;
		}
		else if (tra.direction == 1) {	//内回り（池袋→渋谷）
			if (sta[2].Xs - 200.0 <= tra.x && tra.x < sta[2].Xs + 200.0) {	//目白
				// if (tra.t_stop >= tra.T_delay + 2.0 * 60.0 && tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;　　発車時刻を設定する場合
				// if (tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
				if (tra.Type == 0) {
					if (tra.t_stop - tra.t_stop_old >= (stoptime + wait_time[0])) tra.accelflag = 1;
				}

			}
			else if (sta[3].Xs - 200.0 <= tra.x && tra.x < sta[3].Xs + 200.0) {	//高田馬場
				// if (tra.t_stop >= tra.T_delay + 4.0 * 60.0 + 20 && tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
				// if (tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
				if (tra.Type == 0) {
					if (tra.t_stop - tra.t_stop_old >= (stoptime - wait_time[0] + wait_time[1])) tra.accelflag = 1;
				}

			}
			else if (sta[4].Xs - 200.0 <= tra.x && tra.x < sta[4].Xs + 200.0) {	//新大久保
				// if (tra.t_stop >= tra.T_delay + 6.0 * 60.0 && tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
				// if (tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
				if (tra.Type == 0) {
					if (tra.t_stop - tra.t_stop_old >= (stoptime - wait_time[1] + wait_time[2])) tra.accelflag = 1;
				}

			}
			else if (sta[5].Xs - 200.0 <= tra.x && tra.x < sta[5].Xs + 200.0) {	//新宿
				// if (tra.t_stop >= tra.T_delay + 9.0 * 60.0 && tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
				if (tra.Type == 0) {
					if (tra.t_stop >= tra.T_delay + 9.0 * 60.0 && tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
				}

			}
			else if (sta[6].Xs - 200.0 <= tra.x && tra.x < sta[6].Xs + 200.0) {	//代々木
				// if (tra.t_stop >= tra.T_delay + 10.0 * 60.0 && tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
				// if (tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
				if (tra.Type == 0) {
					if (tra.t_stop - tra.t_stop_old >= (stoptime + wait_time[3])) tra.accelflag = 1;
				}

			}
			else if (sta[7].Xs - 200.0 <= tra.x && tra.x < sta[7].Xs + 200.0) {	//原宿
				// if (tra.t_stop >= tra.T_delay + 13.0 * 60.0 && tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1; 
				// if (tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
				if (tra.Type == 0) {
					if (tra.t_stop - tra.t_stop_old >= (stoptime - wait_time[3] + wait_time[4])) tra.accelflag = 1;
				}

			}
			else if (sta[8].Xs - 200.0 <= tra.x && tra.x < sta[8].Xs + 200.0) {	//渋谷
				if ((tra.t_stop - tra.t_stop_old) >= 7200.0) tra.accelflag = 1;
			}
			else {
				if ((tra.t_stop - tra.t_stop_old) >= 20.0) tra.accelflag = 1;
			}
		}


		else
		{
			tra.accelflag = 4;
		}

	}
	else if (tra.accelflag == 2)		/*惰行中の動作を決定*/
	{
		if (tra.direction == 1.0)		//大宮⇒吹上方面
		{
			if (tra.x + tra.v * tra.v / 2.0 / tra.Bref < tra.X_nextstop)
			{
				tra.accelflag = 2;
			}
			/*	else if(tra.v_c < 50.0)
				{
					tra.accelflag = 1;
				}
				*/
			else
			{
				tra.accelflag = 3;
				tra.brakeflag = 1;
			}
		}
		else							//大宮⇒吹上方面
		{
			if (tra.x - tra.v * tra.v / 2.0 / tra.Bref > tra.X_nextstop)
			{
				tra.accelflag = 2;
			}

			/*else if(tra.v_c < 60.0)
			{
				tra.accelflag = 1;
			}
			*/
			else
			{
				tra.accelflag = 3;
				tra.brakeflag = 1;
			}
		}
	}
	else		/*加速中or減速中の動作を決定*/
	{
		if (tra.direction == 1.0)								//大宮⇒吹上方面
		{
			if (tra.brakeflag == 1)		/*減速中の動作*/
			{
				if (tra.v_c >= 0.0)
				{
					tra.accelflag = 3;		//減速区間
				}
				else
				{
					tra.accelflag = 4;		//停車
					tra.brakeflag = 0;
					tra.t_stop_old = tra.t_stop;
				}
			}
			else
			{
				if (tra.laststopflag == 1)	//駅の手前に停車した場合の救済
				{
					tra.accelflag = 3;
					tra.brakeflag = 1;
				}
				else if (tra.x + tra.v * tra.v / 2.0 / tra.Bref > tra.X_nextstop && tra.v_c > speed_min - 20.0)
				{
					tra.accelflag = 3;
				}
				else if ((tra.v_c <= tra.Speedlimit) && (tra.a >= 0))
				{
					tra.accelflag = 1;		//加速区間
				}
				else
				{
					tra.accelflag = 2;		//惰行区間
				}
			}
		}
		else													//吹上⇒大宮方面
		{
			if (tra.brakeflag == 1)		/*減速中の動作*/
			{
				if (tra.v_c >= 0.0)
				{
					tra.accelflag = 3;		//減速区間
				}
				else
				{
					tra.accelflag = 4;		//停車
					tra.brakeflag = 0;
					tra.t_stop_old = tra.t_stop;
				}
			}
			else
			{
				if (tra.laststopflag == 1)	//駅の手前に停車した場合の救済
				{
					tra.accelflag = 3;
					tra.brakeflag = 1;
				}
				else if (tra.x - tra.v * tra.v / 2.0 / tra.Bref < tra.X_nextstop && tra.v_c > speed_min - 20.0)
				{
					tra.accelflag = 3;
				}
				else if ((tra.v_c <= tra.Speedlimit) && (tra.a >= 0))
				{
					tra.accelflag = 1;		//加速区間
				}
				else
				{
					tra.accelflag = 2;		//惰行区間
				}
			}
		}
	}

}

void SimulateEngine::Calculation_traction_force_brake_force(TRAIN& tra) {

	/*** 各走行モードに応じた動作を指定 ***/
	if (tra.accelflag == 1)				//加速中
	{
		/** ブレーキ力計算(加速中なので0) **/
		tra.b_reg = 0.0;
		tra.b_air = 0.0;
		tra.b_reg_loss = 0.0;

		/** 車両引張力計算 **/
		Traction_force(tra);

		/**力行時d軸電流計算**/
		d_current_power(tra);

		/**力行時q軸電流計算**/
		q_current_power(tra);

		tra.iq_reg = 0.0;

		tra.ws = tra.R2 / tra.L2 * tra.iqtot / tra.idtot;
		tra.we = tra.ws + tra.wr;

		/*d軸電圧*/
		tra.v1d = tra.R1 * tra.idtot - (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.iqtotF;

		/*q軸電圧*/
		tra.v1q = tra.R1 * tra.iqtotF + (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.idtot + (tra.M_IM / tra.L2) * tra.M_IM * tra.we * tra.idtot;
		/*インバータ入力電圧*/
		tra.Em = sqrt(tra.v1d * tra.v1d + tra.v1q * tra.v1q);
	}
	else if (tra.accelflag == 2)			//惰行中
	{
		/** ブレーキ力計算(惰行中なので0) **/
		tra.b_reg = 0.0;
		tra.b_air = 0.0;
		tra.b_reg_loss = 0.0;

		/** 車両引張力計算 **/
		tra.Fmot = 0.0;
		tra.Ftot = 0.0;

		/*電流計算*/
		tra.idtot = 0.0;
		tra.iqtot = 0.0;
		tra.iqtotF = 0.0;
		tra.iq_reg = 0.0;

		//角周波数計算
		tra.ws = 0.0;
		tra.we = tra.ws + tra.wr;

		/*d軸電圧*/
		tra.v1d = 0.0;

		/*q軸電圧*/
		tra.v1q = 0.0;
		tra.Em = 0.0;

	}
	else if (tra.accelflag == 3)			//減速中
	{
		/** 車両ブレーキ力計算 **/
		Regenerative_force(tra);
		d_current_regenerative(tra);
		q_current_regenerative(tra);

		/********** 軽負荷回生中のトルク絞り込み模擬部（軽負荷回生制御下の電気ブレーキ力を計算） **********/
		if (tra.vfc > Vclim && tra.vfc < Vcmax)						//回生絞り込み中
		{
			tra.b_reg = tra.Fmot * (Vcmax - tra.vfc) / (Vcmax - Vclim);
			tra.iq_reg = (tra.b_reg * tra.rd / tra.Gr / tra.num_mot) * tra.L2 / (tra.POLE * tra.M_IM * tra.M_IM * tra.idtot);
			tra.iqtot = tra.iq_reg;
		}
		else if (tra.vfc >= Vcmax)										//回生絞り込み終了
		{
			tra.b_reg = 0.0;
			tra.iq_reg = 0.0;
			tra.iqtot = 0.0;
		}
		else															//回生絞り込みなし
		{
			tra.b_reg = tra.Fmot;
			tra.iq_reg = tra.iq_mot_R;
			tra.iqtot = tra.iq_reg;
		}

		tra.iq_regF = funcdelay(tra.Ts, tra.Tf, tra.iq_regF, tra.iq_reg);
		tra.iqtotF = tra.iq_regF;

		tra.b_reg_loss = tra.Fmot - tra.b_reg;		//回生絞り込み量[N]

		tra.b_air = -tra.BEmax - tra.b_reg;			//足りないブレーキ力は空気ブレーキで補う

		tra.ws = tra.R2 / tra.L2 * tra.iqtot / tra.idtot;
		tra.we = tra.ws + tra.wr;

		/*d軸電圧*/
		tra.v1d = tra.R1 * tra.idtot - (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.iqtotF;

		/*q軸電圧*/
		tra.v1q = tra.R1 * tra.iqtotF + (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.idtot + (tra.M_IM / tra.L2) * tra.M_IM * tra.we * tra.idtot;
		tra.Em = sqrt(tra.v1q * tra.v1q + tra.v1d * tra.v1d);
	}
	else if (tra.accelflag == 4)						//停車中
	{
		tra.Fmot = 0.0;
		tra.Ftot = 0.0;
		tra.idtot = 0.0;
		tra.iqtot = 0.0;
		tra.iqtotF = 0.0;
		tra.v = 0.0;
		tra.v_c = 0.0;
		tra.a = 0.0;
		tra.b_reg = 0.0;
		tra.b_air = 0.0;
		tra.b_reg_loss = 0.0;
		tra.iq_reg = 0.0;

		tra.ws = 0.0;
		tra.we = 0.0;
		tra.wr = 0.0;
		/*d軸電圧*/
		tra.v1d = 0.0;

		/*q軸電圧*/
		tra.v1q = 0.0;
		tra.Em = sqrt(tra.v1d * tra.v1d + tra.v1q * tra.v1q);
	}
	else if (tra.accelflag == 5) {
		if (tra.v_c > tra.con_speed) {
			tra.Fmot = -tra.Bmax * tra.num_mot * (tra.v_c - tra.con_speed) / 5.0;
			tra.Ftot = tra.Fmot;

			d_current_regenerative(tra);
			q_current_regenerative(tra);

			/********** 軽負荷回生中のトルク絞り込み模擬部（軽負荷回生制御下の電気ブレーキ力を計算） **********/
			if (tra.vfc > Vclim && tra.vfc < Vcmax)						//回生絞り込み中
			{
				tra.b_reg = tra.Fmot * (Vcmax - tra.vfc) / (Vcmax - Vclim);
				tra.iq_reg = (tra.b_reg * tra.rd / tra.Gr / tra.num_mot) * tra.L2 / (tra.POLE * tra.M_IM * tra.M_IM * tra.idtot);
				tra.iqtot = tra.iq_reg;
			}
			else if (tra.vfc >= Vcmax)										//回生絞り込み終了
			{
				tra.b_reg = 0.0;
				tra.iq_reg = 0.0;
				tra.iqtot = 0.0;
			}
			else															//回生絞り込みなし
			{
				tra.b_reg = tra.Fmot;
				tra.iq_reg = tra.iq_mot_R;
				tra.iqtot = tra.iq_reg;
			}

			tra.iq_regF = funcdelay(tra.Ts, tra.Tf, tra.iq_regF, tra.iq_reg);
			tra.iqtotF = tra.iq_regF;

			tra.b_reg_loss = tra.Fmot - tra.b_reg;		//回生絞り込み量[N]

			tra.b_air = tra.Ftot - tra.b_reg;			//足りないブレーキ力は空気ブレーキで補う
		}
		else {
			tra.Fmot = tra.Fmax * tra.num_mot * (tra.con_speed - tra.v_c) / 5.0;
			tra.Ftot = tra.Fmot;

			d_current_power(tra);
			q_current_power(tra);

			tra.iq_reg = 0.0;
		}

		tra.ws = tra.R2 / tra.L2 * tra.iqtot / tra.idtot;
		tra.we = tra.ws + tra.wr;

		/*d軸電圧*/
		tra.v1d = tra.R1 * tra.idtot - (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.iqtotF;

		/*q軸電圧*/
		tra.v1q = tra.R1 * tra.iqtotF + (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.idtot + (tra.M_IM / tra.L2) * tra.M_IM * tra.we * tra.idtot;
		tra.Em = sqrt(tra.v1q * tra.v1q + tra.v1d * tra.v1d);
	}

	tra.Fmot_c = tra.Fmot / 1000.0;
	tra.b_reg_c = tra.b_reg / 1000.0;
	tra.Ftot_c = tra.Ftot / 1000.0;
	tra.b_air_c = tra.b_air / 1000.0;
	tra.b_reg_loss_c = tra.b_reg_loss / 1000.0;
}

void SimulateEngine::Calculation_traction_circuit(TRAIN& tra)
{
	tra.isiv = P_SIV * tra.num / tra.vfc;

	if (tra.Type == 3) tra.isiv = P_KARAKIDA / tra.vfc;
	if (tra.Type == 4) tra.isiv = P_SHINYURIGAOKA / 1500.0;

	tra.ifl = tra.isiv + (tra.Pmot + tra.Ploss_mot + tra.Ploss_fe) / tra.vfc;

	tra.Ploss_mot = (tra.R1 * (tra.idtot * tra.idtot + tra.iqtotF * tra.iqtotF) + tra.R2 * (tra.M_IM / tra.L2 * tra.iqtotF) * (tra.M_IM / tra.L2 * tra.iqtotF)) * tra.num_mot;//モータの銅損
	tra.Ploss_fe = tra.g0 * (tra.v1d * tra.v1d + tra.v1q * tra.v1q) * tra.num_mot;
	tra.Ploss_fl = tra.Rfl * tra.ifl * tra.ifl;					//フィルタリアクトルのジュール損失[W]
	tra.Ploss_fl_c = tra.Ploss_fl / 1000.0;					//フィルタリアクトルのジュール損失[kW]

	tra.Ploss_mot_c = (tra.Ploss_mot + tra.Ploss_fe) / 1000.0;				//モータの銅損[kW]

	tra.Ploss_all = tra.Ploss_mot + tra.Ploss_fe + tra.Ploss_fl;
	tra.Ploss_all_c = tra.Ploss_all / 1000.0;
}

void SimulateEngine::Calculation_power_motor(TRAIN& tra)
{
	tra.Ptot = tra.v * tra.Ftot;

	if (tra.accelflag == 1)				/*力行中*/
	{
		tra.Pmot = (tra.v1d * tra.idtot + tra.v1q * tra.iqtotF) * tra.num_mot;
		tra.Pmot_air = 0.0;
		tra.Pmot_reg_loss = 0.0;

	}
	else if (tra.accelflag == 3)		/*回生中*/
	{
		tra.Pmot = (tra.v1d * tra.idtot + tra.v1q * tra.iqtotF) * tra.num_mot;
		tra.Pmot_air = tra.v * tra.b_air;
		tra.Pmot_reg_loss = tra.v * tra.b_reg_loss;
	}
	else if (tra.accelflag == 5) {
		if (tra.v_c > tra.con_speed) {
			tra.Pmot = (tra.v1d * tra.idtot + tra.v1q * tra.iqtotF) * tra.num_mot;
			tra.Pmot_air = tra.v * tra.b_air;
			tra.Pmot_reg_loss = tra.v * tra.b_reg_loss;
		}
		else {
			tra.Pmot = (tra.v1d * tra.idtot + tra.v1q * tra.iqtotF) * tra.num_mot;
			tra.Pmot_air = 0.0;
			tra.Pmot_reg_loss = 0.0;
		}
	}
	else
	{
		tra.Pmot = 0.0;
		tra.Pmot_air = 0.0;
		tra.Pmot_reg_loss = 0.0;
	}

	tra.Pres = tra.v * tra.R_total;		//走行抵抗でのロス

	tra.Ptot_c = tra.Ptot / 1000.0;
	tra.Pmot_c = tra.Pmot / 1000.0;
	tra.Pmot_air_c = tra.Pmot_air / 1000.0;
	tra.Pres_c = tra.Pres / 1000.0;
	tra.Pmot_reg_loss_c = tra.Pmot_reg_loss / 1000.0;

}

void SimulateEngine::Calculation_power_train(TRAIN& tra)
{
	tra.Pvh = tra.vp * tra.ifl;

	if (tra.accelflag == 1 || tra.accelflag == 2)				//力行中,惰行中
	{
		tra.Pvh_in = tra.vp * tra.ifl;
		tra.Pvh_out = 0.0;
		tra.Pvh_st = 0.0;
	}
	else if (tra.accelflag == 3)		//回生中
	{
		tra.Pvh_in = 0.0;
		tra.Pvh_out = tra.vp * tra.ifl;
		tra.Pvh_st = 0.0;
	}
	else if (tra.accelflag == 5) {
		if (tra.v_c > tra.con_speed) {
			tra.Pvh_in = 0.0;
			tra.Pvh_out = tra.vp * tra.ifl;
			tra.Pvh_st = 0.0;
		}
		else {
			tra.Pvh_in = tra.vp * tra.ifl;
			tra.Pvh_out = 0.0;
			tra.Pvh_st = 0.0;
		}
	}
	else
	{
		tra.Pvh_in = 0.0;
		tra.Pvh_out = 0.0;
		tra.Pvh_st = tra.vp * tra.ifl;
	}
	tra.Pvh_c = tra.Pvh / 1000.0;
	tra.Pvh_in_c = tra.Pvh_in / 1000.0;
	tra.Pvh_out_c = tra.Pvh_out / 1000.0;
	tra.Pvh_st_c = tra.Pvh_st / 1000;
}

void SimulateEngine::Calculation_energy_train(TRAIN& tra)
{
	if (tra.accelflag == 1)				//力行中
	{
		tra.Emot_pow += tra.Pmot * dt;

		tra.Eloss_inv_pow += tra.Ploss_inv * dt;
		tra.Eloss_mot_pow += (tra.Ploss_mot + tra.Ploss_fe) * dt;
		tra.Eloss_fl_pow += tra.Ploss_fl * dt;
		tra.Eloss_all_pow += tra.Ploss_all * dt;
		tra.Esiv_pow += P_SIV * dt;

		tra.Eres_pow += tra.Pres * dt;
		//		tra.Emot_pow += (tra.Pvh_in - tra.Ploss_all - P_SIV - tra.Pres)* dt;
	}
	else if (tra.accelflag == 3)		//回生中
	{
		tra.Emot_reg += tra.Pmot * dt;
		tra.Emot_reg_loss += tra.Pmot_reg_loss * dt;

		tra.Eloss_inv_reg += tra.Ploss_inv * dt;
		tra.Eloss_mot_reg += (tra.Ploss_mot + tra.Ploss_fe) * dt;
		tra.Eloss_fl_reg += tra.Ploss_fl * dt;
		tra.Eloss_all_reg += tra.Ploss_all * dt;
		tra.Esiv_reg += P_SIV * dt;

		tra.Eres_reg += tra.Pres * dt;
		//		tra.Emot_reg += (tra.Pvh_out - tra.Ploss_all - P_SIV - tra.Pres)* dt;
	}
	else if (tra.accelflag == 2)
	{
		tra.Eres_coa += tra.Pres * dt;
		tra.Esiv_coa += P_SIV * dt;
	}
	else if (tra.accelflag == 4)
	{
		tra.Esiv_stp += P_SIV * dt;
		if (tra.Type == 3) tra.Esiv_stp += P_KARAKIDA * dt;
		if (tra.Type == 4) tra.Esiv_stp += P_SHINYURIGAOKA * dt;
	}
	else if (tra.accelflag == 5) {
		if (tra.v_c > tra.con_speed) {
			tra.Emot_reg += tra.Pmot * dt;
			tra.Emot_reg_loss += tra.Pmot_reg_loss * dt;

			tra.Eloss_inv_reg += tra.Ploss_inv * dt;
			tra.Eloss_mot_reg += (tra.Ploss_mot + tra.Ploss_fe) * dt;
			tra.Eloss_fl_reg += tra.Ploss_fl * dt;
			tra.Eloss_all_reg += tra.Ploss_all * dt;
			tra.Esiv_reg += P_SIV * dt;

			tra.Eres_reg += tra.Pres * dt;
		}
		else {
			tra.Emot_pow += tra.Pmot * dt;

			tra.Eloss_inv_pow += tra.Ploss_inv * dt;
			tra.Eloss_mot_pow += (tra.Ploss_mot + tra.Ploss_fe) * dt;
			tra.Eloss_fl_pow += tra.Ploss_fl * dt;
			tra.Eloss_all_pow += tra.Ploss_all * dt;
			tra.Esiv_pow += P_SIV * dt;

			tra.Eres_pow += tra.Pres * dt;
		}
	}

	tra.Emot_air += tra.Pmot_air * dt;


	tra.Emot_pow_c = tra.Emot_pow / 1000.0 / 3600.0;
	tra.Emot_reg_c = tra.Emot_reg / 1000.0 / 3600.0;
	tra.Emot_air_c = tra.Emot_air / 1000.0 / 3600.0;
	tra.Emot_reg_loss_c = tra.Emot_reg_loss / 1000.0 / 3600.0;

	tra.Eloss_inv_pow_c = tra.Eloss_inv_pow / 1000.0 / 3600.0;
	tra.Eloss_mot_pow_c = tra.Eloss_mot_pow / 1000.0 / 3600.0;
	tra.Eloss_fl_pow_c = tra.Eloss_fl_pow / 1000.0 / 3600.0;
	tra.Eloss_all_pow_c = tra.Eloss_all_pow / 1000.0 / 3600.0;

	tra.Eloss_inv_reg_c = tra.Eloss_inv_reg / 1000.0 / 3600.0;
	tra.Eloss_mot_reg_c = tra.Eloss_mot_reg / 1000.0 / 3600.0;
	tra.Eloss_fl_reg_c = tra.Eloss_fl_reg / 1000.0 / 3600.0;
	tra.Eloss_all_reg_c = tra.Eloss_all_reg / 1000.0 / 3600.0;

	tra.Eres_pow_c = tra.Eres_pow / 1000.0 / 3600.0;
	tra.Eres_reg_c = tra.Eres_reg / 1000.0 / 3600.0;
	tra.Eres_coa_c = tra.Eres_coa / 1000.0 / 3600.0;

	tra.Esiv_pow_c = tra.Esiv_pow / 1000.0 / 3600.0;
	tra.Esiv_reg_c = tra.Esiv_reg / 1000.0 / 3600.0;
	tra.Esiv_coa_c = tra.Esiv_coa / 1000.0 / 3600.0;
	tra.Esiv_stp_c = tra.Esiv_stp / 1000.0 / 3600.0;

	tra.Evh += tra.Pvh * dt;
	tra.Evh_in += tra.Pvh_in * dt;
	tra.Evh_out += tra.Pvh_out * dt;
	tra.Evh_st += tra.Pvh_st * dt;

	tra.Evh_c = tra.Evh / 1000.0 / 3600.0;
	tra.Evh_in_c = tra.Evh_in / 1000.0 / 3600.0;
	tra.Evh_out_c = tra.Evh_out / 1000.0 / 3600.0;
	tra.Evh_st_c = tra.Evh_st / 1000.0 / 3600.0;
}

void SimulateEngine::Calculation_power_sub(SUB& sub)
{
	sub.Pss = sub.iss * sub.vout;
	sub.Pss_c = sub.Pss / 1000.0;
}

void SimulateEngine::Calculation_energy_sub(SUB& sub)
{
	if (t >= 15.0 * 60.0 && t < 45.0 * 60.0) {
		sub.Ess += sub.Pss * dt;						//[J]
		sub.Ess_c = sub.Ess / 1000.0 / 3600.0;		//[kWh]
		if (sub.diode == 1)
		{
			sub.Wss += sub.Rss * sub.iss * sub.iss * dt;		//[J]
			sub.Wss_c = sub.Wss / 1000 / 3600;				//[kWh]
		}
	}
}

/*** xおよびyが指す要素を交換 ***/
void SimulateEngine::sort_s(const std::unique_ptr<NODE[]>& data, int n)
{

	std::vector<NODE> temp(N_node);
	for (size_t i = 0; i < n; i++)temp.at(i) = data[i];

	std::sort(temp.begin(), temp.end(), [](auto const& ll, auto const& rr) {return ll.X < rr.X; });

	for (size_t i = 0; i < n; i++)data[i] = temp[i];

	/*int k = n - 1;
	while (k >= 0)
	{
		int i, j;
		for (i = 1, j = -1; i <= k; i++)
			if (data[i - 1].X > data[i].X) {
				j = i - 1;
				swap_node(data[i], data[j]);
			}
		k = j;
	}*/

	return;
}


/*** 配列data[]の先頭n個の要素を距離の昇順にソート ***/
void SimulateEngine::sort(const std::unique_ptr<NODE[]>& data, int n)
{
	int k = n - 1;
	while (k >= 0)
	{
		size_t i = 0, j = 0;
		for (i = 1, j = -1; i <= k; i++)
			if (data[i - 1].X > data[i].X) {
				j = i - 1;
				NODE tempI, tempJ;
				tempI = data[i];
				tempJ = data[j];
				data[j] = tempI;
				data[i] = tempJ;
				//swap(data[i], data[j]);
			}
		k = j;
	}
}

void SimulateEngine::Initialize_BRANCH(BRANCH& branch)
{
	//*branch = (BRANCH*)malloc(sizeof(BRANCH));
	//if (*branch == NULL)
	//{
	//	//エラー処理
	//}
}

int SimulateEngine::Count_Trains_Direction(const std::unique_ptr<TRAIN[]>& tra, int direction) {
	int count_d = 0;
	for (int i = 0; i < NUM_tra; i++) {
		if (tra[i].direction == direction) count_d++;
	}
	return count_d;
}

void SimulateEngine::Make_branch(const std::unique_ptr<BRANCH[]>& branch, const std::unique_ptr<NODE[]>& node_order, const std::unique_ptr<TRAIN[]>& tra)
{
	size_t i = 0, j = 0;						//j：正方向に接続されているノード探索用，k：負方向に接続されているノード探索用
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
		//printf("[%d]flag = %d\n", temp[i].Number, temp[i].flag);
		if (node_order[i].flag == 0)				//In case of SS( i Starts with 0 to NUM_sub-1 )
		{
			//puts("SS");
			branch[count_SS].X_start = node_order[i].X;		//位置[m]（基準点からの距離）
			branch[count_SS].X_end = node_order[i].X;		//位置[m]（基準点からの距離）
			branch[count_SS].r = node_order[i].r;			//Resistance [Ohm]
			branch[count_SS].flag = node_order[i].flag;		//ブランチの分類（1：UP, -1：DOWN, 0：SS）
			branch[count_SS].Node_pos = -1;			//正方向側に接続されているノード№
			branch[count_SS].Node_neg = -1;			//負方向側に接続されているノード№
			count_SS = count_SS + 1;
		}
	}

	/*②上りのブランチを作成*/
	//printf("-----------UP--------------\n");
	for (i = 0; i < N_node; i++)
	{
		if (node_order[i].flag >= 0)			//In case of +1 direction
		{
			//printf("[%d]flag = %d\n", temp[i].Number, temp[i].flag);
			if (count_UP < NUM_tra_UP + 2 * NUM_sub - 1)		//終端以外
			{

				for (size_t j = 0; node_order[i + (j + 1)].flag == -1; j++) {}				//temp[i]に対して正方向側に接続されているノードを探索
				//printf("\nj = %d\n", j);

				/*temp[i]に対して正方向に接続されているブランチBを作成*/
				branch[count_UP].Node_pos = node_order[i + (j + 1)].Number;			//ブランチBについて，正方向側に接続されているノード№
				branch[count_UP].Node_neg = node_order[i].Number;					//ブランチBについて，負方向側に接続されているノード№

				branch[count_UP].X_start = node_order[i].X;							//[m]
				branch[count_UP].X_end = node_order[i + (j + 1)].X;					//[m]

				if (node_order[i + (j + 1)].X <= 7500.0) {
					branch[count_UP].r = R_StoN_UP * std::abs(node_order[i + (j + 1)].X - node_order[i].X) / 1000.0;		//Resistance [Ohm]
				}
				else {
					branch[count_UP].r = R_NtoK_UP * std::abs(node_order[i + (j + 1)].X - node_order[i].X) / 1000.0;		//Resistance [Ohm]
				}

				branch[count_UP].flag = 1;										//ブランチの分類（1：UP, -1：DOWN, 0：SS）

				/*次のループのための準備*/
				j = 0;
				count_UP = count_UP + 1;
			}
			else {}		//終端のノードでは何もしない
		}
	}


	/*③下りのブランチを作成*/
	//printf("-----------DOWN--------------\n");
	for (i = 0; i < N_node; i++)
	{
		if (node_order[i].flag <= 0)			//In case of +1 direction
		{
			//printf("[%d]flag = %d\n", temp[i].Number, temp[i].flag);
			if (count_DOWN < N_branch)		//終端以外
			{
				for (size_t j = 0; node_order[i + (j + 1)].flag == 1; j++) {}				//temp[i]に対して正方向側に接続されているノードを探索
				//printf("\nj = %d\n", j);

				/*正方向に接続されているブランチを作成*/
				branch[count_DOWN].Node_pos = node_order[i + (j + 1)].Number;		//正方向側に接続されているノード№
				branch[count_DOWN].Node_neg = node_order[i].Number;					//負方向側に接続されているノード№

				branch[count_DOWN].X_start = node_order[i].X;						//[m]
				branch[count_DOWN].X_end = node_order[i + (j + 1)].X;				//[m]

				if (node_order[i + (j + 1)].X <= 7500.0) {
					branch[count_DOWN].r = R_StoN_DOWN * std::abs(node_order[i + (j + 1)].X - node_order[i].X) / 1000.0;		//Resistance [Ohm]
				}
				else {
					branch[count_DOWN].r = R_NtoK_DOWN * std::abs(node_order[i + (j + 1)].X - node_order[i].X) / 1000.0;		//Resistance [Ohm]
				}

				branch[count_DOWN].flag = -1;									//ブランチの分類（1：UP, -1：DOWN, 0：SS）

				/*次のループのための準備*/
				j = 0;
				count_DOWN = count_DOWN + 1;
			}
			else {}		//終端のノードでは何もしない
		}
	}
}

void SimulateEngine::Initialize_matrix1(const std::unique_ptr<double[]>& Mtr, size_t M, size_t N)		/*M行N列の２次元配列を初期化する関数*/
{

	for (size_t i = 0; i < M; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			Mtr[i * N + j] = 0.0000;
		}
	}
}

void SimulateEngine::Make_H_matrix(const std::unique_ptr<BRANCH[]>& branch, const std::unique_ptr<double[]>& Hm)		/*H matrix作成用関数*/
{
	//int i;	/*行*/
	//int j;	/*列*/

	for (size_t i = 0; i < N_node; i++)
	{
		for (size_t j = 0; j < NUM_sub; j++)
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
	for (size_t i = 0; i < N_node; i++)
	{
		for (size_t j = NUM_sub; j < N_branch; j++)
		{
			if (i == branch[j].Node_pos)
			{
				Hm[i + j * N_node] = 1.0;
			}
			else if (i == branch[j].Node_neg)
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

void SimulateEngine::Make_Y_matrix(const std::unique_ptr<BRANCH[]>& branch, const std::unique_ptr<double[]>& Y)		/*Y matrix作成用関数*/
{
	//int i;	/*行*/
	//int j;	/*列*/

	for (size_t i = 0; i < N_branch; i++)
	{
		for (size_t j = 0; j < N_branch; j++)
		{
			//printf("\ni = %d\nj = %d\n",i,j);
			if (i == j)
			{
				//printf("\ni = %d\nj = %d\n",i,j);
				if (branch[j].r < 0.0001)
				{
					Y[i + j * N_branch] = 10000.0;
				}
				else
				{
					Y[i + j * N_branch] = 1.0 / branch[j].r;		/*ここでコンダクタンス[S]に変換*/
				}
			}
			else
			{
				Y[i + j * N_branch] = 0.0;
			}
		}
	}
}

void SimulateEngine::Make_In_vector(const std::unique_ptr<NODE[]>& data, const std::unique_ptr<double[]>& In)		/*In vector作成用関数*/
{

	for (size_t i = 0; i < N_node; i++)
	{
		In[i] = data[i].I;
	}
}

void SimulateEngine::Transpose(const std::unique_ptr<double[]>& trans, const std::unique_ptr<double[]>& X, size_t row, size_t column)		/*転置行列計算用関数*/
{
	//int i, j;
	double temp = 0.0;

	for (size_t i = 0; i < row; i++)
	{
		for (size_t j = 0; j < column; j++)
		{
			temp = X[i + j * row];
			trans[j + i * column] = temp;
		}
	}
}

int SimulateEngine::mySimulate(std::vector<double> wait_time, CsvWriter& _cw, std::mutex& _mtx)
{
	//Init();

	//H.reset(); H = std::make_unique<double[]>(N_node * N_branch);		//Connection Matrix
	//H_tra.reset(); H_tra = std::make_unique<double[]>(N_branch * N_node);	//Transpose matrix of H
	//Y.reset(); Y = std::make_unique<double[]>(N_branch * N_branch);		//Conductance Matrix
	//A.reset(); A = std::make_unique<double[]>(N_node * N_node);			//A matrix(H*Y*H_tra)
	//TEMP.reset(); TEMP = std::make_unique<double[]>(N_node * N_branch);	//Test
	//In.reset(); In = std::make_unique<double[]>(N_node);					//Node Current Matrix
	//In_cpy.reset(); In_cpy = std::make_unique<double[]>(N_node);				//Copy of Node Current Matrix
	//Vn.reset(); Vn = std::make_unique<double[]>(N_node);					//Node Voltage Matrix
	//Ib.reset(); Ib = std::make_unique<double[]>(N_branch);				//Branch Current Matrix
	//Vb.reset(); Vb = std::make_unique<double[]>(N_branch);
	//ipiv_Vn.reset(); ipiv_Vn = std::make_unique<int[]>(N_node);
	//tra.reset(); tra = std::make_unique<TRAIN[]>(NUM_tra);	//ポインタの配列として定義
	//sub.reset(); sub = std::make_unique<SUB[]>(NUM_sub);
	//node.reset(); node = std::make_unique<NODE[]>(N_node);
	//node_order.reset(); node_order = std::make_unique<NODE[]>(N_node);
	//branch.reset(); branch = std::make_unique<BRANCH[]>(N_branch);

	ERROR = 0;						//回路計算の再計算判定用フラグ（flagを全て加算してNUM_sub以上なら回路計算ループを抜ける。NUM_sub以下なら再計算。）
	ERROR_REG = 0;
	count_loop = 0;
	count_loop_rec = 0;

	t = 0.0;

	sub_total = 0.0;

	/***  列車構造体・変電所構造体・ESD構造体初期化  ***/
	for (size_t i = 0; i < NUM_sub; i++)
	{
		Make_substation(sub[i], ini_sub[i], i);				//ナンバリングは通し番号
	}
	for (size_t i = 0; i < NUM_tra; i++)
	{
		Make_train(tra[i], ini_tra[i], i + NUM_sub);							//ナンバリングは通し番号
	}

	for (size_t i = 0; i < NUM_tra; i++)
	{
		Calculate_BEmax(tra[i]);
	}

	////////////////////////// 数値計算部 /////////////////////////
	for (t = 0.0; t <= 60.0 * 45.0; t += dt)	/***** 1時間シミュレーション *****/
	{
		minute = t / 60.0;
		count_loop_rec = 0;

		////////////////////////// 車両運動方程式計算部 /////////////////////
			/*** 車両走行パターン模擬 & 運動方程式計算 ***/
			/*** 車両パワー計算 ***/

		for (size_t i = 0; i < NUM_tra; i++)
		{
			Calculate_BEmax(tra[i]);

			Run_pattern(tra[i], sta, t, wait_time);

		}

		/////////////////////////////////////////////////////////////////////

		while (count_loop <= 3000002)
		{
			// start_while:				//ダイオード判定が間違っていた場合，continueでここに戻ってくる

			if (count_loop > 3000000)
			{
				printf("\n");
				printf("-----------t = %f--------------\n", t);
				printf("-----------accelflag & position vehicle--------------\n");
				for (size_t i = 0; i < NUM_tra; i++)
				{
					printf("%d, %f\n", tra[i].accelflag, tra[i].x);
				}
				/*				printf("-----------Fmot vehicle--------------\n");
								for (size_t i = 0; i < NUM_tra; i++)
								{
									printf("%f\n", tra[i].Fmot);
								}*/
				printf("-----------Pmot vehicle--------------\n");
				for (size_t i = 0; i < NUM_tra; i++)
				{
					printf("%f\n", tra[i].Pmot);
				}
				printf("-----------Vp vehicle--------------\n");
				for (size_t i = 0; i < NUM_tra; i++)
				{
					printf("%f\n", tra[i].vp);
				}
				printf("-----------Vfc vehicle--------------\n");
				for (size_t i = 0; i < NUM_tra; i++)
				{
					printf("%f\n", tra[i].vfc);
				}
				printf("-----------ifl vehicle--------------\n");
				for (size_t i = 0; i < NUM_tra; i++)
				{
					printf("%f\n", tra[i].ifl);
				}

				printf("\n");
			}

			if (count_loop > 3000001) {
				printf("Holy shit!!\n");
				exit(EXIT_FAILURE);
				break;
			}

			/***  node[]の初期化  ***/
			for (size_t i = 0; i < NUM_sub; i++)
			{
				Make_NODE_SS(node[i], sub[i]);
				Make_NODE_SS(node_order[i], sub[i]);
			}

			for (size_t i = 0; i < NUM_tra; i++)
			{
				Make_NODE_TRAIN(node[i + NUM_sub], tra[i]);
				Make_NODE_TRAIN(node_order[i + NUM_sub], tra[i]);
			}

			for (size_t i = 0; i < N_branch; i++)
			{
				Initialize_BRANCH(branch[i]);
			}

			/*** 各Matrixの初期化 ***/
			Initialize_matrix1(H, N_node, N_branch);
			Initialize_matrix1(H_tra, N_branch, N_node);
			Initialize_matrix1(Y, N_branch, N_branch);
			Initialize_matrix1(A, N_node, N_node);
			Initialize_matrix1(TEMP, N_node, N_branch);
			Initialize_matrix1(In, N_node, 1);
			Initialize_matrix1(In_cpy, N_node, 1);
			Initialize_matrix1(Vn, N_node, 1);
			Initialize_matrix1(Ib, N_branch, 1);
			Initialize_matrix1(Vb, N_branch, 1);

			/***  node_order[]を距離順にソート  ***/
			//sort_s(node_order, N_node);
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
			dgemm(&trans, &trans, &K_node, &K_branch, &K_branch, &alpha, H.get(), &ld_node, Y.get(), &ld_branch, &beta, TEMP.get(), &ld_node);

			//上記行列×転置H行列 (dgemmはintel MKLを用いて行列の掛け算を計算)
			dgemm(&trans, &trans, &K_node, &K_node, &K_branch, &alpha, TEMP.get(), &ld_node, H_tra.get(), &ld_branch, &beta, A.get(), &ld_node);

			//連立1次方程式を解くためにH行列×Y行列×転置H行列の解A(N行N列)をLU分解
			dgetrf(&ld_node, &K_branch, A.get(), &ld_node, ipiv_Vn.get(), &info_Vn);

			//dgetrsで連立1次方程式「A*Vn=In」を解き，Inに解(ノード電圧行列)を上書き
			dgetrs(&trans, &K_node, &nrhs, A.get(), &ld_node, ipiv_Vn.get(), In_cpy.get(), &ld_node, &info_Vn);	//In_cpyに方程式の解Vnが格納されていることに注意

			Make_Y_matrix(branch, Y);		//なぜかVnを求める際の"dgetrf"でY-matrixが書き換えられているため，ここで再計算させてる（実際は無駄な計算なのでしたくないが…）

			for (size_t i = 0; i < N_node; i++)	/*Vnベクトルはマイナスをつける*/
			{
				Vn[i] = -In_cpy[i];
			}

			//転置H行列×Vnベクトルでブランチ電圧ベクトルを計算
			dgemv(&trans, &K_branch, &K_node, &alpha, H_tra.get(), &ld_branch, Vn.get(), &incx, &beta, Vb.get(), &incy);

			//Y行列×ブランチ電圧ベクトルでブランチ電流ベクトルを計算
			dgemv(&trans, &K_branch, &K_branch, &alpha, Y.get(), &ld_branch, Vb.get(), &incx, &beta, Ib.get(), &incy);


			for (size_t i = 0; i < N_node; i++)			//ノード電圧の計算結果をノード構造体へ返す（結果的には，Vn配列はノード構造体の順番通りに並んでいる？）
			{
				node[i].V = Vn[i];
			}

			for (size_t i = 0; i < NUM_sub; i++)			//ノード電圧・電流とブランチ電流から，変電所出力電圧・電流を計算
			{
				sub[i].vout = node[i].V;
				sub[i].iss = -In[i] + Ib[i];
				Error_Detection_SS(sub[i]);
				ERROR += sub[i].flag;
			}

			for (size_t i = 0; i < NUM_tra; i++)
			{
				tra[i].vfc_old = tra[i].vfc;

				for (size_t j = 0; j < N_node; j++)
				{
					if (node[j].Number == tra[i].name_train)	//ノード番号と列車識別番号が一致していたらノード電圧を列車パンタ点電圧に返す
					{
						tra[i].vp = node[j].V;
						tra[i].vfc = tra[i].vp - tra[i].Rfl * tra[i].ifl;

						//引張力・ブレーキ力・モータ電流・モータ電圧を計算
						Calculation_traction_force_brake_force(tra[i]);

						//dq軸電流・電圧からモータパワーを計算
						Calculation_power_motor(tra[i]);

						//主回路計算
						Calculation_traction_circuit(tra[i]);
					}
				}
			}

			for (size_t i = 0; i < NUM_tra; i++)
			{
				tra[i].judgement = std::abs(tra[i].vfc - tra[i].vfc_old) / tra[i].vfc_old;

				if (tra[i].judgement < 0.001) {
					tra[i].reg_flag = 1;
				}
				else {
					tra[i].reg_flag = 0;
				}

				if (tra[i].vfc > Vcmax) tra[i].reg_flag = 0;

				ERROR_REG += tra[i].reg_flag;
			}

			/*		for (size_t i = 0; i < NUM_tra; i++)
					{
						if (tra[i].accelflag == 3 || (tra[i].accelflag == 5 && tra[i].v_c > tra[i].con_speed)) {
							Error_Detection_REG(tra[i]);
						}
						else {
							tra[i].reg_flag = 1;
						}

						ERROR_REG += tra[i].reg_flag;
					}*/

					/*		if (t > 1020.05 && t < 1020.09) {
								printf("-----------t = %f--------------\n", t);
								printf("-----------accelflag vehicle--------------\n");
								for (size_t i = 0; i < NUM_tra; i++)
								{
									printf("%d\n", tra[i].accelflag);
								}
								printf("-----------Fmot vehicle--------------\n");
								for (size_t i = 0; i < NUM_tra; i++)
								{
									printf("%f\n", tra[i].Fmot);
								}
								printf("-----------Pmot vehicle--------------\n");
								for (size_t i = 0; i < NUM_tra; i++)
								{
									printf("%f\n", tra[i].Pmot);
								}
								printf("-----------Vp vehicle--------------\n");
								for (size_t i = 0; i < NUM_tra; i++)
								{
									printf("%f\n", tra[i].vp);
								}
								printf("-----------ifl vehicle--------------\n");
								for (size_t i = 0; i < NUM_tra; i++)
								{
									printf("%f\n", tra[i].ifl);
								}
								printf("-----------vss substation--------------\n");
								for (size_t i = 0; i < NUM_sub; i++)
								{
									printf("%f\n", sub[i].vout);
								}
								printf("-----------iss substation--------------\n");
								for (size_t i = 0; i < NUM_sub; i++)
								{
									printf("%f\n", sub[i].iss);
								}
								printf("\n");
								printf("\n");
								printf("\n");
							}*/

			if (ERROR < NUM_sub || ERROR_REG < NUM_tra)
			{
				/*for (size_t i = 0; i < N_node; i++)
				{
					free(node[i]);
					free(node_order[i]);
				}
				for (size_t i = 0; i < N_branch; i++)
				{
					free(branch[i]);
				}*/


				ERROR = 0;
				ERROR_REG = 0;
				count_loop++;
				//goto start_while;					//while文の先頭に戻る
				continue;
			}
			else
			{
				ERROR = 0;
				ERROR_REG = 0;
				count_loop_rec = count_loop;
				count_loop = 0;
				break;
			}
		}

		/***  次の状態を計算  ***/
		for (size_t i = 0; i < NUM_tra; i++)
		{
			if (tra[i].accelflag != 4) Solve_next_state(tra[i]);
		}

		/***  車両エネルギー計算  ***/
		for (size_t i = 0; i < NUM_tra; i++)
		{
			Calculation_power_train(tra[i]);
			Calculation_energy_train(tra[i]);
		}

		for (size_t i = 0; i < NUM_sub; i++)
		{
			Calculation_power_sub(sub[i]);
			Calculation_energy_sub(sub[i]);
		}

		/*for (size_t i = 0; i < N_node; i++)
		{
			free(node[i]);
			free(node_order[i]);
		}

		for (size_t i = 0; i < N_branch; i++)
		{
			free(branch[i]);
		}*/

		/*node = std::make_unique<NODE[]>(N_node);
		node_order = std::make_unique<NODE[]>(N_node);
		branch = std::make_unique<BRANCH[]>(N_branch);*/

		/////////////////////////////  cmd表示部  ////////////////////////////

		//////////////////////////////////////////////////////////////////////
	}
	/********** 数値計算ループ終了 **********/

	for (size_t i = 0; i < NUM_sub; i++)
	{
		sub_total += sub[i].Ess_c;
	}

	//fprintf(file, "%.1f, %f, %f, %f, %f, %f\n", wait_time[0], wait_time[1], wait_time[2], wait_time[3], wait_time[4], sub_total);

	std::string wl = "";
	for (size_t i = 0; i < wait_time.size(); i++) {
		wl += std::to_string(wait_time.at(i)) + ",";
	}

	wl += std::to_string(sub_total) + ",";
	//std::lock_guard<std::mutex> _lock(_mtx);
	_cw.PushTask(wl, _mtx);


	/*****  メモリ開放  *****/
	// スマートポインタなので不要
	/*for (size_t i = 0; i < NUM_tra; i++)
	{
		free(tra[i]);
	}
	for (size_t i = 0; i < NUM_sub; i++)
	{
		free(sub[i]);
	}*/

	//std::unique_ptr<TRAIN[]> tra = std::make_unique<TRAIN[]>(NUM_tra);	//ポインタの配列として定義
	//std::unique_ptr<SUB[]> sub = std::make_unique<SUB[]>(NUM_sub);

	//std::this_thread::sleep_for(std::chrono::milliseconds(1000));
	int sleeper = 0;
	while (sleeper < 10000000)sleeper++;


	return  0;

}


SimulateEngine::~SimulateEngine()
{

	//std::this_thread::sleep_for(std::chrono::milliseconds(2000));

	int sleeper = 0;
	while (sleeper < 10000000)sleeper++;

}

bool ExeSimulate(SimulateEngine& SE, std::vector<double> _wt, CsvWriter& _cw, std::mutex& _mtx)
{

	SE.mySimulate(_wt, std::ref(_cw), std::ref(_mtx));


	return true;
}
