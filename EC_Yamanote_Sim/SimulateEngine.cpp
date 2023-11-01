#include "SimulateEngine.hpp"

SimulateEngine::SimulateEngine()
{
	Init();
}

bool SimulateEngine::Init()
{
	/***  �e�s��錾��  ***/
	/***  �e�s��錾��  ***/
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

	ERROR = 0;						//��H�v�Z�̍Čv�Z����p�t���O�iflag��S�ĉ��Z����NUM_sub�ȏ�Ȃ��H�v�Z���[�v�𔲂���BNUM_sub�ȉ��Ȃ�Čv�Z�B�j
	ERROR_REG = 0;
	count_loop = 0;
	count_loop_rec = 0;
	//	test_flag = 0;

		/***  �}�g���b�N�X�v�Z�֌W  ***/
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
	/***  �\���̐錾��  ***/
	tra = std::make_unique<TRAIN[]>(NUM_tra);	//�|�C���^�̔z��Ƃ��Ē�`
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
	//	//�G���[����
	//}

	tra.name_train = i;
	tra.Type = ini.Type;
	tra.notch = 3.0;
	tra.brake = 5.0;

	tra.mass_M = ini.mass_M * ini.num_M;			//[kg]�����p(100%���)
	tra.mass_T = ini.mass_T * ini.num_T;			//[kg]�����p(100%���)
	tra.mass_P = (ini.num_M * 152.0 + (ini.num_T - 2.0) * 152.0 + 2.0 * 144.0) * 55.0 / 2.0;

	if (tra.Type == 2) tra.mass_P = 0.0;

	tra.mass = tra.mass_M + tra.mass_T + tra.mass_P + 15.8 * 1000.0 * (ini.num_M + ini.num_T) / 10.0;				//[kg](250%���)
	tra.mass_c = tra.mass / 1000.0;										//[t](250%���)

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

	tra.Fmax = 18.818 * 1000.0;         //[N/MM]  �����p(100%��Ԏ�)(�Q�l�F�N�������x0.639[m/s/s])
	tra.Bmax = 15.556 * 1000.0;		//[N/MM] (100%��Ԏ�)
	tra.Bref = 1.111 * 5.0 / 7.0;
	tra.BEmax = 0.0;

	tra.i1d_P = 90.42780281;  //�͍s���C��g���N�̈�ł�d���d��
	tra.i1q_P = 237.1149987; //�͍s���C��g���N�̈�ł�q���d��
	tra.i1d_R = 91.09231795;  //�񐶎��C��g���N�̈�ł�d���d��
	tra.i1q_R = -211.1226486; //�񐶎��C��g���N�̈�ł�q���d��

	if (tra.Type == 2) {
		tra.notch = 4.0;
		tra.Fmax = 14.776 * 1000.0;         //[N/MM]  �����p(100%��Ԏ�)(�Q�l�F�N�������x0.639[m/s/s])
		tra.Bmax = 13.910 * 1000.0;		//[N/MM] (100%��Ԏ�)
		tra.Bref = 1.111 * 5.0 / 7.0;
		tra.BEmax = 0.0;

		tra.i1d_P = 91.877;  //�͍s���C��g���N�̈�ł�d���d��
		tra.i1q_P = 183.246; //�͍s���C��g���N�̈�ł�q���d��
		tra.i1d_R = 91.60374;  //�񐶎��C��g���N�̈�ł�d���d��
		tra.i1q_R = -173.0142; //�񐶎��C��g���N�̈�ł�q���d��
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

	tra.Gr = 6.31;  //�M�A��
	tra.Ts = 0.00001;  //��������
	tra.Tf = 0.0001;  //�����萔
	tra.POLE = 2.0;  //�ɑΐ�
	tra.rd = 0.41;   //�ԗ֔��a[m]
	tra.R1 = 0.0970;
	tra.R2 = 0.07327;  //2����R[��]
	tra.L1 = 0.030549;
	tra.L2 = 0.030549;  //2�����ȃC���_�N�^���X[H]
	tra.M_IM = 0.029514;
	tra.g0 = 0.001614;

	tra.Cfc = C_FC * ini.num_M;
	tra.Lfl = L_FL / ini.num_M;
	tra.Rfl = R_FL / ini.num_M;

	tra.theta = 1.0;
	tra.theta_old = 1.0;


	tra.Emot_pow = 0.0;			//���[�^�͍s�G�l���M�[[J]
	tra.Emot_reg = 0.0;			//���[�^�񐶃G�l���M�[[J]
	tra.Emot_air = 0.0;			//�@�B�u���[�L�G�l���M�[[J]

	tra.Emot_pow_c = 0.0;		//���[�^�͍s�G�l���M�[[kWh]
	tra.Emot_reg_c = 0.0;		//���[�^�񐶃G�l���M�[[kWh]
	tra.Emot_air_c = 0.0;		//�@�B�u���[�L�G�l���M�[[kWh]

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



	tra.Evh = 0.0;				//��ԓ��o�̓G�l���M�[[J]
	tra.Evh_in = 0.0;			//��ԓ��̓G�l���M�[[J]
	tra.Evh_out = 0.0;			//��ԏo�̓G�l���M�[[J]
	tra.Evh_st = 0.0;

	tra.Evh_c = 0.0;			//��ԓ��o�̓G�l���M�[[kWh]
	tra.Evh_in_c = 0.0;			//��ԓ��̓G�l���M�[[kWh]
	tra.Evh_out_c = 0.0;		//��ԏo�̓G�l���M�[[kWh]
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
	//	//�G���[����
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
	sub.Jss_ini = -ini.vss / ((ini.e_ss * ini.vss) / (100.0 * ini.Iss));	//Jss�̏����l
	sub.Jss = -ini.vss / ((ini.e_ss * ini.vss) / (100.0 * ini.Iss));		//�X�V�p��Jss

	sub.Ess = 0.0;			//�ϓd���o�̓G�l���M�[[J]
	sub.Ess_c = 0.0;		//�ϓd���o�̓G�l���M�[[kWh]
	sub.Wss = 0.0;
	sub.Wss_c = 0.0;
}

void SimulateEngine::Make_NODE_TRAIN(NODE& node, TRAIN& tra)
{
	//*node = (NODE*)malloc(sizeof(NODE));
	//if (*node == NULL)
	//{
	//	printf("MEMORY?? What is that????\n");
	//	exit(EXIT_FAILURE);//�G���[����
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
	//	exit(EXIT_FAILURE);//�G���[����
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

	if (sub.iss < -0.1)		//�ϓd���ɓd�����������Ă����flag=0�i�G���[�j
	{
		sub.flag = 0;
		sub.Jss = -sub.vout / sub.Rss;
	}
	else
	{
		sub.flag = 1;			//�ϓd���ɓd�����������Ă��Ȃ���flag=1�iOK�j
	}

	/*	if (sub.vout >(Vcmax))	//�d���̐���
		{
			sub.flag = 0;
			sub.Jss = -(Vcmax) / sub.Rss;
		}*/
}

void SimulateEngine::Error_Detection_REG(TRAIN& tra) {

	/********** �y���׉񐶐���p�^�[���𖞂������ǂ����̔��� **********/
	if (tra.vfc > Vclim && tra.vfc < Vcmax)						//�񐶍i�荞�ݒ�
	{
		tra.b_reg_judge = tra.Fmot * (Vcmax - tra.vfc) / (Vcmax - Vclim);
	}
	else if (tra.vfc >= Vcmax)										//�񐶍i�荞�ݏI��
	{
		tra.b_reg_judge = 0.0;
	}
	else															//�񐶍i�荞�݂Ȃ�
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

void SimulateEngine::Calculate_R_run(TRAIN& tra)		//���s��R�̎Z�o
{
	tra.R_run = (tra.mass_M / 1000.0 * (1.65 + 0.0247 * tra.v_c) + tra.mass_T / 1000.0 * (0.78 + 0.0028 * tra.v_c) + (0.028 + 0.0078 * (tra.num - 1)) * tra.v_c * tra.v_c) * g;		//rr=g�~W�~(a+bv+cv^2/W)
	//	tra.R_run = g*(1.32+0.0614 * tra.v_c)*tra.mass_c+(0.028 + 0.0078 * (tra.num - 1)) * tra.v_c * tra.v_c;
}

void SimulateEngine::Calculate_R_total(TRAIN& tra)		//���v��R�̎Z�o
{
	Calculate_R_run(tra);

	tra.R_total = tra.R_run;
}

void SimulateEngine::Traction_force(TRAIN& tra)		//�͍s�������͂̎Z�o
{

	double v1_pp = tra.v1_P * tra.vfc / 1350.0;
	double v2_pp = tra.v2_P * tra.vfc / 1350.0;

	if (tra.Type == 2) tra.Fmax = 14.776 * 1000.0 * (tra.notch / 4.0);

	if (tra.v < v1_pp)		//��g���N�̈�
	{
		tra.Fmot = tra.Fmax * tra.num_mot;
	}
	else if (tra.v < v2_pp)		//��d�͗̈�
	{
		tra.Fmot = tra.Fmax * tra.num_mot * v1_pp / tra.v;
	}
	else     //�����̈�
	{
		tra.Fmot = tra.Fmax * tra.num_mot * v1_pp * v2_pp / (tra.v * tra.v);	//�����̈�
	}

	if (tra.Fmot > tra.Fmax * tra.num_mot) tra.Fmot = tra.Fmax * tra.num_mot;

	tra.Ftot = tra.Fmot;
}

void SimulateEngine::Regenerative_force(TRAIN& tra)		//�񐶎������͂̎Z�o
{
	double v1_rr = tra.v1_R * tra.vfc / 1650.0;
	double v2_rr = tra.v2_R * tra.vfc / 1650.0;

	if (tra.Type == 2) {
		tra.Bmax = 13.910 * 1000.0 * tra.brake / 5.0;
		if (tra.brake > 5.0) {
			tra.Bmax = 15.092 * 1000.0;
		}
	}

	if (tra.v < v1_rr)		//��g���N�̈�
	{
		tra.Fmot_regM = -tra.Bmax * tra.num_mot;
	}
	else if (tra.v < v2_rr)		//��d�͗̈�
	{
		tra.Fmot_regM = -tra.Bmax * tra.num_mot * v1_rr / tra.v;
	}
	else	//�����̈�
	{
		tra.Fmot_regM = -tra.Bmax * tra.num_mot * v1_rr * v2_rr / (tra.v * tra.v);
	}

	tra.Fmot = tra.Fmot_regM;
	if (tra.Fmot < -tra.BEmax) tra.Fmot = -tra.BEmax;
	tra.Ftot = -tra.BEmax;
}

/*** �ꎟ�x��t�B���^ ***/
double SimulateEngine::funcdelay(double aTs, double aTf, double ayold, double au)
{
	double fudely;
	fudely = (aTs * au + aTf * ayold) / (aTf + aTs);
	return(fudely);
}

void SimulateEngine::d_current_power(TRAIN& tra)//�͍s����d���d���Z�o
{
	double v0_PP = tra.v0_P * tra.vfc / 1350.0;
	tra.id_mot_P = tra.i1d_P;

	if (tra.v > v0_PP)		//��ߎ����̈�
	{
		tra.id_mot_P = tra.i1d_P * v0_PP / tra.v;
	}

	tra.idtot = tra.id_mot_P;
}

void SimulateEngine::d_current_regenerative(TRAIN& tra)//�񐶎���d���d���Z�o
{
	double v0_RR = tra.v0_R * tra.vfc / 1650.0;
	tra.id_mot_R = tra.i1d_R;

	if (tra.v > v0_RR)		//��ߎ����̈�
	{
		tra.id_mot_R = tra.i1d_R * v0_RR / tra.v;
	}
	tra.idtot = tra.id_mot_R;
}

void SimulateEngine::q_current_power(TRAIN& tra)//�͍s����q���d���Z�o
{
	tra.iq_mot_P = (tra.Fmot * tra.rd / tra.Gr) * tra.L2 / tra.POLE / (tra.M_IM * tra.M_IM) / tra.idtot / tra.num_mot;
	tra.iqtot = tra.iq_mot_P;
	tra.iqtotF = funcdelay(tra.Ts, tra.Tf, tra.iqtotF, tra.iqtot);
}

void SimulateEngine::q_current_regenerative(TRAIN& tra)//�񐶎���q���d���Z�o
{
	tra.iq_mot_R = (tra.Fmot * tra.rd / tra.Gr) * tra.L2 / tra.POLE / (tra.M_IM * tra.M_IM) / tra.idtot / tra.num_mot;
	tra.iqtot = tra.iq_mot_R;
	tra.iqtotF = funcdelay(tra.Ts, tra.Tf, tra.iqtotF, tra.iqtot);
}

void SimulateEngine::d_q_voltage_calculation(TRAIN& tra) {
	if (tra.accelflag == 2) {
		//�p���g���v�Z
		tra.ws = 0.0;
		tra.we = tra.ws + tra.wr;

		/*d���d��*/
		tra.v1d = 0.0;

		/*q���d��*/
		tra.v1q = 0.0;
		tra.Em = 0.0;
	}
	else if (tra.accelflag == 4) {

		tra.ws = 0.0;
		tra.we = 0.0;
		/*d���d��*/
		tra.v1d = 0.0;

		/*q���d��*/
		tra.v1q = 0.0;
		tra.Em = 0.0;
	}
	else {
		tra.ws = tra.R2 / tra.L2 * tra.iqtot / tra.idtot;
		tra.we = tra.ws + tra.wr;

		/*d���d��*/
		tra.v1d = tra.R1 * tra.idtot - (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.iqtotF;
		//		tra.v1d = - (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.iqtot;

				/*q���d��*/
		tra.v1q = tra.R1 * tra.iqtotF + (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.idtot + (tra.M_IM / tra.L2) * tra.M_IM * tra.we * tra.idtot;
		//		tra.v1q = (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.idtot + (tra.M_IM / tra.L2) * tra.M_IM * tra.we * tra.idtot;

				/*�C���o�[�^���͓d��*/
		tra.Em = sqrt(tra.v1d * tra.v1d + tra.v1q * tra.v1q);
	}
}

void SimulateEngine::Solve_next_state(TRAIN& tra)
{
	tra.a = (tra.Ftot - tra.R_total) / (tra.mass);		//�^��������ma=F
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
		if (tra.direction == 1.0)			//�r�܁ˏa�J
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
			if (i == 0)						//�a�J�˒r��
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
		if (tra.direction == 1.0)			//��{�ː������
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
			if (i == 0)						//����ˑ�{����
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

//������ς�����by Kasai1012
void SimulateEngine::Run_pattern(TRAIN& tra, const STATION* sta, double t, std::vector<double> wait_time)
{
	int i = 0;

	tra.t_stop = t;

	/*** ���̉w��T�� ***/
	tra.X_nextstop_old = tra.X_nextstop;
	if (tra.Type == 0 || tra.Type == 9 || tra.Type == 10 || tra.Type == 11 || tra.Type == 12 || tra.Type == 13 || tra.Type == 14) {
		for (size_t i = 0; i < NUM_station + NUM_final_station; i++)
		{
			if (tra.direction == 1.0)			//�r�܁ˏa�J
			{
				if (sta[i].Xs <= tra.x && sta[i + 1].Xs >= tra.x)
				{
					tra.X_nextstop = sta[i + 1].Xs;
					break;
				}
			}
			else								//�a�J�˒r��
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

	/*���x�����ݒ�*/
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


	/**** ���s��R�v�Z ****/
	Calculate_R_total(tra);

	tra.Bref = 1.111 * tra.brake / 7.0;



	/*** ���s���[�h������ ***/
	if (tra.accelflag == 4)		/*��Ԓ��̓��������*/
	{
		if (tra.laststopflag == 1)
		{
			tra.accelflag = 4;
		}
		else if (t < tra.T_delay)
		{
			tra.accelflag = 4;
		}
		else if (tra.direction == 1) {	//�����i�r�܁��a�J�j
			if (sta[2].Xs - 200.0 <= tra.x && tra.x < sta[2].Xs + 200.0) {	//�ڔ�
				// if (tra.t_stop >= tra.T_delay + 2.0 * 60.0 && tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;�@�@���Ԏ�����ݒ肷��ꍇ
				// if (tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
				if (tra.Type == 0) {
					if (tra.t_stop - tra.t_stop_old >= (stoptime + wait_time[0])) tra.accelflag = 1;
				}

			}
			else if (sta[3].Xs - 200.0 <= tra.x && tra.x < sta[3].Xs + 200.0) {	//���c�n��
				// if (tra.t_stop >= tra.T_delay + 4.0 * 60.0 + 20 && tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
				// if (tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
				if (tra.Type == 0) {
					if (tra.t_stop - tra.t_stop_old >= (stoptime - wait_time[0] + wait_time[1])) tra.accelflag = 1;
				}

			}
			else if (sta[4].Xs - 200.0 <= tra.x && tra.x < sta[4].Xs + 200.0) {	//�V��v��
				// if (tra.t_stop >= tra.T_delay + 6.0 * 60.0 && tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
				// if (tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
				if (tra.Type == 0) {
					if (tra.t_stop - tra.t_stop_old >= (stoptime - wait_time[1] + wait_time[2])) tra.accelflag = 1;
				}

			}
			else if (sta[5].Xs - 200.0 <= tra.x && tra.x < sta[5].Xs + 200.0) {	//�V�h
				// if (tra.t_stop >= tra.T_delay + 9.0 * 60.0 && tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
				if (tra.Type == 0) {
					if (tra.t_stop >= tra.T_delay + 9.0 * 60.0 && tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
				}

			}
			else if (sta[6].Xs - 200.0 <= tra.x && tra.x < sta[6].Xs + 200.0) {	//��X��
				// if (tra.t_stop >= tra.T_delay + 10.0 * 60.0 && tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
				// if (tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
				if (tra.Type == 0) {
					if (tra.t_stop - tra.t_stop_old >= (stoptime + wait_time[3])) tra.accelflag = 1;
				}

			}
			else if (sta[7].Xs - 200.0 <= tra.x && tra.x < sta[7].Xs + 200.0) {	//���h
				// if (tra.t_stop >= tra.T_delay + 13.0 * 60.0 && tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1; 
				// if (tra.t_stop - tra.t_stop_old >= stoptime) tra.accelflag = 1;
				if (tra.Type == 0) {
					if (tra.t_stop - tra.t_stop_old >= (stoptime - wait_time[3] + wait_time[4])) tra.accelflag = 1;
				}

			}
			else if (sta[8].Xs - 200.0 <= tra.x && tra.x < sta[8].Xs + 200.0) {	//�a�J
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
	else if (tra.accelflag == 2)		/*�čs���̓��������*/
	{
		if (tra.direction == 1.0)		//��{�ː������
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
		else							//��{�ː������
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
	else		/*������or�������̓��������*/
	{
		if (tra.direction == 1.0)								//��{�ː������
		{
			if (tra.brakeflag == 1)		/*�������̓���*/
			{
				if (tra.v_c >= 0.0)
				{
					tra.accelflag = 3;		//�������
				}
				else
				{
					tra.accelflag = 4;		//���
					tra.brakeflag = 0;
					tra.t_stop_old = tra.t_stop;
				}
			}
			else
			{
				if (tra.laststopflag == 1)	//�w�̎�O�ɒ�Ԃ����ꍇ�̋~��
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
					tra.accelflag = 1;		//�������
				}
				else
				{
					tra.accelflag = 2;		//�čs���
				}
			}
		}
		else													//����ˑ�{����
		{
			if (tra.brakeflag == 1)		/*�������̓���*/
			{
				if (tra.v_c >= 0.0)
				{
					tra.accelflag = 3;		//�������
				}
				else
				{
					tra.accelflag = 4;		//���
					tra.brakeflag = 0;
					tra.t_stop_old = tra.t_stop;
				}
			}
			else
			{
				if (tra.laststopflag == 1)	//�w�̎�O�ɒ�Ԃ����ꍇ�̋~��
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
					tra.accelflag = 1;		//�������
				}
				else
				{
					tra.accelflag = 2;		//�čs���
				}
			}
		}
	}

}

void SimulateEngine::Calculation_traction_force_brake_force(TRAIN& tra) {

	/*** �e���s���[�h�ɉ�����������w�� ***/
	if (tra.accelflag == 1)				//������
	{
		/** �u���[�L�͌v�Z(�������Ȃ̂�0) **/
		tra.b_reg = 0.0;
		tra.b_air = 0.0;
		tra.b_reg_loss = 0.0;

		/** �ԗ������͌v�Z **/
		Traction_force(tra);

		/**�͍s��d���d���v�Z**/
		d_current_power(tra);

		/**�͍s��q���d���v�Z**/
		q_current_power(tra);

		tra.iq_reg = 0.0;

		tra.ws = tra.R2 / tra.L2 * tra.iqtot / tra.idtot;
		tra.we = tra.ws + tra.wr;

		/*d���d��*/
		tra.v1d = tra.R1 * tra.idtot - (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.iqtotF;

		/*q���d��*/
		tra.v1q = tra.R1 * tra.iqtotF + (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.idtot + (tra.M_IM / tra.L2) * tra.M_IM * tra.we * tra.idtot;
		/*�C���o�[�^���͓d��*/
		tra.Em = sqrt(tra.v1d * tra.v1d + tra.v1q * tra.v1q);
	}
	else if (tra.accelflag == 2)			//�čs��
	{
		/** �u���[�L�͌v�Z(�čs���Ȃ̂�0) **/
		tra.b_reg = 0.0;
		tra.b_air = 0.0;
		tra.b_reg_loss = 0.0;

		/** �ԗ������͌v�Z **/
		tra.Fmot = 0.0;
		tra.Ftot = 0.0;

		/*�d���v�Z*/
		tra.idtot = 0.0;
		tra.iqtot = 0.0;
		tra.iqtotF = 0.0;
		tra.iq_reg = 0.0;

		//�p���g���v�Z
		tra.ws = 0.0;
		tra.we = tra.ws + tra.wr;

		/*d���d��*/
		tra.v1d = 0.0;

		/*q���d��*/
		tra.v1q = 0.0;
		tra.Em = 0.0;

	}
	else if (tra.accelflag == 3)			//������
	{
		/** �ԗ��u���[�L�͌v�Z **/
		Regenerative_force(tra);
		d_current_regenerative(tra);
		q_current_regenerative(tra);

		/********** �y���׉񐶒��̃g���N�i�荞�ݖ͋[���i�y���׉񐶐��䉺�̓d�C�u���[�L�͂��v�Z�j **********/
		if (tra.vfc > Vclim && tra.vfc < Vcmax)						//�񐶍i�荞�ݒ�
		{
			tra.b_reg = tra.Fmot * (Vcmax - tra.vfc) / (Vcmax - Vclim);
			tra.iq_reg = (tra.b_reg * tra.rd / tra.Gr / tra.num_mot) * tra.L2 / (tra.POLE * tra.M_IM * tra.M_IM * tra.idtot);
			tra.iqtot = tra.iq_reg;
		}
		else if (tra.vfc >= Vcmax)										//�񐶍i�荞�ݏI��
		{
			tra.b_reg = 0.0;
			tra.iq_reg = 0.0;
			tra.iqtot = 0.0;
		}
		else															//�񐶍i�荞�݂Ȃ�
		{
			tra.b_reg = tra.Fmot;
			tra.iq_reg = tra.iq_mot_R;
			tra.iqtot = tra.iq_reg;
		}

		tra.iq_regF = funcdelay(tra.Ts, tra.Tf, tra.iq_regF, tra.iq_reg);
		tra.iqtotF = tra.iq_regF;

		tra.b_reg_loss = tra.Fmot - tra.b_reg;		//�񐶍i�荞�ݗ�[N]

		tra.b_air = -tra.BEmax - tra.b_reg;			//����Ȃ��u���[�L�͂͋�C�u���[�L�ŕ₤

		tra.ws = tra.R2 / tra.L2 * tra.iqtot / tra.idtot;
		tra.we = tra.ws + tra.wr;

		/*d���d��*/
		tra.v1d = tra.R1 * tra.idtot - (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.iqtotF;

		/*q���d��*/
		tra.v1q = tra.R1 * tra.iqtotF + (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.idtot + (tra.M_IM / tra.L2) * tra.M_IM * tra.we * tra.idtot;
		tra.Em = sqrt(tra.v1q * tra.v1q + tra.v1d * tra.v1d);
	}
	else if (tra.accelflag == 4)						//��Ԓ�
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
		/*d���d��*/
		tra.v1d = 0.0;

		/*q���d��*/
		tra.v1q = 0.0;
		tra.Em = sqrt(tra.v1d * tra.v1d + tra.v1q * tra.v1q);
	}
	else if (tra.accelflag == 5) {
		if (tra.v_c > tra.con_speed) {
			tra.Fmot = -tra.Bmax * tra.num_mot * (tra.v_c - tra.con_speed) / 5.0;
			tra.Ftot = tra.Fmot;

			d_current_regenerative(tra);
			q_current_regenerative(tra);

			/********** �y���׉񐶒��̃g���N�i�荞�ݖ͋[���i�y���׉񐶐��䉺�̓d�C�u���[�L�͂��v�Z�j **********/
			if (tra.vfc > Vclim && tra.vfc < Vcmax)						//�񐶍i�荞�ݒ�
			{
				tra.b_reg = tra.Fmot * (Vcmax - tra.vfc) / (Vcmax - Vclim);
				tra.iq_reg = (tra.b_reg * tra.rd / tra.Gr / tra.num_mot) * tra.L2 / (tra.POLE * tra.M_IM * tra.M_IM * tra.idtot);
				tra.iqtot = tra.iq_reg;
			}
			else if (tra.vfc >= Vcmax)										//�񐶍i�荞�ݏI��
			{
				tra.b_reg = 0.0;
				tra.iq_reg = 0.0;
				tra.iqtot = 0.0;
			}
			else															//�񐶍i�荞�݂Ȃ�
			{
				tra.b_reg = tra.Fmot;
				tra.iq_reg = tra.iq_mot_R;
				tra.iqtot = tra.iq_reg;
			}

			tra.iq_regF = funcdelay(tra.Ts, tra.Tf, tra.iq_regF, tra.iq_reg);
			tra.iqtotF = tra.iq_regF;

			tra.b_reg_loss = tra.Fmot - tra.b_reg;		//�񐶍i�荞�ݗ�[N]

			tra.b_air = tra.Ftot - tra.b_reg;			//����Ȃ��u���[�L�͂͋�C�u���[�L�ŕ₤
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

		/*d���d��*/
		tra.v1d = tra.R1 * tra.idtot - (1.0 - tra.M_IM * tra.M_IM / tra.L1 / tra.L2) * tra.L1 * tra.we * tra.iqtotF;

		/*q���d��*/
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

	tra.Ploss_mot = (tra.R1 * (tra.idtot * tra.idtot + tra.iqtotF * tra.iqtotF) + tra.R2 * (tra.M_IM / tra.L2 * tra.iqtotF) * (tra.M_IM / tra.L2 * tra.iqtotF)) * tra.num_mot;//���[�^�̓���
	tra.Ploss_fe = tra.g0 * (tra.v1d * tra.v1d + tra.v1q * tra.v1q) * tra.num_mot;
	tra.Ploss_fl = tra.Rfl * tra.ifl * tra.ifl;					//�t�B���^���A�N�g���̃W���[������[W]
	tra.Ploss_fl_c = tra.Ploss_fl / 1000.0;					//�t�B���^���A�N�g���̃W���[������[kW]

	tra.Ploss_mot_c = (tra.Ploss_mot + tra.Ploss_fe) / 1000.0;				//���[�^�̓���[kW]

	tra.Ploss_all = tra.Ploss_mot + tra.Ploss_fe + tra.Ploss_fl;
	tra.Ploss_all_c = tra.Ploss_all / 1000.0;
}

void SimulateEngine::Calculation_power_motor(TRAIN& tra)
{
	tra.Ptot = tra.v * tra.Ftot;

	if (tra.accelflag == 1)				/*�͍s��*/
	{
		tra.Pmot = (tra.v1d * tra.idtot + tra.v1q * tra.iqtotF) * tra.num_mot;
		tra.Pmot_air = 0.0;
		tra.Pmot_reg_loss = 0.0;

	}
	else if (tra.accelflag == 3)		/*�񐶒�*/
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

	tra.Pres = tra.v * tra.R_total;		//���s��R�ł̃��X

	tra.Ptot_c = tra.Ptot / 1000.0;
	tra.Pmot_c = tra.Pmot / 1000.0;
	tra.Pmot_air_c = tra.Pmot_air / 1000.0;
	tra.Pres_c = tra.Pres / 1000.0;
	tra.Pmot_reg_loss_c = tra.Pmot_reg_loss / 1000.0;

}

void SimulateEngine::Calculation_power_train(TRAIN& tra)
{
	tra.Pvh = tra.vp * tra.ifl;

	if (tra.accelflag == 1 || tra.accelflag == 2)				//�͍s��,�čs��
	{
		tra.Pvh_in = tra.vp * tra.ifl;
		tra.Pvh_out = 0.0;
		tra.Pvh_st = 0.0;
	}
	else if (tra.accelflag == 3)		//�񐶒�
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
	if (tra.accelflag == 1)				//�͍s��
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
	else if (tra.accelflag == 3)		//�񐶒�
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

/*** x�����y���w���v�f������ ***/
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


/*** �z��data[]�̐擪n�̗v�f�������̏����Ƀ\�[�g ***/
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
	//	//�G���[����
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
	size_t i = 0, j = 0;						//j�F�������ɐڑ�����Ă���m�[�h�T���p�Ck�F�������ɐڑ�����Ă���m�[�h�T���p
	int count_SS, count_UP, count_DOWN, NUM_tra_UP;	//�u�����`�쐬�p�J�E���g�ϐ�

	count_SS = 0;
	count_UP = NUM_sub;
	NUM_tra_UP = Count_Trains_Direction(tra, 1);

	//printf("%d\n", NUM_tra_UP);

	count_DOWN = 2 * NUM_sub + NUM_tra_UP - 1;

	/*�@SS�̃u�����`���쐬*/
	//printf("-----------SS--------------\n");
	for (i = 0; i < N_node; i++)
	{
		//printf("[%d]flag = %d\n", temp[i].Number, temp[i].flag);
		if (node_order[i].flag == 0)				//In case of SS( i Starts with 0 to NUM_sub-1 )
		{
			//puts("SS");
			branch[count_SS].X_start = node_order[i].X;		//�ʒu[m]�i��_����̋����j
			branch[count_SS].X_end = node_order[i].X;		//�ʒu[m]�i��_����̋����j
			branch[count_SS].r = node_order[i].r;			//Resistance [Ohm]
			branch[count_SS].flag = node_order[i].flag;		//�u�����`�̕��ށi1�FUP, -1�FDOWN, 0�FSS�j
			branch[count_SS].Node_pos = -1;			//���������ɐڑ�����Ă���m�[�h��
			branch[count_SS].Node_neg = -1;			//���������ɐڑ�����Ă���m�[�h��
			count_SS = count_SS + 1;
		}
	}

	/*�A���̃u�����`���쐬*/
	//printf("-----------UP--------------\n");
	for (i = 0; i < N_node; i++)
	{
		if (node_order[i].flag >= 0)			//In case of +1 direction
		{
			//printf("[%d]flag = %d\n", temp[i].Number, temp[i].flag);
			if (count_UP < NUM_tra_UP + 2 * NUM_sub - 1)		//�I�[�ȊO
			{

				for (size_t j = 0; node_order[i + (j + 1)].flag == -1; j++) {}				//temp[i]�ɑ΂��Đ��������ɐڑ�����Ă���m�[�h��T��
				//printf("\nj = %d\n", j);

				/*temp[i]�ɑ΂��Đ������ɐڑ�����Ă���u�����`B���쐬*/
				branch[count_UP].Node_pos = node_order[i + (j + 1)].Number;			//�u�����`B�ɂ��āC���������ɐڑ�����Ă���m�[�h��
				branch[count_UP].Node_neg = node_order[i].Number;					//�u�����`B�ɂ��āC���������ɐڑ�����Ă���m�[�h��

				branch[count_UP].X_start = node_order[i].X;							//[m]
				branch[count_UP].X_end = node_order[i + (j + 1)].X;					//[m]

				if (node_order[i + (j + 1)].X <= 7500.0) {
					branch[count_UP].r = R_StoN_UP * std::abs(node_order[i + (j + 1)].X - node_order[i].X) / 1000.0;		//Resistance [Ohm]
				}
				else {
					branch[count_UP].r = R_NtoK_UP * std::abs(node_order[i + (j + 1)].X - node_order[i].X) / 1000.0;		//Resistance [Ohm]
				}

				branch[count_UP].flag = 1;										//�u�����`�̕��ށi1�FUP, -1�FDOWN, 0�FSS�j

				/*���̃��[�v�̂��߂̏���*/
				j = 0;
				count_UP = count_UP + 1;
			}
			else {}		//�I�[�̃m�[�h�ł͉������Ȃ�
		}
	}


	/*�B����̃u�����`���쐬*/
	//printf("-----------DOWN--------------\n");
	for (i = 0; i < N_node; i++)
	{
		if (node_order[i].flag <= 0)			//In case of +1 direction
		{
			//printf("[%d]flag = %d\n", temp[i].Number, temp[i].flag);
			if (count_DOWN < N_branch)		//�I�[�ȊO
			{
				for (size_t j = 0; node_order[i + (j + 1)].flag == 1; j++) {}				//temp[i]�ɑ΂��Đ��������ɐڑ�����Ă���m�[�h��T��
				//printf("\nj = %d\n", j);

				/*�������ɐڑ�����Ă���u�����`���쐬*/
				branch[count_DOWN].Node_pos = node_order[i + (j + 1)].Number;		//���������ɐڑ�����Ă���m�[�h��
				branch[count_DOWN].Node_neg = node_order[i].Number;					//���������ɐڑ�����Ă���m�[�h��

				branch[count_DOWN].X_start = node_order[i].X;						//[m]
				branch[count_DOWN].X_end = node_order[i + (j + 1)].X;				//[m]

				if (node_order[i + (j + 1)].X <= 7500.0) {
					branch[count_DOWN].r = R_StoN_DOWN * std::abs(node_order[i + (j + 1)].X - node_order[i].X) / 1000.0;		//Resistance [Ohm]
				}
				else {
					branch[count_DOWN].r = R_NtoK_DOWN * std::abs(node_order[i + (j + 1)].X - node_order[i].X) / 1000.0;		//Resistance [Ohm]
				}

				branch[count_DOWN].flag = -1;									//�u�����`�̕��ށi1�FUP, -1�FDOWN, 0�FSS�j

				/*���̃��[�v�̂��߂̏���*/
				j = 0;
				count_DOWN = count_DOWN + 1;
			}
			else {}		//�I�[�̃m�[�h�ł͉������Ȃ�
		}
	}
}

void SimulateEngine::Initialize_matrix1(const std::unique_ptr<double[]>& Mtr, size_t M, size_t N)		/*M�sN��̂Q�����z�������������֐�*/
{

	for (size_t i = 0; i < M; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			Mtr[i * N + j] = 0.0000;
		}
	}
}

void SimulateEngine::Make_H_matrix(const std::unique_ptr<BRANCH[]>& branch, const std::unique_ptr<double[]>& Hm)		/*H matrix�쐬�p�֐�*/
{
	//int i;	/*�s*/
	//int j;	/*��*/

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

void SimulateEngine::Make_Y_matrix(const std::unique_ptr<BRANCH[]>& branch, const std::unique_ptr<double[]>& Y)		/*Y matrix�쐬�p�֐�*/
{
	//int i;	/*�s*/
	//int j;	/*��*/

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
					Y[i + j * N_branch] = 1.0 / branch[j].r;		/*�����ŃR���_�N�^���X[S]�ɕϊ�*/
				}
			}
			else
			{
				Y[i + j * N_branch] = 0.0;
			}
		}
	}
}

void SimulateEngine::Make_In_vector(const std::unique_ptr<NODE[]>& data, const std::unique_ptr<double[]>& In)		/*In vector�쐬�p�֐�*/
{

	for (size_t i = 0; i < N_node; i++)
	{
		In[i] = data[i].I;
	}
}

void SimulateEngine::Transpose(const std::unique_ptr<double[]>& trans, const std::unique_ptr<double[]>& X, size_t row, size_t column)		/*�]�u�s��v�Z�p�֐�*/
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
	//tra.reset(); tra = std::make_unique<TRAIN[]>(NUM_tra);	//�|�C���^�̔z��Ƃ��Ē�`
	//sub.reset(); sub = std::make_unique<SUB[]>(NUM_sub);
	//node.reset(); node = std::make_unique<NODE[]>(N_node);
	//node_order.reset(); node_order = std::make_unique<NODE[]>(N_node);
	//branch.reset(); branch = std::make_unique<BRANCH[]>(N_branch);

	ERROR = 0;						//��H�v�Z�̍Čv�Z����p�t���O�iflag��S�ĉ��Z����NUM_sub�ȏ�Ȃ��H�v�Z���[�v�𔲂���BNUM_sub�ȉ��Ȃ�Čv�Z�B�j
	ERROR_REG = 0;
	count_loop = 0;
	count_loop_rec = 0;

	t = 0.0;

	sub_total = 0.0;

	/***  ��ԍ\���́E�ϓd���\���́EESD�\���̏�����  ***/
	for (size_t i = 0; i < NUM_sub; i++)
	{
		Make_substation(sub[i], ini_sub[i], i);				//�i���o�����O�͒ʂ��ԍ�
	}
	for (size_t i = 0; i < NUM_tra; i++)
	{
		Make_train(tra[i], ini_tra[i], i + NUM_sub);							//�i���o�����O�͒ʂ��ԍ�
	}

	for (size_t i = 0; i < NUM_tra; i++)
	{
		Calculate_BEmax(tra[i]);
	}

	////////////////////////// ���l�v�Z�� /////////////////////////
	for (t = 0.0; t <= 60.0 * 45.0; t += dt)	/***** 1���ԃV�~�����[�V���� *****/
	{
		minute = t / 60.0;
		count_loop_rec = 0;

		////////////////////////// �ԗ��^���������v�Z�� /////////////////////
			/*** �ԗ����s�p�^�[���͋[ & �^���������v�Z ***/
			/*** �ԗ��p���[�v�Z ***/

		for (size_t i = 0; i < NUM_tra; i++)
		{
			Calculate_BEmax(tra[i]);

			Run_pattern(tra[i], sta, t, wait_time);

		}

		/////////////////////////////////////////////////////////////////////

		while (count_loop <= 3000002)
		{
			// start_while:				//�_�C�I�[�h���肪�Ԉ���Ă����ꍇ�Ccontinue�ł����ɖ߂��Ă���

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

			/***  node[]�̏�����  ***/
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

			/*** �eMatrix�̏����� ***/
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

			/***  node_order[]���������Ƀ\�[�g  ***/
			//sort_s(node_order, N_node);
			sort(node_order, N_node);

			Make_branch(branch, node_order, tra);

			//H�s��̐���(�u�񐶎Ԃ��܂ނ��d�V�X�e���̌���Ƃ�����vp.19�Q��)
			Make_H_matrix(branch, H);

			//�A�h�~�^���X�s��̐���
			Make_Y_matrix(branch, Y);

			//�m�[�h�d���s��̐���
			Make_In_vector(node, In);
			Make_In_vector(node, In_cpy);		//"dgetrs"�ł͉E�Ӄx�N�g�����������̉��ŏ㏑������Ă��܂��̂ŁC�v�Z�p��In���R�s�[����In_cpy���쐬

			//H�s���]�u����
			Transpose(H_tra, H, N_node, N_branch);

			//H�s��~Y�s�� (dgemm��intel MKL��p���čs��̊|���Z���v�Z)
			dgemm(&trans, &trans, &K_node, &K_branch, &K_branch, &alpha, H.get(), &ld_node, Y.get(), &ld_branch, &beta, TEMP.get(), &ld_node);

			//��L�s��~�]�uH�s�� (dgemm��intel MKL��p���čs��̊|���Z���v�Z)
			dgemm(&trans, &trans, &K_node, &K_node, &K_branch, &alpha, TEMP.get(), &ld_node, H_tra.get(), &ld_branch, &beta, A.get(), &ld_node);

			//�A��1�����������������߂�H�s��~Y�s��~�]�uH�s��̉�A(N�sN��)��LU����
			dgetrf(&ld_node, &K_branch, A.get(), &ld_node, ipiv_Vn.get(), &info_Vn);

			//dgetrs�ŘA��1���������uA*Vn=In�v�������CIn�ɉ�(�m�[�h�d���s��)���㏑��
			dgetrs(&trans, &K_node, &nrhs, A.get(), &ld_node, ipiv_Vn.get(), In_cpy.get(), &ld_node, &info_Vn);	//In_cpy�ɕ������̉�Vn���i�[����Ă��邱�Ƃɒ���

			Make_Y_matrix(branch, Y);		//�Ȃ���Vn�����߂�ۂ�"dgetrf"��Y-matrix�������������Ă��邽�߁C�����ōČv�Z�����Ă�i���ۂ͖��ʂȌv�Z�Ȃ̂ł������Ȃ����c�j

			for (size_t i = 0; i < N_node; i++)	/*Vn�x�N�g���̓}�C�i�X������*/
			{
				Vn[i] = -In_cpy[i];
			}

			//�]�uH�s��~Vn�x�N�g���Ńu�����`�d���x�N�g�����v�Z
			dgemv(&trans, &K_branch, &K_node, &alpha, H_tra.get(), &ld_branch, Vn.get(), &incx, &beta, Vb.get(), &incy);

			//Y�s��~�u�����`�d���x�N�g���Ńu�����`�d���x�N�g�����v�Z
			dgemv(&trans, &K_branch, &K_branch, &alpha, Y.get(), &ld_branch, Vb.get(), &incx, &beta, Ib.get(), &incy);


			for (size_t i = 0; i < N_node; i++)			//�m�[�h�d���̌v�Z���ʂ��m�[�h�\���̂֕Ԃ��i���ʓI�ɂ́CVn�z��̓m�[�h�\���̂̏��Ԓʂ�ɕ���ł���H�j
			{
				node[i].V = Vn[i];
			}

			for (size_t i = 0; i < NUM_sub; i++)			//�m�[�h�d���E�d���ƃu�����`�d������C�ϓd���o�͓d���E�d�����v�Z
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
					if (node[j].Number == tra[i].name_train)	//�m�[�h�ԍ��Ɨ�Ԏ��ʔԍ�����v���Ă�����m�[�h�d�����ԃp���^�_�d���ɕԂ�
					{
						tra[i].vp = node[j].V;
						tra[i].vfc = tra[i].vp - tra[i].Rfl * tra[i].ifl;

						//�����́E�u���[�L�́E���[�^�d���E���[�^�d�����v�Z
						Calculation_traction_force_brake_force(tra[i]);

						//dq���d���E�d�����烂�[�^�p���[���v�Z
						Calculation_power_motor(tra[i]);

						//���H�v�Z
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
				//goto start_while;					//while���̐擪�ɖ߂�
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

		/***  ���̏�Ԃ��v�Z  ***/
		for (size_t i = 0; i < NUM_tra; i++)
		{
			if (tra[i].accelflag != 4) Solve_next_state(tra[i]);
		}

		/***  �ԗ��G�l���M�[�v�Z  ***/
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

		/////////////////////////////  cmd�\����  ////////////////////////////

		//////////////////////////////////////////////////////////////////////
	}
	/********** ���l�v�Z���[�v�I�� **********/

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


	/*****  �������J��  *****/
	// �X�}�[�g�|�C���^�Ȃ̂ŕs�v
	/*for (size_t i = 0; i < NUM_tra; i++)
	{
		free(tra[i]);
	}
	for (size_t i = 0; i < NUM_sub; i++)
	{
		free(sub[i]);
	}*/

	//std::unique_ptr<TRAIN[]> tra = std::make_unique<TRAIN[]>(NUM_tra);	//�|�C���^�̔z��Ƃ��Ē�`
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
